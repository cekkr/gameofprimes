// optimized_simulation_gpu.js

// Cache for GPU resources to avoid recreating them
let gpuCache = {
    device: null,
    adapter: null,
    computePipeline: null,
    bindGroupLayout: null,
    buffers: new Map(), // Track buffers by size/purpose
    inUse: false
};

async function cleanupGPUResources() {
    if (gpuCache.device) {
        // Destroy all cached buffers
        for (const [key, buffer] of gpuCache.buffers.entries()) {
            buffer.destroy();
        }
        gpuCache.buffers.clear();
        
        // Explicitly mark device for destruction
        if (gpuCache.device.destroy) {
            gpuCache.device.destroy();
        }
        
        // Reset cache
        gpuCache = {
            device: null,
            adapter: null,
            computePipeline: null,
            bindGroupLayout: null,
            buffers: new Map(),
            inUse: false
        };
    }
}

async function calculateForcesGPU(positions, masses, charges, rules) {
    if (gpuCache.inUse) {
        console.warn("GPU calculation already in progress. Skipping this frame.");
        return null;
    }
    
    gpuCache.inUse = true;
    
    try {
        if (!navigator.gpu) {
            console.warn("WebGPU is not supported. Falling back to CPU.");
            gpuCache.inUse = false;
            return null;
        }

        // Initialize GPU device if not already initialized
        if (!gpuCache.device) {
            const adapter = await navigator.gpu.requestAdapter();
            if (!adapter) {
                console.warn("Failed to get GPU adapter. Falling back to CPU.");
                gpuCache.inUse = false;
                return null;
            }
            
            gpuCache.adapter = adapter;
            gpuCache.device = await adapter.requestDevice();
            
            // Create bind group layout - only once
            gpuCache.bindGroupLayout = gpuCache.device.createBindGroupLayout({
                entries: [
                    { binding: 0, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
                    { binding: 1, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
                    { binding: 2, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
                    { binding: 3, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } },
                ],
            });
            
            // Create compute pipeline - only once
            gpuCache.computePipeline = gpuCache.device.createComputePipeline({
                layout: gpuCache.device.createPipelineLayout({ 
                    bindGroupLayouts: [gpuCache.bindGroupLayout] 
                }),
                compute: {
                    module: gpuCache.device.createShaderModule({
                        code: `
                            @group(0) @binding(0) var<storage, read> positions: array<vec3<f32>>;
                            @group(0) @binding(1) var<storage, read> masses: array<f32>;
                            @group(0) @binding(2) var<storage, read> charges: array<f32>;
                            @group(0) @binding(3) var<storage, read_write> forces: array<vec3<f32>>;

                            @compute @workgroup_size(64)
                            fn main(@builtin(global_invocation_id) global_id: vec3<u32>) {
                                let num_molecules = arrayLength(&positions);
                                if (global_id.x >= num_molecules) {
                                    return;
                                }

                                var force = vec3<f32>(0.0, 0.0, 0.0);
                                for (var i = 0u; i < num_molecules; i++) {
                                    if (i == global_id.x) {
                                        continue;
                                    }

                                    let direction = positions[i] - positions[global_id.x];
                                    let distance = length(direction);

                                    if (distance < 0.0001) {
                                        continue;
                                    }

                                    let direction_normalized = normalize(direction);

                                    // *Simplified* force calculation (example: gravity only)
                                    let G = 0.1; // You'd need to pass constants via a uniform buffer
                                    let force_magnitude = G * masses[global_id.x] * masses[i] / (distance * distance);
                                    force = force + force_magnitude * direction_normalized;
                                }
                                forces[global_id.x] = force;
                            }
                        `,
                    }),
                    entryPoint: 'main',
                },
            });
        }

        // Reuse or create buffers
        const getOrCreateBuffer = (key, size, usage, data = null) => {
            const bufferKey = `${key}_${size}`;
            
            // Check if we already have a buffer of this size
            if (gpuCache.buffers.has(bufferKey)) {
                const buffer = gpuCache.buffers.get(bufferKey);
                
                // Update buffer data if provided
                if (data) {
                    gpuCache.device.queue.writeBuffer(buffer, 0, data);
                }
                
                return buffer;
            }
            
            // Create new buffer
            const buffer = gpuCache.device.createBuffer({
                size,
                usage
            });
            
            // Store in cache
            gpuCache.buffers.set(bufferKey, buffer);
            
            // Write data if provided
            if (data) {
                gpuCache.device.queue.writeBuffer(buffer, 0, data);
            }
            
            return buffer;
        };

        // Get or create buffers with data
        const positionBuffer = getOrCreateBuffer(
            'position', 
            positions.byteLength,
            GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
            positions
        );
        
        const massBuffer = getOrCreateBuffer(
            'mass',
            masses.byteLength,
            GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
            masses
        );
        
        const chargeBuffer = getOrCreateBuffer(
            'charge',
            charges.byteLength,
            GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
            charges
        );
        
        const forceBuffer = getOrCreateBuffer(
            'force',
            positions.byteLength,
            GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC
        );
        
        const stagingBuffer = getOrCreateBuffer(
            'staging',
            positions.byteLength,
            GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST
        );

        // Create bind group
        const bindGroup = gpuCache.device.createBindGroup({
            layout: gpuCache.bindGroupLayout,
            entries: [
                { binding: 0, resource: { buffer: positionBuffer } },
                { binding: 1, resource: { buffer: massBuffer } },
                { binding: 2, resource: { buffer: chargeBuffer } },
                { binding: 3, resource: { buffer: forceBuffer } },
            ],
        });

        // Run compute pass
        const commandEncoder = gpuCache.device.createCommandEncoder();
        const passEncoder = commandEncoder.beginComputePass();
        passEncoder.setPipeline(gpuCache.computePipeline);
        passEncoder.setBindGroup(0, bindGroup);
        passEncoder.dispatchWorkgroups(Math.ceil(positions.length / (3 * 64)));
        passEncoder.end();

        // Copy results to staging buffer
        commandEncoder.copyBufferToBuffer(
            forceBuffer, 0, stagingBuffer, 0, stagingBuffer.size
        );

        // Submit commands
        gpuCache.device.queue.submit([commandEncoder.finish()]);

        // Read results
        await stagingBuffer.mapAsync(GPUMapMode.READ);
        const forces = new Float32Array(stagingBuffer.getMappedRange().slice());
        stagingBuffer.unmap();
        
        gpuCache.inUse = false;
        return forces;
    } catch (error) {
        console.error("GPU calculation error:", error);
        gpuCache.inUse = false;
        
        // Try to clean up on error
        try {
            await cleanupGPUResources();
        } catch (cleanupError) {
            console.warn("Error during GPU cleanup:", cleanupError);
        }
        
        return null;
    }
}

// Export functions
export { calculateForcesGPU, cleanupGPUResources };