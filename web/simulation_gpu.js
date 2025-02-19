// simulation_gpu.js

async function calculateForcesGPU(positions, masses, charges, rules) {
    if (!navigator.gpu) {
        console.warn("WebGPU is not supported. Falling back to CPU.");
        return null; // Indicate failure to use GPU
    }

    const adapter = await navigator.gpu.requestAdapter();
    if (!adapter) {
        console.warn("Failed to get GPU adapter. Falling back to CPU.");
        return null;
    }
    const device = await adapter.requestDevice();

    // --- Create Buffers ---

    const positionBuffer = device.createBuffer({
        size: positions.byteLength,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
        mappedAtCreation: true,
    });
    new Float32Array(positionBuffer.getMappedRange()).set(positions);
    positionBuffer.unmap();

    const massBuffer = device.createBuffer({
        size: masses.byteLength,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
        mappedAtCreation: true
    });

    new Float32Array(massBuffer.getMappedRange()).set(masses);
    massBuffer.unmap();

    const chargeBuffer = device.createBuffer({
      size: charges.byteLength,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
      mappedAtCreation: true,
    });
    new Float32Array(chargeBuffer.getMappedRange()).set(charges);
    chargeBuffer.unmap();



    const forceBuffer = device.createBuffer({
        size: positions.byteLength, // Same size as positions (3 floats per molecule)
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC,
    });

    // --- Create Bind Group Layout and Pipeline ---
    const bindGroupLayout = device.createBindGroupLayout({
        entries: [
            { binding: 0, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } }, // Positions
            { binding: 1, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } }, // Masses
            { binding: 2, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },  // Charges
            { binding: 3, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } },          // Forces (read-write)
        ],
    });

    const computePipeline = device.createComputePipeline({
        layout: device.createPipelineLayout({ bindGroupLayouts: [bindGroupLayout] }),
        compute: {
            module: device.createShaderModule({
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

    // --- Create Bind Group ---

    const bindGroup = device.createBindGroup({
        layout: bindGroupLayout,
        entries: [
            { binding: 0, resource: { buffer: positionBuffer } },
            { binding: 1, resource: { buffer: massBuffer } },
            { binding: 2, resource: { buffer: chargeBuffer } },
            { binding: 3, resource: { buffer: forceBuffer } },
        ],
    });

    // --- Run Compute Pass ---

    const commandEncoder = device.createCommandEncoder();
    const passEncoder = commandEncoder.beginComputePass();
    passEncoder.setPipeline(computePipeline);
    passEncoder.setBindGroup(0, bindGroup);
    passEncoder.dispatchWorkgroups(Math.ceil(positions.length / (3*64))); // Dispatch enough workgroups
    passEncoder.end();

    // --- Copy Results ---

    const stagingBuffer = device.createBuffer({
        size: forceBuffer.size,
        usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
    });
    commandEncoder.copyBufferToBuffer(forceBuffer, 0, stagingBuffer, 0, stagingBuffer.size);

    device.queue.submit([commandEncoder.finish()]);

    // --- Read Results ---

    await stagingBuffer.mapAsync(GPUMapMode.READ);
    const forces = new Float32Array(stagingBuffer.getMappedRange().slice()); // Copy the data
    stagingBuffer.unmap();

    return forces;
}

export { calculateForcesGPU };