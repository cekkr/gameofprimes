// simulation_gpu.js

let gpuCache = {
    device: null,
    adapter: null,
    computePipelines: new Map(), // Cache multiple pipelines
    bindGroupLayouts: new Map(), // Cache multiple layouts
    buffers: new Map(),
    inUse: false,
};

async function cleanupGPUResources() {
    if (gpuCache.device) {
        for (const buffer of gpuCache.buffers.values()) {
            buffer.destroy();
        }
        gpuCache.buffers.clear();

        // No need to explicitly destroy the device in WebGPU
        gpuCache.device = null;
        gpuCache.adapter = null;
        gpuCache.computePipelines.clear();
        gpuCache.bindGroupLayouts.clear();
        gpuCache.inUse = false;
        console.log("GPU resources cleaned up.");
    }
}


async function calculateForcesGPU(positions, masses, charges, rules, moleculeNumbers) {
    if (gpuCache.inUse) {
        console.warn("GPU calculation already in progress. Skipping.");
        return null;
    }
    gpuCache.inUse = true;

    try {
        if (!navigator.gpu) {
            throw new Error("WebGPU not supported.");
        }

        const adapter = gpuCache.adapter || (gpuCache.adapter = await navigator.gpu.requestAdapter());
        if (!adapter) {
            throw new Error("No GPU adapter found.");
        }
        const device = gpuCache.device || (gpuCache.device = await adapter.requestDevice());

        // --- Constants and Uniforms ---
        const constants = new Float32Array([
            rules.getConstant('G'),
            rules.getConstant('k'),
            rules.getConstant('min_distance'),
            rules.getConstant('max_force'),
            rules.getConstant('damping'),
            rules.getConstant('family_repulsion'),
            //Aggiungi altre costanti qui
        ]);

        const constantsBuffer = getOrCreateBuffer(device, 'constants', constants.byteLength, GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST, constants);

        // --- Shader Code (Dynamic) ---
        const shaderCode = createShaderCode(rules); // Genera lo shader dinamicamente

        // --- Pipeline and Layout Caching ---
        const pipelineKey = `forces_${shaderCode}`; // Unique key for this shader
        let computePipeline = gpuCache.computePipelines.get(pipelineKey);
        let bindGroupLayout = gpuCache.bindGroupLayouts.get(pipelineKey);

        if (!computePipeline) {
            bindGroupLayout = device.createBindGroupLayout({
                entries: [
                    { binding: 0, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
                    { binding: 1, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
                    { binding: 2, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
                    { binding: 3, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } },
                    { binding: 4, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
                    { binding: 5, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } }, // Add this line for moleculeNumbers
                ],
            });
            gpuCache.bindGroupLayouts.set(pipelineKey, bindGroupLayout);


            computePipeline = device.createComputePipeline({
                layout: device.createPipelineLayout({
                    bindGroupLayouts: [bindGroupLayout]
                }),
                compute: {
                    module: device.createShaderModule({ code: shaderCode }),
                    entryPoint: 'main',
                },
            });
            gpuCache.computePipelines.set(pipelineKey, computePipeline);
            console.log("New compute pipeline created and cached.");
        }

        // --- Buffers ---
        const moleculeCount = positions.length / 3;
        const positionBuffer = getOrCreateBuffer(device, 'position', positions.byteLength, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST, positions);
        const massBuffer = getOrCreateBuffer(device, 'mass', masses.byteLength, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST, masses);
        const chargeBuffer = getOrCreateBuffer(device, 'charge', charges.byteLength, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST, charges);
        const forceBuffer = getOrCreateBuffer(device, 'force', positions.byteLength, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC);
        const moleculeNumbersBuffer = getOrCreateBuffer(device, 'moleculeNumbers', moleculeNumbers.byteLength, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST, moleculeNumbers);
        const stagingBuffer = getOrCreateBuffer(device, 'staging', positions.byteLength, GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST);

        // --- Bind Group ---
        const bindGroup = device.createBindGroup({
            layout: bindGroupLayout,
            entries: [
                { binding: 0, resource: { buffer: positionBuffer } },
                { binding: 1, resource: { buffer: massBuffer } },
                { binding: 2, resource: { buffer: chargeBuffer } },
                { binding: 3, resource: { buffer: forceBuffer } },
                { binding: 4, resource: { buffer: constantsBuffer } },
                { binding: 5, resource: { buffer: moleculeNumbersBuffer } }, // Add this line
            ],
        });

        // --- Command Encoder ---
        const commandEncoder = device.createCommandEncoder();
        const passEncoder = commandEncoder.beginComputePass();
        passEncoder.setPipeline(computePipeline);
        passEncoder.setBindGroup(0, bindGroup);
        passEncoder.dispatchWorkgroups(Math.ceil(moleculeCount / 64));
        passEncoder.end();

        commandEncoder.copyBufferToBuffer(forceBuffer, 0, stagingBuffer, 0, stagingBuffer.size);
        device.queue.submit([commandEncoder.finish()]);

        // --- Read Results ---
        await stagingBuffer.mapAsync(GPUMapMode.READ);
        const forces = new Float32Array(stagingBuffer.getMappedRange().slice());
        stagingBuffer.unmap();

        gpuCache.inUse = false;
        return forces;

    } catch (error) {
        console.error("GPU calculation error:", error);
        gpuCache.inUse = false;
        // Consider *not* cleaning up on every error, only on explicit cleanup
        return null;
    }
}

// --- Helper Functions ---

function getOrCreateBuffer(device, key, size, usage, data = null) {
    const bufferKey = `${key}_${size}`;
    let buffer = gpuCache.buffers.get(bufferKey);

    if (!buffer) {
        buffer = device.createBuffer({ size, usage });
        gpuCache.buffers.set(bufferKey, buffer);
        console.log(`New buffer created: ${bufferKey}`);
    }

    if (data) {
        device.queue.writeBuffer(buffer, 0, data);
    }
    return buffer;
}

// Dynamically generate the shader code based on the simulation rules
function createShaderCode(rules) {

    // Helper to get prime factors as a string (up to a reasonable limit)
    const getPrimeFactorsString = (n) => {
        const factors = {};
        let d = 2;
        while (n > 1) {
            while (n % d === 0) {
                factors[d] = (factors[d] || 0) + 1;
                n /= d;
            }
            d++;
            if (d * d > n && n > 1) {
                factors[n] = (factors[n] || 0) + 1;
                break;
            }
        }
        // Convert to a string representation suitable for WGSL.  Example: "1u, 2u, 0u, 1u, 0u, ..."
        const factorArray = [];
        for (let i = 1; i <= 10; i++) { // Limit to first 10 primes for simplicity
            factorArray.push((factors[i] || 0) + "u");
        }

        return factorArray.join(", ");
    };

    // Now we build the shader string *dynamically*, incorporating the rules.
    const shaderCode = `
        struct Constants {
            G: f32,
            k: f32,
            min_distance: f32,
            max_force: f32,
            damping: f32,
            family_repulsion: f32,
        };
        @group(0) @binding(4) var<uniform> constants: Constants;

        @group(0) @binding(0) var<storage, read> positions: array<vec3<f32>>;
        @group(0) @binding(1) var<storage, read> masses: array<f32>;
        @group(0) @binding(2) var<storage, read> charges: array<f32>;
        @group(0) @binding(3) var<storage, read_write> forces: array<vec3<f32>>;
        @group(0) @binding(5) var<storage, read> moleculeNumbers: array<u32>;


        fn get_prime_factors(n: u32) -> array<u32, 10> {
            // Replace with your prime factorization logic, limited to 10 factors.
            // This is just an example, you'll need a more robust method.
            var factors = array<u32, 10>(0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u);
            var num = n;
            var divisor = 2u;

            var index = 0u;
            while(num > 1u && index < 10u){
                if(num % divisor == 0u){
                    factors[index] = divisor;
                    num = num / divisor;
                    index = index + 1u;
                } else {
                    divisor = divisor + 1u;
                }
            }

            return factors;
        }

        fn calculate_charge(factors: array<u32, 10>) -> f32 {
            var charge = 0.0;
            //Simplified charge calculation for demonstration
            for(var i = 0u; i < 10u; i++){
                let prime = factors[i];
                if(prime == 2u){
                    charge += 2.0;
                } else if (prime % 4u == 1u){
                    charge += 1.0;
                } else if (prime > 0u){
                    charge -= 1.0;
                }
            }
            return charge;
        }

       fn resonance_force(direction_normalized: vec3<f32>, distance: f32, mass1: f32, mass2: f32) -> vec3<f32> {
            let orbital = cross(direction_normalized, vec3<f32>(0.0, 1.0, 0.0));
            let norm = length(orbital);

            if (norm > 0.0) {
                let massEffect = sqrt(mass1 * mass2) / (mass1 + mass2);
                return orbital * massEffect / (distance * distance * distance + 0.1); //Prevent too large forces
            }
            return vec3<f32>(0.0, 0.0, 0.0);
        }

        fn prime_resonance_condition(factors1: array<u32, 10>, factors2: array<u32, 10>, charge1: f32, charge2: f32) -> bool {
            // Rename 'shared' to 'common' since 'shared' is a reserved keyword
            const common = count_common_factors(factors1, factors2);
            
            // Rest of your implementation
            // (You didn't provide the full function body)
            return false; // Replace with your actual implementation
        }

        fn prime_resonance_force(direction_normalized: vec3<f32>, distance: f32, mass1: f32, mass2: f32, charge1: f32, charge2: f32) -> vec3<f32> {
            // Your implementation goes here
            // (You didn't provide the full function body)
            return vec3<f32>(0.0, 0.0, 0.0); // Replace with your actual implementation
        }

        fn cross(a: vec3<f32>, b: vec3<f32>) -> vec3<f32> {
            return vec3<f32>(
                a.y * b.z - a.z * b.y,
                a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x
            );
        }

        ${rules.interaction_rules.map(rule => `
            fn ${rule.name}_condition(factors1: array<u32, 10>, factors2: array<u32, 10>, charge1: f32, charge2: f32) -> bool {
                ${rule.condition.toString().replace(/this\.rules\.getConstant/g, 'constants.')
                  .replace(/factors1/g, 'factors1')
                  .replace(/factors2/g, 'factors2')
                  .replace(/calculateCharge\(\s*factors1\s*\)/g, 'charge1') // Replace calculateCharge calls
                  .replace(/calculateCharge\(\s*factors2\s*\)/g, 'charge2')
                  .replace(/Object\.keys\(\s*factors1\s*\)\.filter\(key => factors2\.hasOwnProperty\(key\)\)/g, 'count_common_factors(factors1, factors2)') //Object.keys
                  .replace(/Object\.keys\(\s*factors1\s*\)/g, 'get_present_factors(factors1)')
                  .replace(/\.hasOwnProperty\(key\)/g, '[key] > 0u') //hasOwnProperty
                  .replace(/^function[^{]+\{/i, '') // Remove function signature
                  .replace(/}$/i, '')} // Remove closing brace

            }

            fn ${rule.name}_force(direction_normalized: vec3<f32>, distance: f32, mass1: f32, mass2: f32, charge1: f32, charge2: f32) -> vec3<f32> {
                ${rule.force_function.toString().replace(/this\.rules\.getConstant/g, 'constants.')
                  .replace(/direction/g, 'direction_normalized')
                  .replace(/crossProduct/g, 'cross')
                  .replace(/^function[^{]+\{/i, '') // Remove function signature
                  .replace(/}$/i, '')} // Remove closing brace
            }
        `).join('\n')}

        // Helper to count common factors
        fn count_common_factors(factors1: array<u32, 10>, factors2: array<u32, 10>) -> u32 {
            var count = 0u;
            for (var i = 0u; i < 10u; i++) {
                if (factors1[i] > 0u && factors2[i] > 0u) {
                    count += 1u;
                }
            }
            return count;
        }
        //Helper to get present factors
        fn get_present_factors(factors: array<u32, 10>) -> array<u32, 10> {
             return factors;
        }


        @compute @workgroup_size(64)
        fn main(@builtin(global_invocation_id) global_id: vec3<u32>) {
            let num_molecules = arrayLength(&positions);
            if (global_id.x >= num_molecules) {
                return;
            }

            var force = vec3<f32>(0.0, 0.0, 0.0);
            let factors1 = get_prime_factors(moleculeNumbers[global_id.x]);
            let charge1 = calculate_charge(factors1);

            for (var i = 0u; i < num_molecules; i++) {
                if (i == global_id.x) {
                    continue;
                }

                let direction = positions[i] - positions[global_id.x];
                let distance = length(direction);

                if (distance < constants.min_distance) {
                    continue;
                }

                let direction_normalized = normalize(direction);
                let factors2 = get_prime_factors(moleculeNumbers[i]);
                let charge2 = calculate_charge(factors2);


                ${rules.interaction_rules.map(rule => `
                    if (${rule.name}_condition(factors1, factors2, charge1, charge2)) {
                        force += ${rule.name}_force(direction_normalized, distance, masses[global_id.x], masses[i], charge1, charge2);
                    }
                `).join('\n')}
            }

            forces[global_id.x] = force;
        }
    `;
    return shaderCode;
}

export { calculateForcesGPU, cleanupGPUResources };