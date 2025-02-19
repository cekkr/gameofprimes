import { createCustomRules } from './rules.js';
import { PrimeMolecule } from './molecule.js'; // Import from molecule.js

let simulation;

onmessage = function(event) {
    if (event.data.type === 'init') {
        //console.log("Worker received init:", event.data); // Debugging
        // Initialize simulation (similar to Python's PrimeChemistry)
        const { size, moleculeCount, maxNumber, timeScale } = event.data;
        const rules = createCustomRules(); // Use the rules from rules.js
        rules.setConstant('time_scale', timeScale)  // Apply time scale
        simulation = new PrimeChemistry(rules, size, moleculeCount, maxNumber);
        //console.log("Simulation initialized in worker:", simulation); // Debugging

        // Initial step to populate data
        simulation.step();
        postMessage({ //SEND INITIAL DATA
            type: 'update',
            molecules: simulation.molecules.map(mol => ({
                number: mol.number,
                position: mol.position,
                velocity: mol.velocity,
                prime_factors: mol.prime_factors,
                mass: mol.mass,
                charge: mol.charge,
                color: mol.color,
            })),
            temperature: simulation.temperature,
        });

    } else if (event.data.type === 'step') {
        // Perform a simulation step
        //console.log("Worker received step command"); // Debugging
        simulation.step();

        // Send updated data back to the main thread
        postMessage({
            type: 'update',
            molecules: simulation.molecules.map(mol => ({
                number: mol.number,
                position: mol.position,
                velocity: mol.velocity,
                prime_factors: mol.prime_factors,
                mass: mol.mass,
                charge: mol.charge,
                color: mol.color,
            })),
            temperature: simulation.temperature, // Send temperature
        });
    }
};

// --- PrimeChemistry Class (Adapted for JS) ---



class PrimeChemistry {
   constructor(rules, size, moleculeCount, maxNumber) {
        this.rules = rules;
        this.size = size;
        this.maxNumber = maxNumber;
        this.molecules = [];
        this.temperature = 1.0;
        this.accumulatedTime = 0.0; // Track accumulated time


        // Initialize random molecules
        for (let i = 0; i < moleculeCount; i++) {
            const pos = [
                Math.random() * size - size / 2,
                Math.random() * size - size / 2,
                Math.random() * size - size / 2
            ];
            const number = Math.floor(Math.random() * 99) + 2; // Numbers 2-100
            this.molecules.push(new PrimeMolecule(number, pos));
        }
    }

    step() {
        const newMolecules = [];
        const removedMolecules = new Set();

        // Apply temperature variation
        this.temperature = 1.0 + 0.1 * Math.sin(performance.now() / 1000);

        const forces = {}; // Use object for easier key access
        const timeScale = this.rules.getConstant('time_scale');
        this.accumulatedTime += timeScale; // accumulate time


        for (let i = 0; i < this.molecules.length; i++) {
            if (removedMolecules.has(i)) continue;
            const mol1 = this.molecules[i];

             // Thermal motion
            mol1.velocity[0] += (Math.random() - 0.5) * 0.02 * this.temperature * timeScale;
            mol1.velocity[1] += (Math.random() - 0.5) * 0.02 * this.temperature * timeScale;
            mol1.velocity[2] += (Math.random() - 0.5) * 0.02 * this.temperature * timeScale;

            forces[i] = forces[i] || [0, 0, 0]; // Initialize if not already present

            for (let j = i + 1; j < this.molecules.length; j++) {
                if (removedMolecules.has(j)) continue;
                const mol2 = this.molecules[j];

                 // Apply interaction rules - forces are now scaled by timeScale
                const [force1, force2] = this.applyRules(mol1, mol2);
                forces[i][0] += force1[0] * timeScale;
                forces[i][1] += force1[1] * timeScale;
                forces[i][2] += force1[2] * timeScale;

                forces[j] = forces[j] || [0, 0, 0];
                forces[j][0] += force2[0] * timeScale;
                forces[j][1] += force2[1] * timeScale;
                forces[j][2] += force2[2] * timeScale;


                // Check for reactions
                const distance = Math.sqrt(
                    (mol2.position[0] - mol1.position[0]) ** 2 +
                    (mol2.position[1] - mol1.position[1]) ** 2 +
                    (mol2.position[2] - mol1.position[2]) ** 2
                );

                if (distance < 2.0 && this.accumulatedTime >= 1.0) {
                    const reactionProducts = this.attemptReactions(mol1, mol2);
                    if (reactionProducts.length > 0) {
                        newMolecules.push(...reactionProducts);
                        removedMolecules.add(i);
                        removedMolecules.add(j);
                        break; // Exit inner loop after reaction
                    }

                }
            }
        }

        // Reset accumulated time if it exceeds 1.0
        if (this.accumulatedTime >= 1.0) {
             this.accumulatedTime = 0.0;
        }


        const damping = this.rules.getConstant('damping');
        for (let i = 0; i < this.molecules.length; i++) {
             if (!removedMolecules.has(i)) {
                const molecule = this.molecules[i];
                // Apply forces - scaled by timeScale
                molecule.velocity[0] += forces[i][0] * 0.1;
                molecule.velocity[1] += forces[i][1] * 0.1;
                molecule.velocity[2] += forces[i][2] * 0.1;


                molecule.position[0] += molecule.velocity[0] * timeScale;
                molecule.position[1] += molecule.velocity[1] * timeScale;
                molecule.position[2] += molecule.velocity[2] * timeScale;

                // Damping
                molecule.velocity[0] *= damping;
                molecule.velocity[1] *= damping;
                molecule.velocity[2] *= damping;

                // Boundary conditions
                for (let axis = 0; axis < 3; axis++) {
                    if (Math.abs(molecule.position[axis]) > this.size / 2) {
                        molecule.position[axis] = Math.sign(molecule.position[axis]) * this.size / 2;
                        molecule.velocity[axis] *= -0.8;
                    }
                }
            }
        }

         // Update molecule list (efficient filtering)
        this.molecules = this.molecules.filter((_, i) => !removedMolecules.has(i)).concat(newMolecules);
    }
    applyRules(mol1, mol2) {
        const direction = [
            mol2.position[0] - mol1.position[0],
            mol2.position[1] - mol1.position[1],
            mol2.position[2] - mol1.position[2]
        ];
        let distance = Math.sqrt(direction[0] ** 2 + direction[1] ** 2 + direction[2] ** 2);

        if (distance < 0.0001) {
            return [ [0, 0, 0], [0, 0, 0] ];
        }

        const directionNormalized = direction.map(x => x / distance);
        let totalForce1 = [0, 0, 0];
        let totalForce2 = [0, 0, 0];

        for (const rule of this.rules.interaction_rules) {
            if (rule.condition(mol1.prime_factors, mol2.prime_factors)) {
                let force = rule.force_function(
                    directionNormalized,
                    distance,
                    rule.force_function.length === 4 && rule.force_function.name.includes('mass') ? mol1.mass : mol1.charge,
                    rule.force_function.length === 4 && rule.force_function.name.includes('mass') ? mol2.mass : mol2.charge,
                );


                force = force.map (x=> x * rule.strength);

                totalForce1[0] += force[0];
                totalForce1[1] += force[1];
                totalForce1[2] += force[2];

                totalForce2[0] -= force[0];
                totalForce2[1] -= force[1];
                totalForce2[2] -= force[2];
            }
        }
        // Apply force limits
        const maxForce = this.rules.getConstant('max_force');
        const clipForce = (f) => Math.max(-maxForce, Math.min(maxForce, f));

        totalForce1 = totalForce1.map(clipForce);
        totalForce2 = totalForce2.map(clipForce);

        return [totalForce1, totalForce2];

    }

    attemptReactions(mol1, mol2) {
        const distance = Math.sqrt(
            (mol2.position[0] - mol1.position[0]) ** 2 +
            (mol2.position[1] - mol1.position[1]) ** 2 +
            (mol2.position[2] - mol1.position[2]) ** 2
        );
        const newMolecules = [];

        for (const rule of this.rules.reaction_rules) {
            if (rule.condition(mol1.prime_factors, mol2.prime_factors)) {
                if (Math.random() < rule.probability) {
                   return rule.effect(mol1, mol2);
                }
            }
        }
         return [];
    }
}