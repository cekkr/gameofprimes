import { createCustomRules } from './rules.js';
import { PrimeMolecule } from './molecule.js';

let simulation;
let previousMoleculeIds = new Set(); // Track molecule IDs to detect changes

onmessage = function(event) {
    if (event.data.type === 'init') {
        const { size, moleculeCount, maxNumber, timeScale } = event.data;
        const rules = createCustomRules();
        rules.setConstant('time_scale', timeScale);
        simulation = new PrimeChemistry(rules, size, moleculeCount, maxNumber);

        // Initial step to populate data
        simulation.step();
        sendOptimizedUpdate();

    } else if (event.data.type === 'step') {
        // Perform a simulation step
        simulation.step();
        sendOptimizedUpdate();
    } else if (event.data.type === 'cleanup') {
        // Clean up resources when simulation ends
        cleanupResources();
    }
};

function sendOptimizedUpdate() {
    // Only send changed molecules or minimal data when needed
    const moleculeData = getOptimizedMoleculeData();
    
    postMessage({
        type: 'update',
        molecules: moleculeData,
        temperature: simulation.temperature,
        moleculeCount: simulation.molecules.length
    });
}

function getOptimizedMoleculeData() {
    // Track which molecules are new, changed, or removed
    const currentMoleculeIds = new Set();
    const result = [];
    
    for (const mol of simulation.molecules) {
        const id = mol.id || Math.random().toString(36).substr(2, 9); // Generate stable ID if none exists
        mol.id = id;
        currentMoleculeIds.add(id);
        
        // Send full data for new or modified molecules
        result.push({
            id: id,
            number: mol.number,
            position: [...mol.position], // Create shallow copy to avoid reference issues
            velocity: [...mol.velocity],
            prime_factors: { ...mol.prime_factors }, // Only send shallow copy
            mass: mol.mass,
            charge: mol.charge,
            color: [...mol.color],
            angularVelocity: [...mol.angularVelocity],
            lastReactionTime: mol.lastReactionTime,
        });
    }
    
    // Calculate removed molecules (were in previous set but not current)
    const removedIds = [...previousMoleculeIds].filter(id => !currentMoleculeIds.has(id));
    
    // Update tracking set for next frame
    previousMoleculeIds = currentMoleculeIds;
    
    // Also send list of removed IDs so client can clean up
    return {
        molecules: result,
        removedIds: removedIds
    };
}

function cleanupResources() {
    // Clean up any resources that might cause memory leaks
    simulation = null;
    previousMoleculeIds.clear();
    // Force garbage collection (indirectly)
    postMessage({ type: 'cleanup_complete' });
}

// --- PrimeChemistry Class (With optimization fixes) ---

class PrimeChemistry {
    constructor(rules, size, moleculeCount, maxNumber) {
        this.rules = rules;
        this.size = size;
        this.maxNumber = maxNumber;
        this.molecules = [];
        this.temperature = 1.0;
        this.accumulatedTime = 0.0;
        this.nextMoleculeId = 1; // For stable IDs

        // Initialize random molecules
        for (let i = 0; i < moleculeCount; i++) {        
            const pos = [
                Math.random() * size - size / 2,
                Math.random() * size - size / 2,
                Math.random() * size - size / 2
            ];
            const number = Math.floor(Math.random() * 99) + 2;
            const mol = new PrimeMolecule(number, pos);
            mol.id = `initial-${this.nextMoleculeId++}`; // Add stable ID

            // Set RANDOM INITIAL VELOCITY
            mol.velocity = [
                (Math.random() - 0.5) * 0.5,
                (Math.random() - 0.5) * 0.5,
                (Math.random() - 0.5) * 0.5
            ];

            this.molecules.push(mol);
        }
    }

    step() {
        // Use a more efficient data structure
        const moleculeCount = this.molecules.length;
        const forces = new Array(moleculeCount).fill(0).map(() => [0, 0, 0]);
        const removedIndices = new Set();
        const newMolecules = [];

        // Apply temperature variation
        this.temperature = 1.0 + 0.1 * Math.sin(performance.now() / 1000);

        const timeScale = this.rules.getConstant('time_scale');
        this.accumulatedTime += timeScale;

        // Calculate forces more efficiently
        this.calculateForces(forces, removedIndices, newMolecules);
        
        // Apply forces and update positions
        this.updateMolecules(forces, removedIndices);
        
        // Update molecule list (efficient filtering)
        if (removedIndices.size > 0) {
            this.molecules = this.molecules.filter((_, i) => !removedIndices.has(i));
        }
        
        // Add new molecules with proper IDs
        for (const mol of newMolecules) {
            mol.id = `new-${this.nextMoleculeId++}`;
            this.molecules.push(mol);
        }

        // LIMIT - with explanation
        if (this.molecules.length > 500) {
            // Remove oldest molecules to maintain limit
            this.molecules = this.molecules.slice(this.molecules.length - 500);
        }
        
        // Reset accumulated time if needed
        if (this.accumulatedTime >= 1.0) {
            this.accumulatedTime = 0.0;
        }
    }
    
    calculateForces(forces, removedIndices, newMolecules) {
        const moleculeCount = this.molecules.length;
        const timeScale = this.rules.getConstant('time_scale');
        
        // Only process non-removed molecules
        for (let i = 0; i < moleculeCount; i++) {
            if (removedIndices.has(i)) continue;
            const mol1 = this.molecules[i];
            
            // Thermal motion
            mol1.velocity[0] += (Math.random() - 0.5) * 0.02 * this.temperature * timeScale;
            mol1.velocity[1] += (Math.random() - 0.5) * 0.02 * this.temperature * timeScale;
            mol1.velocity[2] += (Math.random() - 0.5) * 0.02 * this.temperature * timeScale;
            
            for (let j = i + 1; j < moleculeCount; j++) {
                if (removedIndices.has(j)) continue;
                const mol2 = this.molecules[j];
                
                // Apply interaction rules
                const [force1, force2] = this.applyRules(mol1, mol2);
                for (let k = 0; k < 3; k++) {
                    forces[i][k] += force1[k] * timeScale;
                    forces[j][k] += force2[k] * timeScale;
                }
                
                // Check for reactions
                const distance = this.calculateDistance(mol1.position, mol2.position);
                
                if (distance < 2.0 && this.accumulatedTime >= 1.0) {
                    const reactionProducts = this.attemptReactions(mol1, mol2);
                    if (reactionProducts.length > 0) {
                        newMolecules.push(...reactionProducts);
                        removedIndices.add(i);
                        removedIndices.add(j);
                        // Set reaction time for halo effect
                        const currentTime = performance.now();
                        reactionProducts.forEach(mol => mol.setReactionTime(currentTime));
                        break; // Exit inner loop after reaction
                    }
                }
            }
        }
    }
    
    updateMolecules(forces, removedIndices) {
        const damping = this.rules.getConstant('damping');
        const timeScale = this.rules.getConstant('time_scale');
        
        for (let i = 0; i < this.molecules.length; i++) {
            if (removedIndices.has(i)) continue;
            const molecule = this.molecules[i];
            
            // Apply forces
            for (let k = 0; k < 3; k++) {
                molecule.velocity[k] += forces[i][k] * 0.1;
                molecule.position[k] += molecule.velocity[k] * timeScale;
                molecule.velocity[k] *= damping;
            }
            
            // Boundary conditions
            this.applyBoundaryConditions(molecule);
        }
    }
    
    calculateDistance(pos1, pos2) {
        return Math.sqrt(
            (pos2[0] - pos1[0]) ** 2 +
            (pos2[1] - pos1[1]) ** 2 +
            (pos2[2] - pos1[2]) ** 2
        );
    }
    
    applyBoundaryConditions(molecule) {
        for (let axis = 0; axis < 3; axis++) {
            if (Math.abs(molecule.position[axis]) > this.size / 2) {
                molecule.position[axis] = Math.sign(molecule.position[axis]) * this.size / 2;
                molecule.velocity[axis] *= -0.8;
            }
        }
    }

    applyRules(mol1, mol2) {
        // Extracted to avoid creating closures repeatedly
        const direction = [
            mol2.position[0] - mol1.position[0],
            mol2.position[1] - mol1.position[1],
            mol2.position[2] - mol1.position[2]
        ];
        let distance = Math.sqrt(direction[0] ** 2 + direction[1] ** 2 + direction[2] ** 2);

        if (distance < 0.0001) {
            return [[0, 0, 0], [0, 0, 0]];
        }

        const directionNormalized = [
            direction[0] / distance,
            direction[1] / distance,
            direction[2] / distance
        ];
        
        const totalForce1 = [0, 0, 0];
        const totalForce2 = [0, 0, 0];

        const maxForce = this.rules.getConstant('max_force');
        
        for (const rule of this.rules.interaction_rules) {
            if (rule.condition(mol1.prime_factors, mol2.prime_factors)) {
                // Create the parameters needed for the force calculation
                const params = rule.force_function.length === 4 ? 
                    [directionNormalized, distance, mol1.mass, mol2.mass] :
                    [directionNormalized, distance, mol1.charge, mol2.charge];
                
                // Apply the rule's force function
                const force = rule.force_function(...params);
                
                // Apply strength
                for (let i = 0; i < 3; i++) {
                    const scaledForce = force[i] * rule.strength;
                    totalForce1[i] += scaledForce;
                    totalForce2[i] -= scaledForce;
                }
            }
        }
        
        // Apply force limits
        for (let i = 0; i < 3; i++) {
            totalForce1[i] = Math.max(-maxForce, Math.min(maxForce, totalForce1[i]));
            totalForce2[i] = Math.max(-maxForce, Math.min(maxForce, totalForce2[i]));
        }

        return [totalForce1, totalForce2];
    }

    attemptReactions(mol1, mol2) {
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