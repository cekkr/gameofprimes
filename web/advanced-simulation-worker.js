import { createCustomRules } from './rules.js';
import { PrimeMolecule } from './molecule.js';

/**
 * Aggiornamento della classe MoleculeSerializer per gestire le proprietà
 * delle reazioni di fissione, scissione ed emissione
 */
class MoleculeSerializer {
    /**
     * Creates a serializable object from a PrimeMolecule
     * Includes all properties needed for visualization and UI
     */
    static serialize(molecule) {
      // Properties that should be included (both computed and settable)
      return {
        // Basic properties
        id: molecule.id,
        number: molecule.number,
        position: [...molecule.position],
        velocity: [...molecule.velocity],
        
        // Computed properties needed for visualization
        prime_factors: {...molecule.prime_factors},
        mass: molecule.mass,
        charge: molecule.charge,
        color: [...molecule.color],
        
        // Optional properties
        ...(molecule.angularVelocity ? {angularVelocity: [...molecule.angularVelocity]} : {}),
        ...(molecule.lastReactionTime ? {lastReactionTime: molecule.lastReactionTime} : {}),
        
        // Proprietà di relazione (parentela tra molecole)
        ...(molecule.parentIds ? {parentIds: [...molecule.parentIds]} : {}),
        ...(molecule.reactionType ? {reactionType: molecule.reactionType} : {}),
        
        // Quantum field properties if present
        ...(molecule._isQuantumField ? {
          _isQuantumField: true,
          _creationTime: molecule._creationTime,
          _lifetime: molecule._lifetime,
          _fieldType: molecule._fieldType || 'standard',
          _fieldStrength: molecule._fieldStrength || 1.0
        } : {})
      };
    }
    
    /**
     * Creates a PrimeMolecule from serialized data
     * Only sets properties that are safely settable
     */
    static deserialize(molData, PrimeMolecule) {
      // Create molecule with the fundamental properties
      // This ensures prime_factors will be calculated correctly
      const mol = new PrimeMolecule(molData.number, molData.position);
      
      // Set only properties known to be safely settable
      if (molData.id) {
        mol.id = molData.id;
      }
      
      // Handle velocity with safety checks
      if (molData.velocity && Array.isArray(molData.velocity)) {
        mol.velocity = [...molData.velocity];
      }
      
      // Set other optional properties
      if (molData.angularVelocity && Array.isArray(molData.angularVelocity)) {
        mol.angularVelocity = [...molData.angularVelocity];
      }
      
      if (typeof molData.lastReactionTime === 'number') {
        mol.lastReactionTime = molData.lastReactionTime;
      }
      
      // Imposta proprietà di relazione
      if (molData.parentIds && Array.isArray(molData.parentIds)) {
        mol.parentIds = [...molData.parentIds];
      }
      
      if (molData.reactionType) {
        mol.reactionType = molData.reactionType;
      }
      
      // Set quantum field properties if present
      if (molData._isQuantumField) {
        mol._isQuantumField = true;
        mol._creationTime = molData._creationTime;
        mol._lifetime = molData._lifetime;
        mol._fieldType = molData._fieldType || 'standard';
        mol._fieldStrength = molData._fieldStrength || 1.0;
      }
      
      return mol;
    }
    
    /**
     * Gets properties that are safe to update on an existing molecule
     * Used when merging results from sub-workers
     */
    static getUpdatableProperties(molecule) {
      return {
        id: molecule.id,
        position: [...molecule.position],
        velocity: [...molecule.velocity],
        ...(molecule.angularVelocity ? {angularVelocity: [...molecule.angularVelocity]} : {}),
        ...(molecule.lastReactionTime ? {lastReactionTime: molecule.lastReactionTime} : {}),
        ...(molecule.parentIds ? {parentIds: [...molecule.parentIds]} : {}),
        ...(molecule.reactionType ? {reactionType: molecule.reactionType} : {}),
        ...(molecule._isQuantumField ? {
          _isQuantumField: true,
          _creationTime: molecule._creationTime,
          _lifetime: molecule._lifetime,
          _fieldType: molecule._fieldType,
          _fieldStrength: molecule._fieldStrength
        } : {})
      };
    }
    
    /**
     * Updates an existing molecule with new properties
     * Only updates properties that are safely settable
     */
    static updateMoleculeProperties(molecule, updateData) {
      // Update position if provided
      if (updateData.position && Array.isArray(updateData.position)) {
        molecule.position = [...updateData.position]; 
      }
      
      // Update velocity if provided
      if (updateData.velocity && Array.isArray(updateData.velocity)) {
        molecule.velocity = [...updateData.velocity];
      }
      
      // Update angular velocity if provided
      if (updateData.angularVelocity && Array.isArray(updateData.angularVelocity)) {
        molecule.angularVelocity = [...updateData.angularVelocity];
      }
      
      // Update lastReactionTime if provided
      if (typeof updateData.lastReactionTime === 'number') {
        molecule.lastReactionTime = updateData.lastReactionTime;
      }
      
      // Aggiorna proprietà di relazione
      if (updateData.parentIds && Array.isArray(updateData.parentIds)) {
        molecule.parentIds = [...updateData.parentIds];
      }
      
      if (updateData.reactionType) {
        molecule.reactionType = updateData.reactionType;
      }
      
      // Update quantum field properties
      if (updateData._isQuantumField) {
        molecule._isQuantumField = true;
        molecule._creationTime = updateData._creationTime;
        molecule._lifetime = updateData._lifetime;
        molecule._fieldType = updateData._fieldType || molecule._fieldType || 'standard';
        molecule._fieldStrength = updateData._fieldStrength || molecule._fieldStrength || 1.0;
      }
      
      return molecule;
    }
}

////
////
////

// Configuration for multi-worker setup
const NUM_SUB_WORKERS = navigator.hardwareConcurrency || 4; // Use available cores
const INTERACTION_CACHE_LIFETIME = 50; // Frames to keep interaction cache
//const SPATIAL_GRID_SIZE = 50; // Size of spatial partitioning grid cells

class EnhancedChemistry {
constructor(rules, size, moleculeCount, maxNumber) {
    this.rules = rules;
    this.size = size;
    this.maxNumber = maxNumber;
    this.molecules = [];
    this.temperature = 1.0;
    this.accumulatedTime = 0.0;
    this.nextMoleculeId = 1;
    this.reactionCount = 0;

    const SPATIAL_GRID_SIZE = Math.floor(size / 5); // instead of a fixed size
    
    // Initialize sub-workers
    this.subWorkers = [];
    this.isMainWorker = true;
    this.workerId = 'main';
    this.workerBusy = false;
    
    // Spatial partitioning grid for faster neighbor finding
    this.spatialGrid = new SpatialGrid(size, SPATIAL_GRID_SIZE);
    
    // Interaction cache system
    this.interactionCache = new InteractionCache(INTERACTION_CACHE_LIFETIME);
    
    // Initialize molecules with interesting distribution
    this.initializeMolecules(moleculeCount);
}

async initializeSubWorkers() {
    if (!this.isMainWorker) return; // Only main worker creates sub-workers
    
    const workerUrl = self.location.href;
    console.log(`Initializing ${NUM_SUB_WORKERS} sub-workers from ${workerUrl}`);
    
    for (let i = 0; i < NUM_SUB_WORKERS; i++) {
        try {
            const worker = new Worker(workerUrl, {type: 'module'});
            const workerIndex = i + 1;
            
            worker.onmessage = (event) => this.handleSubWorkerMessage(workerIndex, event);
            
            // Initialize sub-worker with subset of configuration
            worker.postMessage({
                type: 'init_sub',
                workerId: `sub-${workerIndex}`,
                size: this.size,
                maxNumber: this.maxNumber,
                timeScale: this.rules.getConstant('time_scale')
            });
            
            this.subWorkers.push({
                worker,
                id: `sub-${workerIndex}`,
                busy: false,
                lastProcessedMolecules: 0
            });
            
            console.log(`Initialized sub-worker ${workerIndex}`);
        } catch (error) {
            console.error(`Failed to initialize sub-worker ${i + 1}:`, error);
        }
    }
}

handleSubWorkerMessage(workerIndex, event) {
    const subWorker = this.subWorkers[workerIndex - 1];
    
    switch (event.data.type) {
        case 'sub_ready':
            subWorker.busy = false;
            console.log(`Sub-worker ${workerIndex} ready`);
            break;
            
        case 'chunk_processed':
            // Merge results from sub-worker
            this.mergeProcessedChunk(event.data.results);
            subWorker.busy = false;
            subWorker.lastProcessedMolecules = event.data.processedCount;
            break;
            
        case 'reaction_occurred':
            // Handle reaction result from sub-worker
            this.handleRemoteReaction(event.data.reaction);
            break;
            
        case 'error':
            console.error(`Error in sub-worker ${workerIndex}:`, event.data.message);
            subWorker.busy = false;
            break;
    }
}

// Merge results from sub-worker
mergeProcessedChunk(results) {
    if (!results) return;
    
    // Update molecule states
    for (const update of results.moleculeUpdates) {
        const molecule = this.molecules.find(m => m.id === update.id);
        if (molecule) {
            // Use the serialization helper to safely update properties
            MoleculeSerializer.updateMoleculeProperties(molecule, update);
        }
    }
    
    // Add any new molecules from reactions
    if (results.newMolecules && results.newMolecules.length > 0) {
        for (const molData of results.newMolecules) {
            // Create new molecule using the serialization helper
            const newMol = MoleculeSerializer.deserialize(molData, PrimeMolecule);
            newMol.id = `main-${this.nextMoleculeId++}`;
            this.molecules.push(newMol);
        }
    }
    
    // Update reaction count
    if (results.reactionCount) {
        this.reactionCount += results.reactionCount;
    }
}    

// Method to extract the minimal necessary properties for a molecule
handleRemoteReaction(reaction) {
    // Handle a reaction that occurred in a sub-worker
    if (!reaction) return;
    
    // Find and remove reacted molecules
    const mol1Index = this.molecules.findIndex(m => m.id === reaction.reactant1Id);
    const mol2Index = this.molecules.findIndex(m => m.id === reaction.reactant2Id);
    
    if (mol1Index >= 0 && mol2Index >= 0) {
        // Remove the reactants
        const removed = [
            this.molecules[mol1Index],
            this.molecules[mol2Index]
        ];
        
        this.molecules = this.molecules.filter((_, i) => 
            i !== mol1Index && i !== mol2Index);
            
        // Add product molecules
        for (const productData of reaction.products) {
            // Create new molecule using the serialization helper
            const newMol = MoleculeSerializer.deserialize(productData, PrimeMolecule);
            newMol.id = `main-${this.nextMoleculeId++}`;
            this.molecules.push(newMol);
        }
        
        // Update reaction count
        this.reactionCount++;
        
        // Update interaction cache
        this.interactionCache.invalidateForMolecules(removed);
    }
}

initializeMolecules(count) {
    // Calcola i numeri primi fino a this.maxNumber con il Crivello di Eratostene
    const calculatePrimes = (max) => {
        const sieve = Array(max + 1).fill(true);
        sieve[0] = sieve[1] = false;
        
        for (let i = 2; i * i <= max; i++) {
            if (sieve[i]) {
                for (let j = i * i; j <= max; j += i) {
                    sieve[j] = false;
                }
            }
        }
        
        return Array.from({ length: max + 1 }, (_, i) => i)
            .filter(num => sieve[num]);
    };
    
    // Genera i numeri primi e composti in base a this.maxNumber
    const primes = calculatePrimes(this.maxNumber);
    
    // Genera i composti (numeri non primi) fino a this.maxNumber
    const compounds = Array.from(
        { length: this.maxNumber - 1 }, 
        (_, i) => i + 2
    ).filter(num => !primes.includes(num));
    
    // Distribuzione di numeri per le molecole
    let numberChoices = [];
    
    // Aggiungi numeri primi
    numberChoices.push(...primes);
    
    // Aggiungi composti
    numberChoices.push(...compounds);
    
    // Aggiungi alcuni numeri casuali più grandi se necessario
    if (numberChoices.length < count * 2) {
        const additionalCount = Math.min(10, count - numberChoices.length);
        for (let i = 0; i < additionalCount; i++) {
            const largeNumber = Math.floor(Math.random() * (this.maxNumber * 0.5)) + this.maxNumber * 0.5;
            if (!numberChoices.includes(largeNumber)) {
                numberChoices.push(largeNumber);
            }
        }
    }
    
    // Crea le molecole
    for (let i = 0; i < count; i++) {
        // Posizione casuale all'interno dei limiti
        const pos = [
            (Math.random() - 0.5) * this.size,
            (Math.random() - 0.5) * this.size,
            (Math.random() - 0.5) * this.size
        ];
        
        // Scegli un numero dalla distribuzione
        const number = numberChoices[Math.floor(Math.random() * numberChoices.length)];
    
        // Crea molecola con ID stabile
        const mol = new PrimeMolecule(number, pos);
        mol.id = `initial-${this.nextMoleculeId++}`;
        
        // Imposta velocità iniziale casuale
        mol.velocity = [
            (Math.random() - 0.5) * 0.5,
            (Math.random() - 0.5) * 0.5,
            (Math.random() - 0.5) * 0.5
        ];
        
        this.molecules.push(mol);
    }
    
    console.log(`Inizializzate ${this.molecules.length} molecole (Primi: ${primes.length}, Composti: ${compounds.length})`);
    
    // Initialize spatial grid with molecules
    this.updateSpatialGrid();
}

async step() {
    const timeScale = this.rules.getConstant('time_scale');
    this.accumulatedTime += timeScale;
    
    // Update temperature
    this.temperature = 1.0 + 0.2 * Math.sin(performance.now() / 5000);
    
    // Update spatial grid
    this.updateSpatialGrid();
    
    // Arrays for efficient processing
    const newMolecules = [];
    const removedIndices = new Set();
    
    // Process physics using sub-workers if available
    if (this.isMainWorker && this.subWorkers.length > 0) {
        await this.distributeWorkToSubWorkers();
    } else {
        // Fallback to single-threaded processing
        this.updatePhysics(removedIndices, newMolecules);
    }
    
    // Filter removed molecules
    if (removedIndices.size > 0) {
        this.molecules = this.molecules.filter((_, i) => !removedIndices.has(i));
    }
    
    // Add new molecules
    for (const mol of newMolecules) {
        mol.id = `new-${this.nextMoleculeId++}`;
        this.molecules.push(mol);
    }
    
    // Reset accumulated time
    if (this.accumulatedTime >= 1.0) {
        this.accumulatedTime = 0.0;
        
        // Periodically clean up caches
        this.interactionCache.ageEntries();
    }
    
    // Limit molecule count for performance
    this.manageMoleculeCount();
    
    // Update quantum fields
    updateQuantumFields(this);
    
    // Periodically clean up obsolete relationships
    if (Math.random() < 0.05) {
        this.rules.cleanupOldRelations && this.rules.cleanupOldRelations();
    }
}

async mainWorkerStep() {
    // This method is now incorporated into step()
}

async distributeWorkToSubWorkers() {
    if (this.subWorkers.length === 0) return;
    
    // Group molecules by spatial partition for better distribution
    const partitions = this.groupMoleculesByPartition();
    
    // Wait for all sub-workers to be ready
    await this.waitForAvailableWorkers();
    
    // Distribute partitions to available workers
    let partitionIndex = 0;
    const processingPromises = [];
    
    for (const subWorker of this.subWorkers) {
        if (partitionIndex >= partitions.length) break;
        
        const partition = partitions[partitionIndex++];
        if (!partition || partition.length === 0) continue;
        
        subWorker.busy = true;
        
        const promise = new Promise(resolve => {
            const workerMessageHandler = (event) => {
                if (event.data.type === 'chunk_processed') {
                    subWorker.worker.removeEventListener('message', workerMessageHandler);
                    resolve();
                }
            };
            
            subWorker.worker.addEventListener('message', workerMessageHandler);
            
            // Send molecules to process
            subWorker.worker.postMessage({
                type: 'process_chunk',
                molecules: this.serializeMolecules(partition),
                temperature: this.temperature,
                timeScale: this.rules.getConstant('time_scale'),
                damping: this.rules.getConstant('damping'),
                cachedInteractions: this.interactionCache.getSerializableCache()
            });
        });
        
        processingPromises.push(promise);
    }
    
    // Process remaining partitions in main thread if needed
    if (partitionIndex < partitions.length) {
        const remainingMolecules = [];
        for (let i = partitionIndex; i < partitions.length; i++) {
            remainingMolecules.push(...partitions[i]);
        }
        
        const removedIndices = new Set();
        const newMolecules = [];
        
        // Use self interaction cache
        this.updatePhysicsWithCache(remainingMolecules, removedIndices, newMolecules);
        
        // Remove molecules that reacted
        if (removedIndices.size > 0) {
            this.molecules = this.molecules.filter((_, i) => !removedIndices.has(i));
        }
        
        // Add new molecules from reactions
        for (const mol of newMolecules) {
            mol.id = `new-${this.nextMoleculeId++}`;
            this.molecules.push(mol);
        }
    }
    
    // Wait for all workers to complete
    await Promise.all(processingPromises);
}

groupMoleculesByPartition() {
    // Group molecules by spatial grid cells for efficient parallel processing
    const partitions = [];
    const partitionMap = new Map();
    
    // Gather molecules by spatial grid cells
    for (let i = 0; i < this.molecules.length; i++) {
        const mol = this.molecules[i];
        const cellKey = this.spatialGrid.getCellKeyForPosition(mol.position);
        
        if (!partitionMap.has(cellKey)) {
            partitionMap.set(cellKey, []);
        }
        
        partitionMap.get(cellKey).push(mol);
    }
    
    // Merge small partitions to balance workload
    const minPartitionSize = Math.max(10, Math.ceil(this.molecules.length / (NUM_SUB_WORKERS * 2)));
    let currentPartition = [];
    
    for (const molecules of partitionMap.values()) {
        if (currentPartition.length + molecules.length > minPartitionSize && 
            currentPartition.length >= minPartitionSize) {
            // Current partition is full, start a new one
            partitions.push(currentPartition);
            currentPartition = [...molecules];
        } else {
            // Add to current partition
            currentPartition.push(...molecules);
        }
    }
    
    if (currentPartition.length > 0) {
        partitions.push(currentPartition);
    }
    
    // Balance partitions to have similar sizes
    this.balancePartitions(partitions);
    
    return partitions;
}

balancePartitions(partitions) {
    if (partitions.length <= 1) return;
    
    // Calculate average size
    const totalMolecules = partitions.reduce((sum, p) => sum + p.length, 0);
    const targetSize = Math.ceil(totalMolecules / partitions.length);
    
    // Sort partitions by size
    partitions.sort((a, b) => b.length - a.length);
    
    // Balance by moving molecules from larger to smaller partitions
    for (let i = 0; i < partitions.length - 1; i++) {
        if (partitions[i].length <= targetSize) continue;
        
        for (let j = partitions.length - 1; j > i; j--) {
            if (partitions[i].length <= targetSize) break;
            
            while (partitions[i].length > targetSize && 
                   partitions[j].length < targetSize) {
                // Move one molecule to balance
                partitions[j].push(partitions[i].pop());
            }
        }
    }
}

waitForAvailableWorkers() {
    if (this.subWorkers.length === 0) return Promise.resolve();
    
    return new Promise(resolve => {
        const checkWorkers = () => {
            // Check if at least one worker is available
            const anyAvailable = this.subWorkers.some(w => !w.busy);
            
            if (anyAvailable) {
                resolve();
            } else {
                // Check again soon
                setTimeout(checkWorkers, 5);
            }
        };
        
        checkWorkers();
    });
}

serializeMolecules(molecules) {
    return molecules.map(mol => MoleculeSerializer.serialize(mol));
}

getSerializableMolecule(mol) {
    return MoleculeSerializer.serialize(mol);
}

subWorkerStep() {
    // Sub-workers only process chunks they're given
    if (this.workerBusy) return;
    
    // Just signal that we're ready for work
    self.postMessage({ type: 'sub_ready', workerId: this.workerId });
}

updateSpatialGrid() {
    // Clear the grid
    this.spatialGrid.clear();
    
    // Add all molecules to grid
    for (const mol of this.molecules) {
        this.spatialGrid.addMolecule(mol);
    }
}

updatePhysicsWithCache(moleculesToProcess, removedIndices, newMolecules) {
    const timeScale = this.rules.getConstant('time_scale');
    const damping = this.rules.getConstant('damping');
    const now = performance.now();
    
    // Pre-calculate forces for each molecule
    const forces = new Map();
    moleculesToProcess.forEach(mol => forces.set(mol.id, [0, 0, 0]));
    
    // For each molecule
    for (let i = 0; i < moleculesToProcess.length; i++) {
        const mol1 = moleculesToProcess[i];
        if (removedIndices.has(i)) continue;
        
        // Apply random thermal motion
        for (let axis = 0; axis < 3; axis++) {
            mol1.velocity[axis] += (Math.random() - 0.5) * 0.01 * this.temperature;
        }
        
        // Get potential interaction partners from spatial grid
        const neighbors = this.spatialGrid.getPotentialInteractions(mol1);
        
        // Process interactions with neighbors
        for (const mol2 of neighbors) {
            // Skip self and already processed pairs
            if (mol1.id === mol2.id) continue;
            
            // Check if this interaction is cached
            const interactionKey = this.interactionCache.getInteractionKey(mol1.id, mol2.id);
            let force1, force2, distance;
            
            if (this.interactionCache.hasInteraction(interactionKey)) {
                // Use cached interaction data
                const cachedData = this.interactionCache.getInteraction(interactionKey);
                force1 = cachedData.force1;
                force2 = cachedData.force2;
                distance = this.calculateUpdatedDistance(mol1, mol2, cachedData);
            } else {
                // Calculate distance
                distance = this.calculateDistance(mol1, mol2);
                
                // Apply forces
                [force1, force2] = this.applyForces(mol1, mol2);
                
                // Store in cache
                this.interactionCache.storeInteraction(interactionKey, {
                    force1: [...force1],
                    force2: [...force2],
                    distance,
                    mol1Position: [...mol1.position],
                    mol2Position: [...mol2.position],
                    timestamp: now
                });
            }
            
            // Apply calculated forces
            const mol1Force = forces.get(mol1.id);
            for (let k = 0; k < 3; k++) {
                mol1Force[k] += force1[k] * timeScale;
            }
            
            // Apply to mol2 only if it's in our processing list
            if (forces.has(mol2.id)) {
                const mol2Force = forces.get(mol2.id);
                for (let k = 0; k < 3; k++) {
                    mol2Force[k] += force2[k] * timeScale;
                }
            }
            
            // Check for reactions
            if (distance < 1.0 + (mol1.mass + mol2.mass) * 0.1) {
                if (this.shouldReact(mol1, mol2)) {
                    // Find indices in the original array
                    const mol1Index = this.molecules.findIndex(m => m.id === mol1.id);
                    const mol2Index = this.molecules.findIndex(m => m.id === mol2.id);
                    
                    if (mol1Index >= 0 && mol2Index >= 0) {
                        const products = this.processReaction(mol1, mol2);
                        if (products.length > 0) {
                            newMolecules.push(...products);
                            removedIndices.add(mol1Index);
                            removedIndices.add(mol2Index);
                            this.reactionCount++;
                            
                            // Invalidate cache for these molecules
                            this.interactionCache.invalidateForMolecule(mol1.id);
                            this.interactionCache.invalidateForMolecule(mol2.id);
                            break;
                        }
                    }
                }
            }
        }
    }
    
    // Apply forces and update positions
    for (const mol of moleculesToProcess) {
        const molForce = forces.get(mol.id);
        if (!molForce) continue;
        
        // Update velocity with forces
        for (let k = 0; k < 3; k++) {
            mol.velocity[k] += molForce[k] * 0.1;
            mol.position[k] += mol.velocity[k] * timeScale;
            mol.velocity[k] *= damping;
        }
        
        // Check boundaries
        this.enforceBoundaries(mol);
    }
}

calculateUpdatedDistance(mol1, mol2, cachedData) {
    // Estimate current distance based on cached data and current positions
    let sumSquared = 0;
    for (let i = 0; i < 3; i++) {
        const diff = mol2.position[i] - mol1.position[i];
        sumSquared += diff * diff;
    }
    return Math.sqrt(sumSquared);
}

calculateDistance(mol1, mol2) {
    let sumSquared = 0;
    for (let i = 0; i < 3; i++) {
        const diff = mol2.position[i] - mol1.position[i];
        sumSquared += diff * diff;
    }
    return Math.sqrt(sumSquared);
}

applyForces(mol1, mol2) {
    const now = performance.now();
    
    // Cooling period check
    const inCoolingPeriod = this.rules.isInCoolingPeriod && 
                           (this.rules.isInCoolingPeriod(mol1, now) || 
                           this.rules.isInCoolingPeriod(mol2, now));
    
    // Calculate direction vector
    const direction = [
        mol2.position[0] - mol1.position[0],
        mol2.position[1] - mol1.position[1],
        mol2.position[2] - mol1.position[2]
    ];
    
    const distance = Math.sqrt(
        direction[0] * direction[0] + 
        direction[1] * direction[1] + 
        direction[2] * direction[2]
    );
    
    if (distance < 0.001) {
        return [[0, 0, 0], [0, 0, 0]];
    }
    
    // Normalize direction
    const dirNorm = direction.map(d => d / distance);
    
    // Initialize forces
    const force1 = [0, 0, 0];
    const force2 = [0, 0, 0];
    
    // Apply physics rules
    for (const rule of this.rules.interaction_rules) {
        if (rule.condition(mol1.prime_factors, mol2.prime_factors)) {
            let f;
            if (rule.force_function.length === 4) {
                // Force with mass
                f = rule.force_function(dirNorm, distance, mol1.mass, mol2.mass);
            } else {
                // Force with charge
                f = rule.force_function(dirNorm, distance, mol1.charge, mol2.charge);
            }
            
            // Apply rule strength
            for (let i = 0; i < 3; i++) {
                const scaledForce = f[i] * rule.strength;
                force1[i] += scaledForce;
                force2[i] -= scaledForce;
            }
        }
    }
    
    // Add family repulsion if applicable
    if (this.rules.areRelated && this.rules.areRelated(mol1.id, mol2.id)) {
        const repulsionFactor = this.rules.getFamilyRepulsionFactor && 
                               this.rules.getFamilyRepulsionFactor(mol1.id, mol2.id, now);
        
        if (repulsionFactor > 0) {
            const minDistance = this.rules.getConstant('min_distance') || 0.1;
            const effectiveDistance = Math.max(distance, minDistance);
            const repulsiveForce = dirNorm.map(
                x => -x * repulsionFactor / (effectiveDistance ** 1.5)
            );
            
            // Apply to resulting force
            for (let i = 0; i < 3; i++) {
                force1[i] += repulsiveForce[i];
                force2[i] -= repulsiveForce[i];
            }
        }
    }
    
    // Limit maximum force
    const maxForce = this.rules.getConstant('max_force');
    for (let i = 0; i < 3; i++) {
        force1[i] = Math.max(-maxForce, Math.min(maxForce, force1[i]));
        force2[i] = Math.max(-maxForce, Math.min(maxForce, force2[i]));
    }
    
    return [force1, force2];
}

shouldReact(mol1, mol2) {
    // Check cooling period
    const now = performance.now();
    if (this.rules.isInCoolingPeriod && 
       (this.rules.isInCoolingPeriod(mol1, now) || 
        this.rules.isInCoolingPeriod(mol2, now))) {
        return false;
    }
    
    // Calculate relative speed
    const relativeSpeed = Math.sqrt(
        Math.pow(mol1.velocity[0] - mol2.velocity[0], 2) +
        Math.pow(mol1.velocity[1] - mol2.velocity[1], 2) +
        Math.pow(mol1.velocity[2] - mol2.velocity[2], 2)
    );
    
    // Calculate reaction probability
    const baseReactionProb = 0.1 * this.temperature * relativeSpeed;
    const sizeFactor = 1.0 / (1.0 + Math.log(mol1.number + mol2.number));
    
    return Math.random() < baseReactionProb * sizeFactor;
}

// Migliora la funzione processReaction nella classe EnhancedChemistry
processReaction(mol1, mol2) {
    // Look for applicable reaction rules
    for (const rule of this.rules.reaction_rules) {
        // Pass molecules to condition function
        if (rule.condition(mol1.prime_factors, mol2.prime_factors, mol1, mol2)) {
            if (Math.random() < rule.probability * this.temperature) {
                // Determina il tipo di reazione in base alle proprietà delle molecole
                const reactionType = this.determineReactionType(mol1, mol2);
                let products = [];
                
                switch(reactionType) {
                    case 'fusion':
                        // Fusione: crea una molecola combinata
                        products = this.handleFusion(mol1, mol2, rule);
                        break;
                    case 'fission':
                        // Fissione: divide in molecole più piccole
                        products = this.handleFission(mol1, mol2, rule);
                        break;
                    case 'emission':
                        // Emissione: mantiene molecole originali ma emette una particella
                        products = this.handleEmission(mol1, mol2, rule);
                        break;
                    case 'standard':
                    default:
                        // Reazione standard come definita nella regola
                        products = rule.effect(mol1, mol2);
                }
                
                // Set reaction time for visual effect
                const now = performance.now();
                products.forEach(p => {
                    p.setReactionTime(now);
                    
                    // Imposta le relazioni di parentela se necessario
                    if (this.rules.establishRelationship) {
                        this.rules.establishRelationship(mol1.id, p.id, 'parent', now);
                        this.rules.establishRelationship(mol2.id, p.id, 'parent', now);
                    }
                });
                
                // Crea campo quantico per reazioni energetiche se necessario
                if (reactionType === 'fission' || products.length > 2) {
                    this.createQuantumField(mol1, mol2, products);
                }
                
                return products;
            }
        }
    }
    
    return [];
}

determineReactionType(mol1, mol2) {
    // Determina il tipo di reazione in base alle proprietà delle molecole
    const totalMass = mol1.mass + mol2.mass;
    const relativeSpeed = Math.sqrt(
        Math.pow(mol1.velocity[0] - mol2.velocity[0], 2) +
        Math.pow(mol1.velocity[1] - mol2.velocity[1], 2) +
        Math.pow(mol1.velocity[2] - mol2.velocity[2], 2)
    );
    
    // Determina il tipo di reazione in base a massa, carica e velocità
    if (totalMass > 30 && relativeSpeed > 0.5) {
        return 'fission'; // Scissione per molecole grandi ad alta energia
    } else if (mol1.charge * mol2.charge < 0 && Math.abs(mol1.charge) + Math.abs(mol2.charge) > 3) {
        return 'emission'; // Emissione quando cariche opposte forti interagiscono
    } else if (mol1.number <= 10 && mol2.number <= 10 && this.temperature > 1.2) {
        return 'fusion'; // Fusione più probabile per numeri piccoli ad alta temperatura
    } else {
        return 'standard'; // Reazione standard in altri casi
    }
}

handleFusion(mol1, mol2, rule) {
    // Combina le molecole in una più grande
    const midpoint = [
        (mol1.position[0] + mol2.position[0]) / 2,
        (mol1.position[1] + mol2.position[1]) / 2,
        (mol1.position[2] + mol2.position[2]) / 2
    ];
    
    // Combina le velocità proporzionalmente alla massa
    const totalMass = mol1.mass + mol2.mass;
    const combinedVelocity = [
        (mol1.velocity[0] * mol1.mass + mol2.velocity[0] * mol2.mass) / totalMass,
        (mol1.velocity[1] * mol1.mass + mol2.velocity[1] * mol2.mass) / totalMass,
        (mol1.velocity[2] * mol1.mass + mol2.velocity[2] * mol2.mass) / totalMass
    ];
    
    // Usa la regola per determinare il numero della nuova molecola
    let products;
    if (rule.effect) {
        products = rule.effect(mol1, mol2);
    } else {
        // Fallback: crea una molecola con il prodotto dei numeri
        const fusedNumber = mol1.number * mol2.number;
        const fusedMolecule = new PrimeMolecule(
            Math.min(fusedNumber, this.maxNumber),
            midpoint
        );
        fusedMolecule.velocity = combinedVelocity;
        products = [fusedMolecule];
    }
    
    return products;
}

handleFission(mol1, mol2, rule) {
    // Divide le molecole in frammenti più piccoli
    const midpoint = [
        (mol1.position[0] + mol2.position[0]) / 2,
        (mol1.position[1] + mol2.position[1]) / 2,
        (mol1.position[2] + mol2.position[2]) / 2
    ];
    
    // Se abbiamo una regola di effetto, usala
    if (rule.effect) {
        const baseProducts = rule.effect(mol1, mol2);
        if (baseProducts.length > 0) return baseProducts;
    }
    
    // Altrimenti crea una fissione basata sui fattori primi
    const products = [];
    const factorsMol1 = Object.entries(mol1.prime_factors);
    const factorsMol2 = Object.entries(mol2.prime_factors);
    
    // Crea frammenti dai fattori primi di mol1
    for (const [prime, exponent] of factorsMol1) {
        if (exponent > 0) {
            const fragmentNumber = parseInt(prime);
            const fragmentMol = new PrimeMolecule(fragmentNumber, [...midpoint]);
            
            // Aggiungi velocità casuale in direzione opposta al centro
            const direction = [
                Math.random() - 0.5,
                Math.random() - 0.5,
                Math.random() - 0.5
            ];
            const norm = Math.sqrt(direction[0]**2 + direction[1]**2 + direction[2]**2);
            fragmentMol.velocity = [
                direction[0]/norm * 0.5 * this.temperature,
                direction[1]/norm * 0.5 * this.temperature,
                direction[2]/norm * 0.5 * this.temperature
            ];
            
            products.push(fragmentMol);
        }
    }
    
    // Crea frammenti dai fattori primi di mol2
    for (const [prime, exponent] of factorsMol2) {
        if (exponent > 0) {
            const fragmentNumber = parseInt(prime);
            const fragmentMol = new PrimeMolecule(fragmentNumber, [...midpoint]);
            
            // Aggiungi velocità casuale in direzione opposta al centro
            const direction = [
                Math.random() - 0.5,
                Math.random() - 0.5,
                Math.random() - 0.5
            ];
            const norm = Math.sqrt(direction[0]**2 + direction[1]**2 + direction[2]**2);
            fragmentMol.velocity = [
                direction[0]/norm * 0.5 * this.temperature,
                direction[1]/norm * 0.5 * this.temperature,
                direction[2]/norm * 0.5 * this.temperature
            ];
            
            products.push(fragmentMol);
        }
    }
    
    // Limita il numero di prodotti se necessario
    if (products.length > 5) {
        return products.slice(0, 5);
    }
    
    return products.length > 0 ? products : rule.effect(mol1, mol2);
}

handleEmission(mol1, mol2, rule) {
    // Emette una particella ma mantiene le molecole originali modificate
    const products = [];
    
    // Posizione di emissione (media delle posizioni)
    const emissionPoint = [
        (mol1.position[0] + mol2.position[0]) / 2,
        (mol1.position[1] + mol2.position[1]) / 2,
        (mol1.position[2] + mol2.position[2]) / 2
    ];
    
    // Se la regola fornisce un effetto, usalo come base
    if (rule.effect) {
        products.push(...rule.effect(mol1, mol2));
    }
    
    // Se non abbiamo prodotti dalla regola, creiamo versioni modificate delle molecole originali
    if (products.length === 0) {
        // Scegli la molecola più piccola per modificarla
        const smallerMol = mol1.number <= mol2.number ? mol1 : mol2;
        const largerMol = mol1.number > mol2.number ? mol1 : mol2;
        
        // Crea copie modificate delle molecole originali
        const modifiedSmaller = new PrimeMolecule(
            Math.max(2, smallerMol.number - 1),
            [...smallerMol.position]
        );
        modifiedSmaller.velocity = [...smallerMol.velocity];
        
        const modifiedLarger = new PrimeMolecule(
            Math.max(2, largerMol.number - 1),
            [...largerMol.position]
        );
        modifiedLarger.velocity = [...largerMol.velocity];
        
        products.push(modifiedSmaller, modifiedLarger);
    }
    
    // Aggiungi una particella emessa (numero primo piccolo)
    const emittedParticle = new PrimeMolecule(
        this.getSmallPrimeNumber(),
        [...emissionPoint]
    );
    
    // Imposta velocità di emissione in direzione casuale
    const emissionDirection = [
        Math.random() - 0.5,
        Math.random() - 0.5,
        Math.random() - 0.5
    ];
    const norm = Math.sqrt(
        emissionDirection[0]**2 + 
        emissionDirection[1]**2 + 
        emissionDirection[2]**2
    );
    
    const emissionSpeed = 0.8 * this.temperature;
    emittedParticle.velocity = [
        emissionDirection[0]/norm * emissionSpeed,
        emissionDirection[1]/norm * emissionSpeed,
        emissionDirection[2]/norm * emissionSpeed
    ];
    
    // Imposta proprietà speciali per la particella emessa
    emittedParticle._isQuantumField = true;
    emittedParticle._creationTime = performance.now();
    emittedParticle._lifetime = 3000 + Math.random() * 2000; // 3-5 secondi
    
    products.push(emittedParticle);
    
    return products;
}

getSmallPrimeNumber() {
    // Restituisce un numero primo casuale piccolo per particelle emesse
    const smallPrimes = [2, 3, 5, 7, 11, 13];
    return smallPrimes[Math.floor(Math.random() * smallPrimes.length)];
}

createQuantumField(mol1, mol2, products) {
    // Crea un campo quantico temporaneo che influenza le molecole circostanti
    const reactionCenter = [
        (mol1.position[0] + mol2.position[0]) / 2,
        (mol1.position[1] + mol2.position[1]) / 2,
        (mol1.position[2] + mol2.position[2]) / 2
    ];
    
    // Crea una molecola speciale che rappresenta un campo quantico
    const field = new PrimeMolecule(1, [...reactionCenter]);
    field._isQuantumField = true;
    field._creationTime = performance.now();
    field._lifetime = 1500 + Math.random() * 1000; // 1.5-2.5 secondi di vita
    field.velocity = [0, 0, 0]; // Campo stazionario
    
    // Imposta colore e proprietà speciali
    field.color = [0.2, 0.8, 1];
    
    // Aggiungi il campo alla lista delle molecole
    field.id = `field-${this.nextMoleculeId++}`;
    this.molecules.push(field);
}

enforceBoundaries(molecule) {
    const margin = 0.1;
    const bounceElasticity = 0.8;
    
    for (let axis = 0; axis < 3; axis++) {
        const halfSize = this.size/2 - margin;
        
        if (molecule.position[axis] > halfSize) {
            molecule.position[axis] = halfSize;
            molecule.velocity[axis] *= -bounceElasticity;
        } 
        else if (molecule.position[axis] < -halfSize) {
            molecule.position[axis] = -halfSize;
            molecule.velocity[axis] *= -bounceElasticity;
        }
    }
}

manageMoleculeCount() {
    // Limit molecule count for performance
    if (this.molecules.length > 300) {
        // Attention! This is disabled
        // Remove excess molecules, preferring older ones
        //this.molecules = this.molecules.slice(this.molecules.length - 300);
    } else if (this.molecules.length < 50) {
        // Add molecules if there are too few
        this.addRandomMolecules(50 - this.molecules.length);
    }
}

addRandomMolecules(count) {
    for (let i = 0; i < count; i++) {
        // Random position near edges
        const side = Math.floor(Math.random() * 6); // 6 faces of cube
        const pos = [0, 0, 0];
        
        for (let j = 0; j < 3; j++) {
            if (Math.floor(side/2) === j) {
                // On this face
                pos[j] = (side % 2 === 0 ? -1 : 1) * (this.size/2 - 0.2);
            } else {
                // Random position on other axes
                pos[j] = (Math.random() - 0.5) * this.size;
            }
        }
        
        // Create simple molecule
        const number = Math.floor(Math.random() * 48) + 2;
        const mol = new PrimeMolecule(number, pos);
        mol.id = `spawn-${this.nextMoleculeId++}`;
        
        // Velocity toward center
        const dirToCenter = pos.map(p => -p);
        const norm = Math.sqrt(dirToCenter.reduce((s, v) => s + v*v, 0));
        mol.velocity = dirToCenter.map(v => (v/norm) * (0.1 + Math.random() * 0.2));
        
        this.molecules.push(mol);
    }
}

// Method to set temperature
setTemperature(value) {
    this.temperature = value;
}

// Method to set timeScale
setTimeScale(value) {
    if (this.rules) {
        this.rules.setConstant('time_scale', value);
        
        // Propagate to sub-workers
        if (this.isMainWorker && this.subWorkers.length > 0) {
            for (const subWorker of this.subWorkers) {
                subWorker.worker.postMessage({
                    type: 'set_timescale',
                    value: value
                });
            }
        }
    }
}

getRelationshipCount() {
    // Count relationships between molecules
    let count = 0;
    if (this.rules.familyRegistry) {
        for (const [molId, relations] of this.rules.familyRegistry.entries()) {
            count += relations.size;
        }
        return count / 2; // Divide by 2 because each relation is counted twice (bidirectional)
    }
    return 0;
}

cleanup() {
    // Cleanup resources
    if (this.isMainWorker && this.subWorkers.length > 0) {
        // Terminate all sub-workers
        for (const subWorker of this.subWorkers) {
            subWorker.worker.postMessage({type: 'cleanup'});
            subWorker.worker.terminate();
        }
        this.subWorkers = [];
    }
    
    this.molecules = [];
    this.interactionCache.clear();
    this.spatialGrid.clear();
}
}

/**
* InteractionCache - Caches molecule interactions to avoid redundant calculations
*/
class InteractionCache {
constructor(lifetime) {
    this.cache = new Map();
    this.lifetime = lifetime || 10;
    this.invalidationSets = new Map(); // Maps molecule IDs to sets of interaction keys
}

getInteractionKey(mol1Id, mol2Id) {
    // Create consistent key regardless of molecule order
    return mol1Id < mol2Id 
        ? `${mol1Id}:${mol2Id}` 
        : `${mol2Id}:${mol1Id}`;
}

hasInteraction(key) {
    return this.cache.has(key);
}

getInteraction(key) {
    const entry = this.cache.get(key);
    if (entry) {
        entry.age = 0; // Reset age when accessed
        return entry.data;
    }
    return null;
}

storeInteraction(key, data) {
    // Store interaction data with age tracking
    this.cache.set(key, {
        data: data,
        age: 0
    });
    
    // Track which molecules are involved in this interaction
    const [mol1Id, mol2Id] = key.split(':');
    
    // Add to invalidation map for mol1
    if (!this.invalidationSets.has(mol1Id)) {
        this.invalidationSets.set(mol1Id, new Set());
    }
    this.invalidationSets.get(mol1Id).add(key);
    
    // Add to invalidation map for mol2
    if (!this.invalidationSets.has(mol2Id)) {
        this.invalidationSets.set(mol2Id, new Set());
    }
    this.invalidationSets.get(mol2Id).add(key);
}

invalidateForMolecule(molId) {
    // Invalidate all cached interactions involving this molecule
    const keys = this.invalidationSets.get(molId);
    if (keys) {
        for (const key of keys) {
            this.cache.delete(key);
        }
        this.invalidationSets.delete(molId);
    }
}

invalidateForMolecules(molecules) {
    // Invalidate cache for multiple molecules
    for (const mol of molecules) {
        this.invalidateForMolecule(mol.id);
    }
}

ageEntries() {
    // Increase age of all cache entries
    for (const [key, entry] of this.cache.entries()) {
        entry.age++;
        
        // Remove entries that exceed lifetime
        if (entry.age > this.lifetime) {
            this.cache.delete(key);
            
            // Clean up invalidation sets
            const [mol1Id, mol2Id] = key.split(':');
            const set1 = this.invalidationSets.get(mol1Id);
            if (set1) set1.delete(key);
            
            const set2 = this.invalidationSets.get(mol2Id);
            if (set2) set2.delete(key);
        }
    }
}

clear() {
    this.cache.clear();
    this.invalidationSets.clear();
}

getSerializableCache() {
    // Create a serializable version of the cache for sending to sub-workers
    const result = {};
    for (const [key, entry] of this.cache.entries()) {
        if (entry.age <= this.lifetime / 2) { // Only share fresh entries
            result[key] = entry.data;
        }
    }
    return result;
}

importCache(serializedCache) {
    // Import cache from serialized format (from main worker)
    if (!serializedCache) return;
    
    for (const [key, data] of Object.entries(serializedCache)) {
        this.storeInteraction(key, data);
    }
}
}

/**
* SpatialGrid - Spatial partitioning system for efficient neighbor finding
*/
class SpatialGrid {
constructor(worldSize, cellSize) {
    this.worldSize = worldSize;
    this.cellSize = cellSize || 5;
    this.grid = new Map();
    
    // Calculate grid dimensions
    this.gridDimension = Math.ceil(worldSize / this.cellSize);
    this.halfGridDim = Math.floor(this.gridDimension / 2);
}

getCellKeyForPosition(position) {
    // Convert world position to grid cell key
    const x = Math.floor((position[0] + this.worldSize/2) / this.cellSize);
    const y = Math.floor((position[1] + this.worldSize/2) / this.cellSize);
    const z = Math.floor((position[2] + this.worldSize/2) / this.cellSize);
    
    return `${x}:${y}:${z}`;
}

addMolecule(molecule) {
    const cellKey = this.getCellKeyForPosition(molecule.position);
    
    if (!this.grid.has(cellKey)) {
        this.grid.set(cellKey, []);
    }
    
    this.grid.get(cellKey).push(molecule);
}

getPotentialInteractions(molecule) {
    const results = [];
    const cellKey = this.getCellKeyForPosition(molecule.position);
    const [x, y, z] = cellKey.split(':').map(Number);
    
    // Check molecule's cell and all adjacent cells
    for (let dx = -1; dx <= 1; dx++) {
        for (let dy = -1; dy <= 1; dy++) {
            for (let dz = -1; dz <= 1; dz++) {
                const neighborKey = `${x+dx}:${y+dy}:${z+dz}`;
                
                if (this.grid.has(neighborKey)) {
                    const cellMolecules = this.grid.get(neighborKey);
                    results.push(...cellMolecules);
                }
            }
        }
    }
    
    return results;
}

clear() {
    this.grid.clear();
}
}

/**
* Quantum fields update function
*/
/**
 * Funzione migliorata per aggiornare i campi quantici,
 * supportando effetti specifici per fissione, scissione ed emissione
 */
function updateQuantumFields(simulation) {
    const now = performance.now();

    // Filter quantum field molecules
    const quantumFields = simulation.molecules.filter(mol => mol._isQuantumField);
    const fieldsToRemove = [];

    for (const field of quantumFields) {
        // Check if field has expired
        if (now - field._creationTime > field._lifetime) {
            fieldsToRemove.push(field);
            continue;
        }
        
        // Determine field type and strength
        const fieldType = field._fieldType || 'standard';
        const fieldStrength = field._fieldStrength || 1.0;
        
        // Calculate field effect on other molecules
        for (const mol of simulation.molecules) {
            // Skip the field itself and other fields
            if (mol._isQuantumField || mol === field) continue;
            
            // Calculate distance from field
            const distance = simulation.calculateDistance(field, mol);
            const fieldRadius = fieldType === 'fission' ? 8.0 : 5.0;
            
            // Field only affects molecules within certain radius
            if (distance < fieldRadius) {
                // Calculate force direction (from field to molecule)
                const direction = [
                    mol.position[0] - field.position[0],
                    mol.position[1] - field.position[1],
                    mol.position[2] - field.position[2]
                ];
                
                // Normalize direction
                const dirNorm = Math.sqrt(
                    direction[0]**2 + 
                    direction[1]**2 + 
                    direction[2]**2
                );
                
                if (dirNorm > 0.001) {
                    const normalizedDir = direction.map(d => d / dirNorm);
                    
                    // Effetti differenziati in base al tipo di campo
                    switch(fieldType) {
                        case 'fission':
                            // Campo di fissione: forza repulsiva pulsante
                            const fissionPulse = Math.sin((now - field._creationTime) / 100) * 0.5 + 0.5;
                            const fissionStrength = 0.05 * fieldStrength * fissionPulse * (1 - distance/fieldRadius);
                            
                            // Applica forza repulsiva
                            for (let i = 0; i < 3; i++) {
                                mol.velocity[i] += normalizedDir[i] * fissionStrength;
                            }
                            
                            // Occasionalmente causa rotazione
                            if (Math.random() < 0.05 && mol.angularVelocity) {
                                mol.angularVelocity = mol.angularVelocity.map(v => v + (Math.random() - 0.5) * 0.1);
                            }
                            break;
                            
                        case 'emission':
                            // Campo di emissione: oscillazione perpendicolare
                            const emissionFreq = (now - field._creationTime) / 150;
                            const emissionWave = Math.sin(emissionFreq) * Math.cos(emissionFreq * 0.7);
                            const emissionStrength = 0.03 * fieldStrength * (1 - distance/fieldRadius);
                            
                            // Crea un vettore perpendicolare alla direzione
                            const perpVector = [
                                normalizedDir[1] - normalizedDir[2],
                                normalizedDir[2] - normalizedDir[0],
                                normalizedDir[0] - normalizedDir[1]
                            ];
                            
                            // Normalizza il vettore perpendicolare
                            const perpNorm = Math.sqrt(
                                perpVector[0]**2 + 
                                perpVector[1]**2 + 
                                perpVector[2]**2
                            );
                            
                            if (perpNorm > 0.001) {
                                // Applica forza oscillante perpendicolare
                                for (let i = 0; i < 3; i++) {
                                    mol.velocity[i] += (perpVector[i] / perpNorm) * emissionWave * emissionStrength;
                                }
                            }
                            break;
                            
                        case 'standard':
                        default:
                            // Campo standard: oscillazione radiale
                            const oscillation = Math.sin((now - field._creationTime) / 200);
                            const fieldEffect = 0.02 * fieldStrength * (1 - distance/fieldRadius) * oscillation;
                            
                            // Applica forza oscillatoria radiale
                            for (let i = 0; i < 3; i++) {
                                mol.velocity[i] += normalizedDir[i] * fieldEffect;
                            }
                    }
                    
                    // Effetto visivo: le molecole nei campi quantici possono cambiare colore temporaneamente
                    if (Math.random() < 0.02 && distance < fieldRadius * 0.5) {
                        // Leggero shift del colore verso il colore del campo
                        for (let i = 0; i < 3; i++) {
                            mol.color[i] = mol.color[i] * 0.95 + field.color[i] * 0.05;
                        }
                    }
                }
            }
        }
        
        // Aggiorna il campo stesso
        if (fieldType === 'fission') {
            // I campi di fissione si espandono gradualmente
            const expansionFactor = 1.0 + (now - field._creationTime) / field._lifetime * 0.5;
            field._fieldStrength = fieldStrength * (1.0 - (now - field._creationTime) / field._lifetime);
        } else if (fieldType === 'emission') {
            // I campi di emissione pulsano
            field._fieldStrength = fieldStrength * (0.5 + 0.5 * Math.sin((now - field._creationTime) / 300));
        }
    }

    // Remove expired fields
    if (fieldsToRemove.length > 0) {
        simulation.molecules = simulation.molecules.filter(
            mol => !fieldsToRemove.includes(mol));
            
        // Se un campo scompare, a volte può lasciare una particella residua
        for (const field of fieldsToRemove) {
            if (Math.random() < 0.3) {
                const residue = new PrimeMolecule(
                    Math.max(2, Math.floor(Math.random() * 5)),
                    [...field.position]
                );
                
                residue.id = `residue-${simulation.nextMoleculeId++}`;
                residue.velocity = [
                    (Math.random() - 0.5) * 0.2,
                    (Math.random() - 0.5) * 0.2,
                    (Math.random() - 0.5) * 0.2
                ];
                
                simulation.molecules.push(residue);
            }
        }
    }
}

/*
// Main worker code
let simulation;
let previousMoleculeIds = new Set();
let isSubWorker = false;
let workerId = 'main';

// Process messages from main thread or parent worker
onmessage = async function(event) {
try {
    const message = event.data;
    
    if (message.type === 'init_sub') {
        // Initialize as a sub-worker
        handleSubWorkerInitialization(message);
        return;
    }
    
    // For most commands, check if simulation is initialized
    if (message.type !== 'init' && !simulation) {
        console.error(`Worker ${workerId}: simulation not initialized for ${message.type}`);
        postMessage({
            type: 'error',
            message: `Simulation not initialized for ${message.type}`,
            workerId: workerId
        });
        return;
    }
    
    switch (message.type) {
        case 'init':
            await handleInitialization(message);
            break;
            
        case 'process_chunk':
            // Sub-worker specific: process molecule chunk
            handleProcessChunk(message);
            break;
            
        case 'step':
            // Process simulation step
            await simulation.step();
            if (!isSubWorker) {
                sendUpdate();
            }
            break;
            
        case 'cleanup':
            cleanupResources();
            break;
            
        case 'set_temperature':
            // Set temperature
            if (typeof message.value === 'number') {
                simulation.temperature = message.value;
                console.log(`Worker ${workerId}: temperature set to ${message.value}`);
                if (!isSubWorker) {
                    sendUpdate();
                }
            } else {
                console.warn(`Worker ${workerId}: invalid temperature value: ${message.value}`);
            }
            break;
            
        case 'set_timescale':
            // Set timeScale
            if (typeof message.value === 'number') {
                simulation.setTimeScale(message.value);
                console.log(`Worker ${workerId}: timeScale set to ${message.value}`);
            }
            break;
            
        case 'add_molecules':
            // Add molecules
            const count = message.count || 20;
            simulation.addRandomMolecules(count);
            console.log(`Worker ${workerId}: added ${count} new molecules`);
            if (!isSubWorker) {
                sendUpdate();
            }
            break;
            
        case 'set_visualization':
            // Visual mode change (handled by main thread)
            console.log(`Worker ${workerId}: visualization mode set to ${message.mode}`);
            break;
            
        default:
            console.warn(`Worker ${workerId}: unknown message type '${message.type}'`);
    }
} catch (error) {
    console.error(`Worker ${workerId}: error processing message '${event.data?.type}'`, error);
    postMessage({
        type: 'error',
        message: `Error in worker ${workerId}: ${error.message}`,
        workerId: workerId,
        stack: error.stack
    });
}
};

async function handleInitialization(message) {
const { size, moleculeCount, maxNumber, timeScale } = message;
const rules = createCustomRules();
rules.setConstant('time_scale', timeScale || 0.1);

// Initialize simulation
simulation = new EnhancedChemistry(rules, size, moleculeCount, maxNumber);
simulation.isMainWorker = !isSubWorker;
simulation.workerId = workerId;

console.log(`Worker ${workerId}: simulation initialized with ${moleculeCount} molecules, timeScale=${timeScale}`);

// Initialize sub-workers
if (!isSubWorker) {
    try {
        await simulation.initializeSubWorkers();
        console.log(`Initialized ${simulation.subWorkers.length} sub-workers`);
    } catch (error) {
        console.error("Failed to initialize sub-workers:", error);
        // Continue with single-worker mode
    }
}

// Initial step
await simulation.step();

if (!isSubWorker) {
    sendUpdate();
}
}

function handleSubWorkerInitialization(message) {
// Configure as sub-worker
isSubWorker = true;
workerId = message.workerId;
console.log(`Started sub-worker ${workerId}`);

// Create rules instance
const rules = createCustomRules();
rules.setConstant('time_scale', message.timeScale || 0.1);

// Initialize with empty molecule set - we'll receive chunks to process
simulation = new EnhancedChemistry(rules, message.size, 0, message.maxNumber);
simulation.isMainWorker = false;
simulation.workerId = workerId;

// Signal readiness
postMessage({
    type: 'sub_ready',
    workerId: workerId
});
}

function handleProcessChunk(message) {
if (!simulation || !isSubWorker) return;

try {
    // Get molecules to process
    const molecules = message.molecules || [];
    if (molecules.length === 0) {
        postMessage({
            type: 'chunk_processed',
            workerId: workerId,
            processedCount: 0,
            results: {
                moleculeUpdates: [],
                newMolecules: [],
                reactionCount: 0
            }
        });
        return;
    }
    
    // Import cached interactions
    if (message.cachedInteractions) {
        simulation.interactionCache.importCache(message.cachedInteractions);
    }
    
    // Set simulation parameters
    simulation.temperature = message.temperature || 1.0;
    simulation.rules.setConstant('time_scale', message.timeScale || 0.1);
    simulation.rules.setConstant('damping', message.damping || 0.99);
    
    // Create PrimeMolecule instances from serialized data
    const moleculeObjs = molecules.map(molData => 
        MoleculeSerializer.deserialize(molData, PrimeMolecule)
    );
    
    // Process chunk
    const removedIndices = new Set();
    const newMolecules = [];
    const reactionCount = simulation.reactionCount;
    
    // Update molecules in-place
    simulation.updatePhysicsWithCache(moleculeObjs, removedIndices, newMolecules);
    
    // Prepare results
    const results = {
        // Updated molecule states
        moleculeUpdates: moleculeObjs
            .filter((_, i) => !removedIndices.has(i))
            .map(mol => ({
                id: mol.id,
                position: [...mol.position],
                velocity: [...mol.velocity],
                lastReactionTime: mol.lastReactionTime,
                // Include number to help with identification
                number: mol.number
            })),
            
        // Newly created molecules - use serialization helper
        newMolecules: newMolecules.map(mol => MoleculeSerializer.serialize(mol)),
        
        // Number of reactions that occurred
        reactionCount: simulation.reactionCount - reactionCount
    };
    
    // Send results back
    postMessage({
        type: 'chunk_processed',
        workerId: workerId,
        processedCount: molecules.length,
        results: results
    });
    
} catch (error) {
    console.error(`Sub-worker ${workerId}: error processing chunk`, error);
    postMessage({
        type: 'error',
        workerId: workerId,
        message: `Error processing chunk: ${error.message}`
    });
}
}

function sendUpdate() {
try {
    // Prepare data for sending
    const moleculeData = getOptimizedMoleculeData();
    
    // Debug
    console.log(`Worker ${workerId}: sending update with ${moleculeData.molecules.length} molecules, ${moleculeData.removedIds.length} removed`);
    
    // Send message with performance statistics
    postMessage({
        type: 'update',
        molecules: moleculeData,
        temperature: simulation.temperature,
        reactionCount: simulation.reactionCount,
        cacheStats: {
            size: simulation.interactionCache.cache.size,
            hitRate: simulation.interactionCache.hitRate
        },
        workerCount: simulation.subWorkers.length,
        workerStatus: simulation.subWorkers.map(worker => ({
            id: worker.id,
            busy: worker.busy,
            lastProcessed: worker.lastProcessedMolecules
        }))
    });
} catch (error) {
    console.error(`Worker ${workerId}: error sending update`, error);
    postMessage({
        type: 'error',
        message: `Error sending update: ${error.message}`,
        stack: error.stack
    });
}
}

function getOptimizedMoleculeData() {
    const currentMoleculeIds = new Set();
    const result = [];

    // Serialize each molecule with all required properties
    for (const mol of simulation.molecules) {
        // Ensure each molecule has an ID
        const id = mol.id || `mol-${Math.random().toString(36).substring(2, 11)}`;
        mol.id = id;
        currentMoleculeIds.add(id);
        
        // Use serialization helper to ensure all required properties are included
        result.push(MoleculeSerializer.serialize(mol));
    }

    // Calculate removed molecules
    const removedIds = [...previousMoleculeIds].filter(id => !currentMoleculeIds.has(id));

    // Update set for next frame
    previousMoleculeIds = currentMoleculeIds;

    return {
        molecules: result,
        removedIds: removedIds
    };
}

function cleanupResources() {
if (simulation) {
    simulation.cleanup();
}

simulation = null;
previousMoleculeIds.clear();

postMessage({ 
    type: 'cleanup_complete',
    workerId: workerId
});
}
*/

// Extension to PrimeMolecule prototype
PrimeMolecule.prototype.isRelatedTo = function(otherMolecule, rules) {
if (!this.id || !otherMolecule.id) return false;
return rules.areRelated(this.id, otherMolecule.id);
};

///
/// Improved chemistry
///

/**
 * Enhanced Worker Control System
 * Fixes pause functionality and improves caching for molecule interactions
 */

// Extend the EnhancedChemistry class with improved worker control
class ImprovedChemistry extends EnhancedChemistry {
    constructor(rules, size, moleculeCount, maxNumber) {
        super(rules, size, moleculeCount, maxNumber);
        this.isPaused = false;
        
        // Inizializza correttamente la moleculeCache
        this.spatialGrid = null; // Rimuovi la grid-based partitioning
        this.moleculeCache = new MoleculeCentricCache(); // Inizializza la nuova cache
        
        // Altre inizializzazioni
        this.stablePositions = new Map();
        this.positionAccumulators = new Map();
        this.positionSampleCount = new Map();
        this.updateInProgress = false;
        this.pendingPositionUpdate = false;
      }
      
      // Override del metodo step
      async step() {
        if (this.isPaused) {
          if (!this.isSubWorker) {
            sendUpdate();
          }
          return;
        }
    
        // Imposta flag di aggiornamento in corso
        this.updateInProgress = true;
        
        const timeScale = this.rules.getConstant('time_scale');
        this.accumulatedTime += timeScale;
        
        // Aggiorna temperatura
        this.temperature = 1.0 + 0.2 * Math.sin(performance.now() / 5000);
        
        // Arrays per elaborazione efficiente
        const newMolecules = [];
        const removedIndices = new Set();
        
        // Elabora la fisica usando sub-worker se disponibili
        if (this.isMainWorker && this.subWorkers.length > 0) {
          await this.distributeWorkToSubWorkers();
        } else {
          // Fallback a elaborazione single-thread
          this.updatePhysics(removedIndices, newMolecules);
        }
        
        // Filtra molecole rimosse
        if (removedIndices.size > 0) {
          this.molecules = this.molecules.filter((_, i) => !removedIndices.has(i));
        }
        
        // Aggiungi nuove molecole
        for (const mol of newMolecules) {
          mol.id = `new-${this.nextMoleculeId++}`;
          this.molecules.push(mol);
        }
        
        // Reset tempo accumulato
        if (this.accumulatedTime >= 1.0) {
          this.accumulatedTime = 0.0;
          this.moleculeCache.cleanupStaleRelationships();
        }
        
        // Limita numero di molecole per performance
        this.manageMoleculeCount();
        
        // Aggiorna campi quantici
        updateQuantumFields(this);
        
        // Pulizia periodica relazioni obsolete
        if (Math.random() < 0.05) {
          this.rules.cleanupOldRelations && this.rules.cleanupOldRelations();
        }
        
        // Accumula posizioni per tutte le molecole
        for (const mol of this.molecules) {
          this.accumulatePosition(mol);
        }
        
        // Aggiornamento completato
        this.updateInProgress = false;
        
        // Invia aggiornamento solo quando tutte le posizioni sono state accumulate
        if (this.isMainWorker) {
          if (this.pendingPositionUpdate) {
            // Stabilizza le posizioni prima di inviare
            this.stabilizePositions();
            sendUpdate();
            
            // Reset flag e buffer
            this.pendingPositionUpdate = false;
            this.positionAccumulators.clear();
            this.positionSampleCount.clear();
          } else {
            // Imposta flag per aggiornamento posizione al prossimo ciclo
            this.pendingPositionUpdate = true;
          }
        }
      }
      
      // Accumula posizione di una molecola per una media stabile
      accumulatePosition(molecule) {
        const id = molecule.id;
        
        // Inizializza accumulatore se necessario
        if (!this.positionAccumulators.has(id)) {
          this.positionAccumulators.set(id, [0, 0, 0]);
          this.positionSampleCount.set(id, 0);
        }
        
        // Aggiungi posizione attuale all'accumulatore
        const accumulator = this.positionAccumulators.get(id);
        const count = this.positionSampleCount.get(id);
        
        for (let i = 0; i < 3; i++) {
          accumulator[i] += molecule.position[i];
        }
        
        this.positionSampleCount.set(id, count + 1);
      }
      
      // Calcola posizioni stabili come media delle posizioni accumulate
      stabilizePositions() {
        this.stablePositions.clear();
        
        for (const mol of this.molecules) {
          const id = mol.id;
          
          if (this.positionAccumulators.has(id) && this.positionSampleCount.has(id)) {
            const accumulator = this.positionAccumulators.get(id);
            const count = this.positionSampleCount.get(id);
            
            if (count > 0) {
              // Calcola posizione media
              const stablePosition = [
                accumulator[0] / count,
                accumulator[1] / count,
                accumulator[2] / count
              ];
              
              // Memorizza posizione stabile
              this.stablePositions.set(id, stablePosition);
              
              // Aggiorna posizione della molecola con quella stabile
              // (solo per la visualizzazione, non per i calcoli fisici)
              const originalPosition = [...mol.position];
              mol.position = stablePosition;
              
              // Preserva la vera posizione fisica in una proprietà separata
              mol._physicsPosition = originalPosition;
            }
          }
        }
      }
      
      // Override di calculateDistance per usare posizioni fisiche reali
      calculateDistance(mol1, mol2) {
        // Usa _physicsPosition se disponibile, altrimenti usa position
        const pos1 = mol1._physicsPosition || mol1.position;
        const pos2 = mol2._physicsPosition || mol2.position;
        
        let sumSquared = 0;
        for (let i = 0; i < 3; i++) {
          const diff = pos2[i] - pos1[i];
          sumSquared += diff * diff;
        }
        return Math.sqrt(sumSquared);
      }
      
      // Override setTemperature per rispettare flag updateInProgress
      setTemperature(value) {
        this.temperature = value;
        
        // Se aggiornamento in corso, attendi
        if (!this.isMainWorker || !this.updateInProgress) {
          propagateTemperature(value);
        }
      }
      
      // Override setTimeScale per rispettare flag updateInProgress
      setTimeScale(value) {
        if (this.rules) {
          this.rules.setConstant('time_scale', value);
          
          // Se aggiornamento in corso, attendi
          if (!this.isMainWorker || !this.updateInProgress) {
            propagateTimeScale(value);
          }
        }
      }
      
      // Nuovo metodo per ripristinare le posizioni fisiche reali
      restorePhysicsPositions() {
        for (const mol of this.molecules) {
          if (mol._physicsPosition) {
            mol.position = [...mol._physicsPosition];
            delete mol._physicsPosition;
          }
        }
      }
      
      // Override di cleanup
      cleanup() {
        this.isPaused = false;
        this.updateInProgress = false;
        this.pendingPositionUpdate = false;
        this.stablePositions.clear();
        this.positionAccumulators.clear();
        this.positionSampleCount.clear();
        
        // Cleanup resources
        if (this.isMainWorker && this.subWorkers.length > 0) {
          for (const subWorker of this.subWorkers) {
            subWorker.worker.postMessage({type: 'cleanup'});
            subWorker.worker.terminate();
          }
          this.subWorkers = [];
        }
        
        this.molecules = [];
        this.moleculeCache.clear();
      }
  
    // Set pause state
    setPauseState(isPaused) {
      this.isPaused = isPaused;
      
      // Propagate pause state to sub-workers
      if (this.isMainWorker && this.subWorkers.length > 0) {
        for (const subWorker of this.subWorkers) {
          subWorker.worker.postMessage({
            type: 'set_pause',
            isPaused: isPaused
          });
        }
      }
    }
  
    // Override distributeWorkToSubWorkers to include pause state
    async distributeWorkToSubWorkers() {
      if (this.subWorkers.length === 0 || this.isPaused) return;
      
      // Group molecules based on natural clustering instead of spatial grid
      const moleculeClusters = this.groupMoleculesByProximity();
      
      // Wait for all sub-workers to be ready
      await this.waitForAvailableWorkers();
      
      // Distribute clusters to available workers
      let clusterIndex = 0;
      const processingPromises = [];
      
      for (const subWorker of this.subWorkers) {
        if (clusterIndex >= moleculeClusters.length) break;
        
        const cluster = moleculeClusters[clusterIndex++];
        if (!cluster || cluster.length === 0) continue;
        
        subWorker.busy = true;
        
        const promise = new Promise(resolve => {
          const workerMessageHandler = (event) => {
            if (event.data.type === 'chunk_processed') {
              subWorker.worker.removeEventListener('message', workerMessageHandler);
              resolve();
            }
          };
          
          subWorker.worker.addEventListener('message', workerMessageHandler);
          
          // Send molecules to process with cached relationships
          const serializedMolecules = this.serializeMolecules(cluster);
          const moleculeRelationships = this.moleculeCache.getRelationshipsForMolecules(cluster);
          
          subWorker.worker.postMessage({
            type: 'process_chunk',
            molecules: serializedMolecules,
            temperature: this.temperature,
            timeScale: this.rules.getConstant('time_scale'),
            damping: this.rules.getConstant('damping'),
            cachedRelationships: moleculeRelationships,
            isPaused: this.isPaused
          });
        });
        
        processingPromises.push(promise);
      }
      
      // Process remaining clusters in main thread if needed
      if (clusterIndex < moleculeClusters.length) {
        const remainingMolecules = [];
        for (let i = clusterIndex; i < moleculeClusters.length; i++) {
          remainingMolecules.push(...moleculeClusters[i]);
        }
        
        const removedIndices = new Set();
        const newMolecules = [];
        
        // Use molecule-centric caching
        this.updatePhysicsWithMoleculeCache(remainingMolecules, removedIndices, newMolecules);
        
        // Remove molecules that reacted
        if (removedIndices.size > 0) {
          this.molecules = this.molecules.filter((_, i) => !removedIndices.has(i));
        }
        
        // Add new molecules from reactions
        for (const mol of newMolecules) {
          mol.id = `new-${this.nextMoleculeId++}`;
          this.molecules.push(mol);
        }
      }
      
      // Wait for all workers to complete
      await Promise.all(processingPromises);
    }
  
    // New method to group molecules by proximity instead of grid
    groupMoleculesByProximity() {
      if (this.molecules.length <= 1) {
        return [this.molecules];
      }
      
      const clusters = [];
      const visited = new Set();
      const proximityThreshold = this.size * 0.15; // Adaptive threshold based on world size
      
      // Helper function to find molecules in proximity
      const findCluster = (startMol) => {
        const cluster = [startMol];
        const queue = [startMol];
        visited.add(startMol.id);
        
        while (queue.length > 0) {
          const current = queue.shift();
          
          // Find nearby molecules
          for (const mol of this.molecules) {
            if (visited.has(mol.id)) continue;
            
            const distance = this.calculateDistance(current, mol);
            if (distance < proximityThreshold) {
              cluster.push(mol);
              queue.push(mol);
              visited.add(mol.id);
            }
          }
        }
        
        return cluster;
      };
      
      // Form natural clusters based on molecule proximity
      for (const mol of this.molecules) {
        if (!visited.has(mol.id)) {
          const cluster = findCluster(mol);
          clusters.push(cluster);
        }
      }
      
      // Balance cluster sizes
      this.balanceClusters(clusters);
      
      return clusters;
    }
  
    // Modified to balance clusters more effectively
    balanceClusters(clusters) {
      if (clusters.length <= 1) return;
      
      // Calculate ideal cluster size
      const totalMolecules = clusters.reduce((sum, c) => sum + c.length, 0);
      const idealSize = Math.max(
        10,
        Math.min(50, Math.ceil(totalMolecules / (this.subWorkers.length || 1)))
      );
      
      // Sort clusters by size (largest first)
      clusters.sort((a, b) => b.length - a.length);
      
      // Split overly large clusters
      let i = 0;
      while (i < clusters.length) {
        if (clusters[i].length > idealSize * 1.5) {
          // Find natural split point using molecule proximity
          const newClusters = this.splitCluster(clusters[i], idealSize);
          
          // Replace original cluster with first part
          clusters[i] = newClusters[0];
          
          // Add remaining parts as new clusters
          clusters.push(...newClusters.slice(1));
          
          // Resort clusters
          clusters.sort((a, b) => b.length - a.length);
        } else {
          i++;
        }
      }
      
      // Merge very small clusters
      while (clusters.length > 1) {
        const smallest = clusters.pop();
        if (smallest.length < idealSize * 0.5) {
          // Find closest cluster to merge with
          let closestClusterIndex = -1;
          let minDistance = Infinity;
          
          for (let j = 0; j < clusters.length; j++) {
            const distance = this.calculateInterClusterDistance(smallest, clusters[j]);
            if (distance < minDistance) {
              minDistance = distance;
              closestClusterIndex = j;
            }
          }
          
          if (closestClusterIndex >= 0) {
            clusters[closestClusterIndex].push(...smallest);
          } else {
            // Put it back if no good merge candidate
            clusters.push(smallest);
            break;
          }
        } else {
          // Put it back if not too small
          clusters.push(smallest);
          break;
        }
      }
    }
  
    // Helper to split a large cluster
    splitCluster(cluster, targetSize) {
      // Use K-means like approach to split clusters
      const result = [];
      let remaining = [...cluster];
      
      while (remaining.length > targetSize) {
        // Find molecule farthest from center of mass
        const centerOfMass = this.calculateCenterOfMass(remaining);
        let farthestIndex = 0;
        let maxDistance = 0;
        
        for (let i = 0; i < remaining.length; i++) {
          const distance = this.calculateDistanceFromPoint(remaining[i], centerOfMass);
          if (distance > maxDistance) {
            maxDistance = distance;
            farthestIndex = i;
          }
        }
        
        // Use the farthest molecule as a new center
        const newCenter = remaining[farthestIndex];
        const newCluster = [newCenter];
        const newRemaining = [];
        
        // Add nearest neighbors until we reach target size
        remaining.splice(farthestIndex, 1);
        
        // Sort remaining by distance to new center
        const distanceMap = remaining.map((mol, index) => ({
          index,
          distance: this.calculateDistance(mol, newCenter)
        }));
        
        distanceMap.sort((a, b) => a.distance - b.distance);
        
        // Take closest neighbors up to target size
        const neighborsToTake = Math.min(targetSize - 1, remaining.length);
        for (let i = 0; i < remaining.length; i++) {
          if (i < neighborsToTake) {
            newCluster.push(remaining[distanceMap[i].index]);
          } else {
            newRemaining.push(remaining[distanceMap[i].index]);
          }
        }
        
        result.push(newCluster);
        remaining = newRemaining;
      }
      
      // Add any remaining molecules as the last cluster
      if (remaining.length > 0) {
        result.push(remaining);
      }
      
      return result;
    }
  
    // Helper to calculate center of mass
    calculateCenterOfMass(molecules) {
      if (molecules.length === 0) return [0, 0, 0];
      
      const center = [0, 0, 0];
      let totalMass = 0;
      
      for (const mol of molecules) {
        const mass = mol.mass || 1;
        totalMass += mass;
        
        for (let i = 0; i < 3; i++) {
          center[i] += mol.position[i] * mass;
        }
      }
      
      if (totalMass > 0) {
        for (let i = 0; i < 3; i++) {
          center[i] /= totalMass;
        }
      }
      
      return center;
    }
  
    // Helper to calculate distance from a point
    calculateDistanceFromPoint(molecule, point) {
      let sumSquared = 0;
      for (let i = 0; i < 3; i++) {
        const diff = molecule.position[i] - point[i];
        sumSquared += diff * diff;
      }
      return Math.sqrt(sumSquared);
    }
  
    // Helper to calculate minimum distance between clusters
    calculateInterClusterDistance(cluster1, cluster2) {
      let minDistance = Infinity;
      
      for (const mol1 of cluster1) {
        for (const mol2 of cluster2) {
          const distance = this.calculateDistance(mol1, mol2);
          if (distance < minDistance) {
            minDistance = distance;
          }
        }
      }
      
      return minDistance;
    }
  
    // New physics update method using molecule-centric cache
    updatePhysicsWithMoleculeCache(moleculesToProcess, removedIndices, newMolecules) {
      const timeScale = this.rules.getConstant('time_scale');
      const damping = this.rules.getConstant('damping');
      const now = performance.now();
      
      // Pre-calculate forces for each molecule
      const forces = new Map();
      moleculesToProcess.forEach(mol => forces.set(mol.id, [0, 0, 0]));
      
      // Calculate which molecules need relationship calculation
      const uncachedMolecules = this.moleculeCache.getUncachedMolecules(moleculesToProcess);
      
      // Process uncached molecules first (establish new relationships)
      for (const mol1 of uncachedMolecules) {
        for (const mol2 of moleculesToProcess) {
          if (mol1.id === mol2.id) continue;
          
          // Check if they're close enough to establish a relationship
          const distance = this.calculateDistance(mol1, mol2);
          const interactionThreshold = this.calculateInteractionThreshold(mol1, mol2);
          
          if (distance < interactionThreshold) {
            // Create new relationship in cache
            this.moleculeCache.createRelationship(mol1, mol2, distance, now);
          }
        }
      }
      
      // Now process all molecules with their cached relationships
      for (let i = 0; i < moleculesToProcess.length; i++) {
        const mol1 = moleculesToProcess[i];
        if (removedIndices.has(i)) continue;
        
        // Apply random thermal motion based on temperature
        for (let axis = 0; axis < 3; axis++) {
          mol1.velocity[axis] += (Math.random() - 0.5) * 0.01 * this.temperature;
        }
        
        // Get cached relationships for this molecule
        const relationships = this.moleculeCache.getRelationshipsForMolecule(mol1);
        
        // Process each relationship
        for (const rel of relationships) {
          // Find the other molecule in our processing list
          const mol2Index = moleculesToProcess.findIndex(m => m.id === rel.otherId);
          if (mol2Index < 0) continue; // Not in our batch
          
          const mol2 = moleculesToProcess[mol2Index];
          if (removedIndices.has(mol2Index)) continue;
          
          // Update relationship distance
          const currentDistance = this.calculateDistance(mol1, mol2);
          rel.distance = currentDistance;
          rel.lastUpdated = now;
          
          // Check if relationship should be maintained
          const maxDistance = this.calculateMaxRelationshipDistance(mol1, mol2);
          
          if (currentDistance > maxDistance) {
            // Molecules moved too far apart, remove relationship
            this.moleculeCache.removeRelationship(mol1.id, mol2.id);
            continue;
          }
          
          // Calculate and apply forces
          const [force1, force2] = this.calculateForces(mol1, mol2, currentDistance, now);
          
          // Apply calculated forces
          const mol1Force = forces.get(mol1.id);
          const mol2Force = forces.get(mol2.id);
          
          for (let k = 0; k < 3; k++) {
            mol1Force[k] += force1[k] * timeScale;
            mol2Force[k] += force2[k] * timeScale;
          }
          
          // Check for reactions
          if (currentDistance < this.calculateReactionDistance(mol1, mol2)) {
            if (this.shouldReact(mol1, mol2)) {
              // Find indices in the original array
              const mol1Index = this.molecules.findIndex(m => m.id === mol1.id);
              const mol2Index = this.molecules.findIndex(m => m.id === mol2.id);
              
              if (mol1Index >= 0 && mol2Index >= 0) {
                const products = this.processReaction(mol1, mol2);
                if (products.length > 0) {
                  newMolecules.push(...products);
                  removedIndices.add(mol1Index);
                  removedIndices.add(mol2Index);
                  this.reactionCount++;
                  
                  // Remove relationships for reacting molecules
                  this.moleculeCache.removeAllRelationshipsForMolecule(mol1.id);
                  this.moleculeCache.removeAllRelationshipsForMolecule(mol2.id);
                  break;
                }
              }
            }
          }
        }
      }
      
      // Apply forces and update positions
      for (const mol of moleculesToProcess) {
        const molForce = forces.get(mol.id);
        if (!molForce) continue;
        
        // Update velocity with forces
        for (let k = 0; k < 3; k++) {
          mol.velocity[k] += molForce[k] * 0.1;
          mol.position[k] += mol.velocity[k] * timeScale;
          mol.velocity[k] *= damping;
        }
        
        // Check boundaries
        this.enforceBoundaries(mol);
      }
    }
  
    // Helper methods for molecule relationships
    calculateInteractionThreshold(mol1, mol2) {
      // Dynamic threshold based on molecule properties
      return Math.max(
        5.0, // Minimum threshold
        (mol1.mass + mol2.mass) * 1.2
      );
    }
  
    calculateMaxRelationshipDistance(mol1, mol2) {
      // Maximum distance to maintain a relationship
      return Math.max(
        10.0, // Minimum threshold
        (mol1.mass + mol2.mass) * 2.0
      );
    }
  
    calculateReactionDistance(mol1, mol2) {
      // Distance at which molecules can react
      return 1.0 + (mol1.mass + mol2.mass) * 0.1;
    }
  
    // Calculate forces between molecules
    calculateForces(mol1, mol2, distance, now) {
      // Calculate direction vector
      const direction = [
        mol2.position[0] - mol1.position[0],
        mol2.position[1] - mol1.position[1],
        mol2.position[2] - mol1.position[2]
      ];
      
      if (distance < 0.001) {
        return [[0, 0, 0], [0, 0, 0]];
      }
      
      // Normalize direction
      const dirNorm = direction.map(d => d / distance);
      
      // Initialize forces
      const force1 = [0, 0, 0];
      const force2 = [0, 0, 0];
      
      // Apply physics rules
      for (const rule of this.rules.interaction_rules) {
        if (rule.condition(mol1.prime_factors, mol2.prime_factors)) {
          let f;
          if (rule.force_function.length === 4) {
            // Force with mass
            f = rule.force_function(dirNorm, distance, mol1.mass, mol2.mass);
          } else {
            // Force with charge
            f = rule.force_function(dirNorm, distance, mol1.charge, mol2.charge);
          }
          
          // Apply rule strength
          for (let i = 0; i < 3; i++) {
            const scaledForce = f[i] * rule.strength;
            force1[i] += scaledForce;
            force2[i] -= scaledForce;
          }
        }
      }
      
      // Add family repulsion if applicable
      if (this.rules.areRelated && this.rules.areRelated(mol1.id, mol2.id)) {
        const repulsionFactor = this.rules.getFamilyRepulsionFactor && 
                               this.rules.getFamilyRepulsionFactor(mol1.id, mol2.id, now);
        
        if (repulsionFactor > 0) {
          const minDistance = this.rules.getConstant('min_distance') || 0.1;
          const effectiveDistance = Math.max(distance, minDistance);
          const repulsiveForce = dirNorm.map(
            x => -x * repulsionFactor / (effectiveDistance ** 1.5)
          );
          
          // Apply to resulting force
          for (let i = 0; i < 3; i++) {
            force1[i] += repulsiveForce[i];
            force2[i] -= repulsiveForce[i];
          }
        }
      }
      
      // Limit maximum force
      const maxForce = this.rules.getConstant('max_force');
      for (let i = 0; i < 3; i++) {
        force1[i] = Math.max(-maxForce, Math.min(maxForce, force1[i]));
        force2[i] = Math.max(-maxForce, Math.min(maxForce, force2[i]));
      }
      
      return [force1, force2];
    }
  
    // Override cleanup to handle pause state
    cleanup() {
      // Reset pause state
      this.isPaused = false;
      
      // Cleanup resources
      if (this.isMainWorker && this.subWorkers.length > 0) {
        // Terminate all sub-workers
        for (const subWorker of this.subWorkers) {
          subWorker.worker.postMessage({type: 'cleanup'});
          subWorker.worker.terminate();
        }
        this.subWorkers = [];
      }
      
      this.molecules = [];
      this.moleculeCache.clear();
    }
  }
  
/**
 * MoleculeCentricCache - New molecule-centric caching system
 * Maintains relationships between molecules based on their proximity
 * rather than spatial grid partitioning
 */
class MoleculeCentricCache {
    constructor() {
      // Map di ID molecola alle sue relazioni
      this.moleculeRelationships = new Map();
      
      // Soglia di obsolescenza relazioni (ms)
      this.stalenessThreshold = 5000;
      
      // Massimo numero di relazioni per molecola
      this.maxRelationshipsPerMolecule = 20;
    }
    
    // Restituisce tutte le relazioni per una molecola
    getRelationshipsForMolecule(molecule) {
      if (!molecule || !molecule.id) return [];
      return this.moleculeRelationships.get(molecule.id) || [];
    }
    
    // Restituisce tutte le relazioni per una lista di molecole
    getRelationshipsForMolecules(molecules) {
      const result = {};
      
      for (const mol of molecules) {
        if (!mol.id) continue;
        
        const relationships = this.getRelationshipsForMolecule(mol);
        if (relationships.length > 0) {
          result[mol.id] = relationships.map(rel => ({
            otherId: rel.otherId,
            distance: rel.distance,
            lastUpdated: rel.lastUpdated
          }));
        }
      }
      
      return result;
    }
    
    // Crea una nuova relazione tra molecole
    createRelationship(mol1, mol2, distance, timestamp) {
      if (!mol1.id || !mol2.id || mol1.id === mol2.id) return;
      
      // Verifica se la relazione esiste già
      if (this.hasRelationship(mol1.id, mol2.id)) {
        // Aggiorna relazione esistente
        this.updateRelationship(mol1.id, mol2.id, distance, timestamp);
        return;
      }
      
      // Crea nuove relazioni in entrambe le direzioni
      this._ensureRelationshipsList(mol1.id);
      this._ensureRelationshipsList(mol2.id);
      
      // Verifica se è necessario eliminare relazioni
      this._pruneRelationshipsIfNeeded(mol1.id);
      this._pruneRelationshipsIfNeeded(mol2.id);
      
      // Crea relazione bidirezionale
      const rel1 = {
        otherId: mol2.id,
        distance: distance,
        lastUpdated: timestamp
      };
      
      const rel2 = {
        otherId: mol1.id,
        distance: distance,
        lastUpdated: timestamp
      };
      
      this.moleculeRelationships.get(mol1.id).push(rel1);
      this.moleculeRelationships.get(mol2.id).push(rel2);
    }
    
    // Verifica se esiste una relazione
    hasRelationship(mol1Id, mol2Id) {
      if (!this.moleculeRelationships.has(mol1Id)) return false;
      
      return this.moleculeRelationships.get(mol1Id).some(rel => 
        rel.otherId === mol2Id
      );
    }
    
    // Aggiorna una relazione esistente
    updateRelationship(mol1Id, mol2Id, distance, timestamp) {
      if (!this.hasRelationship(mol1Id, mol2Id)) return;
      
      // Aggiorna da mol1 a mol2
      const rels1 = this.moleculeRelationships.get(mol1Id);
      const rel1Index = rels1.findIndex(rel => rel.otherId === mol2Id);
      
      if (rel1Index >= 0) {
        rels1[rel1Index].distance = distance;
        rels1[rel1Index].lastUpdated = timestamp;
      }
      
      // Aggiorna da mol2 a mol1
      if (this.moleculeRelationships.has(mol2Id)) {
        const rels2 = this.moleculeRelationships.get(mol2Id);
        const rel2Index = rels2.findIndex(rel => rel.otherId === mol1Id);
        
        if (rel2Index >= 0) {
          rels2[rel2Index].distance = distance;
          rels2[rel2Index].lastUpdated = timestamp;
        }
      }
    }
    
    // Rimuove una relazione
    removeRelationship(mol1Id, mol2Id) {
      // Rimuove da mol1 a mol2
      if (this.moleculeRelationships.has(mol1Id)) {
        const rels1 = this.moleculeRelationships.get(mol1Id);
        const filteredRels1 = rels1.filter(rel => rel.otherId !== mol2Id);
        this.moleculeRelationships.set(mol1Id, filteredRels1);
      }
      
      // Rimuove da mol2 a mol1
      if (this.moleculeRelationships.has(mol2Id)) {
        const rels2 = this.moleculeRelationships.get(mol2Id);
        const filteredRels2 = rels2.filter(rel => rel.otherId !== mol1Id);
        this.moleculeRelationships.set(mol2Id, filteredRels2);
      }
    }
    
    // Rimuove tutte le relazioni per una molecola
    removeAllRelationshipsForMolecule(molId) {
      if (!this.moleculeRelationships.has(molId)) return;
      
      // Ottieni lista delle altre molecole
      const relationships = this.moleculeRelationships.get(molId);
      const otherIds = relationships.map(rel => rel.otherId);
      
      // Rimuovi questa molecola dalle relazioni di tutte le altre molecole
      for (const otherId of otherIds) {
        if (this.moleculeRelationships.has(otherId)) {
          const otherRels = this.moleculeRelationships.get(otherId);
          const filteredOtherRels = otherRels.filter(rel => rel.otherId !== molId);
          this.moleculeRelationships.set(otherId, filteredOtherRels);
        }
      }
      
      // Rimuovi le relazioni di questa molecola
      this.moleculeRelationships.delete(molId);
    }
    
    // Ottieni molecole senza cache da una lista
    getUncachedMolecules(molecules) {
      return molecules.filter(mol => 
        !mol.id || !this.moleculeRelationships.has(mol.id) || 
        this.moleculeRelationships.get(mol.id).length === 0
      );
    }
    
    // Pulisci relazioni obsolete
    cleanupStaleRelationships() {
      const now = performance.now();
      
      for (const [molId, relationships] of this.moleculeRelationships.entries()) {
        // Filtra le relazioni obsolete
        const freshRelationships = relationships.filter(rel => 
          now - rel.lastUpdated < this.stalenessThreshold
        );
        
        if (freshRelationships.length !== relationships.length) {
          // Aggiorna con le sole relazioni fresche
          this.moleculeRelationships.set(molId, freshRelationships);
        }
      }
      
      // Rimuovi liste di relazioni vuote
      for (const [molId, relationships] of this.moleculeRelationships.entries()) {
        if (relationships.length === 0) {
          this.moleculeRelationships.delete(molId);
        }
      }
    }
    
    // Cancella tutte le relazioni
    clear() {
      this.moleculeRelationships.clear();
    }
    
    // Helper per assicurarsi che una molecola abbia una lista di relazioni
    _ensureRelationshipsList(molId) {
      if (!this.moleculeRelationships.has(molId)) {
        this.moleculeRelationships.set(molId, []);
      }
    }
    
    // Helper per potare relazioni in eccesso
    _pruneRelationshipsIfNeeded(molId) {
      if (!this.moleculeRelationships.has(molId)) return;
      
      const relationships = this.moleculeRelationships.get(molId);
      if (relationships.length >= this.maxRelationshipsPerMolecule) {
        // Ordina per recenza (più recenti prima)
        relationships.sort((a, b) => b.lastUpdated - a.lastUpdated);
        
        // Mantieni solo le relazioni più recenti
        this.moleculeRelationships.set(
          molId, 
          relationships.slice(0, this.maxRelationshipsPerMolecule - 1)
        );
      }
    }
  }
  
  /**
   * Modified Main Worker Code
   * Includes support for pausing/resuming and improved caching
   */
  
  let simulation;
  let previousMoleculeIds = new Set();
  let isSubWorker = false;
  let workerId = 'main';
  
  // Process messages from main thread or parent worker
  onmessage = async function(event) {
    try {
      const message = event.data;
      
      if (message.type === 'init_sub') {
        // Initialize as a sub-worker
        handleSubWorkerInitialization(message);
        return;
      }
      
      // For most commands, check if simulation is initialized
      if (message.type !== 'init' && !simulation) {
        console.error(`Worker ${workerId}: simulation not initialized for ${message.type}`);
        postMessage({
          type: 'error',
          message: `Simulation not initialized for ${message.type}`,
          workerId: workerId
        });
        return;
      }
      
      switch (message.type) {
        case 'init':
          await handleInitialization(message);
          break;
          
        case 'process_chunk':
          // Sub-worker specific: process molecule chunk
          handleProcessChunk(message);
          break;
          
        case 'step':
          // Process simulation step
          await simulation.step();
          if (!isSubWorker) {
            sendUpdate();
          }
          break;
          
        case 'cleanup':
          cleanupResources();
          break;
          
        case 'set_temperature':
          // Set temperature
          if (typeof message.value === 'number') {
            simulation.temperature = message.value;
            console.log(`Worker ${workerId}: temperature set to ${message.value}`);
            if (!isSubWorker) {
              sendUpdate();
            }
          } else {
            console.warn(`Worker ${workerId}: invalid temperature value: ${message.value}`);
          }
          break;
          
        case 'set_timescale':
          // Set timeScale
          if (typeof message.value === 'number') {
            simulation.setTimeScale(message.value);
            console.log(`Worker ${workerId}: timeScale set to ${message.value}`);
          }
          break;
          
        case 'set_pause':
          // NEW: Set pause state
          if (typeof message.isPaused === 'boolean') {
            simulation.setPauseState(message.isPaused);
            console.log(`Worker ${workerId}: pause state set to ${message.isPaused}`);
            
            // If main worker, still send updates when paused
            if (!isSubWorker && message.isPaused) {
              sendUpdate();
            }
          }
          break;
          
        case 'add_molecules':
          // Add molecules
          const count = message.count || 20;
          simulation.addRandomMolecules(count);
          console.log(`Worker ${workerId}: added ${count} new molecules`);
          if (!isSubWorker) {
            sendUpdate();
          }
          break;
          
        case 'set_visualization':
          // Visual mode change (handled by main thread)
          console.log(`Worker ${workerId}: visualization mode set to ${message.mode}`);
          break;
          
        default:
          console.warn(`Worker ${workerId}: unknown message type '${message.type}'`);
      }
    } catch (error) {
      console.error(`Worker ${workerId}: error processing message '${event.data?.type}'`, error);
      postMessage({
        type: 'error',
        message: `Error in worker ${workerId}: ${error.message}`,
        workerId: workerId,
        stack: error.stack
      });
    }
  };
  
  async function handleInitialization(message) {
    const { size, moleculeCount, maxNumber, timeScale } = message;
    const rules = createCustomRules();
    rules.setConstant('time_scale', timeScale || 0.1);
    
    // Initialize simulation with the improved version
    simulation = new ImprovedChemistry(rules, size, moleculeCount, maxNumber);
    simulation.isMainWorker = !isSubWorker;
    simulation.workerId = workerId;
    
    console.log(`Worker ${workerId}: simulation initialized with ${moleculeCount} molecules, timeScale=${timeScale}`);
    
    // Initialize sub-workers
    if (!isSubWorker) {
      try {
        await simulation.initializeSubWorkers();
        console.log(`Initialized ${simulation.subWorkers.length} sub-workers`);
      } catch (error) {
        console.error("Failed to initialize sub-workers:", error);
        // Continue with single-worker mode
      }
    }
    
    // Initial step
    await simulation.step();
    
    if (!isSubWorker) {
      sendUpdate();
    }
  }
  
  function handleSubWorkerInitialization(message) {
    // Configure as sub-worker
    isSubWorker = true;
    workerId = message.workerId;
    console.log(`Started sub-worker ${workerId}`);
    
    // Create rules instance
    const rules = createCustomRules();
    rules.setConstant('time_scale', message.timeScale || 0.1);
    
    // Initialize with empty molecule set - we'll receive chunks to process
    simulation = new ImprovedChemistry(rules, message.size, 0, message.maxNumber);
    simulation.isMainWorker = false;
    simulation.workerId = workerId;
    
    // Signal readiness
    postMessage({
      type: 'sub_ready',
      workerId: workerId
    });
  }
  
  function handleProcessChunk(message) {
    if (!simulation || !isSubWorker) return;
    
    try {
      // Check pause state first
      if (message.isPaused) {
        // Send empty response if paused
        postMessage({
          type: 'chunk_processed',
          workerId: workerId,
          processedCount: 0,
          results: {
            moleculeUpdates: [],
            newMolecules: [],
            reactionCount: 0
          }
        });
        return;
      }
      
      // Get molecules to process
      const molecules = message.molecules || [];
      if (molecules.length === 0) {
        postMessage({
          type: 'chunk_processed',
          workerId: workerId,
          processedCount: 0,
          results: {
            moleculeUpdates: [],
            newMolecules: [],
            reactionCount: 0
          }
        });
        return;
      }
      
      // Import cached relationships
      if (message.cachedRelationships) {
        for (const [molId, relationships] of Object.entries(message.cachedRelationships)) {
          for (const rel of relationships) {
            simulation.moleculeCache._ensureRelationshipsList(molId);
            simulation.moleculeCache._ensureRelationshipsList(rel.otherId);
            
            // Check if relationship already exists
            if (!simulation.moleculeCache.hasRelationship(molId, rel.otherId)) {
              const relEntry = {
                otherId: rel.otherId,
                distance: rel.distance,
                lastUpdated: rel.lastUpdated || performance.now()
              };
              
              simulation.moleculeCache.moleculeRelationships.get(molId).push(relEntry);
              
              // Add reverse relationship if needed
              if (!simulation.moleculeCache.hasRelationship(rel.otherId, molId)) {
                const reverseEntry = {
                  otherId: molId,
                  distance: rel.distance,
                  lastUpdated: rel.lastUpdated || performance.now()
                };
                
                simulation.moleculeCache.moleculeRelationships.get(rel.otherId).push(reverseEntry);
              }
            }
          }
        }
      }
      
      // Set simulation parameters
      simulation.temperature = message.temperature || 1.0;
      simulation.rules.setConstant('time_scale', message.timeScale || 0.1);
      simulation.rules.setConstant('damping', message.damping || 0.99);
      
      // Create PrimeMolecule instances from serialized data
      const moleculeObjs = molecules.map(molData => 
        MoleculeSerializer.deserialize(molData, PrimeMolecule)
      );
      
      // Process chunk
      const removedIndices = new Set();
      const newMolecules = [];
      const reactionCount = simulation.reactionCount;
      
      // Update molecules using new molecule-centric caching
      simulation.updatePhysicsWithMoleculeCache(moleculeObjs, removedIndices, newMolecules);
      
      // Prepare results
      const results = {
        // Updated molecule states
        moleculeUpdates: moleculeObjs
          .filter((_, i) => !removedIndices.has(i))
          .map(mol => MoleculeSerializer.serialize(mol)),
        
        // Newly created molecules - use serialization helper
        newMolecules: newMolecules.map(mol => MoleculeSerializer.serialize(mol)),
        
        // Number of reactions that occurred
        reactionCount: simulation.reactionCount - reactionCount,
        
        // Cache updates
        updatedRelationships: simulation.moleculeCache.getRelationshipsForMolecules(
          moleculeObjs.filter((_, i) => !removedIndices.has(i))
        )
      };
      
      // Send results back
      postMessage({
        type: 'chunk_processed',
        workerId: workerId,
        processedCount: molecules.length,
        results: results
      });
      
    } catch (error) {
      console.error(`Sub-worker ${workerId}: error processing chunk`, error);
      postMessage({
        type: 'error',
        workerId: workerId,
        message: `Error processing chunk: ${error.message}`
      });
    }
  }
  
  
// Modifica alla funzione sendUpdate
function sendUpdate() {
    try {
      // Prepara dati usando posizioni stabili quando disponibili
      const moleculeData = getOptimizedMoleculeData(simulation.stablePositions);
      
      postMessage({
        type: 'update',
        molecules: moleculeData,
        temperature: simulation.temperature,
        reactionCount: simulation.reactionCount,
        isPaused: simulation.isPaused,
        cacheStats: {
          moleculeRelationships: simulation.moleculeCache.moleculeRelationships.size,
          totalRelationships: countTotalRelationships()
        },
        workerCount: simulation.subWorkers.length,
        workerStatus: simulation.subWorkers.map(worker => ({
          id: worker.id,
          busy: worker.busy,
          lastProcessed: worker.lastProcessedMolecules
        }))
      });
      
      // Ripristina posizioni fisiche reali dopo l'invio dell'aggiornamento
      if (simulation.stablePositions.size > 0) {
        simulation.restorePhysicsPositions();
      }
    } catch (error) {
      console.error(`Worker ${workerId}: error sending update`, error);
      postMessage({
        type: 'error',
        message: `Error sending update: ${error.message}`,
        stack: error.stack
      });
    }
  }
  
  // Funzione aggiornata per utilizzare posizioni stabili
  function getOptimizedMoleculeData(stablePositions) {
    const currentMoleculeIds = new Set();
    const result = [];
    
    for (const mol of simulation.molecules) {
      const id = mol.id || `mol-${Math.random().toString(36).substring(2, 11)}`;
      mol.id = id;
      currentMoleculeIds.add(id);
      
      // Serializza utilizzando il serializzatore standard
      const serialized = MoleculeSerializer.serialize(mol);
      
      // Se disponibile, usa la posizione stabile
      if (stablePositions && stablePositions.has(id)) {
        serialized.position = [...stablePositions.get(id)];
      }
      
      result.push(serialized);
    }
    
    // Calcola molecole rimosse
    const removedIds = [...previousMoleculeIds].filter(id => !currentMoleculeIds.has(id));
    
    // Aggiorna set per il prossimo frame
    previousMoleculeIds = currentMoleculeIds;
    
    return {
      molecules: result,
      removedIds: removedIds
    };
  }
  
  // Funzioni di supporto
  
  function propagateTemperature(value) {
    // Propaga ai sub-worker
    if (simulation.isMainWorker && simulation.subWorkers.length > 0) {
      for (const subWorker of simulation.subWorkers) {
        subWorker.worker.postMessage({
          type: 'set_temperature',
          value: value
        });
      }
    }
  }
  
  function propagateTimeScale(value) {
    // Propaga ai sub-worker
    if (simulation.isMainWorker && simulation.subWorkers.length > 0) {
      for (const subWorker of simulation.subWorkers) {
        subWorker.worker.postMessage({
          type: 'set_timescale',
          value: value
        });
      }
    }
  }
  
  function countTotalRelationships() {
    if (!simulation || !simulation.moleculeCache) return 0;
    
    let count = 0;
    for (const relationships of simulation.moleculeCache.moleculeRelationships.values()) {
      count += relationships.length;
    }
    return count;
  }
  
  function cleanupResources() {
    if (simulation) {
      simulation.cleanup();
    }
    
    simulation = null;
    previousMoleculeIds.clear();
    
    postMessage({ 
      type: 'cleanup_complete',
      workerId: workerId
    });
  }