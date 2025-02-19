import {
    createCustomRules
} from './rules.js';
import {
    PrimeMolecule
} from './molecule.js';
import {
    MoleculeCentricCache,
    InteractionCache
} from './molecule-cache.js';

/**
 * Classe MoleculeSerializer per gestire serializzazione/deserializzazione di molecole
 * con supporto per il nuovo calcolo della massa basato su addizione
 */
class MoleculeSerializer {
    /**
     * Crea un oggetto serializzabile da una PrimeMolecule
     * Include tutte le proprietà necessarie per visualizzazione e UI
     */
    static serialize(molecule) {
        return {
            // Proprietà di base
            id: molecule.id,
            number: molecule.number,
            position: [...molecule.position],
            velocity: [...molecule.velocity],

            // Proprietà calcolate necessarie per visualizzazione
            prime_factors: {
                ...molecule.prime_factors
            },
            mass: molecule.mass,
            charge: molecule.charge,
            color: [...molecule.color],

            // Proprietà opzionali
            ...(molecule.angularVelocity ? {
                angularVelocity: [...molecule.angularVelocity]
            } : {}),
            ...(molecule.lastReactionTime ? {
                lastReactionTime: molecule.lastReactionTime
            } : {}),

            // Proprietà di relazione (parentela tra molecole)
            ...(molecule.parentIds && molecule.parentIds.length > 0 ? {
                parentIds: [...molecule.parentIds]
            } : {}),
            ...(molecule.reactionType ? {
                reactionType: molecule.reactionType
            } : {}),

            // Tracciamento posizione fisica
            ...(molecule._physicsPosition ? {
                _physicsPosition: [...molecule._physicsPosition]
            } : {}),

            // Proprietà campo quantico se presenti
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
     * Crea una PrimeMolecule da dati serializzati
     * Imposta solo proprietà sicuramente modificabili
     */
    static deserialize(molData, PrimeMolecule) {
        // Crea molecola con proprietà fondamentali
        const mol = new PrimeMolecule(molData.number, molData.position);

        // Imposta solo proprietà che sappiamo essere sicure da modificare
        if (molData.id) {
            mol.id = molData.id;
        }

        // Gestisci velocità con controlli di sicurezza
        if (molData.velocity && Array.isArray(molData.velocity)) {
            mol.velocity = [...molData.velocity];
        }

        // Imposta altre proprietà opzionali
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

        // Imposta tracciamento posizione fisica
        if (molData._physicsPosition && Array.isArray(molData._physicsPosition)) {
            mol._physicsPosition = [...molData._physicsPosition];
        }

        // Imposta proprietà campo quantico se presenti
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
     * Ottiene proprietà sicure da aggiornare su una molecola esistente
     * Usato quando si uniscono risultati da sub-worker
     */
    static getUpdatableProperties(molecule) {
        return {
            id: molecule.id,
            position: [...molecule.position],
            velocity: [...molecule.velocity],
            ...(molecule.angularVelocity ? {
                angularVelocity: [...molecule.angularVelocity]
            } : {}),
            ...(molecule.lastReactionTime ? {
                lastReactionTime: molecule.lastReactionTime
            } : {}),
            ...(molecule.parentIds ? {
                parentIds: [...molecule.parentIds]
            } : {}),
            ...(molecule.reactionType ? {
                reactionType: molecule.reactionType
            } : {}),
            ...(molecule._physicsPosition ? {
                _physicsPosition: [...molecule._physicsPosition]
            } : {}),
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
     * Aggiorna una molecola esistente con nuove proprietà
     * Aggiorna solo proprietà che possono essere modificate in sicurezza
     */
    static updateMoleculeProperties(molecule, updateData) {
        // Aggiorna posizione se fornita
        if (updateData.position && Array.isArray(updateData.position)) {
            molecule.position = [...updateData.position];
        }

        // Aggiorna velocità se fornita
        if (updateData.velocity && Array.isArray(updateData.velocity)) {
            molecule.velocity = [...updateData.velocity];
        }

        // Aggiorna velocità angolare se fornita
        if (updateData.angularVelocity && Array.isArray(updateData.angularVelocity)) {
            molecule.angularVelocity = [...updateData.angularVelocity];
        }

        // Aggiorna tempo ultima reazione se fornito
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

        // Aggiorna posizione fisica se fornita
        if (updateData._physicsPosition && Array.isArray(updateData._physicsPosition)) {
            molecule._physicsPosition = [...updateData._physicsPosition];
        }

        // Aggiorna proprietà campo quantico
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


// Configurazione per setup multi-worker
const NUM_SUB_WORKERS = navigator.hardwareConcurrency || 4; // Usa core disponibili

/**
 * Classe EnhancedChemistry - implementa simulazione di chimica con nuove regole fisiche
 */
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
        this.isPaused = false;

        // Inizializza sub-worker
        this.subWorkers = [];
        this.isMainWorker = true;
        this.workerId = 'main';
        this.workerBusy = false;

        // Usa il nuovo sistema di caching centrato sulle molecole
        this.moleculeCache = new MoleculeCentricCache();

        // Supporto per posizioni stabili
        this.stablePositions = new Map();
        this.positionAccumulators = new Map();
        this.positionSampleCount = new Map();

        // Flag per aggiornamenti
        this.updateInProgress = false;
        this.pendingPositionUpdate = false;

        // Inizializza molecole con distribuzione interessante
        this.initializeMolecules(moleculeCount);
    }

    /**
     * Inizializza sub-worker per elaborazione parallela
     */
    async initializeSubWorkers() {
        if (!this.isMainWorker) return; // Solo worker principale crea sub-worker

        const workerUrl = self.location.href;
        console.log(`Inizializzazione di ${NUM_SUB_WORKERS} sub-worker da ${workerUrl}`);

        for (let i = 0; i < NUM_SUB_WORKERS; i++) {
            try {
                const worker = new Worker(workerUrl, {
                    type: 'module'
                });
                const workerIndex = i + 1;

                worker.onmessage = (event) => this.handleSubWorkerMessage(workerIndex, event);

                // Inizializza sub-worker con configurazione parziale
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

                console.log(`Inizializzato sub-worker ${workerIndex}`);
            } catch (error) {
                console.error(`Errore inizializzazione sub-worker ${i + 1}:`, error);
            }
        }
    }

    /**
     * Gestisce messaggi dai sub-worker
     */
    handleSubWorkerMessage(workerIndex, event) {
        const subWorker = this.subWorkers[workerIndex - 1];

        switch (event.data.type) {
            case 'sub_ready':
                subWorker.busy = false;
                console.log(`Sub-worker ${workerIndex} pronto`);
                break;

            case 'chunk_processed':
                // Unisce risultati dal sub-worker
                this.mergeProcessedChunk(event.data.results);
                subWorker.busy = false;
                subWorker.lastProcessedMolecules = event.data.processedCount;
                break;

            case 'reaction_occurred':
                // Gestisce risultato reazione dal sub-worker
                this.handleRemoteReaction(event.data.reaction);
                break;

            case 'error':
                console.error(`Errore in sub-worker ${workerIndex}:`, event.data.message);
                subWorker.busy = false;
                break;
        }
    }

    /**
     * Unisce risultati elaborati dal sub-worker
     */
    mergeProcessedChunk(results) {
        if (!results) return;

        // Aggiorna stati molecole
        for (const update of results.moleculeUpdates) {
            const molecule = this.molecules.find(m => m.id === update.id);
            if (molecule) {
                // Usa helper serializzazione per aggiornare proprietà in sicurezza
                MoleculeSerializer.updateMoleculeProperties(molecule, update);
            }
        }

        // Aggiungi eventuali nuove molecole da reazioni
        if (results.newMolecules && results.newMolecules.length > 0) {
            for (const molData of results.newMolecules) {
                // Crea nuova molecola usando helper serializzazione
                const newMol = MoleculeSerializer.deserialize(molData, PrimeMolecule);
                newMol.id = `main-${this.nextMoleculeId++}`;
                this.molecules.push(newMol);
            }
        }

        // Aggiorna conteggio reazioni
        if (results.reactionCount) {
            this.reactionCount += results.reactionCount;
        }

        // Integra relazioni aggiornate
        if (results.updatedRelationships) {
            this.importRelationships(results.updatedRelationships);
        }
    }

    /**
     * Importa relazioni aggiornate dai sub-worker
     */
    importRelationships(updatedRelationships) {
        for (const [molId, relationships] of Object.entries(updatedRelationships)) {
            for (const rel of relationships) {
                const mol = this.molecules.find(m => m.id === molId);
                const otherMol = this.molecules.find(m => m.id === rel.otherId);

                if (mol && otherMol) {
                    this.moleculeCache.createRelationship(
                        mol,
                        otherMol,
                        rel.distance,
                        rel.lastUpdated || performance.now()
                    );
                }
            }
        }
    }

    /**
     * Gestisce una reazione avvenuta in un sub-worker
     */
    handleRemoteReaction(reaction) {
        if (!reaction) return;

        // Trova e rimuovi molecole reagite
        const mol1Index = this.molecules.findIndex(m => m.id === reaction.reactant1Id);
        const mol2Index = this.molecules.findIndex(m => m.id === reaction.reactant2Id);

        if (mol1Index >= 0 && mol2Index >= 0) {
            // Rimuovi reagenti
            const removed = [
                this.molecules[mol1Index],
                this.molecules[mol2Index]
            ];

            this.molecules = this.molecules.filter((_, i) =>
                i !== mol1Index && i !== mol2Index);

            // Aggiungi molecole prodotto
            for (const productData of reaction.products) {
                // Crea nuova molecola usando helper serializzazione
                const newMol = MoleculeSerializer.deserialize(productData, PrimeMolecule);
                newMol.id = `main-${this.nextMoleculeId++}`;
                this.molecules.push(newMol);
            }

            // Aggiorna conteggio reazioni
            this.reactionCount++;

            // Rimuovi relazioni per molecole rimosse
            for (const mol of removed) {
                this.moleculeCache.removeAllRelationshipsForMolecule(mol.id);
            }
        }
    }

    /**
     * Inizializza molecole con distribuzione interessante di numeri
     */
    initializeMolecules(count) {
        // Calcola numeri primi fino a this.maxNumber con Crivello di Eratostene
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

            return Array.from({
                    length: max + 1
                }, (_, i) => i)
                .filter(num => sieve[num]);
        };

        // Genera numeri primi e composti in base a this.maxNumber
        const primes = calculatePrimes(this.maxNumber);

        // Genera composti (numeri non primi) fino a this.maxNumber
        const compounds = Array.from({
                length: this.maxNumber - 1
            },
            (_, i) => i + 2
        ).filter(num => !primes.includes(num));

        // Distribuzione di numeri per molecole
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

        // Crea molecole
        for (let i = 0; i < count; i++) {
            // Posizione casuale entro limiti
            const pos = [
                (Math.random() - 0.5) * this.size,
                (Math.random() - 0.5) * this.size,
                (Math.random() - 0.5) * this.size
            ];

            // Scegli numero dalla distribuzione
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
    }

    /**
     * Esegue un passo della simulazione
     */
    async step() {
        // Se in pausa, non aggiornare fisica
        if (this.isPaused) {
            if (!this.isMainWorker) {
                sendUpdate();
            }
            return;
        }

        // Imposta flag aggiornamento in corso
        this.updateInProgress = true;

        const timeScale = this.rules.getConstant('time_scale');
        this.accumulatedTime += timeScale;

        // Aggiorna temperatura
        this.temperature = 1.0 + 0.2 * Math.sin(performance.now() / 5000);

        // Array per elaborazione efficiente
        const newMolecules = [];
        const removedIndices = new Set();

        // Elabora fisica usando sub-worker se disponibili
        if (this.isMainWorker && this.subWorkers.length > 0) {
            await this.distributeWorkToSubWorkers();
        } else {
            // Fallback a elaborazione single-thread
            this.updatePhysicsWithMoleculeCache(this.molecules, removedIndices, newMolecules);
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

            // Periodicamente pulisci cache
            this.moleculeCache.cleanupStaleRelationships();
        }

        // Limita numero molecole per performance
        this.manageMoleculeCount();

        // Aggiorna campi quantici
        updateQuantumFields(this);

        // Accumula posizioni per tutte molecole per ottenere posizioni stabili
        for (const mol of this.molecules) {
            this.accumulatePosition(mol);
        }

        // Periodicamente pulisci relazioni obsolete
        if (Math.random() < 0.05) {
            this.rules.cleanupOldRelations && this.rules.cleanupOldRelations();
        }

        // Aggiornamento completato
        this.updateInProgress = false;

        // Invia aggiornamento solo quando tutte le posizioni sono state accumulate
        if (this.isMainWorker) {
            if (this.pendingPositionUpdate) {
                // Stabilizza posizioni prima di inviare
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

    /**
     * Accumula posizione di una molecola per stabilità
     */
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

    /**
     * Calcola posizioni stabili come media
     */
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

                    // Preserva posizione fisica in proprietà separata
                    mol._physicsPosition = [...mol.position];
                    // Usa posizione stabile per visualizzazione
                    mol.position = [...stablePosition];
                }
            }
        }
    }

    /**
     * Ripristina posizioni fisiche reali dopo visualizzazione
     */
    restorePhysicsPositions() {
        for (const mol of this.molecules) {
            if (mol._physicsPosition) {
                mol.position = [...mol._physicsPosition];
                delete mol._physicsPosition;
            }
        }
    }

    /**
     * Distribuisce lavoro ai sub-worker
     */
    async distributeWorkToSubWorkers() {
        if (this.subWorkers.length === 0) return;

        // Usa nuovo sistema raggruppamento basato su relazioni
        const partitions = this.groupMoleculesByClusters();

        // Attendi che tutti sub-worker siano pronti
        await this.waitForAvailableWorkers();

        // Distribuisci partizioni ai worker disponibili
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

                // Invia molecole con loro relazioni
                const serializedMolecules = this.serializeMolecules(partition);
                const relationshipsData = this.moleculeCache.getRelationshipsForMolecules(partition);

                subWorker.worker.postMessage({
                    type: 'process_chunk',
                    molecules: serializedMolecules,
                    temperature: this.temperature,
                    timeScale: this.rules.getConstant('time_scale'),
                    damping: this.rules.getConstant('damping'),
                    cachedRelationships: relationshipsData,
                    randomInteractionRate: this.rules.getConstant('random_interaction_probability'),
                    isPaused: this.isPaused
                });
            });

            processingPromises.push(promise);
        }

        // Elabora partizioni rimanenti in thread principale se necessario
        if (partitionIndex < partitions.length) {
            const remainingMolecules = [];
            for (let i = partitionIndex; i < partitions.length; i++) {
                remainingMolecules.push(...partitions[i]);
            }

            const removedIndices = new Set();
            const newMolecules = [];

            // Aggiorna fisica usando nuovo sistema caching
            this.updatePhysicsWithMoleculeCache(remainingMolecules, removedIndices, newMolecules);

            // Rimuovi molecole che hanno reagito
            if (removedIndices.size > 0) {
                this.molecules = this.molecules.filter((_, i) => !removedIndices.has(i));
            }

            // Aggiungi nuove molecole da reazioni
            for (const mol of newMolecules) {
                mol.id = `new-${this.nextMoleculeId++}`;
                this.molecules.push(mol);
            }
        }

        // Attendi che tutti worker completino
        await Promise.all(processingPromises);
    }

    /**
     * Raggruppa molecole in cluster basati su interazione
     */
    groupMoleculesByClusters() {
        if (this.molecules.length <= 1) {
            return [this.molecules];
        }

        const clusters = [];
        const visited = new Set();
        const proximityThreshold = this.size * 0.15;

        // Helper per trovare cluster
        const findCluster = (startMol) => {
            const cluster = [startMol];
            const queue = [startMol];
            visited.add(startMol.id);

            while (queue.length > 0) {
                const current = queue.shift();

                // Trova molecole in prossimità usando cache quando possibile
                const relationships = this.moleculeCache.getRelationshipsForMolecule(current);

                for (const rel of relationships) {
                    const neighborMol = this.molecules.find(m => m.id === rel.otherId);
                    if (!neighborMol || visited.has(neighborMol.id)) continue;

                    cluster.push(neighborMol);
                    queue.push(neighborMol);
                    visited.add(neighborMol.id);
                }

                // Per molecole senza relazioni nella cache, cerca per distanza
                for (const mol of this.molecules) {
                    if (visited.has(mol.id) || mol.id === current.id) continue;

                    const distance = this.calculateDistance(current, mol);
                    if (distance < proximityThreshold) {
                        cluster.push(mol);
                        queue.push(mol);
                        visited.add(mol.id);

                        // Crea relazione nella cache per uso futuro
                        this.moleculeCache.createRelationship(
                            current, mol, distance, performance.now()
                        );
                    }
                }
            }

            return cluster;
        };

        // Forma cluster naturali basati su vicinanza
        for (const mol of this.molecules) {
            if (!visited.has(mol.id)) {
                const cluster = findCluster(mol);
                clusters.push(cluster);
            }
        }

        // Bilancia dimensioni cluster
        this.balanceClusters(clusters);

        return clusters;
    }

    /**
     * Bilanciamento migliorato dei cluster
     */
    balanceClusters(clusters) {
        if (clusters.length <= 1) return;

        // Dimensione ideale cluster
        const totalMolecules = clusters.reduce((sum, c) => sum + c.length, 0);
        const idealSize = Math.max(
            10,
            Math.min(50, Math.ceil(totalMolecules / (this.subWorkers.length || 1)))
        );

        // Ordina cluster per dimensione (più grandi prima)
        clusters.sort((a, b) => b.length - a.length);

        // Dividi cluster troppo grandi
        let i = 0;
        while (i < clusters.length) {
            if (clusters[i].length > idealSize * 1.5) {
                // Trova punto divisione naturale usando prossimità
                const newClusters = this.splitCluster(clusters[i], idealSize);

                // Sostituisci cluster originale con prima parte
                clusters[i] = newClusters[0];

                // Aggiungi parti rimanenti come nuovi cluster
                // Aggiungi parti rimanenti come nuovi cluster
                clusters.push(...newClusters.slice(1));

                // Riordina cluster
                clusters.sort((a, b) => b.length - a.length);
            } else {
                i++;
            }
        }

        // Unisci cluster molto piccoli
        while (clusters.length > 1) {
            const smallest = clusters.pop();
            if (smallest.length < idealSize * 0.5) {
                // Trova cluster più vicino per unione
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
                    // Rimetti indietro se non c'è buon candidato
                    clusters.push(smallest);
                    break;
                }
            } else {
                // Rimetti indietro se non troppo piccolo
                clusters.push(smallest);
                break;
            }
        }
    }

    /**
     * Helper per dividere un cluster grande
     */
    splitCluster(cluster, targetSize) {
        const result = [];
        let remaining = [...cluster];

        while (remaining.length > targetSize) {
            // Trova molecola più lontana dal centro di massa
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

            // Usa molecola più lontana come nuovo centro
            const newCenter = remaining[farthestIndex];
            const newCluster = [newCenter];
            remaining.splice(farthestIndex, 1);

            // Ordina rimanenti per distanza dal nuovo centro
            const distanceMap = remaining.map((mol, index) => ({
                index,
                distance: this.calculateDistance(mol, newCenter)
            }));

            distanceMap.sort((a, b) => a.distance - b.distance);

            // Prendi vicini più prossimi fino a targetSize
            const neighborsToTake = Math.min(targetSize - 1, remaining.length);
            const newRemaining = [];

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

        // Aggiungi molecole rimanenti come ultimo cluster
        if (remaining.length > 0) {
            result.push(remaining);
        }

        return result;
    }

    /**
     * Calcola centro di massa di un gruppo di molecole
     */
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

    /**
     * Calcola distanza da un punto
     */
    calculateDistanceFromPoint(molecule, point) {
        let sumSquared = 0;
        for (let i = 0; i < 3; i++) {
            const diff = molecule.position[i] - point[i];
            sumSquared += diff * diff;
        }
        return Math.sqrt(sumSquared);
    }

    /**
     * Calcola distanza minima tra cluster
     */
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

    /**
     * Calcola distanza tra molecole
     */
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

    /**
     * Attende worker disponibili
     */
    waitForAvailableWorkers() {
        if (this.subWorkers.length === 0) return Promise.resolve();

        return new Promise(resolve => {
            const checkWorkers = () => {
                // Verifica se almeno un worker è disponibile
                const anyAvailable = this.subWorkers.some(w => !w.busy);

                if (anyAvailable) {
                    resolve();
                } else {
                    // Controlla di nuovo presto
                    setTimeout(checkWorkers, 5);
                }
            };

            checkWorkers();
        });
    }

    /**
     * Serializza molecole per invio
     */
    serializeMolecules(molecules) {
        return molecules.map(mol => MoleculeSerializer.serialize(mol));
    }

    /**
     * Aggiornamento fisico usando caching centrato su molecole
     */
    updatePhysicsWithMoleculeCache(moleculesToProcess, removedIndices, newMolecules) {
        const timeScale = this.rules.getConstant('time_scale');
        const damping = this.rules.getConstant('damping');
        const now = performance.now();

        // Pre-calcola forze per ogni molecola
        const forces = new Map();
        moleculesToProcess.forEach(mol => forces.set(mol.id, [0, 0, 0]));

        // Determina quali molecole necessitano calcoli relazione
        const uncachedMolecules = this.moleculeCache.getUncachedMolecules(moleculesToProcess);

        // Elabora prima molecole non in cache (stabilisce nuove relazioni)
        for (const mol1 of uncachedMolecules) {
            for (const mol2 of moleculesToProcess) {
                if (mol1.id === mol2.id) continue;

                // Verifica se abbastanza vicine per stabilire relazione
                const distance = this.calculateDistance(mol1, mol2);
                const interactionThreshold = this.calculateInteractionThreshold(mol1, mol2);

                if (distance < interactionThreshold) {
                    // Crea nuova relazione in cache
                    this.moleculeCache.createRelationship(mol1, mol2, distance, now);
                }
            }
        }

        // Considera anche interazioni casuali per garantire copertura completa
        if (this.rules.getConstant('random_interaction_probability') > 0) {
            this.processRandomInteractions(moleculesToProcess, now);
        }

        // Elabora tutte molecole con relazioni in cache
        for (let i = 0; i < moleculesToProcess.length; i++) {
            const mol1 = moleculesToProcess[i];
            if (removedIndices.has(i)) continue;

            // Applica moto termico casuale basato su temperatura
            this.applyThermalMotion(mol1);

            // Ottieni relazioni da cache per questa molecola
            const relationships = this.moleculeCache.getRelationshipsForMolecule(mol1);

            // Elabora ogni relazione
            for (const rel of relationships) {
                // Trova altra molecola nella lista elaborazione
                const mol2Index = moleculesToProcess.findIndex(m => m.id === rel.otherId);
                if (mol2Index < 0) continue; // Non nel nostro batch

                const mol2 = moleculesToProcess[mol2Index];
                if (removedIndices.has(mol2Index)) continue;

                // Aggiorna distanza relazione
                const currentDistance = this.calculateDistance(mol1, mol2);
                rel.distance = currentDistance;
                rel.lastUpdated = now;

                // Verifica se relazione deve essere mantenuta
                const maxDistance = this.calculateMaxRelationshipDistance(mol1, mol2);

                if (currentDistance > maxDistance) {
                    // Molecole troppo distanti, rimuovi relazione
                    this.moleculeCache.removeRelationship(mol1.id, mol2.id);
                    continue;
                }

                // Calcola e applica forze
                const [force1, force2] = this.calculateForces(mol1, mol2, currentDistance, now);

                // Applica forze calcolate
                const mol1Force = forces.get(mol1.id);
                const mol2Force = forces.get(mol2.id);

                if (mol1Force && mol2Force) {
                    for (let k = 0; k < 3; k++) {
                        mol1Force[k] += force1[k] * timeScale;
                        mol2Force[k] += force2[k] * timeScale;
                    }
                }

                // Verifica reazioni
                if (currentDistance < this.calculateReactionDistance(mol1, mol2)) {
                    if (this.shouldReact(mol1, mol2)) {
                        const products = this.processReaction(mol1, mol2);
                        if (products.length > 0) {
                            newMolecules.push(...products);

                            // Trova indici nell'array originale
                            const mol1Index = this.molecules.findIndex(m => m.id === mol1.id);
                            const mol2Index = this.molecules.findIndex(m => m.id === mol2.id);

                            if (mol1Index >= 0 && mol2Index >= 0) {
                                removedIndices.add(mol1Index);
                                removedIndices.add(mol2Index);
                                this.reactionCount++;

                                // Rimuovi relazioni per molecole reagenti
                                this.moleculeCache.removeAllRelationshipsForMolecule(mol1.id);
                                this.moleculeCache.removeAllRelationshipsForMolecule(mol2.id);
                                break;
                            }
                        }
                    }
                }
            }
        }

        // Applica forze e aggiorna posizioni
        for (const mol of moleculesToProcess) {
            const molForce = forces.get(mol.id);
            if (!molForce) continue;

            // Aggiorna velocità con forze
            for (let k = 0; k < 3; k++) {
                // Usa massa nel calcolo accelerazione (F = ma, quindi a = F/m)
                const acceleration = molForce[k] / mol.mass;
                mol.velocity[k] += acceleration * 0.1;
                mol.position[k] += mol.velocity[k] * timeScale;
                mol.velocity[k] *= damping;
            }

            // Verifica confini
            this.enforceBoundaries(mol);
        }
    }

    /**
     * Processa interazioni casuali per garantire copertura completa
     */
    processRandomInteractions(molecules, timestamp) {
        const randomInteractionRate = this.rules.getConstant('random_interaction_probability');
        const moleculeCount = molecules.length;

        // Limita numero interazioni casuali per frame
        const maxRandomInteractions = Math.min(20, Math.ceil(moleculeCount * randomInteractionRate));
        let processedCount = 0;

        // Elabora interazioni casuali
        for (let attempts = 0; attempts < 50 && processedCount < maxRandomInteractions; attempts++) {
            // Scegli due molecole casuali
            const idx1 = Math.floor(Math.random() * moleculeCount);
            let idx2 = Math.floor(Math.random() * moleculeCount);

            // Evita di scegliere stessa molecola
            while (idx1 === idx2 && moleculeCount > 1) {
                idx2 = Math.floor(Math.random() * moleculeCount);
            }

            if (idx1 === idx2) continue;

            const mol1 = molecules[idx1];
            const mol2 = molecules[idx2];

            // Verifica se questa interazione casuale dovrebbe essere considerata
            if (!this.moleculeCache.hasRelationship(mol1.id, mol2.id)) {
                // Calcola distanza
                const distance = this.calculateDistance(mol1, mol2);
                // Crea relazione
                this.moleculeCache.createRelationship(mol1, mol2, distance, timestamp);
                processedCount++;
            }
        }
    }

    /**
     * Applica moto termico casuale
     */
    applyThermalMotion(molecule) {
        // Effetto termico proporzionale a temperatura/massa
        const thermalFactor = this.temperature / (molecule.mass + 1);

        for (let axis = 0; axis < 3; axis++) {
            molecule.velocity[axis] += (Math.random() - 0.5) * 0.01 * thermalFactor;
        }
    }

    /**
     * Calcola soglia interazione tra molecole
     */
    calculateInteractionThreshold(mol1, mol2) {
        // Soglia dinamica basata su massa
        return Math.max(
            5.0, // Minimo
            (mol1.mass + mol2.mass) * 1.2
        );
    }

    /**
     * Calcola distanza massima per mantenere relazione
     */
    calculateMaxRelationshipDistance(mol1, mol2) {
        return Math.max(
            10.0, // Minimo
            (mol1.mass + mol2.mass) * 2.0
        );
    }

    /**
     * Calcola distanza reazione basata su massa
     */
    calculateReactionDistance(mol1, mol2) {
        return 1.0 + (mol1.mass + mol2.mass) * 0.1;
    }

    /**
     * Calcola forze tra molecole
     */
    calculateForces(mol1, mol2, distance, now) {
        // Verifica periodo raffreddamento
        if (this.rules.isInCoolingPeriod &&
            (this.rules.isInCoolingPeriod(mol1, now) ||
                this.rules.isInCoolingPeriod(mol2, now))) {
            return [
                [0, 0, 0],
                [0, 0, 0]
            ];
        }

        // Calcola vettore direzione
        const direction = [
            mol2.position[0] - mol1.position[0],
            mol2.position[1] - mol1.position[1],
            mol2.position[2] - mol1.position[2]
        ];

        if (distance < 0.001) {
            return [
                [0, 0, 0],
                [0, 0, 0]
            ];
        }

        // Normalizza direzione
        const dirNorm = direction.map(d => d / distance);

        // Inizializza forze
        const force1 = [0, 0, 0];
        const force2 = [0, 0, 0];

        // Applica regole fisiche
        for (const rule of this.rules.interaction_rules) {
            if (rule.condition(mol1.prime_factors, mol2.prime_factors)) {
                let f;
                if (rule.force_function.length === 4) {
                    // Forza con massa
                    f = rule.force_function(dirNorm, distance, mol1.mass, mol2.mass);
                } else {
                    // Forza con carica
                    f = rule.force_function(dirNorm, distance, mol1.charge, mol2.charge);
                }

                // Applica intensità regola
                for (let i = 0; i < 3; i++) {
                    const scaledForce = f[i] * rule.strength;
                    force1[i] += scaledForce;
                    force2[i] -= scaledForce;
                }
            }
        }

        // Aggiungi repulsione tra "parenti" se applicabile
        if (this.rules.areRelated && this.rules.areRelated(mol1.id, mol2.id)) {
            const repulsionFactor = this.rules.getFamilyRepulsionFactor &&
                this.rules.getFamilyRepulsionFactor(mol1.id, mol2.id, now);

            if (repulsionFactor > 0) {
                const minDistance = this.rules.getConstant('min_distance') || 0.1;
                const effectiveDistance = Math.max(distance, minDistance);
                const repulsiveForce = dirNorm.map(
                    x => -x * repulsionFactor / (effectiveDistance ** 1.5)
                );

                // Applica a forza risultante
                for (let i = 0; i < 3; i++) {
                    force1[i] += repulsiveForce[i];
                    force2[i] -= repulsiveForce[i];
                }
            }
        }

        // Limita forza massima
        const maxForce = this.rules.getConstant('max_force');
        for (let i = 0; i < 3; i++) {
            force1[i] = Math.max(-maxForce, Math.min(maxForce, force1[i]));
            force2[i] = Math.max(-maxForce, Math.min(maxForce, force2[i]));
        }

        return [force1, force2];
    }

    /**
     * Determina se molecole dovrebbero reagire
     */
    shouldReact(mol1, mol2) {
        // Verifica periodo raffreddamento
        const now = performance.now();
        if (this.rules.isInCoolingPeriod &&
            (this.rules.isInCoolingPeriod(mol1, now) ||
                this.rules.isInCoolingPeriod(mol2, now))) {
            return false;
        }

        // Calcola velocità relativa
        const relativeSpeed = Math.sqrt(
            Math.pow(mol1.velocity[0] - mol2.velocity[0], 2) +
            Math.pow(mol1.velocity[1] - mol2.velocity[1], 2) +
            Math.pow(mol1.velocity[2] - mol2.velocity[2], 2)
        );

        // Probabilità reazione basata su massa
        const baseReactionProb = 0.1 * this.temperature * relativeSpeed;
        // Molecole più pesanti reagiscono più lentamente
        const sizeFactor = 1.0 / (1.0 + Math.log(mol1.mass + mol2.mass));

        return Math.random() < baseReactionProb * sizeFactor;
    }

    /**
     * Elabora reazione tra molecole
     */
    processReaction(mol1, mol2) {
        // Cerca regole reazione applicabili
        for (const rule of this.rules.reaction_rules) {
            // Passa molecole a funzione condizione
            if (rule.condition(mol1.prime_factors, mol2.prime_factors, mol1, mol2)) {
                if (Math.random() < rule.probability * this.temperature) {
                    // Determina tipo reazione
                    const reactionType = this.determineReactionType(mol1, mol2);
                    let products = [];

                    switch (reactionType) {
                        case 'fusion':
                            // Fusione: crea molecola combinata
                            products = this.handleFusion(mol1, mol2, rule);
                            break;
                        case 'fission':
                            // Fissione: divide in molecole più piccole
                            products = this.handleFission(mol1, mol2, rule);
                            break;
                        case 'emission':
                            // Emissione: mantiene originali ma emette particella
                            products = this.handleEmission(mol1, mol2, rule);
                            break;
                        case 'standard':
                        default:
                            // Reazione standard come definita nella regola
                            products = rule.effect(mol1, mol2);
                    }

                    // Imposta tempo reazione per effetto visivo
                    const now = performance.now();
                    products.forEach(p => {
                        p.setReactionTime(now);

                        // Imposta relazioni parentela se necessario
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

    /**
     * Determina tipo reazione basato su proprietà molecole
     */
    determineReactionType(mol1, mol2) {
        // Determina tipo reazione in base a proprietà molecole
        const totalMass = mol1.mass + mol2.mass;
        const relativeSpeed = Math.sqrt(
            Math.pow(mol1.velocity[0] - mol2.velocity[0], 2) +
            Math.pow(mol1.velocity[1] - mol2.velocity[1], 2) +
            Math.pow(mol1.velocity[2] - mol2.velocity[2], 2)
        );

        // Determina tipo reazione in base a massa, carica e velocità
        if (totalMass > 30 && relativeSpeed > 0.5) {
            return 'fission'; // Scissione per molecole grandi ad alta energia
        } else if (mol1.charge * mol2.charge < 0 && Math.abs(mol1.charge) + Math.abs(mol2.charge) > 3) {
            return 'emission'; // Emissione quando cariche opposte forti interagiscono
        } else if (mol1.mass <= 20 && mol2.mass <= 20 && this.temperature > 1.2) {
            return 'fusion'; // Fusione più probabile per masse piccole ad alta temperatura
        } else {
            return 'standard'; // Reazione standard in altri casi
        }
    }

    /**
     * Gestisce fusione di molecole
     */
    handleFusion(mol1, mol2, rule) {
        // Combina molecole in una più grande
        const midpoint = [
            (mol1.position[0] + mol2.position[0]) / 2,
            (mol1.position[1] + mol2.position[1]) / 2,
            (mol1.position[2] + mol2.position[2]) / 2
        ];

        // Combina velocità proporzionalmente alla massa
        const totalMass = mol1.mass + mol2.mass;
        const combinedVelocity = [
            (mol1.velocity[0] * mol1.mass + mol2.velocity[0] * mol2.mass) / totalMass,
            (mol1.velocity[1] * mol1.mass + mol2.velocity[1] * mol2.mass) / totalMass,
            (mol1.velocity[2] * mol1.mass + mol2.velocity[2] * mol2.mass) / totalMass
        ];

        // Usa regola per determinare numero nuova molecola
        let products;
        if (rule.effect) {
            products = rule.effect(mol1, mol2);
        } else {
            // Fallback: crea molecola con prodotto numeri
            const fusedNumber = mol1.number * mol2.number;
            const fusedMolecule = new PrimeMolecule(
                Math.min(fusedNumber, this.maxNumber),
                midpoint
            );
            fusedMolecule.velocity = combinedVelocity;
            fusedMolecule.parentIds = [mol1.id, mol2.id];
            fusedMolecule.reactionType = 'fusion';
            products = [fusedMolecule];
        }

        // Assicura che relazioni parentela siano impostate
        products.forEach(p => {
            if (!p.parentIds || p.parentIds.length === 0) {
                p.parentIds = [mol1.id, mol2.id];
            }
            if (!p.reactionType) {
                p.reactionType = 'fusion';
            }
        });

        return products;
    }

    /**
     * Gestisce fissione di molecole
     */
    handleFission(mol1, mol2, rule) {
        // Divide molecole in frammenti più piccoli
        const midpoint = [
            (mol1.position[0] + mol2.position[0]) / 2,
            (mol1.position[1] + mol2.position[1]) / 2,
            (mol1.position[2] + mol2.position[2]) / 2
        ];

        // Se abbiamo regola di effetto, usala
        if (rule.effect) {
            const baseProducts = rule.effect(mol1, mol2);
            if (baseProducts.length > 0) {
                // Assicura che relazioni parentela siano impostate
                baseProducts.forEach(p => {
                    if (!p.parentIds || p.parentIds.length === 0) {
                        p.parentIds = [mol1.id, mol2.id];
                    }
                    if (!p.reactionType) {
                        p.reactionType = 'fission';
                    }
                });
                return baseProducts;
            }
        }

        // Altrimenti crea fissione basata su fattori primi
        const products = [];
        const largerMol = mol1.mass > mol2.mass ? mol1 : mol2;
        const factorsLarger = Object.entries(largerMol.prime_factors);

        // Crea frammenti dai fattori primi della molecola più grande
        for (const [prime, exponent] of factorsLarger) {
            if (exponent > 0) {
                const fragmentNumber = parseInt(prime);
                const fragmentMol = new PrimeMolecule(fragmentNumber, [...midpoint]);

                // Aggiungi velocità casuale in direzione opposta al centro
                const direction = [
                    Math.random() - 0.5,
                    Math.random() - 0.5,
                    Math.random() - 0.5
                ];
                const norm = Math.sqrt(direction[0] ** 2 + direction[1] ** 2 + direction[2] ** 2);
                if (norm > 0.001) {
                    fragmentMol.velocity = [
                        direction[0] / norm * 0.5 * this.temperature,
                        direction[1] / norm * 0.5 * this.temperature,
                        direction[2] / norm * 0.5 * this.temperature
                    ];
                } else {
                    fragmentMol.velocity = [
                        (Math.random() - 0.5) * 0.5 * this.temperature,
                        (Math.random() - 0.5) * 0.5 * this.temperature,
                        (Math.random() - 0.5) * 0.5 * this.temperature
                    ];
                }

                fragmentMol.parentIds = [largerMol.id];
                fragmentMol.reactionType = 'fission';
                products.push(fragmentMol);
            }
        }

        // Limita numero prodotti se necessario
        if (products.length > 5) {
            return products.slice(0, 5);
        }

        if (products.length > 0) {
            return products;
        } else {
            // Usa regola originale se non abbiamo potuto creare frammenti
            const defaultProducts = rule.effect(mol1, mol2);
            defaultProducts.forEach(p => {
                if (!p.parentIds || p.parentIds.length === 0) {
                    p.parentIds = [mol1.id, mol2.id];
                }
                if (!p.reactionType) {
                    p.reactionType = 'fission';
                }
            });
            return defaultProducts;
        }
    }

    /**
     * Gestisce emissione di particelle
     */
    handleEmission(mol1, mol2, rule) {
        // Emette particella ma mantiene molecole originali modificate
        const products = [];

        // Posizione emissione (media posizioni)
        const emissionPoint = [
            (mol1.position[0] + mol2.position[0]) / 2,
            (mol1.position[1] + mol2.position[1]) / 2,
            (mol1.position[2] + mol2.position[2]) / 2
        ];

        // Se regola fornisce effetto, usalo come base
        if (rule.effect) {
            const baseProducts = rule.effect(mol1, mol2);
            if (baseProducts.length > 0) {
                baseProducts.forEach(p => {
                    if (!p.parentIds || p.parentIds.length === 0) {
                        p.parentIds = [mol1.id, mol2.id];
                    }
                    if (!p.reactionType) {
                        p.reactionType = 'emission';
                    }
                });
                products.push(...baseProducts);
            }
        }

        // Se non abbiamo prodotti da regola, creiamo versioni modificate molecole originali
        if (products.length === 0) {
            // Scegli molecola più piccola per modificarla
            const smallerMol = mol1.mass <= mol2.mass ? mol1 : mol2;
            const largerMol = mol1.mass > mol2.mass ? mol1 : mol2;

            // Crea copie modificate molecole originali
            const modifiedSmaller = new PrimeMolecule(
                Math.max(2, smallerMol.number - 1),
                [...smallerMol.position]
            );
            modifiedSmaller.velocity = [...smallerMol.velocity];
            modifiedSmaller.parentIds = [smallerMol.id];
            modifiedSmaller.reactionType = 'emission';

            const modifiedLarger = new PrimeMolecule(
                Math.max(2, largerMol.number - 1),
                [...largerMol.position]
            );
            modifiedLarger.velocity = [...largerMol.velocity];
            modifiedLarger.parentIds = [largerMol.id];
            modifiedLarger.reactionType = 'emission';

            products.push(modifiedSmaller, modifiedLarger);
        }

        // Aggiungi particella emessa (numero primo piccolo)
        const emittedParticle = new PrimeMolecule(
            this.getSmallPrimeNumber(),
            [...emissionPoint]
        );

        // Imposta velocità emissione in direzione casuale
        const emissionDirection = [
            Math.random() - 0.5,
            Math.random() - 0.5,
            Math.random() - 0.5
        ];
        const norm = Math.sqrt(
            emissionDirection[0] ** 2 +
            emissionDirection[1] ** 2 +
            emissionDirection[2] ** 2
        );

        const emissionSpeed = 0.8 * this.temperature;
        if (norm > 0.001) {
            emittedParticle.velocity = [
                emissionDirection[0] / norm * emissionSpeed,
                emissionDirection[1] / norm * emissionSpeed,
                emissionDirection[2] / norm * emissionSpeed
            ];
        } else {
            emittedParticle.velocity = [
                (Math.random() - 0.5) * emissionSpeed,
                (Math.random() - 0.5) * emissionSpeed,
                (Math.random() - 0.5) * emissionSpeed
            ];
        }

        // Imposta proprietà speciali per particella emessa
        emittedParticle._isQuantumField = true;
        emittedParticle._creationTime = performance.now();
        emittedParticle._lifetime = 3000 + Math.random() * 2000; // 3-5 secondi
        emittedParticle.parentIds = [mol1.id, mol2.id];
        emittedParticle.reactionType = 'emission';

        products.push(emittedParticle);

        return products;
    }

    /**
     * Restituisce numero primo casuale piccolo per particelle emesse
     */
    getSmallPrimeNumber() {
        const smallPrimes = [2, 3, 5, 7, 11, 13];
        return smallPrimes[Math.floor(Math.random() * smallPrimes.length)];
    }

    /**
     * Crea campo quantico temporaneo
     */
    createQuantumField(mol1, mol2, products) {
        // Crea campo quantico temporaneo che influenza molecole circostanti
        const reactionCenter = [
            (mol1.position[0] + mol2.position[0]) / 2,
            (mol1.position[1] + mol2.position[1]) / 2,
            (mol1.position[2] + mol2.position[2]) / 2
        ];

        // Crea molecola speciale che rappresenta campo quantico
        const field = new PrimeMolecule(1, [...reactionCenter]);
        field._isQuantumField = true;
        field._creationTime = performance.now();
        field._lifetime = 1500 + Math.random() * 1000; // 1.5-2.5 secondi vita
        field.velocity = [0, 0, 0]; // Campo stazionario

        // Imposta colore e proprietà speciali
        field.color = [0.2, 0.8, 1];
        field.parentIds = [mol1.id, mol2.id];
        field.reactionType = 'quantum_field';

        // Aggiungi campo a lista molecole
        field.id = `field-${this.nextMoleculeId++}`;
        this.molecules.push(field);
    }

    /**
     * Mantiene molecole entro confini della simulazione
     */
    enforceBoundaries(molecule) {
        const margin = 0.1;
        const bounceElasticity = 0.8;

        for (let axis = 0; axis < 3; axis++) {
            const halfSize = this.size / 2 - margin;

            if (molecule.position[axis] > halfSize) {
                molecule.position[axis] = halfSize;
                molecule.velocity[axis] *= -bounceElasticity;
            } else if (molecule.position[axis] < -halfSize) {
                molecule.position[axis] = -halfSize;
                molecule.velocity[axis] *= -bounceElasticity;
            }
        }
    }

    /**
     * Gestisce numero totale molecole per performance
     */
    manageMoleculeCount() {
        // Limita numero molecole per performance
        if (this.molecules.length > 300) {
            // Rimuovi molecole in eccesso, preferendo quelle più vecchie
            // Ordina per ID (assumendo che ID più bassi siano più vecchi)
            this.molecules.sort((a, b) => {
                const idA = parseInt(a.id.split('-')[1]) || 0;
                const idB = parseInt(b.id.split('-')[1]) || 0;
                return idA - idB;
            });

            // Rimuovi più vecchie, ma conserva energie più alte
            const toKeep = this.molecules.slice(-250);
            const potentialToRemove = this.molecules.slice(0, -250);

            // Mantieni anche alcune molecole ad alta energia
            const highEnergyThreshold = 0.8; // Valore arbitrario per energia alta
            const highEnergyMolecules = potentialToRemove.filter(mol => {
                const speed = Math.sqrt(
                    mol.velocity[0] ** 2 + mol.velocity[1] ** 2 + mol.velocity[2] ** 2
                );
                const kineticEnergy = 0.5 * mol.mass * speed ** 2;
                return kineticEnergy > highEnergyThreshold;
            });

            // Limita a 50 molecole ad alta energia da mantenere
            const highEnergyToKeep = highEnergyMolecules.slice(-50);

            // Combina molecole da mantenere
            this.molecules = [...toKeep, ...highEnergyToKeep];

            // Pulizia relazioni per molecole rimosse
            this.moleculeCache.cleanupStaleRelationships();
        } else if (this.molecules.length < 50) {
            // Aggiungi molecole se troppo poche
            this.addRandomMolecules(50 - this.molecules.length);
        }
    }

    /**
     * Aggiunge molecole casuali al sistema
     */
    addRandomMolecules(count) {
        for (let i = 0; i < count; i++) {
            // Posizione casuale vicino bordi
            const side = Math.floor(Math.random() * 6); // 6 facce del cubo
            const pos = [0, 0, 0];

            for (let j = 0; j < 3; j++) {
                if (Math.floor(side / 2) === j) {
                    // Su questa faccia
                    pos[j] = (side % 2 === 0 ? -1 : 1) * (this.size / 2 - 0.2);
                } else {
                    // Posizione casuale su altri assi
                    pos[j] = (Math.random() - 0.5) * this.size;
                }
            }

            // Crea molecola semplice
            const number = Math.floor(Math.random() * 48) + 2;
            const mol = new PrimeMolecule(number, pos);
            mol.id = `spawn-${this.nextMoleculeId++}`;

            // Velocità verso centro
            const dirToCenter = pos.map(p => -p);
            const norm = Math.sqrt(dirToCenter.reduce((s, v) => s + v * v, 0));
            if (norm > 0.001) {
                mol.velocity = dirToCenter.map(v => (v / norm) * (0.1 + Math.random() * 0.2));
            } else {
                mol.velocity = [
                    (Math.random() - 0.5) * 0.3,
                    (Math.random() - 0.5) * 0.3,
                    (Math.random() - 0.5) * 0.3
                ];
            }

            this.molecules.push(mol);
        }
    }

    /**
     * Imposta temperatura sistema
     */
    setTemperature(value) {
        this.temperature = value;
    }

    /**
     * Imposta fattore di scala tempo
     */
    setTimeScale(value) {
        if (this.rules) {
            this.rules.setConstant('time_scale', value);

            // Propaga a sub-worker
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

    /**
     * Imposta stato pausa
     */
    setPauseState(isPaused) {
        this.isPaused = isPaused;

        // Propaga stato pausa a sub-worker
        if (this.isMainWorker && this.subWorkers.length > 0) {
            for (const subWorker of this.subWorkers) {
                subWorker.worker.postMessage({
                    type: 'set_pause',
                    isPaused: isPaused
                });
            }
        }
    }

    /**
     * Ottieni conteggio relazioni
     */
    getRelationshipCount() {
        // Ottieni statistiche relazioni
        if (this.moleculeCache) {
            return this.moleculeCache.getStats();
        }
        return {
            totalRelationships: 0
        };
    }

    /**
     * Pulisci risorse
     */
    cleanup() {
        // Pulisci risorse
        if (this.isMainWorker && this.subWorkers.length > 0) {
            // Termina tutti sub-worker
            for (const subWorker of this.subWorkers) {
                subWorker.worker.postMessage({
                    type: 'cleanup'
                });
                subWorker.worker.terminate();
            }
            this.subWorkers = [];
        }

        this.molecules = [];
        this.moleculeCache.clear();
        this.stablePositions.clear();
        this.positionAccumulators.clear();
        this.positionSampleCount.clear();
        this.isPaused = false;
        this.updateInProgress = false;
        this.pendingPositionUpdate = false;
    }
}

/**
 * Funzione per aggiornare campi quantici
 * Supporta effetti specifici per fissione, scissione ed emissione
 */
function updateQuantumFields(simulation) {
    const now = performance.now();

    // Filtra molecole campo quantico
    const quantumFields = simulation.molecules.filter(mol => mol._isQuantumField);
    const fieldsToRemove = [];

    for (const field of quantumFields) {
        // Verifica se campo è scaduto
        if (now - field._creationTime > field._lifetime) {
            fieldsToRemove.push(field);
            continue;
        }

        // Determina tipo e intensità campo
        const fieldType = field._fieldType || 'standard';
        const fieldStrength = field._fieldStrength || 1.0;

        // Calcola effetto campo su altre molecole
        for (const mol of simulation.molecules) {
            // Salta campo stesso e altri campi
            if (mol._isQuantumField || mol === field) continue;

            // Calcola distanza da campo
            const distance = simulation.calculateDistance(field, mol);
            const fieldRadius = fieldType === 'fission' ? 8.0 : 5.0;

            // Campo influenza solo molecole entro certo raggio
            if (distance < fieldRadius) {
                // Calcola direzione forza (da campo a molecola)
                const direction = [
                    mol.position[0] - field.position[0],
                    mol.position[1] - field.position[1],
                    mol.position[2] - field.position[2]
                ];

                // Normalizza direzione
                const dirNorm = Math.sqrt(
                    direction[0] ** 2 +
                    direction[1] ** 2 +
                    direction[2] ** 2
                );

                if (dirNorm > 0.001) {
                    const normalizedDir = direction.map(d => d / dirNorm);

                    // Effetti differenziati in base a tipo campo
                    switch (fieldType) {
                        case 'fission':
                            // Campo fissione: forza repulsiva pulsante
                            const fissionPulse = Math.sin((now - field._creationTime) / 100) * 0.5 + 0.5;
                            const fissionStrength = 0.05 * fieldStrength * fissionPulse * (1 - distance / fieldRadius);

                            // Applica forza repulsiva inversamente proporzionale alla massa
                            const massEffect = 1 / Math.sqrt(mol.mass);
                            for (let i = 0; i < 3; i++) {
                                mol.velocity[i] += normalizedDir[i] * fissionStrength * massEffect;
                            }

                            // Occasionalmente causa rotazione
                            if (Math.random() < 0.05 && mol.angularVelocity) {
                                mol.angularVelocity = mol.angularVelocity.map(v => v + (Math.random() - 0.5) * 0.1);
                            }
                            break;

                        case 'emission':
                            // Campo emissione: oscillazione perpendicolare
                            const emissionFreq = (now - field._creationTime) / 150;
                            const emissionWave = Math.sin(emissionFreq) * Math.cos(emissionFreq * 0.7);
                            const emissionStrength = 0.03 * fieldStrength * (1 - distance / fieldRadius);

                            // Crea vettore perpendicolare alla direzione
                            const perpVector = [
                                normalizedDir[1] - normalizedDir[2],
                                normalizedDir[2] - normalizedDir[0],
                                normalizedDir[0] - normalizedDir[1]
                            ];

                            // Normalizza vettore perpendicolare
                            const perpNorm = Math.sqrt(
                                perpVector[0] ** 2 +
                                perpVector[1] ** 2 +
                                perpVector[2] ** 2
                            );

                            if (perpNorm > 0.001) {
                                // Applica forza oscillante perpendicolare, scala in base a massa
                                const perpMassEffect = 1 / (0.5 + Math.sqrt(mol.mass));
                                for (let i = 0; i < 3; i++) {
                                    mol.velocity[i] += (perpVector[i] / perpNorm) * emissionWave * emissionStrength * perpMassEffect;
                                }
                            }
                            break;

                        case 'standard':
                        default:
                            // Campo standard: oscillazione radiale
                            const oscillation = Math.sin((now - field._creationTime) / 200);
                            const fieldEffect = 0.02 * fieldStrength * (1 - distance / fieldRadius) * oscillation;

                            // Applica forza oscillatoria radiale con effetto massa
                            const standardMassEffect = 1 / (1 + Math.log(mol.mass));
                            for (let i = 0; i < 3; i++) {
                                mol.velocity[i] += normalizedDir[i] * fieldEffect * standardMassEffect;
                            }
                    }

                    // Effetto visivo: molecole nei campi quantici possono cambiare colore temporaneamente
                    if (Math.random() < 0.02 && distance < fieldRadius * 0.5) {
                        // Leggero shift colore verso colore campo
                        for (let i = 0; i < 3; i++) {
                            mol.color[i] = mol.color[i] * 0.95 + field.color[i] * 0.05;
                        }
                    }
                }
            }
        }

        // Aggiorna campo stesso
        if (fieldType === 'fission') {
            // Campi fissione si espandono gradualmente
            const expansionFactor = 1.0 + (now - field._creationTime) / field._lifetime * 0.5;
            field._fieldStrength = fieldStrength * (1.0 - (now - field._creationTime) / field._lifetime);
        } else if (fieldType === 'emission') {
            // Campi emissione pulsano
            field._fieldStrength = fieldStrength * (0.5 + 0.5 * Math.sin((now - field._creationTime) / 300));
        }
    }

    // Rimuovi campi scaduti
    if (fieldsToRemove.length > 0) {
        simulation.molecules = simulation.molecules.filter(
            mol => !fieldsToRemove.includes(mol));

        // Se campo scompare, a volte può lasciare particella residua
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

                if (field.parentIds && field.parentIds.length > 0) {
                    residue.parentIds = [...field.parentIds];
                }
                residue.reactionType = 'quantum_residue';

                simulation.molecules.push(residue);
            }
        }
    }
}

// Extension to PrimeMolecule prototype
PrimeMolecule.prototype.isRelatedTo = function(otherMolecule, rules) {
    if (!this.id || !otherMolecule.id) return false;
    return rules.areRelated(this.id, otherMolecule.id);
};

///
/// Improved chemistry
///

let simulation;
let previousMoleculeIds = new Set();
let isSubWorker = false;
let workerId = 'main';

/**
 * MoleculeInteractionManager - gestisce efficientemente le interazioni tra molecole
 * Implementa un sistema di calcolo non ridondante e ottimizzato per la parallelizzazione
 */
class MoleculeInteractionManager {
    constructor() {
        // Mappa delle interazioni già calcolate (evita ridondanza)
        this.calculatedPairs = new Set();
        
        // Cache per interazioni - utilizza ID unici per le coppie
        this.interactionCache = new Map();
        
        // Statistiche sulle prestazioni
        this.stats = {
            calculationsAvoided: 0,
            totalCalculations: 0,
            cacheHits: 0
        };
    }

    /**
     * Genera un ID univoco per ogni coppia di molecole
     * Garantisce che (A,B) e (B,A) abbiano lo stesso ID
     */
    getPairId(mol1Id, mol2Id) {
        // Ordina gli ID per garantire che (A,B) e (B,A) producano lo stesso ID
        const [smaller, larger] = mol1Id < mol2Id 
            ? [mol1Id, mol2Id] 
            : [mol2Id, mol1Id];
        return `${smaller}:${larger}`;
    }

    /**
     * Verifica se una coppia di molecole è già stata elaborata
     */
    isPairProcessed(mol1Id, mol2Id) {
        const pairId = this.getPairId(mol1Id, mol2Id);
        return this.calculatedPairs.has(pairId);
    }

    /**
     * Marca una coppia come elaborata
     */
    markPairProcessed(mol1Id, mol2Id) {
        const pairId = this.getPairId(mol1Id, mol2Id);
        this.calculatedPairs.add(pairId);
    }

    /**
     * Memorizza il risultato dell'interazione nella cache
     */
    cacheInteraction(mol1Id, mol2Id, interactionData) {
        const pairId = this.getPairId(mol1Id, mol2Id);
        this.interactionCache.set(pairId, {
            data: interactionData,
            timestamp: performance.now()
        });
    }

    /**
     * Recupera i dati di interazione dalla cache, se disponibili
     */
    getCachedInteraction(mol1Id, mol2Id) {
        const pairId = this.getPairId(mol1Id, mol2Id);
        const cached = this.interactionCache.get(pairId);
        
        if (cached) {
            this.stats.cacheHits++;
            return cached.data;
        }
        
        return null;
    }

    /**
     * Determina l'ordine ottimale di elaborazione per un gruppo di molecole
     * Costruisce una lista di coppie da elaborare evitando duplicazioni
     */
    buildInteractionWorkPlan(molecules) {
        // Reset del piano di elaborazione
        this.resetWorkSession();
        
        const workPlan = [];
        const moleculeCount = molecules.length;
        
        // Costruisci la lista utilizzando una matrice triangolare superiore
        // In questo modo ogni coppia viene considerata una sola volta
        for (let i = 0; i < moleculeCount; i++) {
            for (let j = i + 1; j < moleculeCount; j++) {
                const mol1 = molecules[i];
                const mol2 = molecules[j];
                
                // Controlla se la coppia è già stata calcolata in precedenza
                if (this.isPairProcessed(mol1.id, mol2.id)) {
                    this.stats.calculationsAvoided++;
                    continue;
                }
                
                // Aggiungi la coppia al piano di lavoro
                workPlan.push({
                    mol1Index: i,
                    mol2Index: j,
                    mol1Id: mol1.id,
                    mol2Id: mol2.id,
                    // Priorità basata sulla distanza se disponibile nella cache
                    priority: this.estimateInteractionPriority(mol1, mol2)
                });
                
                this.stats.totalCalculations++;
            }
        }
        
        // Ordina il piano di lavoro in base alla priorità
        return workPlan.sort((a, b) => b.priority - a.priority);
    }
    
    /**
     * Stima una priorità per l'interazione basata sulla massa e posizione relativa
     */
    estimateInteractionPriority(mol1, mol2) {
        // Calcola distanza approssimativa
        const pos1 = mol1._physicsPosition || mol1.position;
        const pos2 = mol2._physicsPosition || mol2.position;
        
        let distanceSquared = 0;
        for (let i = 0; i < 3; i++) {
            const diff = pos2[i] - pos1[i];
            distanceSquared += diff * diff;
        }
        
        // Molecole più vicine e con massa maggiore hanno priorità più alta
        const combinedMass = (mol1.mass || 1) + (mol2.mass || 1);
        return combinedMass / (distanceSquared + 1);
    }
    
    /**
     * Resetta le strutture dati per una nuova sessione di lavoro
     * ma mantiene la cache delle interazioni per riutilizzarla
     */
    resetWorkSession() {
        this.calculatedPairs.clear();
        this.stats = {
            calculationsAvoided: 0,
            totalCalculations: 0,
            cacheHits: 0
        };
    }
    
    /**
     * Pulisce le interazioni obsolete dalla cache
     */
    cleanupCache(maxAge = 5000) {
        const now = performance.now();
        let removed = 0;
        
        for (const [pairId, entry] of this.interactionCache.entries()) {
            if (now - entry.timestamp > maxAge) {
                this.interactionCache.delete(pairId);
                removed++;
            }
        }
        
        return removed;
    }
    
    /**
     * Divide il piano di lavoro in lotti bilanciati per elaborazione parallela
     * Corretto per garantire che ciascuna interazione sia assegnata a una sola partizione
     */
    partitionWorkPlan(workPlan, workerCount) {
        if (workPlan.length === 0 || workerCount <= 1) {
            return [workPlan];
        }
        
        // Determina il numero ottimale di partizioni
        const partitionCount = Math.min(workerCount, Math.ceil(workPlan.length / 10));
        const partitions = Array(partitionCount).fill().map(() => []);
        
        // Distribuisci in modo che ogni partizione abbia elementi unici
        // Ogni elemento del workPlan va in una sola partizione
        for (let i = 0; i < workPlan.length; i++) {
            // Distribuisci in modo bilanciato tra le partizioni
            const partitionIndex = Math.floor(i / Math.ceil(workPlan.length / partitionCount));
            
            // Assicurati di non superare il numero di partizioni
            if (partitionIndex < partitionCount) {
                partitions[partitionIndex].push(workPlan[i]);
            } else {
                // Nel caso in cui ci siano più elementi di quanti previsti, aggiungi all'ultima partizione
                partitions[partitionCount - 1].push(workPlan[i]);
            }
        }
        
        return partitions;
    }
    
    /**
     * Genera un piano di distribuzione del lavoro completo per i subworker
     */
    generateDistributionPlan(molecules, workerCount) {
        // Crea il piano di lavoro completo
        const workPlan = this.buildInteractionWorkPlan(molecules);
        
        // Dividi il piano in partizioni bilanciate
        const partitions = this.partitionWorkPlan(workPlan, workerCount);
        
        // Prepara il piano di distribuzione
        return partitions.map(partition => ({
            interactions: partition,
            moleculeIndices: this.getRequiredMoleculeIndices(partition)
        }));
    }
    
    /**
     * Determina quali molecole sono necessarie per una partizione
     */
    getRequiredMoleculeIndices(partition) {
        const indices = new Set();
        
        partition.forEach(item => {
            indices.add(item.mol1Index);
            indices.add(item.mol2Index);
        });
        
        return [...indices].sort((a, b) => a - b);
    }

    /**
     * Genera un piano di distribuzione del lavoro basato su interazioni prioritarie
     * @param {Array} molecules - Tutte le molecole nella simulazione
     * @param {number} workerCount - Numero di worker disponibili
     * @param {Map} significantPairs - Mappa delle coppie significative con priorità
     * @returns {Array} Piano di distribuzione ottimizzato
     */
    generateDistributionPlanWithPriorities(molecules, workerCount, significantPairs) {
        // Reset del piano di elaborazione
        this.resetWorkSession();
        
        const workPlan = [];
        const moleculeCount = molecules.length;
        const moleculeIndices = new Map();
        
        // Crea mappa degli indici per accesso veloce
        molecules.forEach((mol, index) => {
            moleculeIndices.set(mol.id, index);
        });
        
        // FASE 1: Aggiungi prima le interazioni significative
        if (significantPairs && significantPairs.size > 0) {
            for (const [pairId, pairInfo] of significantPairs.entries()) {
                const mol1Index = moleculeIndices.get(pairInfo.mol1Id);
                const mol2Index = moleculeIndices.get(pairInfo.mol2Id);
                
                // Verifica che entrambe le molecole esistano
                if (mol1Index === undefined || mol2Index === undefined) continue;
                
                // Controlla se la coppia è già stata calcolata in precedenza
                if (this.isPairProcessed(pairInfo.mol1Id, pairInfo.mol2Id)) {
                    this.stats.calculationsAvoided++;
                    continue;
                }
                
                // Aggiungi al piano con priorità elevata
                workPlan.push({
                    mol1Index: mol1Index,
                    mol2Index: mol2Index,
                    mol1Id: pairInfo.mol1Id,
                    mol2Id: pairInfo.mol2Id,
                    priority: pairInfo.priority || 10.0 // Priorità alta per default
                });
                
                this.stats.totalCalculations++;
            }
        }
        
        // FASE 2: Aggiungi un sottoinsieme di altre interazioni potenziali
        // Utilizza una matrice triangolare superiore ma con limite al numero di interazioni
        
        // Calcola limite interazioni basato sulla complessità
        const maxTotalInteractions = Math.min(
            2000, // Limite assoluto per prestazioni
            Math.ceil(moleculeCount * Math.sqrt(moleculeCount) * 0.2) // Limite scalato
        );
        
        // Se abbiamo già abbastanza interazioni significative, limitiamo quelle aggiuntive
        const remainingSlots = Math.max(0, maxTotalInteractions - workPlan.length);
        const skipFactor = Math.max(1, Math.floor((moleculeCount * (moleculeCount - 1) / 2) / remainingSlots));
        
        let interactionCounter = 0;
        for (let i = 0; i < moleculeCount; i++) {
            for (let j = i + 1; j < moleculeCount; j++) {
                interactionCounter++;
                
                // Salta alcune interazioni per limitare il carico
                if (interactionCounter % skipFactor !== 0 && workPlan.length >= maxTotalInteractions / 2) {
                    continue;
                }
                
                const mol1 = molecules[i];
                const mol2 = molecules[j];
                
                // Salta se già presente nelle significative o già processata
                const pairId = this.getPairId(mol1.id, mol2.id);
                if (significantPairs && significantPairs.has(pairId)) continue;
                if (this.isPairProcessed(mol1.id, mol2.id)) {
                    this.stats.calculationsAvoided++;
                    continue;
                }
                
                // Aggiungi la coppia al piano di lavoro
                workPlan.push({
                    mol1Index: i,
                    mol2Index: j,
                    mol1Id: mol1.id,
                    mol2Id: mol2.id,
                    // Priorità basata sulla distanza
                    priority: this.estimateInteractionPriority(mol1, mol2)
                });
                
                this.stats.totalCalculations++;
                
                // Limita il numero totale di interazioni per mantenere le prestazioni
                if (workPlan.length >= maxTotalInteractions) {
                    i = moleculeCount; // Forza uscita dal loop esterno
                    break;
                }
            }
        }
        
        // Ordina il piano di lavoro in base alla priorità (più alta prima)
        const sortedPlan = workPlan.sort((a, b) => b.priority - a.priority);
        
        // Dividi il piano in partizioni bilanciate
        const partitions = this.partitionWorkPlan(sortedPlan, workerCount);
        
        // Prepara il piano di distribuzione
        return partitions.map(partition => ({
            interactions: partition,
            moleculeIndices: this.getRequiredMoleculeIndices(partition)
        }));
    }
}


/**
 * Estendi EnhancedChemistry con il nuovo sistema di gestione interazioni
 */
class OptimizedChemistry extends EnhancedChemistry {
    constructor(rules, size, moleculeCount, maxNumber) {
        super(rules, size, moleculeCount, maxNumber);
        
        // Inizializza il gestore interazioni
        this.interactionManager = new MoleculeInteractionManager();
    }
    
    /**
     * Distribuisce il lavoro ai sub-worker in modo ottimizzato
     * utilizzando la cache delle interazioni significative
     */
    async distributeWorkToSubWorkers() {
        if (this.subWorkers.length === 0 || this.isPaused) return;
        
        // Attendi che tutti i sub-worker siano disponibili
        await this.waitForAvailableWorkers();
        
        // PASSO 1: Identifica le interazioni significative usando la cache
        const significantPairs = this.identifySignificantInteractions();
        
        // PASSO 2: Genera piano di distribuzione ottimizzato basato sulle interazioni significative
        const distributionPlan = this.interactionManager.generateDistributionPlanWithPriorities(
            this.molecules, 
            this.subWorkers.length,
            significantPairs
        );
        
        // Distribuisci il lavoro ai worker
        const processingPromises = [];
        for (let i = 0; i < Math.min(this.subWorkers.length, distributionPlan.length); i++) {
            const partition = distributionPlan[i];
            const subWorker = this.subWorkers[i];
            
            if (!partition || partition.interactions.length === 0) continue;
            
            subWorker.busy = true;
            
            const promise = this.assignWorkToSubWorker(subWorker, partition);
            processingPromises.push(promise);
        }
        
        // Attendi che tutti i worker completino
        await Promise.all(processingPromises);
        
        // Esegui pulizia periodica della cache
        if (Math.random() < 0.05) {
            const removed = this.interactionManager.cleanupCache();
            const moleculeRemoved = this.moleculeCache.cleanupStaleRelationships();
            console.log(`Rimosse ${removed} interazioni obsolete dalla cache e ${moleculeRemoved} relazioni obsolete`);
        }
    }

    /**
     * Identifica le coppie di molecole con interazioni significative
     * utilizzando la cache delle molecole e le proprietà fisiche
     */
    identifySignificantInteractions() {
        const significantPairs = new Map();
        const interactionThreshold = this.rules.getConstant('interaction_threshold') || 5.0;
        const now = performance.now();
        
        // 1. Considera prima le relazioni già presenti nella cache delle molecole
        for (const [molId, relationships] of this.moleculeCache.moleculeRelationships.entries()) {
            const mol1 = this.molecules.find(m => m.id === molId);
            if (!mol1) continue;
            
            for (const rel of relationships) {
                // Salta relazioni obsolete
                if (now - rel.lastUpdated > this.moleculeCache.stalenessThreshold) continue;
                
                const mol2 = this.molecules.find(m => m.id === rel.otherId);
                if (!mol2) continue;
                
                // Se la distanza è abbastanza piccola, considera l'interazione significativa
                const currentDistance = this.calculateDistance(mol1, mol2);
                if (currentDistance <= this.calculateMaxRelationshipDistance(mol1, mol2)) {
                    const pairId = this.interactionManager.getPairId(mol1.id, mol2.id);
                    significantPairs.set(pairId, {
                        mol1Id: mol1.id,
                        mol2Id: mol2.id,
                        priority: 1.0 / (currentDistance + 0.1) * (mol1.mass + mol2.mass)
                    });
                }
            }
        }
        
        // 2. Aggiungi alcune interazioni casuali per scoprire nuove relazioni potenziali
        // ma limita il numero per evitare di sovraccaricare il sistema
        const randomPairCount = Math.min(50, Math.ceil(this.molecules.length * 0.05));
        const randomInteractionRate = this.rules.getConstant('random_interaction_probability') || 0.05;
        
        if (randomInteractionRate > 0 && this.molecules.length > 10) {
            for (let attempts = 0; attempts < randomPairCount * 2; attempts++) {
                if (significantPairs.size >= significantPairs.size + randomPairCount) break;
                
                // Scegli molecole casuali
                const idx1 = Math.floor(Math.random() * this.molecules.length);
                let idx2 = Math.floor(Math.random() * this.molecules.length);
                
                // Evita coppie identiche
                while (idx1 === idx2 && this.molecules.length > 1) {
                    idx2 = Math.floor(Math.random() * this.molecules.length);
                }
                
                const mol1 = this.molecules[idx1];
                const mol2 = this.molecules[idx2];
                const pairId = this.interactionManager.getPairId(mol1.id, mol2.id);
                
                // Aggiungi solo se non già considerata
                if (!significantPairs.has(pairId)) {
                    const distance = this.calculateDistance(mol1, mol2);
                    significantPairs.set(pairId, {
                        mol1Id: mol1.id,
                        mol2Id: mol2.id,
                        priority: 0.1 * (1.0 / (distance + 1.0) * (mol1.mass + mol2.mass))
                    });
                }
            }
        }
        
        return significantPairs;
    }
    
    /**
     * Assegna una partizione di lavoro a un sub-worker
     */
    assignWorkToSubWorker(subWorker, partition) {
        return new Promise(resolve => {
            const workerMessageHandler = (event) => {
                if (event.data.type === 'chunk_processed') {
                    subWorker.worker.removeEventListener('message', workerMessageHandler);
                    resolve();
                }
            };
            
            subWorker.worker.addEventListener('message', workerMessageHandler);
            
            // Prepara le molecole necessarie per questa partizione
            const requiredMolecules = partition.moleculeIndices.map(
                index => this.molecules[index]
            );
            
            const serializedMolecules = this.serializeMolecules(requiredMolecules);
            const cachedRelationships = this.moleculeCache.getRelationshipsForMolecules(requiredMolecules);
            
            // Invia al worker la partizione e le molecole necessarie
            subWorker.worker.postMessage({
                type: 'process_optimized_chunk',
                molecules: serializedMolecules,
                interactionPlan: partition.interactions,
                moleculeIndices: partition.moleculeIndices,
                cachedRelationships: cachedRelationships,
                temperature: this.temperature,
                timeScale: this.rules.getConstant('time_scale'),
                damping: this.rules.getConstant('damping'),
                isPaused: this.isPaused
            });
        });
    }
    
    /**
     * Versione ottimizzata dell'aggiornamento fisico
     * Utilizza un piano prestabilito per le interazioni
     * e sfrutta la cache per ridurre i calcoli
     */
    updatePhysicsWithInteractionPlan(moleculesToProcess, removedIndices, newMolecules) {
        const timeScale = this.rules.getConstant('time_scale');
        const damping = this.rules.getConstant('damping');
        const now = performance.now();
        
        // Pre-calcola forze per ogni molecola
        const forces = new Map();
        moleculesToProcess.forEach(mol => forces.set(mol.id, [0, 0, 0]));
        
        // 1. Prima usando le relazioni già in cache
        this.updatePhysicsFromCache(moleculesToProcess, removedIndices, newMolecules, forces, timeScale, now);
        
        // 2. Poi usa il piano per altre interazioni potenzialmente significative
        // ma limita il numero per evitare calcoli eccessivi
        const interactionPlan = this.interactionManager.buildInteractionWorkPlan(moleculesToProcess);
        
        // Limita il numero di interazioni da processare in base alla dimensione
        const maxInteractions = Math.min(
            1000,
            Math.ceil(moleculesToProcess.length * Math.sqrt(moleculesToProcess.length) * 0.1)
        );
        
        // Processa solo le interazioni a priorità più alta
        const limitedPlan = interactionPlan.slice(0, maxInteractions);
        
        // Processa interazioni secondo il piano limitato
        for (const interaction of limitedPlan) {
            const mol1 = moleculesToProcess[interaction.mol1Index];
            const mol2 = moleculesToProcess[interaction.mol2Index];
            
            // Salta se una delle molecole è stata rimossa in un'iterazione precedente
            if (removedIndices.has(interaction.mol1Index) || 
                removedIndices.has(interaction.mol2Index)) {
                continue;
            }
            
            // Calcola distanza attuale
            const currentDistance = this.calculateDistance(mol1, mol2);
            
            // Aggiorna o crea relazione nella cache solo se la distanza è significativa
            const interactionThreshold = this.calculateInteractionThreshold(mol1, mol2);
            if (currentDistance <= interactionThreshold) {
                this.updateMoleculeRelationship(mol1, mol2, currentDistance, now);
                
                // Calcola e applica forze solo se abbastanza vicine
                const maxDistance = this.calculateMaxRelationshipDistance(mol1, mol2);
                if (currentDistance <= maxDistance) {
                    this.processInteraction(mol1, mol2, currentDistance, forces, removedIndices, newMolecules, now, timeScale);
                }
            }
            
            // Segna questa interazione come completata
            this.interactionManager.markPairProcessed(mol1.id, mol2.id);
        }
        
        // Applica moto termico e aggiorna posizioni
        this.applyForcesAndUpdatePositions(moleculesToProcess, forces, timeScale, damping);
    }

    /**
     * Elabora interazioni già memorizzate nella cache
     */
    updatePhysicsFromCache(molecules, removedIndices, newMolecules, forces, timeScale, now) {
        // Per ogni molecola
        for (const mol1 of molecules) {
            if (removedIndices.has(mol1.id)) continue;
            
            // Ottieni relazioni dalla cache
            const relationships = this.moleculeCache.getRelationshipsForMolecule(mol1);
            if (!relationships || relationships.length === 0) continue;
            
            // Per ogni relazione in cache
            for (const rel of relationships) {
                if (now - rel.lastUpdated > this.moleculeCache.stalenessThreshold) continue;
                
                // Trova l'altra molecola
                const mol2 = molecules.find(m => m.id === rel.otherId);
                if (!mol2 || removedIndices.has(mol2.id)) continue;
                
                // Ricalcola distanza attuale
                const currentDistance = this.calculateDistance(mol1, mol2);
                
                // Aggiorna timestamp e distanza nella cache
                rel.distance = currentDistance;
                rel.lastUpdated = now;
                
                // Verifica se la relazione è ancora significativa
                const maxDistance = this.calculateMaxRelationshipDistance(mol1, mol2);
                if (currentDistance <= maxDistance) {
                    this.processInteraction(mol1, mol2, currentDistance, forces, removedIndices, newMolecules, now, timeScale);
                }
            }
        }
    }

    /**
     * Processa una singola interazione tra due molecole
     */
    processInteraction(mol1, mol2, distance, forces, removedIndices, newMolecules, now, timeScale) {
        // Calcola forze
        const [force1, force2] = this.calculateForces(mol1, mol2, distance, now);
        
        // Applica forze calcolate
        const mol1Force = forces.get(mol1.id);
        const mol2Force = forces.get(mol2.id);
        
        if (mol1Force && mol2Force) {
            for (let k = 0; k < 3; k++) {
                mol1Force[k] += force1[k] * timeScale;
                mol2Force[k] += force2[k] * timeScale;
            }
        }
        
        // Verifica reazioni solo se abbastanza vicine
        if (distance < this.calculateReactionDistance(mol1, mol2)) {
            if (this.shouldReact(mol1, mol2)) {
                // Trova indici nell'array originale
                const mol1OriginalIndex = this.molecules.findIndex(m => m.id === mol1.id);
                const mol2OriginalIndex = this.molecules.findIndex(m => m.id === mol2.id);
                
                if (mol1OriginalIndex >= 0 && mol2OriginalIndex >= 0) {
                    const products = this.processReaction(mol1, mol2);
                    if (products.length > 0) {
                        newMolecules.push(...products);
                        removedIndices.add(mol1OriginalIndex);
                        removedIndices.add(mol2OriginalIndex);
                        this.reactionCount++;
                        
                        // Rimuovi relazioni per molecole reagenti
                        this.moleculeCache.removeAllRelationshipsForMolecule(mol1.id);
                        this.moleculeCache.removeAllRelationshipsForMolecule(mol2.id);
                    }
                }
            }
        }
    }

    /**
     * Calcola soglia di interazione tra due molecole
     */
    calculateInteractionThreshold(mol1, mol2) {
        // Calcola soglia basata su massa e carica
        const massFactor = Math.sqrt(mol1.mass + mol2.mass);
        const chargeFactor = Math.abs(mol1.charge) + Math.abs(mol2.charge);
        
        return Math.max(
            3.0, // Minima distanza di interazione
            massFactor * 0.8 + chargeFactor * 0.5
        );
    }
    
    /**
     * Aggiorna una relazione tra molecole nella cache
     */
    updateMoleculeRelationship(mol1, mol2, distance, timestamp) {
        // Crea o aggiorna relazione nella cache
        this.moleculeCache.createRelationship(mol1, mol2, distance, timestamp);
        
        // Memorizza anche nell'interaction manager
        this.interactionManager.cacheInteraction(mol1.id, mol2.id, {
            distance: distance,
            timestamp: timestamp
        });
    }
    
    /**
     * Applica forze e aggiorna posizioni delle molecole
     */
    applyForcesAndUpdatePositions(molecules, forces, timeScale, damping) {
        for (const mol of molecules) {
            const molForce = forces.get(mol.id);
            if (!molForce) continue;
            
            // Aggiorna velocità con forze
            for (let k = 0; k < 3; k++) {
                const acceleration = molForce[k] / mol.mass;
                mol.velocity[k] += acceleration * 0.1;
                mol.position[k] += mol.velocity[k] * timeScale;
                mol.velocity[k] *= damping;
            }
            
            // Verifica confini
            this.enforceBoundaries(mol);
        }
    }
    
    /**
     * Override dell'aggiornamento fisico standard per usare l'approccio ottimizzato
     */
    updatePhysicsWithMoleculeCache(moleculesToProcess, removedIndices, newMolecules) {
        // Usa il nuovo metodo ottimizzato
        this.updatePhysicsWithInteractionPlan(moleculesToProcess, removedIndices, newMolecules);
    }
    
    /**
     * Implementazione migliorata del raggruppamento di molecole in cluster
     * usando l'approccio non ridondante
     */
    groupMoleculesByClusters() {
        if (this.molecules.length <= 1) {
            return [this.molecules];
        }
        
        const clusters = [];
        const visited = new Set();
        const proximityThreshold = this.size * 0.15;
        
        // Helper per trovare cluster utilizzando il nuovo sistema
        const findCluster = (startMol) => {
            const cluster = [startMol];
            const queue = [startMol];
            visited.add(startMol.id);
            
            while (queue.length > 0) {
                const current = queue.shift();
                
                // Usa relazioni dalla cache quando possibile
                const relationships = this.moleculeCache.getRelationshipsForMolecule(current);
                
                for (const rel of relationships) {
                    const neighborMol = this.molecules.find(m => m.id === rel.otherId);
                    if (!neighborMol || visited.has(neighborMol.id)) continue;
                    
                    // Verifica se relazione è ancora valida in termini di distanza
                    const distance = rel.distance;
                    if (distance < proximityThreshold) {
                        cluster.push(neighborMol);
                        queue.push(neighborMol);
                        visited.add(neighborMol.id);
                    }
                }
                
                // Per molecole senza relazioni nella cache, crea nuovo piano di interazione
                const unvisitedMolecules = this.molecules.filter(m => 
                    !visited.has(m.id) && m.id !== current.id);
                
                if (unvisitedMolecules.length > 0) {
                    // Costruisci mini piano solo per questa molecola
                    const miniPlan = [];
                    for (const mol of unvisitedMolecules) {
                        if (this.interactionManager.isPairProcessed(current.id, mol.id)) {
                            continue;
                        }
                        
                        // Verifica distanza
                        const distance = this.calculateDistance(current, mol);
                        if (distance < proximityThreshold) {
                            miniPlan.push({
                                mol: mol,
                                distance: distance
                            });
                            
                            // Crea relazione nella cache
                            this.moleculeCache.createRelationship(
                                current, mol, distance, performance.now()
                            );
                            
                            // Marca come processata
                            this.interactionManager.markPairProcessed(current.id, mol.id);
                        }
                    }
                    
                    // Aggiungi molecole vicine al cluster
                    for (const item of miniPlan) {
                        if (!visited.has(item.mol.id)) {
                            cluster.push(item.mol);
                            queue.push(item.mol);
                            visited.add(item.mol.id);
                        }
                    }
                }
            }
            
            return cluster;
        };
        
        // Forma cluster naturali basati su vicinanza
        for (const mol of this.molecules) {
            if (!visited.has(mol.id)) {
                const cluster = findCluster(mol);
                clusters.push(cluster);
            }
        }
        
        // Bilancia dimensioni cluster
        this.balanceClusters(clusters);
        
        return clusters;
    }
}

/**
 * Gestione dei messaggi nei sub-worker che supporta il nuovo approccio
 * ottimizzato per le interazioni
 */
function handleOptimizedProcessChunk(message) {
    if (!simulation || !isSubWorker) return;
    
    try {
        // Verifica stato pausa
        if (message.isPaused) {
            sendEmptyChunkResult();
            return;
        }
        
        // Ottieni molecole e piano interazioni
        const serializedMolecules = message.molecules || [];
        const interactionPlan = message.interactionPlan || [];
        const moleculeIndices = message.moleculeIndices || [];
        
        if (serializedMolecules.length === 0 || interactionPlan.length === 0) {
            sendEmptyChunkResult();
            return;
        }
        
        // Importa cache relazioni
        importCachedRelationships(message.cachedRelationships);
        
        // Imposta parametri simulazione
        setSimulationParameters(message);
        
        // Deserializza molecole
        const moleculeObjs = serializedMolecules.map(
            molData => MoleculeSerializer.deserialize(molData, PrimeMolecule)
        );
        
        // Converti piano interazioni per adattarlo alle molecole deserializzate
        const adaptedPlan = adaptInteractionPlan(interactionPlan, moleculeObjs, moleculeIndices);
        
        // Processa chunk
        const removedIndices = new Set();
        const newMolecules = [];
        const reactionCount = simulation.reactionCount;
        
        // Usa piano ottimizzato per aggiornare molecole
        processInteractionPlan(moleculeObjs, adaptedPlan, removedIndices, newMolecules);
        
        // Prepara risultati
        const results = prepareProcessedResults(
            moleculeObjs, 
            removedIndices, 
            newMolecules, 
            reactionCount
        );
        
        // Invia risultati
        postMessage({
            type: 'chunk_processed',
            workerId: workerId,
            processedCount: interactionPlan.length,
            results: results
        });
        
    } catch (error) {
        reportProcessingError(error);
    }
}

/**
 * Aggiornamento al gestore messaggi worker principale
 * Supporta il nuovo tipo di messaggio 'process_optimized_chunk'
 */
onmessage = async function(event) {
    try {
        const message = event.data;

        if (message.type === 'init_sub') {
            handleSubWorkerInitialization(message);
            return;
        }

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
                // Metodo legacy per compatibilità
                handleProcessChunk(message);
                break;
                
            case 'process_optimized_chunk':
                // Nuovo metodo ottimizzato
                handleOptimizedProcessChunk(message);
                break;

            case 'step':
                await simulation.step();
                if (!isSubWorker) {
                    sendUpdate();
                }
                break;

            case 'cleanup':
                cleanupResources();
                break;

            case 'set_temperature':
                if (typeof message.value === 'number') {
                    simulation.temperature = message.value;
                    console.log(`Worker ${workerId}: temperature set to ${message.value}`);
                    if (!isSubWorker) {
                        sendUpdate();
                    }
                }
                break;

            case 'set_timescale':
                if (typeof message.value === 'number') {
                    simulation.setTimeScale(message.value);
                    console.log(`Worker ${workerId}: timeScale set to ${message.value}`);
                }
                break;

            case 'set_pause':
                if (typeof message.isPaused === 'boolean') {
                    simulation.setPauseState(message.isPaused);
                    console.log(`Worker ${workerId}: pause state set to ${message.isPaused}`);

                    if (!isSubWorker && message.isPaused) {
                        sendUpdate();
                    }
                }
                break;

            case 'add_molecules':
                const count = message.count || 20;
                simulation.addRandomMolecules(count);
                console.log(`Worker ${workerId}: added ${count} new molecules`);
                if (!isSubWorker) {
                    sendUpdate();
                }
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

///
/// Worker and sub-workers
///

/**
 * Inizializza questo worker come sub-worker
 * @param {Object} message - Messaggio di inizializzazione
 */
function handleSubWorkerInitialization(message) {
    // Configura come sub-worker
    isSubWorker = true;
    workerId = message.workerId;
    console.log(`Started sub-worker ${workerId}`);

    // Crea istanza regole
    const rules = createCustomRules();
    rules.setConstant('time_scale', message.timeScale || 0.1);

    // Inizializza con set molecole vuoto - riceveremo chunk da processare
    simulation = new OptimizedChemistry(rules, message.size, 0, message.maxNumber);
    simulation.isMainWorker = false;
    simulation.workerId = workerId;

    // Segnala disponibilità
    postMessage({
        type: 'sub_ready',
        workerId: workerId
    });
}

/**
 * Gestisce un chunk standard (compatibilità legacy)
 * @param {Object} message - Messaggio con chunk da processare
 */
function handleProcessChunk(message) {
    if (!simulation || !isSubWorker) return;

    try {
        // Verifica stato pausa
        if (message.isPaused) {
            sendEmptyChunkResult();
            return;
        }

        // Ottieni molecole da processare
        const molecules = message.molecules || [];
        if (molecules.length === 0) {
            sendEmptyChunkResult();
            return;
        }

        // Importa relazioni dalla cache
        if (message.cachedRelationships) {
            importCachedRelationships(message.cachedRelationships);
        }

        // Imposta parametri simulazione
        setSimulationParameters(message);

        // Deserializza molecole
        const moleculeObjs = molecules.map(molData => 
            MoleculeSerializer.deserialize(molData, PrimeMolecule)
        );

        // Processa chunk
        const removedIndices = new Set();
        const newMolecules = [];
        const reactionCount = simulation.reactionCount;

        // Aggiorna fisica usando caching centrato su molecole
        simulation.updatePhysicsWithMoleculeCache(moleculeObjs, removedIndices, newMolecules);

        // Prepara risultati
        const results = {
            // Aggiornamenti molecole
            moleculeUpdates: moleculeObjs
                .filter((_, i) => !removedIndices.has(i))
                .map(mol => MoleculeSerializer.serialize(mol)),

            // Nuove molecole
            newMolecules: newMolecules.map(mol => 
                MoleculeSerializer.serialize(mol)
            ),

            // Numero reazioni avvenute
            reactionCount: simulation.reactionCount - reactionCount,
            
            // Aggiornamenti relazioni
            updatedRelationships: simulation.moleculeCache.getRelationshipsForMolecules(
                moleculeObjs.filter((_, i) => !removedIndices.has(i))
            )
        };

        // Invia risultati
        postMessage({
            type: 'chunk_processed',
            workerId: workerId,
            processedCount: molecules.length,
            results: results
        });

    } catch (error) {
        reportProcessingError(error);
    }
}

/**
 * Invia risultato vuoto per chunk
 */
function sendEmptyChunkResult() {
    postMessage({
        type: 'chunk_processed',
        workerId: workerId,
        processedCount: 0,
        results: {
            moleculeUpdates: [],
            newMolecules: [],
            reactionCount: 0,
            updatedRelationships: {}
        }
    });
}

/**
 * Importa relazioni dalla cache
 * @param {Object} relationships - Relazioni da importare
 */
function importCachedRelationships(relationships) {
    if (!relationships) return;
    
    for (const [molId, rels] of Object.entries(relationships)) {
        for (const rel of rels) {
            simulation.moleculeCache._ensureRelationshipsList(molId);
            simulation.moleculeCache._ensureRelationshipsList(rel.otherId);
            
            // Verifica se relazione già esiste
            if (!simulation.moleculeCache.hasRelationship(molId, rel.otherId)) {
                const relEntry = {
                    otherId: rel.otherId,
                    distance: rel.distance,
                    lastUpdated: rel.lastUpdated || performance.now()
                };
                
                simulation.moleculeCache.moleculeRelationships.get(molId).push(relEntry);
                
                // Aggiungi relazione inversa se necessario
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

/**
 * Imposta parametri di simulazione
 * @param {Object} message - Messaggio con parametri
 */
function setSimulationParameters(message) {
    simulation.temperature = message.temperature || 1.0;
    simulation.rules.setConstant('time_scale', message.timeScale || 0.1);
    simulation.rules.setConstant('damping', message.damping || 0.99);
}

/**
 * Segnala errore elaborazione
 * @param {Error} error - Errore riscontrato
 */
function reportProcessingError(error) {
    console.error(`Sub-worker ${workerId}: errore elaborazione chunk`, error);
    postMessage({
        type: 'error',
        workerId: workerId,
        message: `Error processing chunk: ${error.message}`
    });
}

/**
 * Adatta il piano di interazione alle molecole deserializzate
 * @param {Array} plan - Piano interazioni originale
 * @param {Array} molecules - Molecole deserializzate 
 * @param {Array} indices - Indici molecole
 * @returns {Array} Piano adattato
 */
function adaptInteractionPlan(plan, molecules, indices) {
    return plan.map(interaction => ({
        mol1: molecules[indices.indexOf(interaction.mol1Index)],
        mol2: molecules[indices.indexOf(interaction.mol2Index)],
        mol1Id: interaction.mol1Id,
        mol2Id: interaction.mol2Id,
        priority: interaction.priority
    }));
}

/**
 * Elabora le interazioni secondo il piano ottimizzato
 * @param {Array} molecules - Molecole da processare
 * @param {Array} plan - Piano interazioni
 * @param {Set} removedIndices - Set per indici rimossi
 * @param {Array} newMolecules - Array per nuove molecole
 */
function processInteractionPlan(molecules, plan, removedIndices, newMolecules) {
    const timeScale = simulation.rules.getConstant('time_scale');
    const damping = simulation.rules.getConstant('damping');
    const now = performance.now();
    
    // Pre-calcola forze per ogni molecola
    const forces = new Map();
    molecules.forEach(mol => forces.set(mol.id, [0, 0, 0]));
    
    // Processa ogni interazione pianificata
    for (const interaction of plan) {
        const mol1 = interaction.mol1;
        const mol2 = interaction.mol2;
        
        // Salta se una delle molecole è stata rimossa
        if (!mol1 || !mol2 || removedIndices.has(mol1.id) || removedIndices.has(mol2.id)) {
            continue;
        }
        
        // Calcola distanza attuale
        const currentDistance = simulation.calculateDistance(mol1, mol2);
        
        // Aggiorna relazione nella cache
        simulation.moleculeCache.createRelationship(mol1, mol2, currentDistance, now);
        
        // Calcola forze se molecole abbastanza vicine
        if (currentDistance <= simulation.calculateMaxRelationshipDistance(mol1, mol2)) {
            const [force1, force2] = simulation.calculateForces(mol1, mol2, currentDistance, now);
            
            // Applica forze
            const mol1Force = forces.get(mol1.id);
            const mol2Force = forces.get(mol2.id);
            
            if (mol1Force && mol2Force) {
                for (let k = 0; k < 3; k++) {
                    mol1Force[k] += force1[k] * timeScale;
                    mol2Force[k] += force2[k] * timeScale;
                }
            }
            
            // Verifica reazioni
            processReactionIfNeeded(mol1, mol2, currentDistance, removedIndices, newMolecules);
        }
    }
    
    // Applica forze e aggiorna posizioni
    applyForcesAndMove(molecules, forces, timeScale, damping);
}

/**
 * Processa una possibile reazione tra molecole
 * @param {Object} mol1 - Prima molecola
 * @param {Object} mol2 - Seconda molecola
 * @param {number} distance - Distanza tra molecole
 * @param {Set} removedIndices - Set per indici rimossi
 * @param {Array} newMolecules - Array per nuove molecole
 */
function processReactionIfNeeded(mol1, mol2, distance, removedIndices, newMolecules) {
    if (distance < simulation.calculateReactionDistance(mol1, mol2)) {
        if (simulation.shouldReact(mol1, mol2)) {
            const products = simulation.processReaction(mol1, mol2);
            if (products.length > 0) {
                newMolecules.push(...products);
                
                // Aggiungi a indici da rimuovere
                const mol1Index = simulation.molecules.findIndex(m => m.id === mol1.id);
                const mol2Index = simulation.molecules.findIndex(m => m.id === mol2.id);
                
                if (mol1Index >= 0) removedIndices.add(mol1.id);
                if (mol2Index >= 0) removedIndices.add(mol2.id);
                
                simulation.reactionCount++;
                
                // Rimuovi relazioni
                simulation.moleculeCache.removeAllRelationshipsForMolecule(mol1.id);
                simulation.moleculeCache.removeAllRelationshipsForMolecule(mol2.id);
            }
        }
    }
}

/**
 * Applica forze e aggiorna posizioni delle molecole
 * @param {Array} molecules - Molecole da aggiornare
 * @param {Map} forces - Mappa delle forze
 * @param {number} timeScale - Scala temporale
 * @param {number} damping - Fattore smorzamento
 */
function applyForcesAndMove(molecules, forces, timeScale, damping) {
    for (const mol of molecules) {
        const molForce = forces.get(mol.id);
        if (!molForce) continue;
        
        // Applica moto termico casuale
        for (let axis = 0; axis < 3; axis++) {
            mol.velocity[axis] += (Math.random() - 0.5) * 0.01 * simulation.temperature;
        }
        
        // Aggiorna velocità con forze
        for (let k = 0; k < 3; k++) {
            const acceleration = molForce[k] / mol.mass;
            mol.velocity[k] += acceleration * 0.1;
            mol.position[k] += mol.velocity[k] * timeScale;
            mol.velocity[k] *= damping;
        }
        
        // Verifica confini
        simulation.enforceBoundaries(mol);
    }
}

/**
 * Prepara i risultati dell'elaborazione
 * @param {Array} molecules - Molecole processate
 * @param {Set} removedIds - ID molecole rimosse
 * @param {Array} newMolecules - Nuove molecole create
 * @param {number} initialReactionCount - Conteggio reazioni iniziale
 * @returns {Object} Risultati elaborazione
 */
function prepareProcessedResults(molecules, removedIds, newMolecules, initialReactionCount) {
    return {
        // Aggiornamenti molecole
        moleculeUpdates: molecules
            .filter(mol => !removedIds.has(mol.id))
            .map(mol => MoleculeSerializer.serialize(mol)),
            
        // Nuove molecole create
        newMolecules: newMolecules.map(mol => 
            MoleculeSerializer.serialize(mol)
        ),
        
        // Numero reazioni avvenute
        reactionCount: simulation.reactionCount - initialReactionCount,
        
        // Aggiornamenti relazioni
        updatedRelationships: simulation.moleculeCache.getRelationshipsForMolecules(
            molecules.filter(mol => !removedIds.has(mol.id))
        )
    };
}

/**
 * Ottiene dati ottimizzati molecole per aggiornamento
 * @param {Map} stablePositions - Posizioni stabilizzate
 * @returns {Object} Dati ottimizzati
 */
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

/**
 * Conta relazioni totali nella cache
 * @returns {number} Numero totale relazioni
 */
function countTotalRelationships() {
    if (!simulation || !simulation.moleculeCache) return 0;

    let count = 0;
    for (const relationships of simulation.moleculeCache.moleculeRelationships.values()) {
        count += relationships.length;
    }
    return count;
}

/**
 * Statistiche avanzate per il monitoraggio delle prestazioni
 * @returns {Object} Statistiche prestazioni
 */
function getPerformanceStats() {
    if (!simulation || !simulation.interactionManager) {
        return {};
    }
    
    return {
        // Statistiche molecole
        moleculeCount: simulation.molecules.length,
        totalReactions: simulation.reactionCount,
        
        // Statistiche cache e calcoli
        relationshipCacheSize: simulation.moleculeCache.moleculeRelationships.size,
        totalCachedRelationships: countTotalRelationships(),
        
        // Statistiche ottimizzazione
        calculationsAvoided: simulation.interactionManager.stats.calculationsAvoided,
        totalCalculations: simulation.interactionManager.stats.totalCalculations,
        cacheHits: simulation.interactionManager.stats.cacheHits,
        
        // Statistiche worker
        workerCount: simulation.subWorkers.length,
        busyWorkers: simulation.subWorkers.filter(w => w.busy).length,
        
        // Stato simulazione
        isPaused: simulation.isPaused,
        temperature: simulation.temperature,
        elapsedTime: simulation.accumulatedTime
    };
}

/**
 * Invia aggiornamento migliorato con statistiche prestazioni
 */
function sendUpdate() {
    try {
        // Prepara dati usando posizioni stabili quando disponibili
        const moleculeData = getOptimizedMoleculeData(simulation.stablePositions);
        const perfStats = getPerformanceStats();

        postMessage({
            type: 'update',
            molecules: moleculeData,
            temperature: simulation.temperature,
            reactionCount: simulation.reactionCount,
            isPaused: simulation.isPaused,
            
            // Statistiche avanzate
            stats: perfStats,
            
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

/**
 * Pulisce risorse
 */
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

/**
 * Funzione migliorata per inizializzazione
 * @param {Object} message - Messaggio di inizializzazione
 */
async function handleInitialization(message) {
    const { size, moleculeCount, maxNumber, timeScale } = message;
    const rules = createCustomRules();
    rules.setConstant('time_scale', timeScale || 0.1);

    // Inizializza simulazione con versione ottimizzata
    simulation = new OptimizedChemistry(rules, size, moleculeCount, maxNumber);
    simulation.isMainWorker = !isSubWorker;
    simulation.workerId = workerId;

    console.log(`Worker ${workerId}: simulation initialized with ${moleculeCount} molecules, timeScale=${timeScale}`);

    // Inizializza sub-worker
    if (!isSubWorker) {
        try {
            await simulation.initializeSubWorkers();
            console.log(`Initialized ${simulation.subWorkers.length} sub-workers`);
        } catch (error) {
            console.error("Failed to initialize sub-workers:", error);
            // Continua in modalità single-worker
        }
    }

    // Passo iniziale
    await simulation.step();

    if (!isSubWorker) {
        sendUpdate();
    }
}