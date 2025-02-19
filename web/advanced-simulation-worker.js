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
const INTERACTION_CACHE_LIFETIME = 50; // Frame per cui mantenere cache interazione

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
                subWorker.worker.postMessage({
                    type: 'cleanup'
                });
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
            return [
                [0, 0, 0],
                [0, 0, 0]
            ];
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
                subWorker.worker.postMessage({
                    type: 'cleanup'
                });
                subWorker.worker.terminate();
            }
            this.subWorkers = [];
        }

        this.molecules = [];
        this.moleculeCache.clear();
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
    const {
        size,
        moleculeCount,
        maxNumber,
        timeScale
    } = message;
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