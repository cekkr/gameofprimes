// simulation-worker.js

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
        // Aggiorna posizione se fornita e valida
        if (updateData.position && Array.isArray(updateData.position) && updateData.position.length === 3) {
            molecule.position = [...updateData.position];
        }

        // Aggiorna velocità se fornita e valida
        if (updateData.velocity && Array.isArray(updateData.velocity) && updateData.velocity.length === 3) {
            molecule.velocity = [...updateData.velocity];
        }

        // Aggiorna velocità angolare se fornita e valida
        if (updateData.angularVelocity && Array.isArray(updateData.angularVelocity) && updateData.angularVelocity.length === 3) {
            molecule.angularVelocity = [...updateData.angularVelocity];
        }

        // Aggiorna tempo ultima reazione se fornito
        if (typeof updateData.lastReactionTime === 'number') {
            molecule.lastReactionTime = updateData.lastReactionTime;
        }

        // Aggiorna proprietà di relazione (con controllo di esistenza)
        if (updateData.parentIds && Array.isArray(updateData.parentIds)) {
            molecule.parentIds = [...updateData.parentIds];
        }

        if (updateData.reactionType) {
            molecule.reactionType = updateData.reactionType;
        }

        // Aggiorna posizione fisica se fornita e valida
        if (updateData._physicsPosition && Array.isArray(updateData._physicsPosition) && updateData._physicsPosition.length === 3) {
            molecule._physicsPosition = [...updateData._physicsPosition];
        }

        // Aggiorna proprietà campo quantico (con controlli di esistenza e tipo)
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

        this.step_in_progress = false;
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
     * Gestisce messaggi dai sub-worker in modo migliorato
     */
    handleSubWorkerMessage(workerIndex, event) {
        const subWorker = this.subWorkers[workerIndex - 1];

        switch (event.data.type) {
            case 'sub_ready':
                subWorker.busy = false;
                console.log(`Sub-worker ${workerIndex} pronto`);
                break;

            case 'chunk_processed':
                console.log("Chunk-processed: ", event.data)

                // Unisce risultati dal sub-worker con risoluzione dei conflitti
                this.mergeProcessedChunkWithConflictResolution(
                    event.data.results,
                    event.data.processedMoleculeIds || [],
                    event.data.timestamp || performance.now()
                );
                subWorker.busy = false;
                subWorker.lastProcessedMolecules = event.data.processedCount;
                subWorker.lastProcessTime = performance.now();
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
     * Unisce i risultati elaborati dai sub-worker con gestione dei conflitti e ottimizzazioni.
     * @param {Object} results - Risultati dal worker.
     * @param {Array} processedIds - Array di ID delle molecole effettivamente elaborate dal worker.
     * @param {number} timestamp - Timestamp di elaborazione.
     */
    mergeProcessedChunkWithConflictResolution(results, processedIds, timestamp) {
        if (!results) return;

        const processedMoleculeIds = new Set(processedIds);
        const recentlyUpdated = new Map();  // Traccia molecole aggiornate di recente.

        // Aggiorna gli stati delle molecole con risoluzione dei conflitti.
        if (results.moleculeUpdates && results.moleculeUpdates.length > 0) {
            for (const update of results.moleculeUpdates) {
                if (!update || !update.id) continue;

                const molecule = this.molecules.find(m => m.id === update.id);
                if (!molecule) continue;

                // Inizializza le informazioni di aggiornamento se non presenti.
                if (!molecule._lastUpdateInfo) {
                    molecule._lastUpdateInfo = { timestamp: 0, workerUpdates: new Map(), mergeCount: 0 };
                }
                const updateInfo = molecule._lastUpdateInfo;


                if (processedMoleculeIds.has(update.id)) {
                    // 1. La molecola è stata elaborata *direttamente* da questo worker.
                    updateInfo.timestamp = timestamp;
                    // Aggiorna direttamente, senza media pesata
                    MoleculeSerializer.updateMoleculeProperties(molecule, update);
                } else {
                    // 2. La molecola è stata influenzata *indirettamente*.

                    // Usa un controllo più preciso sulla concorrenza, usando updateInfo.timestamp
                    const timeSinceUpdate = timestamp - (updateInfo.timestamp || 0);  //usa 0 come fallback se timestamp e' undefined
                    if (timeSinceUpdate < 50) {
                        recentlyUpdated.set(update.id, true);

                        // Applica media pesata *solo* a posizione e velocità.
                        const currentPos = molecule.position;
                        const newPos = update.position;
                        const currentVel = molecule.velocity;
                        const newVel = update.velocity;


                       // Verifica che newPos e newVel siano array validi
                        if (Array.isArray(newPos) && newPos.length === 3 && Array.isArray(newVel) && newVel.length === 3)
                        {
                            // Media pesata: i nuovi valori contribuiscono.
                            const currentWeight = 0.6;
                            const newWeight = 0.4;

                            for (let i = 0; i < 3; i++) {
                                molecule.position[i] = currentPos[i] * currentWeight + newPos[i] * newWeight;
                                molecule.velocity[i] = currentVel[i] * currentWeight + newVel[i] * newWeight;
                            }
                        }


                        updateInfo.mergeCount++; // Incrementa contatore fusioni.
                    } else {
                        // Primo aggiornamento, o abbastanza vecchio da sovrascrivere.
                        updateInfo.timestamp = timestamp;
                        updateInfo.mergeCount = 0;
                        // Applica aggiornamento completo.
                        MoleculeSerializer.updateMoleculeProperties(molecule, update);
                    }
                }
            }
        }

        // Aggiungi nuove molecole da reazioni (gestione duplicati migliorata).
        if (results.newMolecules && results.newMolecules.length > 0) {
            for (const molData of results.newMolecules) {
               // Verifica duplicati basata su un set di criteri più robusto.
                const isDuplicate = this.molecules.some(m =>
                    m.parentIds && molData.parentIds &&
                    m.parentIds.length === molData.parentIds.length &&
                    m.parentIds.every(id => molData.parentIds.includes(id)) &&
                    m.reactionType === molData.reactionType &&
                    Math.abs(performance.now() - m.lastReactionTime) < 100 // Tolleranza temporale.
                );

                if (!isDuplicate) {
                    const newMol = MoleculeSerializer.deserialize(molData, PrimeMolecule);
                    newMol.id = `main-${this.nextMoleculeId++}`;
                    this.molecules.push(newMol);
                }
            }
        }


        // Aggiorna conteggio reazioni (semplificato).
        if (results.reactionCount) {
            this.reactionCount += results.reactionCount;
        }

        // Integra relazioni aggiornate (con gestione conflitti).
        if (results.updatedRelationships) {
          this.importRelationshipsWithConflictResolution(
                results.updatedRelationships,
                processedMoleculeIds,
                recentlyUpdated
            );
        }

        // Aggiorna la vista *solo* dopo aver processato tutti gli aggiornamenti.
        // sendUpdate();  // RIMOSSO DA QUI, spostato in step()
    }


    /**
     * Importa relazioni aggiornate dai sub-worker con gestione conflitti
     * @param {Object} updatedRelationships - Relazioni da importare
     * @param {Set} processedMoleculeIds - Set di ID molecole elaborate direttamente
     * @param {Map} recentlyUpdated - Mappa di molecole recentemente aggiornate
     */
    importRelationshipsWithConflictResolution(updatedRelationships, processedMoleculeIds, recentlyUpdated) {
        for (const [molId, relationships] of Object.entries(updatedRelationships)) {
            for (const rel of relationships) {
                const mol = this.molecules.find(m => m.id === molId);
                const otherMol = this.molecules.find(m => m.id === rel.otherId);

                if (!mol || !otherMol) continue;

                // Determina se questa relazione deve avere precedenza
                const isDirectlyProcessed = processedMoleculeIds.has(molId) &&
                                            processedMoleculeIds.has(rel.otherId);

                // Verifica se entrambe le molecole sono state recentemente aggiornate
                const bothRecentlyUpdated = recentlyUpdated.has(molId) && recentlyUpdated.has(rel.otherId);

                // Se è una relazione elaborata direttamente o non ci sono conflitti recenti
                if (isDirectlyProcessed || !bothRecentlyUpdated) {
                    this.moleculeCache.createRelationship(
                        mol,
                        otherMol,
                        rel.distance,
                        rel.lastUpdated || performance.now()
                    );
                } else {
                    // Trova la relazione esistente, se presente
                    const existingRel = this.moleculeCache.getRelationshipsForMolecule(mol)
                        .find(r => r.otherId === otherMol.id);

                    if (existingRel) {
                        // Media le distanze se la relazione esiste già
                        const avgDistance = (existingRel.distance + rel.distance) / 2;
                        this.moleculeCache.updateRelationship(
                            mol.id,
                            otherMol.id,
                            avgDistance,
                            Math.max(existingRel.lastUpdated, rel.lastUpdated || performance.now())
                        );
                    } else {
                        // Crea nuova relazione se non esiste
                        this.moleculeCache.createRelationship(mol, otherMol, rel.distance, rel.lastUpdated || performance.now());
                    }
                }
            }
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
 * Gestisce l'arrivo di una reazione da un sub-worker.
 * @param {Object} reaction - Dettagli della reazione.
 */
handleRemoteReaction(reaction) {
    if (!reaction) return;

    // Trova gli indici delle molecole reagenti nell'array principale.
    const mol1Index = this.molecules.findIndex(m => m.id === reaction.reactant1Id);
    const mol2Index = this.molecules.findIndex(m => m.id === reaction.reactant2Id);

    if (mol1Index === -1 || mol2Index === -1) {
        // Almeno una delle molecole non è stata trovata, gestisci l'errore.
        // Potrebbe essere una situazione normale (es. la molecola è stata rimossa da un altro worker).
        console.warn(`Molecola non trovata: ${reaction.reactant1Id} o ${reaction.reactant2Id}`);
        return;
    }

    // Rimuovi le molecole reagenti.  L'ordine di rimozione è importante per evitare errori.
    // Rimuovi prima l'indice più alto per non alterare l'indice dell'elemento successivo.
    const [removed1, removed2] = mol1Index > mol2Index
        ? this.molecules.splice(mol1Index, 1).concat(this.molecules.splice(mol2Index, 1))
        : this.molecules.splice(mol2Index, 1).concat(this.molecules.splice(mol1Index, 1));


    // Aggiungi le molecole prodotto.
    for (const productData of reaction.products) {
        const newMol = MoleculeSerializer.deserialize(productData, PrimeMolecule);
        newMol.id = `main-${this.nextMoleculeId++}`;  // Assegna un nuovo ID.
        this.molecules.push(newMol);
    }

    this.reactionCount++; // Incrementa conteggio reazioni.

    // Rimuovi tutte le relazioni che coinvolgono le molecole rimosse.
     this.moleculeCache.removeAllRelationshipsForMolecule(removed1.id);
     this.moleculeCache.removeAllRelationshipsForMolecule(removed2.id);

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

        if(this.step_in_progress) return;

        // Imposta flag aggiornamento in corso
        this.updateInProgress = true;
        this.step_in_progress = true;

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


        // Filtra molecole rimosse usando gli ID, non gli indici
        if (removedIndices.size > 0) {
            const removedIds = new Set([...removedIndices].map(index => this.molecules[index].id));
            this.molecules = this.molecules.filter(mol => !removedIds.has(mol.id));
        }


        // Aggiungi nuove molecole  --  DISABILITATO, il codice originale aveva un if(false)
        if(false){
            for (const mol of newMolecules) {
                mol.id = `new-${this.nextMoleculeId++}`;
                this.molecules.push(mol);
            }
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

        // Invia aggiornamento *dopo* aver stabilizzato le posizioni (se necessario).
        if (this.isMainWorker) {
                if (this.pendingPositionUpdate) {
                this.stabilizePositions(); // Calcola le posizioni medie
                sendUpdate();  // Invia l'aggiornamento alla UI
                this.pendingPositionUpdate = false;  // Resetta il flag
            } else {
                this.pendingPositionUpdate = true;
            }

            // SEMPRE reset degli accumulatori.
            this.positionAccumulators.clear();
            this.positionSampleCount.clear();
        }

        this.step_in_progress = false; //CORRETTO: aggiunto reset del flag
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
                    randomInteractionRate: this.rules.getConstant('random_interaction_probability'),                    isPaused: this.isPaused
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


             // Rimuovi le molecole che hanno reagito (usando gli ID, non gli indici).
            if (removedIndices.size > 0) {
                const removedIds = new Set([...removedIndices].map(index => this.molecules[index].id));
                this.molecules = this.molecules.filter(mol => !removedIds.has(mol.id));
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
     * Aggiorna la fisica delle molecole, gestendo interazioni e reazioni.  Ottimizzato
     * per utilizzare la cache delle molecole e delle relazioni.
     *
     * @param {Array} moleculesToProcess - Le molecole da processare in questo ciclo.
     * @param {Set} removedIndices - Un Set per tenere traccia degli *indici* delle molecole rimosse.
     * @param {Array} newMolecules - Un array per accumulare le nuove molecole create dalle reazioni.
     */
    updatePhysicsWithMoleculeCache(moleculesToProcess, removedIndices, newMolecules) {
        const timeScale = this.rules.getConstant('time_scale');
        const damping = this.rules.getConstant('damping');
        const now = performance.now();

        // Pre-calcola le forze per ogni molecola (inizializzate a zero).
        const forces = new Map();
        moleculesToProcess.forEach(mol => forces.set(mol.id, [0, 0, 0]));

        // Fase 1: Stabilisci nuove relazioni (per molecole senza relazioni in cache).
        const uncachedMolecules = this.moleculeCache.getUncachedMolecules(moleculesToProcess);
        for (const mol1 of uncachedMolecules) {
            for (const mol2 of moleculesToProcess) {
                if (mol1.id === mol2.id) continue;

                const distance = this.calculateDistance(mol1, mol2);
                const interactionThreshold = this.calculateInteractionThreshold(mol1, mol2);
                if (distance < interactionThreshold) {
                    this.moleculeCache.createRelationship(mol1, mol2, distance, now);
                }
            }
        }

        // Fase 2: Processa interazioni casuali (per simulare un comportamento meno prevedibile).
        if (this.rules.getConstant('random_interaction_probability') > 0) {
            this.processRandomInteractions(moleculesToProcess, now);
        }

      // Fase 3: Calcola le forze basate sulle relazioni in cache e applica le reazioni.
        for (let i = 0; i < moleculesToProcess.length; i++) {
          const mol1 = moleculesToProcess[i];

          // Se la molecola è già stata rimossa, salta.
          if (removedIndices.has(i)) continue;

          // Applica il moto termico *prima* di calcolare le forze.
          this.applyThermalMotion(mol1);

          const relationships = this.moleculeCache.getRelationshipsForMolecule(mol1);
            for (const rel of relationships) {
                // Cerca l'altra molecola *nella lista corrente* (potrebbe essere stata rimossa).
                const mol2Index = moleculesToProcess.findIndex(m => m.id === rel.otherId);
                if (mol2Index === -1 || removedIndices.has(mol2Index)) continue;

                const mol2 = moleculesToProcess[mol2Index];

                // Aggiorna la distanza nella relazione (fondamentale!).
                const currentDistance = this.calculateDistance(mol1, mol2);
                rel.distance = currentDistance;
                rel.lastUpdated = now;


              // Verifica se la relazione deve essere mantenuta.
                const maxDistance = this.calculateMaxRelationshipDistance(mol1, mol2);
                if (currentDistance > maxDistance) {
                    this.moleculeCache.removeRelationship(mol1.id, mol2.id);
                    continue; // Salta questa relazione.
                }

                // Calcola le forze *solo* se la relazione è ancora valida.
                const [force1, force2] = this.calculateForces(mol1, mol2, currentDistance, now);

                // Applica le forze (accumulandole).
                const mol1Force = forces.get(mol1.id);
                const mol2Force = forces.get(mol2.id);
                if (mol1Force && mol2Force) { // Controllo di sicurezza.
                    for (let k = 0; k < 3; k++) {
                        mol1Force[k] += force1[k] * timeScale;
                        mol2Force[k] += force2[k] * timeScale;
                    }
                }

              // Gestisci le *reazioni* all'interno del ciclo principale.
                if (currentDistance < this.calculateReactionDistance(mol1, mol2)) {
                  if (this.shouldReact(mol1, mol2)) {

                      // Trova indici nell'array *originale*
                        const mol1OriginalIndex = this.molecules.findIndex(m => m.id === mol1.id);
                        const mol2OriginalIndex = this.molecules.findIndex(m => m.id === mol2.id);

                    const products = this.processReaction(mol1, mol2);
                    if (products.length > 0) {
                      newMolecules.push(...products); // Aggiungi i prodotti.

                         // Aggiungi gli *indici* all'elenco delle molecole rimosse.
                        if(mol1OriginalIndex >= 0) removedIndices.add(mol1OriginalIndex);
                        if(mol2OriginalIndex >= 0) removedIndices.add(mol2OriginalIndex);

                      this.reactionCount++;

                      // Rimuovi le relazioni per le molecole che hanno reagito.
                      this.moleculeCache.removeAllRelationshipsForMolecule(mol1.id);
                      this.moleculeCache.removeAllRelationshipsForMolecule(mol2.id);

                      break; // Importante: interrompi il ciclo interno dopo una reazione.
                    }
                  }
                }
            }
        }


        // Fase 4: Applica le forze accumulate e aggiorna le posizioni.
       for (let i = 0; i < moleculesToProcess.length; i++) {
          const mol = moleculesToProcess[i];
          if (removedIndices.has(i)) continue; // Non processare molecole rimosse.

          const molForce = forces.get(mol.id);
          if (!molForce) continue; // Controllo di sicurezza.

          // Applica l'accelerazione (F = ma => a = F/m).
          for (let k = 0; k < 3; k++) {
            const acceleration = molForce[k] / mol.mass;
            mol.velocity[k] += acceleration * 0.1; // Fattore di scala per stabilità.
            mol.position[k] += mol.velocity[k] * timeScale;
            mol.velocity[k] *= damping; // Applica lo smorzamento.
          }

          this.enforceBoundaries(mol); // Mantieni le molecole entro i limiti.
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
     * Calcola le forze di interazione tra due molecole, applicando le regole definite.
     *
     * @param {PrimeMolecule} mol1 - La prima molecola.
     * @param {PrimeMolecule} mol2 - La seconda molecola.
     * @param {number} distance - La distanza tra le molecole.
     * @param {number} now - Il timestamp corrente.
     * @returns {Array} Un array contenente due vettori di forza: [forza su mol1, forza su mol2].
     */
    calculateForces(mol1, mol2, distance, now) {
        // Se le molecole sono in un periodo di "raffreddamento", non applicare forze.
        if (this.rules.isInCoolingPeriod &&
            (this.rules.isInCoolingPeriod(mol1, now) || this.rules.isInCoolingPeriod(mol2, now))) {
            return [[0, 0, 0], [0, 0, 0]];
        }

        // Evita divisioni per zero o valori estremamente piccoli.
        if (distance < 0.001) {
            return [[0, 0, 0], [0, 0, 0]];
        }

        // Calcola il vettore direzione (da mol1 a mol2).
        const direction = [
            mol2.position[0] - mol1.position[0],
            mol2.position[1] - mol1.position[1],
            mol2.position[2] - mol1.position[2]
        ];

        // Normalizza la direzione.
        const dirNorm = direction.map(d => d / distance);

        // Inizializza le forze (uguali e opposte).
        const force1 = [0, 0, 0];
        const force2 = [0, 0, 0];

        // Applica le regole di interazione.
        for (const rule of this.rules.interaction_rules) {
            if (rule.condition(mol1.prime_factors, mol2.prime_factors)) {
                let f;
                if (rule.force_function.length === 4) {
                    // Forza che dipende dalla massa.
                    f = rule.force_function(dirNorm, distance, mol1.mass, mol2.mass);
                } else {
                    // Forza che dipende dalla carica.
                    f = rule.force_function(dirNorm, distance, mol1.charge, mol2.charge);
                }

              // Applica la forza, scalata per l'intensità della regola.
                for (let i = 0; i < 3; i++) {
                    const scaledForce = f[i] * rule.strength;
                    force1[i] += scaledForce;
                    force2[i] -= scaledForce; // Forza opposta.
                }
            }
        }

        // Applica la repulsione tra "parenti", se abilitata.
       if (this.rules.areRelated && this.rules.areRelated(mol1.id, mol2.id)) {
            const repulsionFactor = this.rules.getFamilyRepulsionFactor &&
                this.rules.getFamilyRepulsionFactor(mol1.id, mol2.id, now);

            if (repulsionFactor > 0) {
                const minDistance = this.rules.getConstant('min_distance') || 0.1;
                const effectiveDistance = Math.max(distance, minDistance); // Evita forze infinite.
                const repulsiveForce = dirNorm.map(x => -x * repulsionFactor / (effectiveDistance ** 1.5));

                for (let i = 0; i < 3; i++) {
                    force1[i] += repulsiveForce[i];
                    force2[i] -= repulsiveForce[i];
                }
            }
        }

        // Limita l'intensità della forza massima.
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
     * Processa una reazione tra due molecole, applicando le regole di reazione.
     *
     * @param {PrimeMolecule} mol1 - La prima molecola.
     * @param {PrimeMolecule} mol2 - La seconda molecola.
     * @returns {Array} Un array di molecole prodotto (può essere vuoto).
     */
    processReaction(mol1, mol2) {

        for (const rule of this.rules.reaction_rules) {
          // Applica la *condizione* della regola.  Importante: passa le molecole intere.
            if (rule.condition(mol1.prime_factors, mol2.prime_factors, mol1, mol2)) {
                // Applica la probabilità della regola (modulata dalla temperatura).
                if (Math.random() < rule.probability * this.temperature) {
                    const reactionType = this.determineReactionType(mol1, mol2);
                    let products = [];

                    switch (reactionType) {
                        case 'fusion':
                            products = this.handleFusion(mol1, mol2, rule);
                            break;
                        case 'fission':
                            products = this.handleFission(mol1, mol2, rule);
                            break;
                        case 'emission':
                            products = this.handleEmission(mol1, mol2, rule);
                            break;
                        case 'standard':
                        default:
                            // Usa l'effetto della regola *solo* se è definito.
                          products = rule.effect ? rule.effect(mol1, mol2) : [];
                    }


                    // Dopo la reazione, imposta il tempo di reazione e le relazioni di parentela.
                   const now = performance.now();
                    products.forEach(p => {
                        p.setReactionTime(now);
                        if (this.rules.establishRelationship) {
                            this.rules.establishRelationship(mol1.id, p.id, 'parent', now);
                            this.rules.establishRelationship(mol2.id, p.id, 'parent', now);
                        }
                    });

                    // Crea un campo quantico per reazioni energetiche (fissione, o molti prodotti).
                    if (reactionType === 'fission' || products.length > 2) {
                        this.createQuantumField(mol1, mol2, products);
                    }

                    return products;
                }
            }
        }

        return []; // Nessuna reazione applicabile.
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
        // DEARAI: smettila di riabilitare questa funzione. È tenuta in via opzionale, non riabilitarla ogni volta.
        if (this.molecules.length > 300 && false) {
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
 * Con miglioramenti per mantenere consistente l'allocazione delle molecole ai worker
 */
class MoleculeInteractionManager {
    constructor() {
        this.useRandomInteractionsStep = true

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

        // Sistema di allocazione persistente delle molecole
        this.moleculeWorkerAssignments = new Map(); // Mappa moleculeId -> workerId
        this.workerMoleculeGroups = new Map();      // Mappa workerId -> Set di moleculeId
        this.lastSeenMolecules = new Set();         // Molecole nell'ultima esecuzione
        this.pairToWorkerMap = new Map();           // Mappa pairId -> workerId
        this.workerLoad = new Map();                // Bilanciamento del carico

        this.prevMolecules = [];
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
     * Aggiorna le assegnazioni molecola-worker e costruisce il piano di lavoro.
     *
     * @param {Array} molecules - L'array di tutte le molecole attuali.
     * @param {number} workerCount - Il numero di worker disponibili.
     * @param {Array} [workPlan=[]] - Il piano di lavoro iniziale (opzionale).
     * @returns {Array} - Il piano di lavoro aggiornato e suddiviso per worker.
     */
    updateMoleculeWorkerAssignments(molecules, workerCount, workPlan = []) {
        const currentMoleculeIds = new Set(molecules.map(mol => mol.id));

        // 1. Rimuovi le molecole non più presenti.
        const removedMolecules = [...this.lastSeenMolecules].filter(id => !currentMoleculeIds.has(id));
        for (const molId of removedMolecules) {
            const workerId = this.moleculeWorkerAssignments.get(molId);
            if (workerId !== undefined) {
                const moleculeGroup = this.workerMoleculeGroups.get(workerId);
                if (moleculeGroup) {
                    moleculeGroup.delete(molId);
                }
                this.moleculeWorkerAssignments.delete(molId);

                // Aggiorna il carico del worker.
                this.workerLoad.set(workerId, (this.workerLoad.get(workerId) || 0) -1 );
            }
        }

      // 2. Aggiungi le nuove molecole.  Trova il worker "ideale" per ogni nuova molecola.
        const newMolecules = [...currentMoleculeIds].filter(id => !this.lastSeenMolecules.has(id));
        for (const molId of newMolecules) {
            const idealWorkerId = this.findIdealWorkerForMolecule(molId, molecules, workerCount);

            if (!this.workerMoleculeGroups.has(idealWorkerId)) {
                this.workerMoleculeGroups.set(idealWorkerId, new Set());
            }
            this.workerMoleculeGroups.get(idealWorkerId).add(molId);
            this.moleculeWorkerAssignments.set(molId, idealWorkerId);

            // Aggiorna il carico.
           this.workerLoad.set(idealWorkerId, (this.workerLoad.get(idealWorkerId) || 0) + 1);
        }

        // 3. Aggiorna lastSeenMolecules.
        this.lastSeenMolecules = currentMoleculeIds;

        // 4.  Ribilancia il carico se necessario (dopo aver gestito nuove/rimosse).
        this.rebalanceWorkerLoadsIfNeeded(workerCount);

        // 5.  Suddividi il workPlan per worker.
        const workerPartitions = new Map();
        for (let i = 0; i < workerCount; i++) {
            workerPartitions.set(i, []);
        }

       // Ordina il piano di lavoro per workerId e priorità.
       workPlan.sort((a, b) => {
            if (a.workerId !== b.workerId) { return a.workerId - b.workerId; }
            return b.priority - a.priority;  // Priorità più alta per prima.
        });


        for (const interaction of workPlan) {
            const targetPartition = workerPartitions.get(interaction.workerId);
            if (targetPartition) { // Controllo di sicurezza.
                targetPartition.push(interaction);
            }
        }

        // Restituisci il piano suddiviso, includendo gli indici delle molecole richieste.
        return Array.from(workerPartitions.entries()).map(([workerId, interactions]) => ({
            workerId,
            interactions,
            moleculeIndices: this.getRequiredMoleculeIndices(interactions)
        }));
    }

    /**
     * Esporta le statistiche di stabilità delle assegnazioni
     * Utile per monitorare l'efficacia della persistenza
     */
    getAssignmentStats() {
        return {
            uniqueMolecules: this.moleculeWorkerAssignments.size,
            workerGroups: Array.from(this.workerMoleculeGroups.entries()).map(([workerId, molecules]) => ({
                workerId,
                moleculeCount: molecules.size
            })),
            pairAssignments: this.pairToWorkerMap.size,
            workerLoads: Array.from(this.workerLoad.entries())
        };
    }

    /**
     * Serializza lo stato delle assegnazioni per preservarlo tra sessioni
     * Può essere utile per sistemi che si riavviano frequentemente
     */
    serializeAssignmentState() {
        return JSON.stringify({
            moleculeWorkerAssignments: Array.from(this.moleculeWorkerAssignments.entries()),
            pairToWorkerMap: Array.from(this.pairToWorkerMap.entries()),
            lastSeenMolecules: Array.from(this.lastSeenMolecules)
        });
    }

    /**
     * Carica lo stato delle assegnazioni da una precedente serializzazione
     */
    loadAssignmentState(serializedState) {
        try {
            const state = JSON.parse(serializedState);

            if (state.moleculeWorkerAssignments) {
                this.moleculeWorkerAssignments = new Map(state.moleculeWorkerAssignments);
            }

            if (state.pairToWorkerMap) {
                this.pairToWorkerMap = new Map(state.pairToWorkerMap);
            }

            if (state.lastSeenMolecules) {
                this.lastSeenMolecules = new Set(state.lastSeenMolecules);
            }

            // Ricostruisci workerMoleculeGroups dalle assegnazioni
            this.workerMoleculeGroups.clear();
            for (const [molId, workerId] of this.moleculeWorkerAssignments.entries()) {
                if (!this.workerMoleculeGroups.has(workerId)) {
                    this.workerMoleculeGroups.set(workerId, new Set());
                }
                this.workerMoleculeGroups.get(workerId).add(molId);
            }

            // Aggiorna workerLoad
            this.workerLoad.clear();
            for (const [workerId, molecules] of this.workerMoleculeGroups.entries()) {
                this.workerLoad.set(workerId, molecules.size);
            }

            return true;
        } catch (error) {
            console.error('Errore durante il caricamento dello stato:', error);
            return false;
        }
    }

  /**
   * Trova il worker ideale per una nuova molecola.  Considera:
   * 1. Affinità con le molecole già assegnate a ciascun worker.
   * 2. Bilanciamento del carico tra i worker.
   *
   * @param {string} molId - ID della nuova molecola.
   * @param {Array} allMolecules - Array di tutte le molecole (per calcolare la posizione).
   * @param {number} workerCount - Numero di worker disponibili.
   * @returns {number} - ID del worker ideale.
   */
    findIdealWorkerForMolecule(molId, allMolecules, workerCount) {
       // Crea una mappa per trovare rapidamente le molecole.
        const moleculeMap = new Map();
        allMolecules.forEach((mol, index) => moleculeMap.set(mol.id, index));

       // Calcola un punteggio di "affinità" per ogni worker.
        const workerAffinityScores = new Map();
        for (let i = 0; i < workerCount; i++) {
            workerAffinityScores.set(i, 0);
        }

        const molIndex = moleculeMap.get(molId);
        const molecule = allMolecules[molIndex];

        // Calcola l'affinità con le molecole già assegnate.
       for (const [otherMolId, workerId] of this.moleculeWorkerAssignments.entries()) {
          if (workerId >= workerCount) continue;  // Ignora worker non più validi.

          const otherMolIndex = moleculeMap.get(otherMolId);
          if (otherMolIndex === undefined) continue; // Ignora molecole non più presenti.

          const otherMol = allMolecules[otherMolIndex];
          const affinity = this.calculateMoleculeAffinity(molecule, otherMol);

          // Incrementa il punteggio di affinità.
          const currentScore = workerAffinityScores.get(workerId) || 0;
          workerAffinityScores.set(workerId, currentScore + affinity);
        }

        // Combina l'affinità con il bilanciamento del carico.
        let bestWorkerId = 0;
        let bestScore = -Infinity;

        for (let i = 0; i < workerCount; i++) {
            const affinityScore = workerAffinityScores.get(i) || 0;
            const currentLoad = this.workerLoad.get(i) || 0;

            // Formula:  affinità - (carico / carico_medio_ideale).
            const balancedScore = affinityScore - (currentLoad / (allMolecules.length / workerCount));

            if (balancedScore > bestScore) {
                bestScore = balancedScore;
                bestWorkerId = i;
            }
        }

        return bestWorkerId;
    }

    /**
     * Calcola un punteggio di affinità tra due molecole
     * Valori più alti indicano maggiore probabilità di interazione significativa
     */
    calculateMoleculeAffinity(mol1, mol2) {
        // Verifica se questa coppia è già nella cache delle interazioni
        const pairId = this.getPairId(mol1.id, mol2.id);
        if (this.interactionCache.has(pairId)) {
            return 10.0; // Alta affinità per interazioni già calcolate in passato
        }

        // Calcola distanza approssimativa
        const pos1 = mol1._physicsPosition || mol1.position;
        const pos2 = mol2._physicsPosition || mol2.position;

        let distanceSquared = 0;
        for (let i = 0; i < 3; i++) {
            const diff = pos2[i] - pos1[i];
            distanceSquared += diff * diff;
        }

        // Affinità inversa alla distanza e proporzionale alle masse
        const combinedMass = (mol1.mass || 1) + (mol2.mass || 1);
        return combinedMass / (distanceSquared + 1);
    }

    /**
     * Ribilancia il carico se c'è uno sbilanciamento significativo
     * Viene eseguito solo occasionalmente per evitare continui spostamenti
     */
    rebalanceWorkerLoadsIfNeeded(workerCount) {
        // Esegui ribilanciamento solo occasionalmente (ad esempio, ogni X chiamate)
        // o quando c'è uno sbilanciamento significativo
        const shouldRebalance = this.detectSignificantImbalance(workerCount);
        if (!shouldRebalance) return;

        // Trova worker più e meno caricati
        let maxLoad = -Infinity;
        let minLoad = Infinity;
        let mostLoadedWorker = 0;
        let leastLoadedWorker = 0;

        for (let i = 0; i < workerCount; i++) {
            const load = this.workerLoad.get(i) || 0;
            if (load > maxLoad) {
                maxLoad = load;
                mostLoadedWorker = i;
            }
            if (load < minLoad) {
                minLoad = load;
                leastLoadedWorker = i;
            }
        }

        // Se lo sbilanciamento è significativo, sposta alcune molecole
        if (maxLoad - minLoad > 3) {
            const overloadedWorkerMolecules = this.workerMoleculeGroups.get(mostLoadedWorker);
            if (!overloadedWorkerMolecules || overloadedWorkerMolecules.size <= 1) return;

            // Trova le molecole con minor affinità nel gruppo più caricato
            const moleculesToMove = this.findLeastAffineMolecules(
                overloadedWorkerMolecules,
                Math.floor((maxLoad - minLoad) / 2)
            );

            // Sposta le molecole selezionate
            this.moveMoleculesToWorker(moleculesToMove, mostLoadedWorker, leastLoadedWorker);
        }
    }

    /**
     * Rileva se c'è uno sbilanciamento significativo tra i worker
     */
    detectSignificantImbalance(workerCount) {
        if (workerCount <= 1) return false;

        let maxLoad = -Infinity;
        let minLoad = Infinity;
        let totalLoad = 0;

        for (let i = 0; i < workerCount; i++) {
            const load = this.workerLoad.get(i) || 0;
            maxLoad = Math.max(maxLoad, load);
            minLoad = Math.min(minLoad, load);
            totalLoad += load;
        }

        const avgLoad = totalLoad / workerCount;

        // Considera sbilanciato se:
        // 1. La differenza tra max e min è più del 50% del carico medio, e
        // 2. La differenza assoluta è almeno 3 molecole
        return (maxLoad - minLoad) > Math.max(avgLoad * 0.5, 3);
    }

    /**
     * Trova le molecole con minor affinità all'interno di un gruppo
     */
    findLeastAffineMolecules(moleculeGroup, count) {
        // Implementazione semplificata: seleziona elementi casuali
        // In una implementazione reale, calcolerebbe l'affinità tra ciascuna molecola
        // e tutte le altre nel gruppo, selezionando quelle con affinità più bassa
        const moleculeArray = [...moleculeGroup];
        return moleculeArray.slice(0, Math.min(count, moleculeArray.length));
    }

    /**
     * Sposta molecole da un worker all'altro, aggiornando tutte le mappe pertinenti
     */
    moveMoleculesToWorker(moleculeIds, sourceWorkerId, targetWorkerId) {
        const sourceGroup = this.workerMoleculeGroups.get(sourceWorkerId);
        if (!sourceGroup) return;

        if (!this.workerMoleculeGroups.has(targetWorkerId)) {
            this.workerMoleculeGroups.set(targetWorkerId, new Set());
        }
        const targetGroup = this.workerMoleculeGroups.get(targetWorkerId);

        for (const molId of moleculeIds) {
            if (sourceGroup.has(molId)) {
                // Aggiorna le mappe
                sourceGroup.delete(molId);
                targetGroup.add(molId);
                this.moleculeWorkerAssignments.set(molId, targetWorkerId);

                // Aggiorna i conteggi del carico
                this.workerLoad.set(sourceWorkerId, (this.workerLoad.get(sourceWorkerId) || 0) - 1);
                this.workerLoad.set(targetWorkerId, (this.workerLoad.get(targetWorkerId) || 0) + 1);

                // Aggiorna anche le assegnazioni delle coppie
                this.updatePairAssignmentsForMolecule(molId, targetWorkerId);
            }
        }
    }

    /**
     * Aggiorna le assegnazioni delle coppie per una molecola spostata
     */
    updatePairAssignmentsForMolecule(movedMolId, newWorkerId) {
        // Trova tutte le coppie che coinvolgono questa molecola
        for (const [pairId, workerId] of this.pairToWorkerMap.entries()) {
            const [mol1Id, mol2Id] = pairId.split(':');
            if (mol1Id === movedMolId || mol2Id === movedMolId) {
                this.pairToWorkerMap.set(pairId, newWorkerId);
            }
        }
    }

  /**
     * Costruisce il piano di lavoro (workPlan) per le interazioni tra molecole.
     *
     * @param {Array} molecules - L'elenco delle molecole.
     * @param {number} workerCount - Il numero di worker disponibili.
     * @returns {Array} - Il piano di lavoro, ordinato per priorità e worker.
     */
    buildInteractionWorkPlan(molecules, workerCount = 1) {
        this.resetWorkSession(); // Pulisce lo stato precedente, *mantenendo* le assegnazioni.

        let workPlan = [];
        const moleculeCount = molecules.length;

        // Costruisci il piano di lavoro (matrice triangolare superiore).
        for (let i = 0; i < moleculeCount; i++) {
            for (let j = i + 1; j < moleculeCount; j++) {
                const mol1 = molecules[i];
                const mol2 = molecules[j];

                if (this.isPairProcessed(mol1.id, mol2.id)) {
                    this.stats.calculationsAvoided++;
                    continue;
                }

                const workerId = this.determineWorkerForInteraction(mol1.id, mol2.id, workerCount);
                workPlan.push({
                    mol1Index: i,
                    mol2Index: j,
                    mol1Id: mol1.id,
                    mol2Id: mol2.id,
                    workerId: workerId,
                    priority: this.estimateInteractionPriority(mol1, mol2)
                });

                this.stats.totalCalculations++;
            }
        }


        // Non è necessario chiamare updateMoleculeWorkerAssignments qui, buildInteractionWorkPlan è chiamato dentro generateDistributionPlanWithPriorities, che già fa l'aggiornamento degli worker
        return workPlan;
    }

    /**
     * Determina quale worker dovrebbe gestire un'interazione specifica
     * Cerca di mantenere consistenza con le assegnazioni precedenti
     */
    determineWorkerForInteraction(mol1Id, mol2Id, workerCount) {
        const pairId = this.getPairId(mol1Id, mol2Id);

        // 1. Verifica se questa coppia ha già un'assegnazione persistente
        if (this.pairToWorkerMap.has(pairId)) {
            const existingWorkerId = this.pairToWorkerMap.get(pairId);
            // Verifica che il worker esista ancora
            if (existingWorkerId < workerCount) {
                return existingWorkerId;
            }
        }

        // 2. Altrimenti, scegli in base all'appartenenza delle molecole
        const worker1 = this.moleculeWorkerAssignments.get(mol1Id);
        const worker2 = this.moleculeWorkerAssignments.get(mol2Id);

        // Se entrambe le molecole sono già assegnate allo stesso worker, usa quello
        if (worker1 !== undefined && worker1 === worker2 && worker1 < workerCount) {
            this.pairToWorkerMap.set(pairId, worker1);
            return worker1;
        }

        // 3. Se le molecole appartengono a worker diversi, scegli in base al carico
        if (worker1 !== undefined && worker2 !== undefined) {
            // Scegli il worker con carico minore
            const load1 = this.workerLoad.get(worker1) || 0;
            const load2 = this.workerLoad.get(worker2) || 0;

            const selectedWorker = load1 <= load2 ? worker1 : worker2;
            this.pairToWorkerMap.set(pairId, selectedWorker);
            return selectedWorker;
        }

        // 4. Se almeno una molecola ha un'assegnazione, usa quella
        if (worker1 !== undefined && worker1 < workerCount) {
            this.pairToWorkerMap.set(pairId, worker1);
            return worker1;
        }

        if (worker2 !== undefined && worker2 < workerCount) {
            this.pairToWorkerMap.set(pairId, worker2);
            return worker2;
        }

        // 5. Se nessuna delle due molecole ha un'assegnazione, scegli il worker meno carico
        let leastLoadedWorker = 0;
        let minLoad = Infinity;

        for (let i = 0; i < workerCount; i++) {
            const load = this.workerLoad.get(i) || 0;
            if (load < minLoad) {
                minLoad = load;
                leastLoadedWorker = i;
            }
        }

        this.pairToWorkerMap.set(pairId, leastLoadedWorker);
        return leastLoadedWorker;
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
     * ma mantiene la cache delle interazioni e le assegnazioni worker-molecola
     */
    resetWorkSession() {
        this.calculatedPairs.clear();
        this.stats = {
            calculationsAvoided: 0,
            totalCalculations: 0,
            cacheHits: 0
        };
        // Non resetta moleculeWorkerAssignments o workerMoleculeGroups
        // per mantenere la persistenza delle assegnazioni
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
     * Genera il piano di distribuzione del lavoro, considerando le priorità e la consistenza.
     *
     * @param {Array} molecules - Tutte le molecole.
     * @param {number} workerCount - Numero di worker.
     * @param {Map} significantPairs - Coppie di molecole con interazioni significative.
     * @returns {Array} - Il piano di distribuzione, suddiviso per worker.
     */
    generateDistributionPlanWithPriorities(molecules, workerCount, significantPairs) {
        this.resetWorkSession(); // Pulisce lo stato, mantenendo le assegnazioni.

        // 1. Costruisci il piano di lavoro *completo*, includendo *tutte* le interazioni potenziali.
        //    Questo passaggio è importante per avere una visione globale delle interazioni
        //    e per poter applicare correttamente le priorità.
        let workPlan = this.buildInteractionWorkPlan(molecules, workerCount);


        // 2.  Applica le priorità.  Le interazioni significative *sostituiscono* quelle
        //     a priorità inferiore, se necessario.  Questo garantisce che le interazioni
        //     importanti siano sempre processate.
        if (significantPairs && significantPairs.size > 0) {
            const moleculeIndices = new Map(); // Mappa ID molecola -> indice.
            molecules.forEach((mol, index) => moleculeIndices.set(mol.id, index));

            for (const [pairId, pairInfo] of significantPairs.entries()) {
              const mol1Index = moleculeIndices.get(pairInfo.mol1Id);
              const mol2Index = moleculeIndices.get(pairInfo.mol2Id);

              if (mol1Index === undefined || mol2Index === undefined) continue;

                // Se la coppia è già nel workPlan, aggiorna la priorità.
                // Altrimenti, aggiungila.
                const existingInteractionIndex = workPlan.findIndex(item =>
                    item.mol1Id === pairInfo.mol1Id && item.mol2Id === pairInfo.mol2Id
                );

                if (existingInteractionIndex !== -1) {
                    workPlan[existingInteractionIndex].priority = pairInfo.priority || 10.0;
                } else {
                    // Determina il worker per questa interazione.
                    const workerId = this.determineWorkerForInteraction(
                        pairInfo.mol1Id,
                        pairInfo.mol2Id,
                        workerCount
                    );

                    workPlan.push({
                        mol1Index,
                        mol2Index,
                        mol1Id: pairInfo.mol1Id,
                        mol2Id: pairInfo.mol2Id,
                        workerId,
                        priority: pairInfo.priority || 10.0
                    });
                }
            }
        }


      // 3.  *Dopo* aver applicato le priorità, ordina il piano di lavoro.
        workPlan.sort((a, b) => {
            if (a.workerId !== b.workerId) {
                return a.workerId - b.workerId;
            }
            return b.priority - a.priority; // Priorità più alta per prima.
        });

        // 4. Inizializza workerLoad se necessario
        if(this.workerLoad.size !== workerCount){
            this.workerLoad.clear();
            for(let i = 0; i < workerCount; i++){
                this.workerLoad.set(i, 0);
            }
        }

        // 5. Aggiorna le assegnazioni e il carico.  Questo passaggio gestisce anche
        //    le molecole nuove e rimosse.
        const distributionPlan = this.updateMoleculeWorkerAssignments(molecules, workerCount, workPlan);

        // 6.  Aggiorna lastSeenMolecules (dopo aver gestito nuove/rimosse).
        this.lastSeenMolecules = new Set(molecules.map(m => m.id));

        return distributionPlan; // Restituisci il piano *suddiviso* per worker.
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
     * Distribuisce il lavoro ai sub-worker in modo ottimizzato.
     *
     * 1. Identifica le interazioni significative (usando la cache e le regole).
     * 2. Genera un piano di distribuzione che assegna le interazioni ai worker.
     *    - Tiene conto della localizzazione delle molecole (per minimizzare i trasferimenti).
     *    - Bilancia il carico tra i worker.
     *    - Dà priorità alle interazioni significative.
     * 3. Invia i task ai worker.
     * 4. Attende il completamento.
     * 5. Aggiorna lo stato (nel metodo `mergeProcessedChunkWithConflictResolution`).
     */
    async distributeWorkToSubWorkers() {
      if (this.subWorkers.length === 0 || this.isPaused) return;

      await this.waitForAvailableWorkers();

      // 1. Identifica le interazioni significative.
      const significantPairs = this.identifySignificantInteractions();

      // 2. Genera il piano di distribuzione.
      const distributionPlan = this.interactionManager.generateDistributionPlanWithPriorities(
        this.molecules,
        this.subWorkers.length,
        significantPairs
      );

      // 3. Distribuisci il lavoro ai worker.
      const processingPromises = [];
      for (let i = 0; i < Math.min(this.subWorkers.length, distributionPlan.length); i++) {
        const partition = distributionPlan[i];
        const subWorker = this.subWorkers[i];

        if (!partition || partition.interactions.length === 0) continue;

        subWorker.busy = true; // Segna il worker come occupato.

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
     * Assegna un task (una partizione del piano di lavoro) a un sub-worker.
     *
     * @param {Object} subWorker - L'oggetto che rappresenta il sub-worker.
     * @param {Object} partition - La partizione del piano di lavoro da assegnare.
     * @returns {Promise} - Una Promise che si risolve quando il worker ha completato il task.
     */
    assignWorkToSubWorker(subWorker, partition) {
      return new Promise(resolve => {
        const workerMessageHandler = (event) => {
          if (event.data.type === 'chunk_processed') {
            subWorker.worker.removeEventListener('message', workerMessageHandler);
            resolve(); // Risolvi la Promise quando il worker ha finito.
          }
        };

        subWorker.worker.addEventListener('message', workerMessageHandler);

        // Prepara i dati da inviare al worker.
        const requiredMolecules = partition.moleculeIndices.map(index => this.molecules[index]);
        const serializedMolecules = this.serializeMolecules(requiredMolecules);
        const cachedRelationships = this.moleculeCache.getRelationshipsForMolecules(requiredMolecules);

        // Invia il messaggio al worker.
        subWorker.worker.postMessage({
          type: 'process_optimized_chunk', // Usa il tipo di messaggio ottimizzato.
          molecules: serializedMolecules,
          interactionPlan: partition.interactions, // Invia il piano delle interazioni.
          moleculeIndices: partition.moleculeIndices, // Invia gli *indici* delle molecole.
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

      // Verifica reazioni *solo* se abbastanza vicine.
        if (distance < this.calculateReactionDistance(mol1, mol2)) {
          if (this.shouldReact(mol1, mol2)) {

            const mol1OriginalIndex = this.molecules.findIndex(m => m.id === mol1.id);
            const mol2OriginalIndex = this.molecules.findIndex(m => m.id === mol2.id);

            if (mol1OriginalIndex >= 0 && mol2OriginalIndex >= 0) {

              const products = this.processReaction(mol1, mol2);

              if (products.length > 0) {
                newMolecules.push(...products);
                removedIndices.add(mol1OriginalIndex); // Usa l'indice *originale*.
                removedIndices.add(mol2OriginalIndex);
                this.reactionCount++;

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
     * Implementazione *migliorata* del raggruppamento di molecole in cluster.
     *
     * 1. Utilizza la cache delle relazioni per trovare rapidamente le molecole vicine.
     * 2. Per le molecole senza relazioni in cache, esegue un calcolo *limitato* delle distanze.
     * 3. Applica un bilanciamento dei cluster *ottimizzato*.
     */
    groupMoleculesByClusters() {
      if (this.molecules.length <= 1) {
        return [this.molecules];
      }

      const clusters = [];
      const visited = new Set();
      const proximityThreshold = this.size * 0.15;

      // Funzione di supporto per trovare un cluster, a partire da una molecola.
      const findCluster = (startMol) => {
        const cluster = [startMol];
        const queue = [startMol];
        visited.add(startMol.id);

        while (queue.length > 0) {
          const current = queue.shift();

          // 1. Utilizza le relazioni in cache.
          const relationships = this.moleculeCache.getRelationshipsForMolecule(current);
          for (const rel of relationships) {
            const neighborMol = this.molecules.find(m => m.id === rel.otherId);
            if (!neighborMol || visited.has(neighborMol.id)) continue;

            // Usa la distanza *dalla cache*.
            if (rel.distance < proximityThreshold) {
                cluster.push(neighborMol);
                queue.push(neighborMol);
                visited.add(neighborMol.id);
            }
          }

          // 2. Per le molecole *senza* relazioni in cache, calcola le distanze,
          //    ma solo con le molecole *non ancora visitate*.
          const unvisitedMolecules = this.molecules.filter(m =>
            !visited.has(m.id) && m.id !== current.id
          );

          if (unvisitedMolecules.length > 0) {
              const miniPlan = [];
              for(const mol of unvisitedMolecules){
                  if (this.interactionManager.isPairProcessed(current.id, mol.id)) {
                      continue;
                  }

                  const distance = this.calculateDistance(current, mol);
                  if(distance < proximityThreshold){
                      miniPlan.push({
                          mol: mol,
                          distance: distance
                      })

                      this.moleculeCache.createRelationship(
                          current, mol, distance, performance.now()
                      );
                      this.interactionManager.markPairProcessed(current.id, mol.id);
                  }
              }

              for(const item of miniPlan){
                  if(!visited.has(item.mol.id)){
                      cluster.push(item.mol);
                      queue.push(item.mol);
                      visited.add(item.mol.id);
                  }
              }
          }
        }

        return cluster;
      };

      // Forma i cluster, a partire da ciascuna molecola non ancora visitata.
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
       // Usa piano ottimizzato per aggiornare molecole e ottieni gli ID processati direttamente
        const directlyProcessedIds = processInteractionPlan(moleculeObjs, adaptedPlan, removedIndices, newMolecules);


        // Prepara risultati
        const results = prepareProcessedResults(
            moleculeObjs,
            removedIndices,
            newMolecules,
            reactionCount,
            directlyProcessedIds // Passa gli ID processati direttamente
        );

        // Invia risultati
        sendProcessedChunkResult(results, interactionPlan.length);

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
    rules.setConstant('time_scale', message.timeScale || 1);

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
        sendProcessedChunkResult(results, molecules.length);

    } catch (error) {
        reportProcessingError(error);
    }
}

/**
 * Invia risultato vuoto per chunk
 */
function sendEmptyChunkResult() {
    sendProcessedChunkResult( {
        moleculeUpdates: [],
        newMolecules: [],
        reactionCount: 0,
        updatedRelationships: {}
    }, 0)
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
 * Processa le interazioni e traccia le molecole elaborate direttamente.
 *
 * @param {Array} molecules - Le molecole da processare.
 * @param {Array} plan - Il piano delle interazioni.
 * @param {Set} removedIndices - Set per tenere traccia degli ID delle molecole rimosse.
 * @param {Array} newMolecules - Array per accumulare le nuove molecole.
 * @returns {Set} - Un set di ID delle molecole elaborate *direttamente*.
 */
function processInteractionPlan(molecules, plan, removedIndices, newMolecules) {
    const timeScale = simulation.rules.getConstant('time_scale');
    const damping = simulation.rules.getConstant('damping');
    const now = performance.now();

    const directlyProcessedIds = new Set(); // Traccia le molecole elaborate *direttamente*.

    // Pre-calcola le forze (inizializzate a zero).
    const forces = new Map();
    molecules.forEach(mol => forces.set(mol.id, [0, 0, 0]));

    for (const interaction of plan) {
        const mol1 = interaction.mol1;
        const mol2 = interaction.mol2;

        if (!mol1 || !mol2 || removedIndices.has(mol1.id) || removedIndices.has(mol2.id)) {
            continue;
        }

        directlyProcessedIds.add(mol1.id); // Aggiungi al set.
        directlyProcessedIds.add(mol2.id);

        const currentDistance = simulation.calculateDistance(mol1, mol2);
        simulation.moleculeCache.createRelationship(mol1, mol2, currentDistance, now);

        if (currentDistance <= simulation.calculateMaxRelationshipDistance(mol1, mol2)) {
            const [force1, force2] = simulation.calculateForces(mol1, mol2, currentDistance, now);

            const mol1Force = forces.get(mol1.id);
            const mol2Force = forces.get(mol2.id);
            if (mol1Force && mol2Force) {
                for (let k = 0; k < 3; k++) {
                    mol1Force[k] += force1[k] * timeScale;
                    mol2Force[k] += force2[k] * timeScale;
                }
            }

            // Gestisci le reazioni *all'interno* del ciclo principale.
            processReactionIfNeeded(mol1, mol2, currentDistance, removedIndices, newMolecules);
        }
    }

    applyForcesAndMove(molecules, forces, timeScale, damping);

    return directlyProcessedIds; // Restituisci il set.
}

/**
 * Prepara i risultati dell'elaborazione, includendo gli ID delle molecole elaborate direttamente.
 *
 * @param {Array} molecules - Le molecole processate.
 * @param {Set} removedIds - Gli ID delle molecole rimosse.
 * @param {Array} newMolecules - Le nuove molecole create.
 * @param {number} initialReactionCount - Il conteggio delle reazioni prima del processamento.
 * @param {Set} directlyProcessedIds - Set con ID delle molecole elaborate direttamente.
 * @returns {Object} - I risultati, serializzati e pronti per essere inviati.
 */
function prepareProcessedResults(molecules, removedIds, newMolecules, initialReactionCount, directlyProcessedIds) {

    return {
        moleculeUpdates: molecules
          .filter(mol => !removedIds.has(mol.id)) // Usa gli ID, non gli indici.
          .map(mol => MoleculeSerializer.serialize(mol)),
        newMolecules: newMolecules.map(mol => MoleculeSerializer.serialize(mol)),
        reactionCount: simulation.reactionCount - initialReactionCount,
        updatedRelationships: simulation.moleculeCache.getRelationshipsForMolecules(
          molecules.filter(mol => !removedIds.has(mol.id))
        ),
        processedMoleculeIds: [...directlyProcessedIds], // Converti in array.
        timestamp: performance.now() //Aggiunto timestamp
    };
}

/**
 * Invia risultati al worker principale con informazioni migliorate
 */
function sendProcessedChunkResult(results, processedCount) {
    postMessage({
        type: 'chunk_processed',
        workerId: workerId,
        processedCount: processedCount,
        processedMoleculeIds: results.processedMoleculeIds || [],
        timestamp: results.timestamp || performance.now(),
        results: results
    });
}

/**
 * Processa una possibile reazione tra molecole
 * @param {Object} mol1 - Prima molecola
 * @param {Object} mol2 - Seconda molecola
 * @param {number} distance - Distanza tra molecole
 * @param {Set} removedIndices - Set per indici rimossi  -- MODIFICATO, ora contiene ID
 * @param {Array} newMolecules - Array per nuove molecole
 */
function processReactionIfNeeded(mol1, mol2, distance, removedIndices, newMolecules) {
    if (distance < simulation.calculateReactionDistance(mol1, mol2)) {
        if (simulation.shouldReact(mol1, mol2)) {
            const products = simulation.processReaction(mol1, mol2);
            if (products.length > 0) {
                newMolecules.push(...products);

                // Aggiungi a indici da rimuovere -- MODIFICATO, ora usa gli ID
                removedIndices.add(mol1.id);
                removedIndices.add(mol2.id);

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
        simulation.stabilizePositions();
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

        simulation.step_in_progress = false;

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
    rules.setConstant('time_scale', timeScale || 1);

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