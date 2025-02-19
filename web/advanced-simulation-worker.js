// enhanced-simulation-worker.js
import { createCustomRules } from './rules.js';
import { PrimeMolecule } from './molecule.js';


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
        
        // Inizializza molecole con distribuzione interessante
        this.initializeMolecules(moleculeCount);
    }
    
    initializeMolecules(count) {
        // Distribuzione di numeri più interessante
        const numberChoices = [];
        
        // Aggiungi numeri primi puri
        const primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];
        numberChoices.push(...primes);
        
        // Aggiungi qualche composto semplice
        const compounds = [4, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 22, 24, 25, 27, 30];
        numberChoices.push(...compounds);
        
        // Aggiungi qualche numero più grande
        for (let i = 0; i < 10; i++) {
            numberChoices.push(Math.floor(Math.random() * 90) + 50);
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
        
        console.log(`Inizializzate ${this.molecules.length} molecole`);
    }
    
    step() {
        const timeScale = this.rules.getConstant('time_scale');
        this.accumulatedTime += timeScale;
        
        // Temperatura variabile
        this.temperature = 1.0 + 0.2 * Math.sin(performance.now() / 5000);
        
        // Arrays per elaborazione efficiente
        const newMolecules = [];
        const removedIndices = new Set();
        
        // Calcola movimenti e forze
        this.updatePhysics(removedIndices, newMolecules);
        
        // Filtra molecole rimosse
        if (removedIndices.size > 0) {
            this.molecules = this.molecules.filter((_, i) => !removedIndices.has(i));
        }
        
        // Aggiungi nuove molecole
        for (const mol of newMolecules) {
            mol.id = `new-${this.nextMoleculeId++}`;
            this.molecules.push(mol);
        }
        
        // Reset accumulated time
        if (this.accumulatedTime >= 1.0) {
            this.accumulatedTime = 0.0;
        }
        
        // Limita il numero di molecole per prestazioni
        if (this.molecules.length > 300 && false) { //todo: parametrized it
            // Rimuovi molecole in eccesso preferendo quelle più vecchie
            this.molecules = this.molecules.slice(this.molecules.length - 300);
        } else if (this.molecules.length < 50) {
            // Aggiungi molecole se ce ne sono troppo poche
            this.addRandomMolecules(50 - this.molecules.length);
        }
    }
    
    updatePhysics(removedIndices, newMolecules) {
        const moleculeCount = this.molecules.length;
        const timeScale = this.rules.getConstant('time_scale');
        const damping = this.rules.getConstant('damping');
        
        // Pre-calcola forze
        const forces = Array(moleculeCount).fill(0).map(() => [0, 0, 0]);
        
        // Calcola interazioni tra molecole
        for (let i = 0; i < moleculeCount; i++) {
            if (removedIndices.has(i)) continue;
            const mol1 = this.molecules[i];
            
            // Movimento termico casuale
            for (let axis = 0; axis < 3; axis++) {
                mol1.velocity[axis] += (Math.random() - 0.5) * 0.01 * this.temperature;
            }
            
            // Interazioni con altre molecole
            for (let j = i + 1; j < moleculeCount; j++) {
                if (removedIndices.has(j)) continue;
                const mol2 = this.molecules[j];
                
                // Calcola distanza
                const distance = this.calculateDistance(mol1, mol2);
                
                // Applica forze
                const [force1, force2] = this.applyForces(mol1, mol2);
                for (let k = 0; k < 3; k++) {
                    forces[i][k] += force1[k] * timeScale;
                    forces[j][k] += force2[k] * timeScale;
                }
                
                // Verifica reazioni chimiche
                if (distance < 1.0 + (mol1.mass + mol2.mass) * 0.1) {
                    if (this.shouldReact(mol1, mol2)) {
                        const products = this.processReaction(mol1, mol2);
                        if (products.length > 0) {
                            newMolecules.push(...products);
                            removedIndices.add(i);
                            removedIndices.add(j);
                            this.reactionCount++;
                            break;
                        }
                    }
                }
            }
        }
        
        // Applica forze e aggiorna posizioni
        for (let i = 0; i < moleculeCount; i++) {
            if (removedIndices.has(i)) continue;
            const mol = this.molecules[i];
            
            // Aggiorna velocità con forze
            for (let k = 0; k < 3; k++) {
                mol.velocity[k] += forces[i][k] * 0.1;
                mol.position[k] += mol.velocity[k] * timeScale;
                mol.velocity[k] *= damping;
            }
            
            // Verifica limiti
            this.enforceBoundaries(mol);
        }
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
        // Calcola vettore direzione
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
        
        // Normalizza direzione
        const dirNorm = direction.map(d => d / distance);
        
        // Calcola forze nette
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
        
        // Limita forza massima
        const maxForce = this.rules.getConstant('max_force');
        for (let i = 0; i < 3; i++) {
            force1[i] = Math.max(-maxForce, Math.min(maxForce, force1[i]));
            force2[i] = Math.max(-maxForce, Math.min(maxForce, force2[i]));
        }
        
        return [force1, force2];
    }
    
    shouldReact(mol1, mol2) {
        // Le reazioni avvengono solo se l'energia è sufficiente
        const relativeSpeed = Math.sqrt(
            Math.pow(mol1.velocity[0] - mol2.velocity[0], 2) +
            Math.pow(mol1.velocity[1] - mol2.velocity[1], 2) +
            Math.pow(mol1.velocity[2] - mol2.velocity[2], 2)
        );
        
        // Maggiore è la temperatura, più probabili sono le reazioni
        const baseReactionProb = 0.1 * this.temperature * relativeSpeed;
        
        // Molecole con numeri più piccoli reagiscono più facilmente
        const sizeFactor = 1.0 / (1.0 + Math.log(mol1.number + mol2.number));
        
        return Math.random() < baseReactionProb * sizeFactor;
    }
    
    processReaction(mol1, mol2) {
        // Cerca regole di reazione applicabili
        for (const rule of this.rules.reaction_rules) {
            if (rule.condition(mol1.prime_factors, mol2.prime_factors)) {
                if (Math.random() < rule.probability * this.temperature) {
                    // Applica la regola di reazione
                    const products = rule.effect(mol1, mol2);
                    
                    // Imposta tempo di reazione per effetto visivo
                    const now = performance.now();
                    products.forEach(p => p.setReactionTime(now));
                    
                    return products;
                }
            }
        }
        
        return [];
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
    
    addRandomMolecules(count) {
        for (let i = 0; i < count; i++) {
            // Posizione casuale vicino ai bordi
            const side = Math.floor(Math.random() * 6); // 6 facce del cubo
            const pos = [0, 0, 0];
            
            for (let j = 0; j < 3; j++) {
                if (Math.floor(side/2) === j) {
                    // Su questa faccia
                    pos[j] = (side % 2 === 0 ? -1 : 1) * (this.size/2 - 0.2);
                } else {
                    // Posizione casuale su altri assi
                    pos[j] = (Math.random() - 0.5) * this.size;
                }
            }
            
            // Crea una molecola semplice
            const number = Math.floor(Math.random() * 48) + 2;
            const mol = new PrimeMolecule(number, pos);
            mol.id = `spawn-${this.nextMoleculeId++}`;
            
            // Velocità verso il centro
            const dirToCenter = pos.map(p => -p);
            const norm = Math.sqrt(dirToCenter.reduce((s, v) => s + v*v, 0));
            mol.velocity = dirToCenter.map(v => (v/norm) * (0.1 + Math.random() * 0.2));
            
            this.molecules.push(mol);
        }
    }
}

///
///
///

let simulation;
let previousMoleculeIds = new Set();

// AGGIORNAMENTO HANDLER ONMESSAGE DEL WORKER
// Sostituire l'attuale handler onmessage nel worker con questa versione

onmessage = function(event) {
    console.log("Worker: ricevuto messaggio", event.data.type);
    
    // Verifica se simulazione è inizializzata per la maggior parte dei comandi
    if (event.data.type !== 'init' && !simulation) {
        console.error("Worker: simulazione non inizializzata per "+event.data.type);
        postMessage({
            type: 'error',
            message: 'Simulazione non inizializzata per '+event.data.type
        });
        return;
    }

    try {
        switch (event.data.type) {
            case 'init':
                const { size, moleculeCount, maxNumber, timeScale } = event.data;
                const rules = createCustomRules();
                rules.setConstant('time_scale', timeScale || 0.1);

                // Inizializza la simulazione
                simulation = new EnhancedChemistry(rules, size, moleculeCount, maxNumber);
                console.log(`Worker: simulazione inizializzata con ${moleculeCount} molecole, timeScale=${timeScale}`);

                // Passo iniziale
                simulation.step();
                sendUpdate();
                break;

            case 'step':
                // Esegui passo simulazione
                simulation.step();
                sendUpdate();
                break;

            case 'cleanup':
                cleanupResources();
                break;

            case 'set_temperature':
                // Imposta temperatura
                if (typeof event.data.value === 'number') {
                    simulation.temperature = event.data.value;
                    console.log(`Worker: temperatura impostata a ${event.data.value}`);
                    sendUpdate();
                } else {
                    console.error("Worker: valore temperatura non valido", event.data.value);
                }
                break;

            case 'set_timescale':
                // Imposta timeScale
                if (typeof event.data.value === 'number') {
                    const newTimeScale = event.data.value;
                    simulation.rules.setConstant('time_scale', newTimeScale);
                    console.log(`Worker: timeScale impostato a ${newTimeScale}`);
                }
                break;

            case 'add_molecules':
                // Aggiungi molecole
                const count = event.data.count || 20;
                simulation.addRandomMolecules(count);
                console.log(`Worker: aggiunte ${count} nuove molecole`);
                sendUpdate();
                break;

            case 'set_visualization':
                // Al momento il worker non necessita di fare nulla per questo messaggio
                // serve solo per la parte grafica nel thread principale
                console.log(`Worker: modalità visualizzazione ${event.data.mode} (ignorata nel worker)`);
                break;

            case 'pause':
                // Il worker non ha bisogno di fare nulla per la pausa,
                // è il thread principale che decide se inviare o meno 'step'
                console.log("Worker: ricevuto comando pausa (ignorato nel worker)");
                break;

            default:
                console.warn(`Worker: messaggio sconosciuto '${event.data.type}'`);
        }
    } catch (error) {
        console.error(`Worker: errore durante l'elaborazione del messaggio '${event.data.type}'`, error);
        postMessage({
            type: 'error',
            message: `Errore nel worker: ${error.message}`
        });
    }
};

function sendUpdate() {
    try {
        // Prepara dati per l'invio
        const moleculeData = getOptimizedMoleculeData();

        // Debug
        console.log(`Worker: invio update con ${moleculeData.molecules.length} molecole, ${moleculeData.removedIds.length} rimosse`);

        // Invia messaggio
        postMessage({
            type: 'update',
            molecules: moleculeData,
            temperature: simulation.temperature,
            reactionCount: simulation.reactionCount
        });
    } catch (error) {
        console.error("Worker: errore durante l'invio dell'update", error);
        postMessage({
            type: 'error',
            message: `Errore nell'invio update: ${error.message}`
        });
    }
}

// MIGLIORAMENTO DELLA CLASSE ENHANCEDCHEMISTRY
// Aggiungi questi metodi alla classe EnhancedChemistry

// Metodo per impostare temperatura
EnhancedChemistry.prototype.setTemperature = function(value) {
    this.temperature = value;
};

// Metodo per impostare timeScale
EnhancedChemistry.prototype.setTimeScale = function(value) {
    if (this.rules) {
        this.rules.setConstant('time_scale', value);
    }
};

function getOptimizedMoleculeData() {
    const currentMoleculeIds = new Set();
    const result = [];
    
    for (const mol of simulation.molecules) {
        // Assicurati che ogni molecola abbia un ID
        const id = mol.id || `mol-${Math.random().toString(36).substring(2, 11)}`;
        mol.id = id;
        currentMoleculeIds.add(id);
        
        // Prepara dati per la molecola
        result.push({
            id: id,
            number: mol.number,
            position: [...mol.position],
            velocity: [...mol.velocity],
            prime_factors: {...mol.prime_factors},
            mass: mol.mass,
            charge: mol.charge,
            color: [...mol.color],
            angularVelocity: [...mol.angularVelocity],
            lastReactionTime: mol.lastReactionTime,
        });
    }
    
    // Calcola molecole rimosse
    const removedIds = [...previousMoleculeIds].filter(id => !currentMoleculeIds.has(id));
    
    // Aggiorna set per prossimo frame
    previousMoleculeIds = currentMoleculeIds;
    
    return {
        molecules: result,
        removedIds: removedIds
    };
}

function cleanupResources() {
    simulation = null;
    previousMoleculeIds.clear();
    postMessage({ type: 'cleanup_complete' });
}

///
/// Extension
///

// Estendi la classe PrimeMolecule per supportare il nuovo sistema di parentela
PrimeMolecule.prototype.isRelatedTo = function(otherMolecule, rules) {
    if (!this.id || !otherMolecule.id) return false;
    return rules.areRelated(this.id, otherMolecule.id);
};

// Funzione per gestire il campo quantistico
function updateQuantumFields(simulation) {
    const now = performance.now();
    
    // Filtra tutte le molecole quantistiche
    const quantumFields = simulation.molecules.filter(mol => mol._isQuantumField);
    const fieldsToRemove = [];
    
    for (const field of quantumFields) {
        // Verifica se il campo è scaduto
        if (now - field._creationTime > field._lifetime) {
            fieldsToRemove.push(field);
            continue;
        }
        
        // Calcola l'effetto del campo sulle altre molecole
        for (const mol of simulation.molecules) {
            // Ignora il campo stesso e altri campi
            if (mol._isQuantumField || mol === field) continue;
            
            // Calcola distanza dal campo
            const distance = simulation.calculateDistance(field, mol);
            
            // Il campo ha effetto solo entro un certo raggio
            if (distance < 5.0) {
                // Calcola direzione dal campo alla molecola
                const direction = [
                    mol.position[0] - field.position[0],
                    mol.position[1] - field.position[1],
                    mol.position[2] - field.position[2]
                ];
                
                // Normalizza
                const dirNorm = Math.sqrt(
                    direction[0]**2 + 
                    direction[1]**2 + 
                    direction[2]**2
                );
                
                if (dirNorm > 0.001) {
                    const normalizedDir = direction.map(d => d / dirNorm);
                    
                    // Calcola effetto oscillatorio
                    const oscillation = Math.sin((now - field._creationTime) / 200);
                    const fieldStrength = 0.02 * (1 - distance/5.0) * oscillation;
                    
                    // Applica forza oscillatoria alla molecola
                    for (let i = 0; i < 3; i++) {
                        mol.velocity[i] += normalizedDir[i] * fieldStrength;
                    }
                }
            }
        }
    }
    
    // Rimuovi campi scaduti
    if (fieldsToRemove.length > 0) {
        simulation.molecules = simulation.molecules.filter(mol => 
            !fieldsToRemove.includes(mol));
    }
}

// Fix for updatePhysics method
EnhancedChemistry.prototype.updatePhysics = function(removedIndices, newMolecules) {
    const moleculeCount = this.molecules.length;
    const timeScale = this.rules.getConstant('time_scale');
    const damping = this.rules.getConstant('damping');
    const now = performance.now(); // Add this line to define 'now'
    
    // Pre-calcola forze
    const forces = Array(moleculeCount).fill(0).map(() => [0, 0, 0]);
    
    // Calcola interazioni tra molecole
    for (let i = 0; i < moleculeCount; i++) {
        if (removedIndices.has(i)) continue;
        const mol1 = this.molecules[i];
        
        // Movimento termico casuale
        for (let axis = 0; axis < 3; axis++) {
            mol1.velocity[axis] += (Math.random() - 0.5) * 0.01 * this.temperature;
        }
        
        // Interazioni con altre molecole
        for (let j = i + 1; j < moleculeCount; j++) {
            if (removedIndices.has(j)) continue;
            const mol2 = this.molecules[j];
            
            // Calcola distanza
            const distance = this.calculateDistance(mol1, mol2);
            
            // Applica forze
            const [force1, force2] = this.applyForces(mol1, mol2);
            for (let k = 0; k < 3; k++) {
                forces[i][k] += force1[k] * timeScale;
                forces[j][k] += force2[k] * timeScale;
            }
            
            // Verifica reazioni chimiche
            if (distance < 1.0 + (mol1.mass + mol2.mass) * 0.1) {
                if (this.shouldReact(mol1, mol2)) {
                    const products = this.processReaction(mol1, mol2);
                    if (products.length > 0) {
                        newMolecules.push(...products);
                        removedIndices.add(i);
                        removedIndices.add(j);
                        this.reactionCount++;
                        break;
                    }
                }
            }
        }
    }
    
    // Applica forze e aggiorna posizioni
    for (let i = 0; i < moleculeCount; i++) {
        if (removedIndices.has(i)) continue;
        const mol = this.molecules[i];
        
        // Aggiorna velocità con forze
        for (let k = 0; k < 3; k++) {
            mol.velocity[k] += forces[i][k] * 0.1;
            mol.position[k] += mol.velocity[k] * timeScale;
            mol.velocity[k] *= damping;
        }
        
        // Verifica limiti
        this.enforceBoundaries(mol);
    }
}

// Fix for applyForces method
EnhancedChemistry.prototype.applyForces = function(mol1, mol2) {
    const now = performance.now(); // Add this line to define 'now'
    
    // Verifica periodo di raffreddamento per reazioni
    const inCoolingPeriod = this.rules.isInCoolingPeriod && 
                           (this.rules.isInCoolingPeriod(mol1, now) || 
                           this.rules.isInCoolingPeriod(mol2, now));
    
    // Calcola vettore direzione
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
    
    // Normalizza direzione
    const dirNorm = direction.map(d => d / distance);
    
    // Calcola forze nette
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
    
    // Aggiungi forza di repulsione familiare se necessario
    if (this.rules.areRelated && this.rules.areRelated(mol1.id, mol2.id)) {
        const repulsionFactor = this.rules.getFamilyRepulsionFactor && 
                               this.rules.getFamilyRepulsionFactor(mol1.id, mol2.id, now);
        
        if (repulsionFactor > 0) {
            const minDistance = this.rules.getConstant('min_distance');
            const effectiveDistance = Math.max(distance, minDistance);
            const repulsiveForce = dirNorm.map(
                x => -x * repulsionFactor / (effectiveDistance ** 1.5)
            );
            
            // Applica alla forza risultante
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
};

// Fix for shouldReact method
EnhancedChemistry.prototype.shouldReact = function(mol1, mol2) {
    // Verifica periodo di raffreddamento
    const now = performance.now(); // Add this line to define 'now'
    if (this.rules.isInCoolingPeriod && 
       (this.rules.isInCoolingPeriod(mol1, now) || 
        this.rules.isInCoolingPeriod(mol2, now))) {
        return false;
    }
    
    // Continua con la logica originale
    const relativeSpeed = Math.sqrt(
        Math.pow(mol1.velocity[0] - mol2.velocity[0], 2) +
        Math.pow(mol1.velocity[1] - mol2.velocity[1], 2) +
        Math.pow(mol1.velocity[2] - mol2.velocity[2], 2)
    );
    
    const baseReactionProb = 0.1 * this.temperature * relativeSpeed;
    const sizeFactor = 1.0 / (1.0 + Math.log(mol1.number + mol2.number));
    
    return Math.random() < baseReactionProb * sizeFactor;
};

// Estensione del metodo step di EnhancedChemistry per supportare campi quantistici
const originalStep = EnhancedChemistry.prototype.step;
EnhancedChemistry.prototype.step = function() {
    // Esegui passo normale
    originalStep.call(this);
    
    // Aggiorna campi quantistici
    updateQuantumFields(this);
    
    // Periodicamente, puliamo le relazioni obsolete
    if (Math.random() < 0.05) {
        this.rules.cleanupOldRelations();
    }
};

// Modifica processReaction per passare i parametri aggiuntivi alle funzioni di condizione
EnhancedChemistry.prototype.processReaction = function(mol1, mol2) {
    // Cerca regole di reazione applicabili
    for (const rule of this.rules.reaction_rules) {
        // Passa anche le molecole complete alla funzione di condizione
        if (rule.condition(mol1.prime_factors, mol2.prime_factors, mol1, mol2)) {
            if (Math.random() < rule.probability * this.temperature) {
                // Applica la regola di reazione
                const products = rule.effect(mol1, mol2);
                
                // Imposta tempo di reazione per effetto visivo
                const now = performance.now();
                products.forEach(p => p.setReactionTime(now));
                
                return products;
            }
        }
    }
    
    return [];
};

// Aggiungi metodo per visualizzare le relazioni tra molecole (utile per debug)
EnhancedChemistry.prototype.getRelationshipCount = function() {
    let count = 0;
    for (const [molId, relations] of this.rules.familyRegistry.entries()) {
        count += relations.size;
    }
    return count / 2; // Dividiamo per 2 perché ogni relazione è contata due volte (bidirezionale)
};