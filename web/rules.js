import { PrimeMolecule } from './molecule.js';

class InteractionRule {
    constructor(name, priority, condition, forceFunction, strength = 1.0, description = "") {
        this.name = name;
        this.priority = priority;
        this.condition = condition;
        this.force_function = forceFunction;
        this.strength = strength;
        this.description = description;
    }
}

class ReactionRule {
    constructor(name, priority, condition, effect, probability = 0.5, description = "") {
        this.name = name;
        this.priority = priority;
        this.condition = condition;
        this.effect = effect;
        this.probability = probability;
        this.description = description;
    }
}

class SimulationRules {
    constructor() {
        this.interaction_rules = [];
        this.reaction_rules = [];
        this.constants = {
            'G': 0.1,
            'k': 1.0,
            'min_distance': 0.5,
            'max_force': 10.0,
            'damping': 0.95,
            'temperature_factor': 1.0,
            'quantum_strength': 0.2,
            'time_scale': 1.0,
            'repulsion_decay_time': 2000, // Tempo in ms per cui le particelle restano ripulsive
            'cooling_period': 500, // Tempo in ms durante il quale le molecole non reagiscono dopo una reazione
            'entropy_factor': 0.3, // Fattore di entropia per prevenire il collasso gravitazionale
            'family_repulsion': 1.5 // Intensità della repulsione tra "parenti"
        };
        
        // Registro di relazioni parentali: mantiene traccia di quali molecole sono "parenti"
        this.familyRegistry = new Map();
    }

    addInteractionRule(rule) {
        this.interaction_rules.push(rule);
        this.interaction_rules.sort((a, b) => b.priority - a.priority);
    }

    addReactionRule(rule) {
        this.reaction_rules.push(rule);
        this.reaction_rules.sort((a, b) => b.priority - a.priority);
    }

    setConstant(name, value) {
        this.constants[name] = value;
    }

    getConstant(name) {
        return this.constants[name] || 0.0;
    }
    
    // Registra una relazione di parentela tra molecole
    registerFamily(parentId, childId, timestamp) {
        if (!this.familyRegistry.has(parentId)) {
            this.familyRegistry.set(parentId, new Map());
        }
        if (!this.familyRegistry.has(childId)) {
            this.familyRegistry.set(childId, new Map());
        }
        
        // Registrazione bidirezionale della relazione
        this.familyRegistry.get(parentId).set(childId, timestamp);
        this.familyRegistry.get(childId).set(parentId, timestamp);
        
        // Pulisci vecchie relazioni ogni tanto
        if (Math.random() < 0.01) {
            this.cleanupOldRelations();
        }
    }
    
    // Verifica se due molecole sono parenti
    areRelated(mol1Id, mol2Id) {
        // Se una delle molecole non ha relazioni registrate
        if (!this.familyRegistry.has(mol1Id) || !this.familyRegistry.has(mol2Id)) {
            return false;
        }
        
        // Verifica relazione diretta
        return this.familyRegistry.get(mol1Id).has(mol2Id);
    }
    
    // Calcola il fattore di repulsione in base al tempo trascorso
    getFamilyRepulsionFactor(mol1Id, mol2Id, currentTime) {
        if (!this.areRelated(mol1Id, mol2Id)) return 0;
        
        const relationTime = this.familyRegistry.get(mol1Id).get(mol2Id);
        const timePassed = currentTime - relationTime;
        const repulsionDecayTime = this.getConstant('repulsion_decay_time');
        
        // Repulsione massima subito dopo la reazione, diminuisce nel tempo
        if (timePassed > repulsionDecayTime) {
            // Dopo il tempo di decadimento, rimuovi la relazione
            this.familyRegistry.get(mol1Id).delete(mol2Id);
            this.familyRegistry.get(mol2Id).delete(mol1Id);
            return 0;
        }
        
        // Repulsione che diminuisce linearmente col tempo
        return this.getConstant('family_repulsion') * (1 - timePassed / repulsionDecayTime);
    }
    
    // Rimuovi relazioni più vecchie del tempo di decadimento
    cleanupOldRelations() {
        const now = performance.now();
        const maxAge = this.getConstant('repulsion_decay_time');
        
        for (const [molId, relations] of this.familyRegistry.entries()) {
            for (const [relatedId, timestamp] of relations.entries()) {
                if (now - timestamp > maxAge) {
                    relations.delete(relatedId);
                }
            }
            
            // Rimuovi completamente le molecole senza più relazioni
            if (relations.size === 0) {
                this.familyRegistry.delete(molId);
            }
        }
    }
    
    // Determina se una molecola è in "periodo di raffreddamento" dopo una reazione
    isInCoolingPeriod(molecule, currentTime) {
        return (currentTime - molecule.lastReactionTime) < this.getConstant('cooling_period');
    }
}

function createCustomRules() {
    const rules = new SimulationRules();

    // ==== REGOLE DI INTERAZIONE FISICA ====

    // 1. Gravità universale modificata con fattore di entropia
    function gravityCondition(factors1, factors2) { 
        return true; // Sempre attiva
    }
    
    function gravityForce(direction, distance, mass1, mass2) {
        const G = rules.getConstant('G');
        const minDistance = rules.getConstant('min_distance');
        const entropyFactor = rules.getConstant('entropy_factor');
        
        // Previeni divisione per zero
        const effectiveDistance = Math.max(distance, minDistance);
        
        // Forza gravitazionale standard
        const baseForce = direction.map(x => x * G * mass1 * mass2 / (effectiveDistance ** 2));
        
        // A grandi distanze la forza rimane attrattiva, ma a distanze ravvicinate
        // aggiungiamo un termine entropico che contrasta il collasso gravitazionale
        if (distance < 5.0) {
            // Termine entropico: repulsione che aumenta a distanze molto ravvicinate
            const entropyTerm = direction.map(x => -x * entropyFactor * mass1 * mass2 / (effectiveDistance ** 4));
            
            // Combina i due effetti
            return baseForce.map((f, i) => f + entropyTerm[i]);
        }
        
        return baseForce;
    }
    
    rules.addInteractionRule(
        new InteractionRule("Modified Gravity", 1.0, gravityCondition, gravityForce, 0.15, 
        "Attrazione gravitazionale con termine entropico per prevenire il collasso")
    );

    // 2. Risonanza di fattori primi con oscillazione
    function resonanceCondition(factors1, factors2) {
        const shared = Object.keys(factors1).filter(key => factors2.hasOwnProperty(key));
        return shared.length >= 1;
    }
    
    function resonanceForce(direction, distance, charge1, charge2) {
        // Crea un movimento orbitale oscillante
        const orbital = crossProduct(direction, [0, 1, 0]);
        const norm = Math.sqrt(orbital[0]**2 + orbital[1]**2 + orbital[2]**2);
        
        if (norm > 0) {
            // Modula la forza con una funzione sinusoidale in base alla distanza
            // Questo crea orbite più stabili e interessanti
            const oscillationFactor = Math.sin(distance * 1.5) * 0.5 + 0.5;
            return orbital.map(x => (x / norm / (distance ** 0.8)) * oscillationFactor); 
        }
        return [0, 0, 0];
    }
    
    rules.addInteractionRule(
        new InteractionRule("Prime Resonance", 2.0, resonanceCondition, resonanceForce, 0.3,
        "Risonanza orbitale oscillante tra molecole con fattori primi comuni")
    );

    // 3. Repulsione elettrostatica tra molecole con carica dello stesso segno
    function chargeRepulsionCondition(factors1, factors2) {
        // Ottieni le cariche (positivo = attrattiva, negativo = repulsiva)
        const charge1 = calculateCharge(factors1);
        const charge2 = calculateCharge(factors2);
        // Attiva solo se le cariche hanno lo stesso segno
        return Math.sign(charge1) === Math.sign(charge2);
    }
    
    function chargeRepulsionForce(direction, distance, charge1, charge2) {
        const k = rules.getConstant('k');
        const minDistance = rules.getConstant('min_distance');
        // Previeni divisione per zero
        const effectiveDistance = Math.max(distance, minDistance);
        // Forza repulsiva
        return direction.map(x => -x * k * Math.abs(charge1) * Math.abs(charge2) / (effectiveDistance ** 2));
    }
    
    rules.addInteractionRule(
        new InteractionRule("Charge Repulsion", 1.5, chargeRepulsionCondition, chargeRepulsionForce, 0.2,
        "Repulsione tra molecole con cariche dello stesso segno")
    );

    // 4. Attrazione elettrostatica tra molecole con cariche opposte
    function chargeAttractionCondition(factors1, factors2) {
        const charge1 = calculateCharge(factors1);
        const charge2 = calculateCharge(factors2);
        // Attiva solo se le cariche hanno segno opposto
        return Math.sign(charge1) !== Math.sign(charge2) && charge1 !== 0 && charge2 !== 0;
    }
    
    function chargeAttractionForce(direction, distance, charge1, charge2) {
        const k = rules.getConstant('k');
        const minDistance = rules.getConstant('min_distance');
        const effectiveDistance = Math.max(distance, minDistance);
        // Forza attrattiva
        return direction.map(x => x * k * Math.abs(charge1) * Math.abs(charge2) / (effectiveDistance ** 2));
    }
    
    rules.addInteractionRule(
        new InteractionRule("Charge Attraction", 1.5, chargeAttractionCondition, chargeAttractionForce, 0.25,
        "Attrazione tra molecole con cariche di segno opposto")
    );
    
    // 5. NUOVA REGOLA: Repulsione tra molecole "parenti"
    function familyRepulsionCondition(factors1, factors2, mol1, mol2) {
        // Verifica se le molecole hanno una relazione registrata
        return rules.areRelated(mol1.id, mol2.id);
    }
    
    function familyRepulsionForce(direction, distance, mass1, mass2, mol1, mol2) {
        const now = performance.now();
        const repulsionFactor = rules.getFamilyRepulsionFactor(mol1.id, mol2.id, now);
        
        if (repulsionFactor <= 0) return [0, 0, 0];
        
        const minDistance = rules.getConstant('min_distance');
        const effectiveDistance = Math.max(distance, minDistance);
        
        // Forza repulsiva che diminuisce con la distanza
        return direction.map(x => -x * repulsionFactor / (effectiveDistance ** 1.5));
    }
    
    rules.addInteractionRule(
        new InteractionRule("Family Repulsion", 3.0, familyRepulsionCondition, familyRepulsionForce, 1.0,
        "Repulsione temporanea tra molecole che hanno recentemente reagito")
    );

    // 6. NUOVA REGOLA: Attrazione/repulsione basata su congruenza modulo
    function modularInteractionCondition(factors1, factors2) {
        // Calcola i numeri totali
        const num1 = Object.entries(factors1).reduce((acc, [prime, count]) => acc * (prime ** count), 1);
        const num2 = Object.entries(factors2).reduce((acc, [prime, count]) => acc * (prime ** count), 1);
        
        // Attiva se almeno uno dei numeri è maggiore di 10
        return num1 > 10 || num2 > 10;
    }
    
    function modularInteractionForce(direction, distance, mass1, mass2, mol1, mol2) {
        // Usa la congruenza modulo per determinare attrazione o repulsione
        const modPrime = 7; // Un numero primo per il calcolo del modulo
        const res1 = mol1.number % modPrime;
        const res2 = mol2.number % modPrime;
        
        // Se i residui sono uguali: repulsione; se diversi: attrazione
        const isRepulsive = res1 === res2;
        
        const minDistance = rules.getConstant('min_distance');
        const effectiveDistance = Math.max(distance, minDistance);
        const forceMagnitude = 0.2 / (effectiveDistance ** 1.8);
        
        // Applica forza attrattiva o repulsiva basata sulla congruenza
        return direction.map(x => x * forceMagnitude * (isRepulsive ? -1 : 1));
    }
    
    rules.addInteractionRule(
        new InteractionRule("Modular Interaction", 1.8, modularInteractionCondition, modularInteractionForce, 0.25,
        "Interazione basata sulla congruenza modulo che crea cluster diversificati")
    );

    // ==== REGOLE DI REAZIONE CHIMICA MIGLIORATE ====

    // 1. Fusione migliorata - con periodo di raffreddamento
    function fusionCondition(factors1, factors2, mol1, mol2) {
        // Verifica periodo di raffreddamento
        const now = performance.now();
        if (rules.isInCoolingPeriod(mol1, now) || rules.isInCoolingPeriod(mol2, now)) {
            return false;
        }
        
        // Devono avere almeno un fattore primo in comune
        const shared = Object.keys(factors1).filter(key => factors2.hasOwnProperty(key));
        // E ciascuna deve avere almeno un fattore non in comune con l'altra
        const diff1 = Object.keys(factors1).filter(key => !factors2.hasOwnProperty(key));
        const diff2 = Object.keys(factors2).filter(key => !factors1.hasOwnProperty(key));
        
        return shared.length >= 1 && diff1.length > 0 && diff2.length > 0;
    }
    
    function fusionEffect(mol1, mol2) {
        // Prodotto delle due molecole
        let newNumber = mol1.number * mol2.number;
        
        // Calcola i fattori primi condivisi
        const sharedPrimes = Object.keys(mol1.prime_factors)
            .filter(key => mol2.prime_factors.hasOwnProperty(key))
            .map(Number);
            
        // Emetti il fattore primo più piccolo tra quelli condivisi
        const emittedPrime = sharedPrimes.length > 0 ? 
            Math.min(...sharedPrimes) : 2;
            
        // Dividi per il fattore emesso
        newNumber = Math.floor(newNumber / emittedPrime);

        // Posizione pesata per massa
        const totalMass = mol1.mass + mol2.mass;
        const newPos = [
            (mol1.position[0] * mol1.mass + mol2.position[0] * mol2.mass) / totalMass,
            (mol1.position[1] * mol1.mass + mol2.position[1] * mol2.mass) / totalMass,
            (mol1.position[2] * mol1.mass + mol2.position[2] * mol2.mass) / totalMass
        ];
        
        // Velocità pesata per massa (conservazione momento)
        const newVel = [
            (mol1.velocity[0] * mol1.mass + mol2.velocity[0] * mol2.mass) / totalMass,
            (mol1.velocity[1] * mol1.mass + mol2.velocity[1] * mol2.mass) / totalMass,
            (mol1.velocity[2] * mol1.mass + mol2.velocity[2] * mol2.mass) / totalMass
        ];
        
        // Crea nuova molecola
        const newMol = new PrimeMolecule(newNumber, newPos);
        newMol.velocity = newVel;
        
        // Puoi anche creare un frammento (il fattore primo emesso)
        const fragmentPos = [
            newPos[0] + (Math.random() - 0.5) * 0.5,
            newPos[1] + (Math.random() - 0.5) * 0.5,
            newPos[2] + (Math.random() - 0.5) * 0.5
        ];
        
        const fragment = new PrimeMolecule(emittedPrime, fragmentPos);
        
        // Velocità opposta per conservazione del momento
        const fragmentVel = newVel.map(v => -v * 1.5);
        fragment.velocity = fragmentVel;
        
        // Registra le relazioni di parentela
        const now = performance.now();
        newMol.setReactionTime(now);
        fragment.setReactionTime(now);
        
        // Genera ID temporanei per poter registrare le relazioni
        const tempIdNewMol = `temp-${Math.random().toString(36).substr(2, 9)}`;
        const tempIdFragment = `temp-${Math.random().toString(36).substr(2, 9)}`;
        newMol.id = tempIdNewMol;
        fragment.id = tempIdFragment;
        
        // Registra relazioni tra tutte le molecole coinvolte
        rules.registerFamily(mol1.id, tempIdNewMol, now);
        rules.registerFamily(mol2.id, tempIdNewMol, now);
        rules.registerFamily(mol1.id, tempIdFragment, now);
        rules.registerFamily(mol2.id, tempIdFragment, now);
        rules.registerFamily(tempIdNewMol, tempIdFragment, now);
        
        return [newMol, fragment];
    }
    
    rules.addReactionRule(
        new ReactionRule("Fusion", 3.0, fusionCondition, fusionEffect, 0.25,
        "Fusione di molecole con fattori comuni con emissione di frammento e periodo di raffreddamento")
    );

    // 2. Fissione migliorata - con periodo di raffreddamento
    function fissionCondition(factors1, factors2, mol1, mol2) {
        // Verifica periodo di raffreddamento
        const now = performance.now();
        if (rules.isInCoolingPeriod(mol1, now) || rules.isInCoolingPeriod(mol2, now)) {
            return false;
        }
        
        // Verifica se almeno una delle molecole è abbastanza grande
        const num1 = Object.entries(factors1)
            .reduce((acc, [prime, count]) => acc * (prime ** count), 1);
        const num2 = Object.entries(factors2)
            .reduce((acc, [prime, count]) => acc * (prime ** count), 1);
            
        // La fissione avviene quando una molecola grande collide con una piccola
        const largerNumber = Math.max(num1, num2);
        const smallerNumber = Math.min(num1, num2);
        
        return largerNumber > 80 && smallerNumber < 40;
    }
    
    function fissionEffect(mol1, mol2) {
        // Identifica molecola più grande
        const bigger = mol1.number > mol2.number ? mol1 : mol2;
        const smaller = mol1.number > mol2.number ? mol2 : mol1;
        
        // Verifica se ha abbastanza fattori primi
        const factors = Object.entries(bigger.prime_factors);
        if (factors.length < 2) {
            return []; // Non può dividersi
        }
        
        // Dividi i fattori in due prodotti approssimativamente uguali
        let product1 = 1;
        let product2 = 1;
        
        for (const [prime, count] of factors) {
            if (product1 <= product2) {
                product1 *= prime ** count;
            } else {
                product2 *= prime ** count;
            }
        }
        
        // Calcola velocità di espulsione
        const directionOut = [
            Math.random() * 2 - 1, 
            Math.random() * 2 - 1, 
            Math.random() * 2 - 1
        ];
        
        // Normalizza il vettore
        const norm = Math.sqrt(
            directionOut[0]**2 + 
            directionOut[1]**2 + 
            directionOut[2]**2
        );
        
        const normalizedDir = directionOut.map(x => x / norm * 2.0);
        
        // Crea i due frammenti
        const mol1New = new PrimeMolecule(product1, [
            bigger.position[0] + (Math.random() * 0.2 - 0.1),
            bigger.position[1] + (Math.random() * 0.2 - 0.1),
            bigger.position[2] + (Math.random() * 0.2 - 0.1)
        ]);
        
        const mol2New = new PrimeMolecule(product2, [
            bigger.position[0] + (Math.random() * 0.2 - 0.1),
            bigger.position[1] + (Math.random() * 0.2 - 0.1),
            bigger.position[2] + (Math.random() * 0.2 - 0.1)
        ]);
        
        // Velocità opposte
        mol1New.velocity = normalizedDir;
        mol2New.velocity = normalizedDir.map(x => -x);
        
        // La molecola piccola ottiene una velocità casuale
        smaller.velocity = [
            Math.random() * 2 - 1, 
            Math.random() * 2 - 1, 
            Math.random() * 2 - 1
        ];
        
        // Imposta il tempo di reazione e registra relazioni
        const now = performance.now();
        mol1New.setReactionTime(now);
        mol2New.setReactionTime(now);
        
        // Genera ID temporanei
        const tempId1 = `temp-${Math.random().toString(36).substr(2, 9)}`;
        const tempId2 = `temp-${Math.random().toString(36).substr(2, 9)}`;
        mol1New.id = tempId1;
        mol2New.id = tempId2;
        
        // Registra relazioni
        rules.registerFamily(bigger.id, tempId1, now);
        rules.registerFamily(bigger.id, tempId2, now);
        rules.registerFamily(smaller.id, tempId1, now);
        rules.registerFamily(smaller.id, tempId2, now);
        rules.registerFamily(tempId1, tempId2, now);
        
        return [mol1New, mol2New, smaller];
    }
    
    rules.addReactionRule(
        new ReactionRule("Fission", 2.0, fissionCondition, fissionEffect, 0.35,
        "Scissione di molecole grandi in frammenti più piccoli con periodo di raffreddamento")
    );

    // 3. Decadimento spontaneo migliorato
    function decayCondition(factors1, factors2, mol1, mol2) {
        // Verifica periodo di raffreddamento
        const now = performance.now();
        if (rules.isInCoolingPeriod(mol1, now) || rules.isInCoolingPeriod(mol2, now)) {
            return false;
        }
        
        // Verifica se la prima molecola è instabile (numero grande)
        const num1 = Object.entries(factors1)
            .reduce((acc, [prime, count]) => acc * (prime ** count), 1);
        const num2 = Object.entries(factors2)
            .reduce((acc, [prime, count]) => acc * (prime ** count), 1);
            
        return num1 > 150 || num2 > 150; // Molecole oltre 150 sono instabili
    }
    
    function decayEffect(mol1, mol2) {
        let mol = (mol1.number > 150) ? mol1 : mol2;
        if (mol.number <= 150) mol = (mol2.number > 150) ? mol2 : mol1;
        if (mol.number <= 150) return [];
        
        const factors = Object.entries(mol.prime_factors);
        if (factors.length < 1) {
            return []; // Niente da decadere
        }
        
        // Decadi in più particelle
        const numFragments = Math.min(factors.length, 3); // Limita frammenti
        const newMolecules = [];
        let remainingNumber = mol.number;
        
        const now = performance.now();
        const tempIds = [];
        
        for (let i = 0; i < numFragments; i++) {
            const [prime, count] = factors[i];
            const fragmentNumber = parseInt(prime);
            
            // Riduci il numero rimanente
            remainingNumber = Math.floor(remainingNumber / fragmentNumber);
            
            // Crea un frammento con questo fattore primo
            const fragment = new PrimeMolecule(fragmentNumber, mol.position.slice());
            fragment.setReactionTime(now);
            
            // Dai al frammento una velocità di espulsione
            const decayDirection = [
                Math.random() * 2 - 1, 
                Math.random() * 2 - 1, 
                Math.random() * 2 - 1
            ];
            
            const norm = Math.sqrt(
                decayDirection[0]**2 + 
                decayDirection[1]**2 + 
                decayDirection[2]**2
            );
            
            const velocity = decayDirection.map(x => x / norm * 2.5);
            fragment.velocity = velocity;
            
            // Genera ID temporaneo
            const tempId = `temp-${Math.random().toString(36).substr(2, 9)}`;
            fragment.id = tempId;
            tempIds.push(tempId);
            
            newMolecules.push(fragment);
        }
        
        // Aggiungi un frammento residuo se c'è ancora qualcosa
        if (remainingNumber > 1) {
            const remainder = new PrimeMolecule(remainingNumber, mol.position.slice());
            remainder.setReactionTime(now);
            
            // Velocità residuo in direzione opposta per conservare momento
            remainder.velocity = newMolecules[0].velocity.map(v => -v * 0.5);
            
            // Genera ID temporaneo
            const tempId = `temp-${Math.random().toString(36).substr(2, 9)}`;
            remainder.id = tempId;
            tempIds.push(tempId);
            
            newMolecules.push(remainder);
        }
        
        // Registra relazioni tra tutti i frammenti e la molecola originale
        for (const id of tempIds) {
            rules.registerFamily(mol.id, id, now);
            
            // Registra anche relazioni tra i frammenti
            for (const otherId of tempIds) {
                if (id !== otherId) {
                    rules.registerFamily(id, otherId, now);
                }
            }
        }
        
        return newMolecules;
    }
    
    rules.addReactionRule(
        new ReactionRule("Spontaneous Decay", 1.8, decayCondition, decayEffect, 0.15,
        "Decadimento spontaneo di molecole instabili con ripulsione tra frammenti")
    );
    
    rules.addReactionRule(
        new ReactionRule("Twin Prime Oscillation", 2.5, twinPrimeOscillationCondition, twinPrimeOscillationEffect, 0.4,
        "Oscillazione quantistica tra numeri primi gemelli che crea un campo di energia temporaneo")
    );
    
    // NUOVA REGOLA: Risonanza Antipolare
    function antipolarResonanceCondition(factors1, factors2, mol1, mol2) {
        // Verifica periodo di raffreddamento
        const now = performance.now();
        if (rules.isInCoolingPeriod(mol1, now) || rules.isInCoolingPeriod(mol2, now)) {
            return false;
        }
        
        // Calcola le cariche
        const charge1 = calculateCharge(factors1);
        const charge2 = calculateCharge(factors2);
        
        // Verifica se hanno cariche opposte e abbastanza forti
        return Math.sign(charge1) !== Math.sign(charge2) && 
               Math.abs(charge1) > 0.5 && 
               Math.abs(charge2) > 0.5;
    }
    
    function antipolarResonanceEffect(mol1, mol2) {
        const now = performance.now();
        
        // Calcola il numero massimo di ciascuna molecola
        const num1 = mol1.number;
        const num2 = mol2.number;
        
        // Calcola nuovo numero tramite operazione XOR sui bit
        const newNumber = num1 ^ num2;
        
        // Calcola posizione centrale
        const midPos = [
            (mol1.position[0] + mol2.position[0]) / 2,
            (mol1.position[1] + mol2.position[1]) / 2,
            (mol1.position[2] + mol2.position[2]) / 2
        ];
        
        // Crea nuova molecola risonante
        const resonantMol = new PrimeMolecule(newNumber, midPos);
        resonantMol.setReactionTime(now);
        
        // Calcola velocità basata sulla conservazione del momento
        const m1 = mol1.mass;
        const m2 = mol2.mass;
        for (let i = 0; i < 3; i++) {
            resonantMol.velocity[i] = (mol1.velocity[i] * m1 + mol2.velocity[i] * m2) / (m1 + m2);
            // Aggiungi una componente casuale per evitare collassi gravitazionali
            resonantMol.velocity[i] += (Math.random() - 0.5) * 0.2;
        }
        
        // Genera ID temporaneo
        const tempId = `resonant-${Math.random().toString(36).substr(2, 9)}`;
        resonantMol.id = tempId;
        
        // Registra relazioni di parentela per ripulsione temporanea
        rules.registerFamily(mol1.id, tempId, now);
        rules.registerFamily(mol2.id, tempId, now);
        
        return [resonantMol];
    }
    
    rules.addReactionRule(
        new ReactionRule("Antipolar Resonance", 2.3, antipolarResonanceCondition, antipolarResonanceEffect, 0.3,
        "Risonanza tra molecole con cariche opposte che crea una nuova molecola con proprietà ibride")
    );

    // NUOVA REGOLA: Oscillazione quantistica di numeri primi gemelli
    function twinPrimeOscillationCondition(factors1, factors2, mol1, mol2) {
        // Verifica periodo di raffreddamento
        const now = performance.now();
        if (rules.isInCoolingPeriod(mol1, now) || rules.isInCoolingPeriod(mol2, now)) {
            return false;
        }
        
        // Controlla se sono entrambi numeri primi puri
        const isPure1 = Object.keys(factors1).length === 1 && 
                       factors1[Object.keys(factors1)[0]] === 1;
        const isPure2 = Object.keys(factors2).length === 1 && 
                       factors2[Object.keys(factors2)[0]] === 1;
                       
        if (!isPure1 || !isPure2) return false;
        
        // Ottieni i numeri primi
        const prime1 = parseInt(Object.keys(factors1)[0]);
        const prime2 = parseInt(Object.keys(factors2)[0]);
        
        // Verifica se sono primi gemelli (differenza di 2)
        return Math.abs(prime1 - prime2) === 2;
    }
    
    function twinPrimeOscillationEffect(mol1, mol2) {
        const now = performance.now();
        
        // Crea un campo di energia quantistica tra i due numeri primi gemelli
        // che influenza le particelle circostanti
        const prime1 = parseInt(Object.keys(mol1.prime_factors)[0]);
        const prime2 = parseInt(Object.keys(mol2.prime_factors)[0]);
        
        // Calcola la posizione media come centro del campo quantistico
        const centerPos = [
            (mol1.position[0] + mol2.position[0]) / 2,
            (mol1.position[1] + mol2.position[1]) / 2,
            (mol1.position[2] + mol2.position[2]) / 2
        ];
        
        // Crea un effetto di entanglement che mantiene i due numeri primi in oscillazione
        // senza fonderli, creando uno stato quantistico speciale
        const quantumField = new PrimeMolecule(prime1 * prime2, centerPos);
        quantumField.setReactionTime(now);
        
        // Imposta proprietà speciali per questo campo
        quantumField._isQuantumField = true;
        quantumField._linkedPrimes = [prime1, prime2];
        quantumField._originIds = [mol1.id, mol2.id];
        quantumField._creationTime = now;
        quantumField._lifetime = 5000; // Durata in ms
        
        // Genera un ID temporaneo
        const tempId = `quantum-${Math.random().toString(36).substr(2, 9)}`;
        quantumField.id = tempId;
        
        // Modifica la velocità dei primi gemelli per creare un effetto orbitale
        const distanceVector = [
            mol2.position[0] - mol1.position[0],
            mol2.position[1] - mol1.position[1],
            mol2.position[2] - mol1.position[2]
        ];
        
        // Normalizza
        const distance = Math.sqrt(
            distanceVector[0]**2 + 
            distanceVector[1]**2 + 
            distanceVector[2]**2
        );
        
        // Crea un vettore perpendicolare per l'orbita
        const orbitVector = crossProduct(distanceVector, [0, 1, 0]);
        const orbitNorm = Math.sqrt(
            orbitVector[0]**2 + 
            orbitVector[1]**2 + 
            orbitVector[2]**2
        );
        
        if (orbitNorm > 0) {
            const orbitVelocity = orbitVector.map(v => v / orbitNorm * 0.5);
            
            // Aggiorna velocità per creare un'orbita
            mol1.velocity = orbitVelocity;
            mol2.velocity = orbitVelocity.map(v => -v);
            
            // Aggiungi una componente di attrazione verso il centro del campo
            for (let i = 0; i < 3; i++) {
                mol1.velocity[i] += (centerPos[i] - mol1.position[i]) * 0.01;
                mol2.velocity[i] += (centerPos[i] - mol2.position[i]) * 0.01;
            }
        }
        
        // Registra relazioni di parentela
        rules.registerFamily(mol1.id, tempId, now);
        rules.registerFamily(mol2.id, tempId, now);
        
        // Teniamo i numeri primi originali e aggiungiamo il campo quantistico
        return [mol1, mol2, quantumField];