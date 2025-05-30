// rules.js 

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
            'G': 0.1,           // Costante gravitazionale (ridotta per renderla più debole)
            'k': 1.0,            // Costante elettrostatica
            'min_distance': 0.05, // Distanza minima per evitare singolarità
            'max_force': 100.0,   // Forza massima applicabile
            'damping': 0.99,     // Smorzamento
            'temperature_factor': 1.0,
            'quantum_strength': 0.5,
            'time_scale': 0.8,   // Rallentato per gravità più graduale
            'repulsion_decay_time': 1000, // Tempo in ms per cui le particelle restano ripulsive
            'cooling_period': 10000, // Tempo in ms durante il quale le molecole non reagiscono dopo una reazione
            'entropy_factor': 0.3, // Fattore di entropia per prevenire il collasso gravitazionale
            'family_repulsion': 2.0, // Intensità della repulsione tra "parenti"
            'random_interaction_probability': 0.1 // Probabilità di considerare molecole casuali nel caching
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
    
    // Nuova funzione per stabilire relazioni tra molecole
    establishRelationship(parentId, childId, relationType, timestamp) {
        this.registerFamily(parentId, childId, timestamp);
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
        return (molecule.lastReactionTime && 
                currentTime - molecule.lastReactionTime < this.getConstant('cooling_period'));
    }
}

// Miglioramenti alle regole di reazione in rules.js
function createCustomRules() {
    const rules = new SimulationRules();

    // ==== REGOLE DI INTERAZIONE FISICA ====

    // 1. Gravità universale (tra tutte le molecole)
    // Ora è più debole e il suo effetto dipende dalla massa
    function gravityCondition(factors1, factors2) { 
        return true; // Sempre attiva
    }
    
    function gravityForce(direction, distance, mass1, mass2) {
        const G = rules.getConstant('G');
        const minDistance = rules.getConstant('min_distance');
        // Previeni divisione per zero
        const effectiveDistance = Math.max(distance, minDistance);
        
        // Effetto gravitazionale proporzionale alle masse
        // Più forte per molecole più massive, più debole per molecole leggere
        return direction.map(x => x * G * mass1 * mass2 / (effectiveDistance ** 2));
    }
    
    rules.addInteractionRule(
        new InteractionRule("Gravity", 1.0, gravityCondition, gravityForce, 0.2, 
        "Attrazione gravitazionale universale - debole ma sempre presente")
    );

    // 2. Risonanza di fattori primi
    function resonanceCondition(factors1, factors2) {
        const shared = Object.keys(factors1).filter(key => factors2.hasOwnProperty(key));
        return shared.length >= 1;
    }
    
    function resonanceForce(direction, distance, mass1, mass2) {
        // Crea un movimento orbitale proporzionale alle masse
        const orbital = crossProduct(direction, [0, 1, 0]);
        const norm = Math.sqrt(orbital[0]**2 + orbital[1]**2 + orbital[2]**2);
        
        if (norm > 0) { 
            // Aggiunge dipendenza dalle masse
            const massEffect = Math.sqrt(mass1 * mass2) / (mass1 + mass2);
            // Forza più forte quando più vicini, ma non troppo da esplodere
            return orbital.map(x => x / norm * massEffect / (distance ** 0.8)); 
        }
        return [0, 0, 0];
    }
    
    rules.addInteractionRule(
        new InteractionRule("Prime Resonance", 2.0, resonanceCondition, resonanceForce, 0.3,
        "Risonanza orbitale tra molecole con fattori primi comuni - influenzata dalla massa")
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
        return direction.map(x => -x * k * Math.abs(charge1 * charge2) / (effectiveDistance ** 2));
    }
    
    rules.addInteractionRule(
        new InteractionRule("Charge Repulsion", 0.5, chargeRepulsionCondition, chargeRepulsionForce, 0.3,
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
        return direction.map(x => x * k * Math.abs(charge1 * charge2) / (effectiveDistance ** 2));
    }
    
    rules.addInteractionRule(
        new InteractionRule("Charge Attraction", 1.0, chargeAttractionCondition, chargeAttractionForce, 0.25,
        "Attrazione tra molecole con cariche di segno opposto")
    );

    // ==== REGOLE DI REAZIONE CHIMICA ====

    // 1. Fusione - quando molecole con fattori primi comuni collidono
    function fusionCondition(factors1, factors2, mol1, mol2) {
        // Devono avere almeno un fattore primo in comune
        const shared = Object.keys(factors1).filter(key => factors2.hasOwnProperty(key));
        // E ciascuna deve avere almeno un fattore non in comune con l'altra
        const diff1 = Object.keys(factors1).filter(key => !factors2.hasOwnProperty(key));
        const diff2 = Object.keys(factors2).filter(key => !factors1.hasOwnProperty(key));
        
        // Verifica anche se la massa combinata non è eccessiva
        const combinedMass = mol1.mass + mol2.mass;
        
        return shared.length >= 1 && diff1.length > 0 && diff2.length > 0 && combinedMass < 100;
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
        newMol.parentIds = [mol1.id, mol2.id];
        newMol.reactionType = 'fusion';
        
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
        fragment.parentIds = [mol1.id, mol2.id];
        fragment.reactionType = 'emission';
        
        return [newMol, fragment];
    }
    
    rules.addReactionRule(
        new ReactionRule("Fusion", 3.0, fusionCondition, fusionEffect, 0.25,
        "Fusione di molecole con fattori comuni con emissione di frammento")
    );

    // 2. Fissione - scissione di molecole grandi
    function fissionCondition(factors1, factors2, mol1, mol2) {
        // La fissione avviene quando una molecola grande collide con una piccola
        // Usa il nuovo sistema di massa invece del numero
        const largerMass = Math.max(mol1.mass, mol2.mass);
        const smallerMass = Math.min(mol1.mass, mol2.mass);
        
        return largerMass > 50 && smallerMass < 25;
    }
    
    function fissionEffect(mol1, mol2) {
        // Identifica molecola più grande
        const bigger = mol1.mass > mol2.mass ? mol1 : mol2;
        const smaller = mol1.mass > mol2.mass ? mol2 : mol1;
        
        // Verifica se ha abbastanza fattori primi
        const factors = Object.entries(bigger.prime_factors);
        if (factors.length < 2) {
            return []; // Non può dividersi
        }
        
        // Dividi i fattori in due prodotti approssimativamente uguali
        let product1 = 1;
        let product2 = 1;
        
        for (const [prime, count] of factors) {
            const primeNum = parseInt(prime);
            if (product1 <= product2) {
                product1 *= primeNum ** count;
            } else {
                product2 *= primeNum ** count;
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
        
        // Imposta relazioni di parentela
        mol1New.parentIds = [bigger.id];
        mol2New.parentIds = [bigger.id];
        mol1New.reactionType = 'fission';
        mol2New.reactionType = 'fission';
        
        // La molecola piccola ottiene una velocità casuale
        const smallerMolecule = new PrimeMolecule(smaller.number, [...smaller.position]);
        smallerMolecule.velocity = [
            Math.random() * 2 - 1, 
            Math.random() * 2 - 1, 
            Math.random() * 2 - 1
        ];
        smallerMolecule.parentIds = [smaller.id];
        smallerMolecule.reactionType = 'catalyst';
        
        return [mol1New, mol2New, smallerMolecule];
    }
    
    rules.addReactionRule(
        new ReactionRule("Fission", 2.0, fissionCondition, fissionEffect, 0.35,
        "Scissione di molecole grandi in frammenti più piccoli")
    );

    // 3. Decadimento spontaneo di molecole molto grandi
    function decayCondition(factors1, factors2, mol1, mol2) {
        // Verifica se una delle molecole è instabile in base alla massa
        return mol1.mass > 80 || mol2.mass > 80;
    }
    
    function decayEffect(mol1, mol2) {
        let mol = (mol1.mass > 80) ? mol1 : mol2;
        if (mol.mass <= 80) return [];
        
        const factors = Object.entries(mol.prime_factors);
        if (factors.length < 1) {
            return []; // Niente da decadere
        }
        
        // Decadi in più particelle
        const numFragments = Math.min(factors.length, 3); // Limita frammenti
        const newMolecules = [];
        let remainingNumber = mol.number;
        
        for (let i = 0; i < numFragments; i++) {
            const [prime, count] = factors[i];
            const fragmentNumber = parseInt(prime);
            
            // Riduci il numero rimanente
            remainingNumber = Math.floor(remainingNumber / fragmentNumber);
            
            // Crea un frammento con questo fattore primo
            const fragment = new PrimeMolecule(fragmentNumber, mol.position.slice());
            
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
            fragment.parentIds = [mol.id];
            fragment.reactionType = 'decay';
            
            newMolecules.push(fragment);
        }
        
        // Aggiungi un frammento residuo se c'è ancora qualcosa
        if (remainingNumber > 1) {
            const remainder = new PrimeMolecule(remainingNumber, mol.position.slice());
            
            // Velocità residuo in direzione opposta per conservare momento
            remainder.velocity = newMolecules[0].velocity.map(v => -v * 0.5);
            remainder.parentIds = [mol.id];
            remainder.reactionType = 'decay';
            
            newMolecules.push(remainder);
        }
        
        return newMolecules;
    }
    
    rules.addReactionRule(
        new ReactionRule("Spontaneous Decay", 1.8, decayCondition, decayEffect, 0.15,
        "Decadimento spontaneo di molecole instabili in frammenti più piccoli")
    );

    // 4. Scambio di esponenti (quando due molecole hanno fattori primi in comune ma con esponenti diversi)
    function exponentExchangeCondition(factors1, factors2) {
        // Trova i fattori primi condivisi
        const sharedPrimes = Object.keys(factors1)
            .filter(key => factors2.hasOwnProperty(key));
            
        // Verifica se hanno esponenti diversi
        for (const prime of sharedPrimes) {
            if (factors1[prime] !== factors2[prime]) {
                return true;
            }
        }
        
        return false;
    }
    
    function exponentExchangeEffect(mol1, mol2) {
        // Trova i fattori primi condivisi
        const sharedPrimes = Object.keys(mol1.prime_factors)
            .filter(key => mol2.prime_factors.hasOwnProperty(key));
            
        // Trova un fattore primo con esponenti diversi
        const transferPrime = sharedPrimes.find(
            prime => mol1.prime_factors[prime] !== mol2.prime_factors[prime]
        );
        
        if (!transferPrime) return [];
        
        let newMol1 = new PrimeMolecule(mol1.number, mol1.position.slice());
        let newMol2 = new PrimeMolecule(mol2.number, mol2.position.slice());
        
        // Trasferisci parte del fattore primo
        if (newMol1.prime_factors[transferPrime] > newMol2.prime_factors[transferPrime]) {
            newMol1.number = Math.floor(newMol1.number / parseInt(transferPrime));
            newMol2.number = newMol2.number * parseInt(transferPrime);
        } else {
            newMol2.number = Math.floor(newMol2.number / parseInt(transferPrime));
            newMol1.number = newMol1.number * parseInt(transferPrime);
        }
        
        // Aggiorna velocità (scambio di momento)
        newMol1.velocity = [...mol1.velocity];
        newMol2.velocity = [...mol2.velocity];
        
        // Relazioni di parentela
        newMol1.parentIds = [mol1.id, mol2.id];
        newMol2.parentIds = [mol1.id, mol2.id];
        newMol1.reactionType = 'exchange';
        newMol2.reactionType = 'exchange';
        
        return [newMol1, newMol2];
    }
    
    rules.addReactionRule(
        new ReactionRule("Exponent Exchange", 2.5, exponentExchangeCondition, exponentExchangeEffect, 0.4,
        "Scambio di esponenti tra molecole con fattori primi comuni")
    );

    // 5. Scambio di fattori primi (quando hanno fattori sia comuni che diversi)
    function primeExchangeCondition(factors1, factors2) {
        const shared = Object.keys(factors1).filter(key => factors2.hasOwnProperty(key));
        const diff1 = Object.keys(factors1).filter(key => !factors2.hasOwnProperty(key));
        const diff2 = Object.keys(factors2).filter(key => !factors1.hasOwnProperty(key));
        
        return shared.length > 0 && diff1.length > 0 && diff2.length > 0;
    }
    
    function primeExchangeEffect(mol1, mol2) {
        // Identifica fattori primi unici di ciascuna molecola
        const diff1 = Object.keys(mol1.prime_factors)
            .filter(key => !mol2.prime_factors.hasOwnProperty(key));
            
        const diff2 = Object.keys(mol2.prime_factors)
            .filter(key => !mol1.prime_factors.hasOwnProperty(key));
            
        if (diff1.length === 0 || diff2.length === 0) return [];
        
        // Scegli un fattore primo da scambiare da ciascuna molecola
        const prime1 = diff1[0];
        const prime2 = diff2[0];
        
        // Crea nuove molecole con i numeri scambiati
        let newMol1 = new PrimeMolecule(mol1.number, mol1.position.slice());
        let newMol2 = new PrimeMolecule(mol2.number, mol2.position.slice());
        
        // Scambia i fattori primi
        newMol1.number = Math.floor(newMol1.number / parseInt(prime1));
        newMol1.number = newMol1.number * parseInt(prime2);
        
        newMol2.number = Math.floor(newMol2.number / parseInt(prime2));
        newMol2.number = newMol2.number * parseInt(prime1);
        
        // Leggera attrazione dopo lo scambio
        const midpoint = [
            (mol1.position[0] + mol2.position[0]) / 2,
            (mol1.position[1] + mol2.position[1]) / 2,
            (mol1.position[2] + mol2.position[2]) / 2
        ];
        
        // Velocità verso il centro con componente originale
        for (let i = 0; i < 3; i++) {
            const toCenter = midpoint[i] - mol1.position[i];
            newMol1.velocity[i] = mol1.velocity[i] * 0.8 + toCenter * 0.05;
            
            const toCenter2 = midpoint[i] - mol2.position[i];
            newMol2.velocity[i] = mol2.velocity[i] * 0.8 + toCenter2 * 0.05;
        }
        
        // Relazioni di parentela
        newMol1.parentIds = [mol1.id, mol2.id];
        newMol2.parentIds = [mol1.id, mol2.id];
        newMol1.reactionType = 'exchange';
        newMol2.reactionType = 'exchange';
        
        return [newMol1, newMol2];
    }
    
    rules.addReactionRule(
        new ReactionRule("Prime Exchange", 2.2, primeExchangeCondition, primeExchangeEffect, 0.3,
        "Scambio di fattori primi tra molecole")
    );
    
    // 6. Catalisi - molecole prime accelerano reazioni vicine
    function catalysisCondition(factors1, factors2, mol1, mol2) {
        // Verifica se una delle molecole è un numero primo puro
        const isPrimePure1 = Object.keys(factors1).length === 1 && 
                            factors1[Object.keys(factors1)[0]] === 1;
                            
        const isPrimePure2 = Object.keys(factors2).length === 1 && 
                            factors2[Object.keys(factors2)[0]] === 1;
                            
        // L'altra deve essere composta
        return (isPrimePure1 && !isPrimePure2) || (isPrimePure2 && !isPrimePure1);
    }
    
    function catalysisEffect(mol1, mol2) {
        // Identifica catalizzatore (molecola prima pura) e substrato
        const isPrimePure1 = Object.keys(mol1.prime_factors).length === 1 && 
                            mol1.prime_factors[Object.keys(mol1.prime_factors)[0]] === 1;
                            
        const catalyst = isPrimePure1 ? mol1 : mol2;
        const substrate = isPrimePure1 ? mol2 : mol1;
        
        // Il catalizzatore rimane invariato, il substrato può trasformarsi
        const catalystPrime = parseInt(Object.keys(catalyst.prime_factors)[0]);
        
        // Ottieni una copia del catalizzatore originale
        const catalystCopy = new PrimeMolecule(catalyst.number, catalyst.position.slice());
        catalystCopy.velocity = [...catalyst.velocity];
        catalystCopy.parentIds = [catalyst.id];
        catalystCopy.reactionType = 'catalyst';
        
        // Possibili trasformazioni:
        // 1. Se il substrato è divisibile per il catalizzatore, si divide
        if (substrate.number % catalystPrime === 0) {
            const quotient = Math.floor(substrate.number / catalystPrime);
            
            // Crea nuova molecola divisa
            const product = new PrimeMolecule(quotient, substrate.position.slice());
            
            // Velocità leggermente modificata
            product.velocity = substrate.velocity.map(v => v * 1.1);
            product.parentIds = [substrate.id, catalyst.id];
            product.reactionType = 'catalyzed';
            
            // Il catalizzatore rimbalza
            catalystCopy.velocity = catalyst.velocity.map(v => -v * 0.9);
            
            return [catalystCopy, product];
        }
        
        // 2. Altrimenti, prova una moltiplicazione condizionale
        const result = substrate.number * catalystPrime;
        if (result < 1000) { // Limita la crescita eccessiva
            // Crea prodotto moltiplicato
            const product = new PrimeMolecule(result, substrate.position.slice());
            
            // Rallenta dopo l'aumento di massa
            product.velocity = substrate.velocity.map(v => v * 0.7);
            product.parentIds = [substrate.id, catalyst.id];
            product.reactionType = 'catalyzed';
            
            // Il catalizzatore rimbalza
            catalystCopy.velocity = catalyst.velocity.map(v => -v * 0.9);
            
            return [catalystCopy, product];
        }
        
        // Se nessuna reazione è possibile, restituisci le molecole originali
        catalystCopy.parentIds = [catalyst.id];
        
        const substrateCopy = new PrimeMolecule(substrate.number, substrate.position.slice());
        substrateCopy.velocity = [...substrate.velocity];
        substrateCopy.parentIds = [substrate.id];
        
        return [catalystCopy, substrateCopy];
    }
    
    rules.addReactionRule(
        new ReactionRule("Catalysis", 2.7, catalysisCondition, catalysisEffect, 0.25,
        "Molecole prime pure agiscono da catalizzatori trasformando altre molecole")
    );

    // Funzioni helper
    function crossProduct(a, b) {
        return [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]
        ];
    }
    
    function calculateCharge(factors) {
        let charge = 0;
        for (const prime in factors) {
            const count = factors[prime];
            const p = parseInt(prime);
            
            if (p === 2) {
                charge += count * 2; // 2 ha carica fortemente positiva
            } else if (p % 4 === 1) {
                charge += count; // Primi 4k+1 hanno carica positiva
            } else {
                charge -= count; // Altri primi hanno carica negativa
            }
        }
        
        // Calcola il totale dei fattori invece di usare il numero stesso
        const totalFactors = Object.entries(factors)
            .reduce((acc, [prime, exponent]) => acc + parseInt(prime) * exponent, 0);
            
        return charge / (1 + Math.sqrt(totalFactors));
    }

    return rules;
}

export { createCustomRules, InteractionRule, ReactionRule, SimulationRules };