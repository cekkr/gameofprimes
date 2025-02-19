// rules.js
import { PrimeMolecule } from './molecule.js'; // Import from molecule.js

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
            'time_scale': 1.0
        };
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
}



function createCustomRules() {
    const rules = new SimulationRules();

    // ==== INTERACTION RULES ====

    // Base gravitational attraction
    function gravityCondition(factors1, factors2) {
        return true;
    }

    function gravityForce(direction, distance, mass1, mass2) {
        const G = rules.getConstant('G');
        return direction.map(x => x * G * mass1 * mass2 / (distance ** 2));
    }

    rules.addInteractionRule(new InteractionRule(
        "Gravity",
        1.0,
        gravityCondition,
        gravityForce,
        0.1
    ));

    // Prime resonance orbital motion
    function resonanceCondition(factors1, factors2) {
        const shared = Object.keys(factors1).filter(key => factors2.hasOwnProperty(key));
        return shared.length >= 1;
    }

   function resonanceForce(direction, distance, charge1, charge2) {
        const orbital = crossProduct(direction, [0, 1, 0]);
        const norm = Math.sqrt(orbital[0]**2 + orbital[1]**2 + orbital[2]**2);
        if (norm > 0) {
            return orbital.map(x => x / norm / (distance ** 0.5));
        }
        return [0, 0, 0];
    }

    rules.addInteractionRule(new InteractionRule(
        "Prime Resonance",
        2.0,
        resonanceCondition,
        resonanceForce,
        0.3
    ));


    // ==== REACTION RULES ====

    // Fusion with emission
   function fusionEmissionCondition(factors1, factors2) {
        const shared = Object.keys(factors1).filter(key => factors2.hasOwnProperty(key));
        const count2in1 = factors1[2] || 0;
        const count2in2 = factors2[2] || 0;

        return shared.length >= 2 && (count2in1 + count2in2 ) >= 1;
    }

    function fusionEmissionEffect(mol1, mol2) {
        let newNumber = mol1.number * mol2.number;
        const sharedPrimes = Object.keys(mol1.prime_factors).filter(key => mol2.prime_factors.hasOwnProperty(key));
        const emittedPrime = Math.min(...sharedPrimes.map(Number));  //Convert to numbers for Math.min
        newNumber = Math.floor(newNumber / emittedPrime);

        const newPos = [
            (mol1.position[0] * mol1.mass + mol2.position[0] * mol2.mass) / (mol1.mass + mol2.mass),
            (mol1.position[1] * mol1.mass + mol2.position[1] * mol2.mass) / (mol1.mass + mol2.mass),
            (mol1.position[2] * mol1.mass + mol2.position[2] * mol2.mass) / (mol1.mass + mol2.mass)
        ];
        const mainProduct = new PrimeMolecule(newNumber, newPos);

        const particles = [];
        for (let i = 0; i < 3; i++) {
            const emissionDirection = [Math.random() * 2 - 1, Math.random() * 2 - 1, Math.random() * 2 - 1];
            const norm = Math.sqrt(emissionDirection[0]**2 + emissionDirection[1]**2 + emissionDirection[2]**2);
            const normalizedDirection = emissionDirection.map(x => x/norm);
            const emissionPos = [
                newPos[0] + normalizedDirection[0] * 0.5,
                newPos[1] + normalizedDirection[1] * 0.5,
                newPos[2] + normalizedDirection[2] * 0.5,
            ]

            const particle = new PrimeMolecule(emittedPrime, emissionPos);
            particle.velocity = normalizedDirection.map(x=> x* 3.0); //high velocity
            particles.push(particle);
        }
        return [mainProduct, ...particles];
    }

    rules.addReactionRule(new ReactionRule(
        "Fusion with Emission",
        2.0,
        fusionEmissionCondition,
        fusionEmissionEffect,
        0.4
    ));


    // Fission reaction for large molecules
    function fissionCondition(factors1, factors2) {
        const num1 = Object.entries(factors1).reduce((acc, [prime, count]) => acc * (prime ** count), 1);
        const num2 = Object.entries(factors2).reduce((acc, [prime, count]) => acc * (prime ** count), 1);

        return Math.max(num1, num2) > 100 && Math.min(num1, num2) < 50;
    }


    function fissionEffect(mol1, mol2) {
      const bigger = mol1.number > mol2.number ? mol1 : mol2;
      const smaller = mol1.number > mol2.number ? mol2 : mol1;

      const factors = Object.entries(bigger.prime_factors);
      if (factors.length < 2) {
        return [];
      }

      let product1 = 1;
      let product2 = 1;

      for (const [prime, count] of factors) {
        if (product1 <= product2) {
            product1 *= prime ** count;
        } else {
            product2 *= prime ** count;
        }
      }

        const vel = [Math.random() * 2 - 1, Math.random() * 2 - 1, Math.random() * 2 - 1];
        const norm = Math.sqrt(vel[0]**2 + vel[1]**2 + vel[2]**2);
        const normalizedVel = vel.map(x => x / norm * 2.0);

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

        mol1New.velocity = normalizedVel;
        mol2New.velocity = normalizedVel.map(x => -x);

        smaller.velocity = [Math.random() * 2-1, Math.random()* 2-1, Math.random()* 2-1];
        return [mol1New, mol2New, smaller];

    }


    rules.addReactionRule(new ReactionRule(
        "Fission",
        1.5,
        fissionCondition,
        fissionEffect,
        0.3
    ));

    // Catalytic decomposition
    function catalyticCondition(factors1, factors2) {
        const isPrime1 = Object.keys(factors1).length === 1;
        const isPrime2 = Object.keys(factors2).length === 1;
        const isComposite1 = Object.keys(factors1).length > 1;
        const isComposite2 = Object.keys(factors2).length > 1;

        return (isPrime1 && isComposite2) || (isPrime2 && isComposite1);

    }


    function catalyticEffect(mol1, mol2) {
        const catalyst = Object.keys(mol1.prime_factors).length === 1 ? mol1: mol2;
        const target   = Object.keys(mol1.prime_factors).length === 1 ? mol2: mol1;

        const newMolecules = [];
        for (const prime in target.prime_factors){
            const count = target.prime_factors[prime];
            for (let i = 0; i < count; i++){
                const pos = [
                    target.position[0] + (Math.random() * 1.0 - 0.5),
                    target.position[1] + (Math.random() * 1.0 - 0.5),
                    target.position[2] + (Math.random() * 1.0 - 0.5)
                ];
                const mol = new PrimeMolecule(parseInt(prime), pos); // Parse prime to integer
                const vel = [Math.random() * 2 - 1, Math.random() * 2 - 1, Math.random() * 2 - 1];
                const norm = Math.sqrt(vel[0]**2 + vel[1]**2 + vel[2]**2);
                mol.velocity = vel.map( x=> x/norm * 2.0);
                newMolecules.push(mol);
            }
        }

        catalyst.velocity[0] += Math.random() * 1.0 - 0.5;
        catalyst.velocity[1] += Math.random() * 1.0 - 0.5;
        catalyst.velocity[2] += Math.random() * 1.0 - 0.5;
        newMolecules.push(catalyst);

        return newMolecules;

    }

    rules.addReactionRule(new ReactionRule(
        "Catalytic Decomposition",
        1.0,
        catalyticCondition,
        catalyticEffect,
        0.25
    ));


    //Helper function
    function crossProduct(a, b) {
        return [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]
        ];
    }

    return rules;
}


export { createCustomRules, InteractionRule, ReactionRule, SimulationRules };