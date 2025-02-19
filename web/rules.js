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

    function gravityCondition(factors1, factors2) { return true; }
    function gravityForce(direction, distance, mass1, mass2) {
        const G = rules.getConstant('G');
        return direction.map(x => x * G * mass1 * mass2 / (distance ** 2));
    }
    rules.addInteractionRule(new InteractionRule("Gravity", 1.0, gravityCondition, gravityForce, 0.1));

    function resonanceCondition(factors1, factors2) {
        const shared = Object.keys(factors1).filter(key => factors2.hasOwnProperty(key));
        return shared.length >= 1;
    }
    function resonanceForce(direction, distance, charge1, charge2) {
        const orbital = crossProduct(direction, [0, 1, 0]);
        const norm = Math.sqrt(orbital[0]**2 + orbital[1]**2 + orbital[2]**2);
        if (norm > 0) { return orbital.map(x => x / norm / (distance ** 0.5)); }
        return [0, 0, 0];
    }
    rules.addInteractionRule(new InteractionRule("Prime Resonance", 2.0, resonanceCondition, resonanceForce, 0.3));


    // ==== REACTION RULES ====

    // 1. Fusion (More Restrictive)
    function fusionCondition(factors1, factors2) {
        const shared = Object.keys(factors1).filter(key => factors2.hasOwnProperty(key));
        const diff1 = Object.keys(factors1).filter(key => !factors2.hasOwnProperty(key));
        const diff2 = Object.keys(factors2).filter(key => !factors1.hasOwnProperty(key));
		return shared.length >= 1 && diff1.length>0 && diff2.length > 0; // Shared and different primes
    }
    function fusionEffect(mol1, mol2) {
        let newNumber = mol1.number * mol2.number;
		const sharedPrimes = Object.keys(mol1.prime_factors).filter(key => mol2.prime_factors.hasOwnProperty(key));
        const emittedPrime = sharedPrimes.length > 0 ? Math.min(...sharedPrimes.map(Number)) : 2;  //Emit 2 if no shared primes
        newNumber = Math.floor(newNumber / emittedPrime);

        const newPos = [
            (mol1.position[0] * mol1.mass + mol2.position[0] * mol2.mass) / (mol1.mass + mol2.mass),
            (mol1.position[1] * mol1.mass + mol2.position[1] * mol2.mass) / (mol1.mass + mol2.mass),
            (mol1.position[2] * mol1.mass + mol2.position[2] * mol2.mass) / (mol1.mass + mol2.mass)
        ];
        const newMol = new PrimeMolecule(newNumber, newPos);
        return [newMol]; // Simplified fusion
    }
    rules.addReactionRule(new ReactionRule("Fusion", 3.0, fusionCondition, fusionEffect, 0.2));

    // 2. Fission (Size and proximity based)
       function fissionCondition(factors1, factors2) {
      const num1 = Object.entries(factors1).reduce((acc, [prime, count]) => acc * (prime ** count), 1);
      const num2 = Object.entries(factors2).reduce((acc, [prime, count]) => acc * (prime ** count), 1);
      return Math.max(num1, num2) > 80 && Math.min(num1, num2) < 40;
    }

    function fissionEffect(mol1, mol2) {
        // ... (Same fission logic) ...
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
    rules.addReactionRule(new ReactionRule("Fission", 1.5, fissionCondition, fissionEffect, 0.3));

    // 3. Spontaneous Decay (Lower probability, affects larger molecules)
   function decayCondition(factors1, factors2) {
        const num1 = Object.entries(factors1).reduce((acc, [prime, count]) => acc * (prime ** count), 1);
        return num1 > 150;  // Larger molecules are unstable
    }
      function decayEffect(mol1, mol2) {
       let mol = (mol1.number > 150) ? mol1 : mol2;
        if (mol.number <= 150) return [];

      const factors = Object.entries(mol.prime_factors);
      if (factors.length < 1) { // Corrected this line
        return [];  // Nothing to decay into
      }

      //Decay in more than one particle
      const numFragments = Math.min(factors.length, 3); // Limit fragments
      const newMolecules = [];
      let remainingNumber = mol.number;

      for (let i = 0; i < numFragments; i++){
        const [prime, count] = factors[i];
        const fragmentNumber = prime;

        remainingNumber = Math.floor(remainingNumber/fragmentNumber);
        const fragment = new PrimeMolecule(fragmentNumber, mol.position.slice());
        const decayDirection = [Math.random() * 2 - 1, Math.random() * 2 - 1, Math.random() * 2 - 1];
        const norm = Math.sqrt(decayDirection[0]**2 + decayDirection[1]**2 + decayDirection[2]**2);
        const velocity = decayDirection.map(x => x/norm*2.5);
        fragment.velocity = velocity;
        newMolecules.push(fragment);
      }

        if(remainingNumber > 1){ //There could be something left
            const remainder = new PrimeMolecule(remainingNumber, mol.position.slice());
            newMolecules.push(remainder);
        }

      return newMolecules;

    }
    rules.addReactionRule(new ReactionRule("Spontaneous Decay", 1.8, decayCondition, decayEffect, 0.02));

    // 4. Exponent-Difference Reaction
    function exponentDiffCondition(factors1, factors2) {
        const sharedPrimes = Object.keys(factors1).filter(key => factors2.hasOwnProperty(key));
        for (const prime of sharedPrimes) {
            if (factors1[prime] !== factors2[prime]) {
                return true; // Different exponents for a shared prime
            }
        }
        return false;
    }
    function exponentDiffEffect(mol1, mol2) {
        const sharedPrimes = Object.keys(mol1.prime_factors).filter(key => mol2.prime_factors.hasOwnProperty(key));
        const transferPrime = sharedPrimes.find(prime => mol1.prime_factors[prime] !== mol2.prime_factors[prime]);

        let newMol1 = new PrimeMolecule(mol1.number, mol1.position.slice());
        let newMol2 = new PrimeMolecule(mol2.number, mol2.position.slice());

        // Transfer *part* of the prime factor
        if (newMol1.prime_factors[transferPrime] > newMol2.prime_factors[transferPrime]){
            newMol1.number = Math.floor(newMol1.number / transferPrime);
            newMol2.number = newMol2.number * transferPrime;
        } else {
            newMol2.number = Math.floor(newMol2.number / transferPrime);
            newMol1.number = newMol1.number * transferPrime;
        }
        return [newMol1, newMol2];
    }
    rules.addReactionRule(new ReactionRule("Exponent Difference", 2.5, exponentDiffCondition, exponentDiffEffect, 0.4));

    // 5. Prime Exchange (If they have some shared, and some different)
    function primeExchangeCondition(factors1, factors2) {
      const shared = Object.keys(factors1).filter(key => factors2.hasOwnProperty(key));
      const diff1 = Object.keys(factors1).filter(key => !factors2.hasOwnProperty(key));
      const diff2 = Object.keys(factors2).filter(key => !factors1.hasOwnProperty(key));
      return shared.length > 0 && diff1.length > 0 && diff2.length > 0;
    }

    function primeExchangeEffect(mol1, mol2){
      const diff1 = Object.keys(mol1.prime_factors).filter(key => !mol2.prime_factors.hasOwnProperty(key));
      const diff2 = Object.keys(mol2.prime_factors).filter(key => !mol1.prime_factors.hasOwnProperty(key));

      const prime1 = diff1[0];
      const prime2 = diff2[0];

      let newMol1 = new PrimeMolecule(mol1.number, mol1.position.slice());
      let newMol2 = new PrimeMolecule(mol2.number, mol2.position.slice());

      newMol1.number = Math.floor(newMol1.number / prime1);
      newMol1.number = newMol1.number * prime2;

      newMol2.number = Math.floor(newMol2.number/prime2);
      newMol2.number = newMol2.number * prime1;

      return [newMol1, newMol2];

    }
    rules.addReactionRule(new ReactionRule("Prime Exchange", 2.2, primeExchangeCondition, primeExchangeEffect, 0.3));

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