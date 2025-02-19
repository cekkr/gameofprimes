// minimal-worker.js
// Una implementazione minimale del worker che dovrebbe funzionare sempre

console.log("Minimal worker starting");

// Classe semplificata PrimeMolecule
class SimpleMolecule {
  constructor(number, position, id) {
    this.id = id;
    this.number = number;
    this.position = position;
    this.velocity = [
      (Math.random() - 0.5) * 0.5,
      (Math.random() - 0.5) * 0.5,
      (Math.random() - 0.5) * 0.5
    ];
    this.prime_factors = this.factorize(number);
    this.mass = Math.log2(number) * 2;
    this.charge = this.calculateCharge();
    this.color = this.generateColor();
    this.angularVelocity = [
      (Math.random() - 0.5) * 0.2,
      (Math.random() - 0.5) * 0.2,
      (Math.random() - 0.5) * 0.2
    ];
    this.lastReactionTime = -Infinity;
  }

  factorize(n) {
    const factors = {};
    let d = 2;
    let remaining = n;
    
    while (remaining > 1) {
      while (remaining % d === 0) {
        factors[d] = (factors[d] || 0) + 1;
        remaining /= d;
      }
      d += 1;
      if (d * d > remaining) {
        if (remaining > 1) {
          factors[remaining] = (factors[remaining] || 0) + 1;
        }
        break;
      }
    }
    return factors;
  }

  calculateCharge() {
    let charge = 0;
    for (const prime in this.prime_factors) {
      const count = this.prime_factors[prime];
      if (prime == 2) {
        charge += count * 2;
      } else if (prime % 4 === 1) {
        charge += count;
      } else {
        charge -= count;
      }
    }
    return charge / (1 + Math.log(this.number));
  }

  generateColor() {
    const h = Math.random();
    const s = 0.7 + Math.random() * 0.3;
    const v = 0.7 + Math.random() * 0.3;
    
    const i = Math.floor(h * 6);
    const f = h * 6 - i;
    const p = v * (1 - s);
    const q = v * (1 - s * f);
    const t = v * (1 - s * (1 - f));

    let r, g, b;
    switch (i % 6) {
      case 0: r = v, g = t, b = p; break;
      case 1: r = q, g = v, b = p; break;
      case 2: r = p, g = v, b = t; break;
      case 3: r = p, g = q, b = v; break;
      case 4: r = t, g = p, b = v; break;
      case 5: r = v, g = p, b = q; break;
    }
    
    return [r, g, b];
  }
}

// Simulazione semplificata
class SimpleSimulation {
  constructor(size, moleculeCount, maxNumber) {
    this.size = size;
    this.maxNumber = maxNumber;
    this.molecules = [];
    this.temperature = 1.0;
    this.nextId = 1;
    
    // Inizializza molecole random
    for (let i = 0; i < moleculeCount; i++) {
      const pos = [
        Math.random() * size - size / 2,
        Math.random() * size - size / 2,
        Math.random() * size - size / 2
      ];
      const number = Math.floor(Math.random() * (maxNumber-2)) + 2;
      const id = `mol-${this.nextId++}`;
      this.molecules.push(new SimpleMolecule(number, pos, id));
    }
  }
  
  step() {
    // Aggiorna posizioni
    this.molecules.forEach(mol => {
      // Aggiungi una piccola randomicità alla velocità
      mol.velocity[0] += (Math.random() - 0.5) * 0.01;
      mol.velocity[1] += (Math.random() - 0.5) * 0.01;
      mol.velocity[2] += (Math.random() - 0.5) * 0.01;
      
      // Aggiorna posizione
      mol.position[0] += mol.velocity[0];
      mol.position[1] += mol.velocity[1];
      mol.position[2] += mol.velocity[2];
      
      // Boundary check
      for (let i = 0; i < 3; i++) {
        if (Math.abs(mol.position[i]) > this.size/2) {
          mol.position[i] = Math.sign(mol.position[i]) * this.size/2;
          mol.velocity[i] *= -0.9;
        }
      }
    });
    
    // Simulazione di temperatura
    this.temperature = 1.0 + 0.1 * Math.sin(performance.now() / 1000);
  }
  
  // Ottiene dati ottimizzati per l'invio
  getMoleculeData() {
    return {
      molecules: this.molecules.map(mol => ({
        id: mol.id,
        number: mol.number,
        position: [...mol.position],
        velocity: [...mol.velocity],
        prime_factors: {...mol.prime_factors},
        mass: mol.mass,
        charge: mol.charge,
        color: [...mol.color],
        angularVelocity: [...mol.angularVelocity],
        lastReactionTime: mol.lastReactionTime
      })),
      removedIds: []
    };
  }
}

// Simulazione
let simulation;

// Gestione messaggi
self.onmessage = function(event) {
  console.log("Minimal worker received message:", event.data);
  
  if (event.data.type === 'init') {
    console.log("Initializing simulation...");
    const { size, moleculeCount, maxNumber } = event.data;
    simulation = new SimpleSimulation(size, moleculeCount, maxNumber);
    
    // Invia dati iniziali
    self.postMessage({
      type: 'update',
      molecules: simulation.getMoleculeData(),
      temperature: simulation.temperature
    });
    
    console.log("Simulation initialized with", moleculeCount, "molecules");
  } 
  else if (event.data.type === 'step') {
    if (!simulation) {
      console.error("Cannot step: simulation not initialized");
      return;
    }
    
    // Esegui un passo di simulazione
    simulation.step();
    
    // Invia dati aggiornati
    self.postMessage({
      type: 'update',
      molecules: simulation.getMoleculeData(),
      temperature: simulation.temperature
    });
  }
  else if (event.data.type === 'cleanup') {
    simulation = null;
    self.postMessage({ type: 'cleanup_complete' });
  }
};

console.log("Minimal worker initialized and ready");