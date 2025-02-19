// molecule.js - Versione Modificata

export class PrimeMolecule {
    constructor(number, position) {
        this.id = null; // Will be assigned by PrimeChemistry
        this.number = number;
        this.position = position;
        this.velocity = [0, 0, 0]; 
        this.acceleration = [0, 0, 0];
        this._primeFactors = null; // Lazy load
        this._mass = null; // Lazy load
        this._charge = null; // Lazy load  
        this._color = null; // Lazy load
        this.lastReactionTime = -Infinity;
        
        // Ensure angular velocity is initialized
        this.angularVelocity = [
            (Math.random() - 0.5) * 0.2,
            (Math.random() - 0.5) * 0.2,
            (Math.random() - 0.5) * 0.2
        ];
        
        // Tracking per le relazioni
        this.parentIds = []; // IDs delle molecole "genitori"
        this.reactionType = null; // Tipo di reazione che ha creato questa molecola
    }

    // Lazy-loaded getters
    get prime_factors() {
        if (!this._primeFactors) {
            this._primeFactors = this.factorize(this.number);
        }
        return this._primeFactors;
    }
    
    get mass() {
        if (this._mass === null) {
            // Nuovo calcolo della massa: somma dei (primi * esponenti)
            this._mass = this.calculateMassByAddition();
        }
        return this._mass;
    }
    
    get charge() {
        if (this._charge === null) {
            this._charge = this.calculate_charge();
        }
        return this._charge;
    }
    
    get color() {
        if (!this._color) {
            this._color = this.generate_color();
        }
        return this._color;
    }

    set color(val) {
        this._color = val;
    }

    // Calcolo della massa con addizione di (primi * esponenti)
    calculateMassByAddition() {
        let totalMass = 0;
        for (const prime in this.prime_factors) {
            const primeValue = parseInt(prime);
            const exponent = this.prime_factors[prime];
            
            // Aggiungi alla massa (primo * esponente)
            totalMass += primeValue * exponent;
        }
        
        // Scala la massa per renderla significativa nella simulazione
        // Evita che sia troppo piccola o troppo grande
        const scaledMass = Math.log(1 + totalMass) * 5;
        
        return Math.max(1, scaledMass); // Massa minima di 1
    }

    // Factorization più efficiente che consente più esponenti dello stesso numero primo
    factorize(n) {
        // Static cache for factorization results
        if (!PrimeMolecule.factorCache) {
            PrimeMolecule.factorCache = new Map();
        }
        
        // Check cache first
        if (PrimeMolecule.factorCache.has(n)) {
            return { ...PrimeMolecule.factorCache.get(n) }; // Return copy to avoid mutations
        }
        
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
        
        // Cache the result (only for smaller numbers to avoid memory bloat)
        if (n < 10000) {
            PrimeMolecule.factorCache.set(n, { ...factors });
        }
        
        return factors;
    }    

    calculate_charge() {
        let charge = 0;
        for (const prime in this.prime_factors) {
            const count = this.prime_factors[prime];
            const primeNum = parseInt(prime);
            
            if (primeNum == 2) {
                charge += count * 2;
            } else if (primeNum % 4 === 1) {
                charge += count;
            } else {
                charge -= count;
            }
        }
        
        // Normalizza la carica in base alla massa invece che al logaritmo del numero
        return charge / (1 + Math.sqrt(this.mass));
    }

   generate_color() {
        if (Object.keys(this.prime_factors).length === 0) {
            return [0.5, 0.5, 0.5]; // Gray for 1
        }

        // Cache HSV to RGB conversions
        if (!PrimeMolecule.hsvCache) {
            PrimeMolecule.hsvCache = new Map();
        }

        let h = 0;
        let s = 0;
        let v = 0;

        // Get the primes sorted for consistent coloring
        const primes = Object.keys(this.prime_factors).map(Number).sort((a, b) => a - b);
        
        for (let i = 0; i < primes.length; i++) {
            const prime = primes[i];
            const count = this.prime_factors[prime];

            h += (0.618033988749895 * prime) + (i * 0.1);
            h %= 1; // Keep hue within 0-1

            // Usa gli esponenti per influenzare la saturazione
            s += count / (1 + Math.sqrt(this.mass));
            v += 1 / (1 + Math.log(prime));
        }
        
        // Ensure values stay between 0 and 1
        s = Math.min(1, Math.max(0.3, s));
        v = Math.min(1, Math.max(0.5, v));

        // Check if we've cached this HSV combination
        const key = `${h.toFixed(3)}_${s.toFixed(3)}_${v.toFixed(3)}`;
        if (PrimeMolecule.hsvCache.has(key)) {
            return [...PrimeMolecule.hsvCache.get(key)]; // Return copy
        }

        // Convert HSV to RGB
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
        
        const result = [r, g, b];
        
        // Cache the result (limit cache size to prevent memory bloat)
        if (PrimeMolecule.hsvCache.size < 1000) {
            PrimeMolecule.hsvCache.set(key, [...result]);
        }
        
        return result;
    }

    setReactionTime(time) {
        this.lastReactionTime = time;
    }

    getHaloColor() {
        // Calculate alpha based on time since last reaction
        const timeSinceReaction = performance.now() - this.lastReactionTime;
        const alpha = Math.max(0, 1 - timeSinceReaction / 1000); // Fade out over 1 second
        return [...this.color, alpha]; // Return color with alpha
    }

    getHaloSize() {
        const timeSinceReaction = performance.now() - this.lastReactionTime;
        // Usa la massa invece del logaritmo del numero
        const size = 0.3 + 0.05 * this.mass;
        const factor = Math.max(0, 1 - timeSinceReaction/500); // Shrink
        return size * (1 + factor);
    }
    
    // Clean up instance for garbage collection
    dispose() {
        this._primeFactors = null;
        this._color = null;
    }
    
    // Static method to clear caches
    static clearCaches() {
        if (PrimeMolecule.factorCache) {
            PrimeMolecule.factorCache.clear();
        }
        if (PrimeMolecule.hsvCache) {
            PrimeMolecule.hsvCache.clear();
        }
    }
}