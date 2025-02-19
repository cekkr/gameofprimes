// molecule.js
export class PrimeMolecule {
    constructor(number, position) {
        this.number = number;
        this.position = position;
        this.velocity = [0, 0, 0]; // Initialized to zero, will be set in PrimeChemistry
        this.acceleration = [0, 0, 0];
        this.prime_factors = this.factorize(number);
        this.mass = Math.log2(number) * 2;
        this.charge = this.calculate_charge();
        this.color = this.generate_color();
        this.lastReactionTime = -Infinity; // For halo effect

        // Ensure angular velocity is initialized here
        this.angularVelocity = [
            (Math.random() - 0.5) * 0.2,  // Adjust range as needed
            (Math.random() - 0.5) * 0.2,
            (Math.random() - 0.5) * 0.2
        ];
    }

    factorize(n) {
        const factors = {};
        let d = 2;
        while (n > 1) {
            while (n % d === 0) {
                factors[d] = (factors[d] || 0) + 1;
                n /= d;
            }
            d += 1;
            if (d * d > n) {
                if (n > 1) {
                    factors[n] = (factors[n] || 0) + 1;
                }
                break;
            }
        }
        return factors;
    }

    calculate_charge() {
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

   generate_color() {
        if (Object.keys(this.prime_factors).length === 0) {
            return [0.5, 0.5, 0.5]; // Gray for 1
        }

        let h = 0;
        let s = 0;
        let v = 0;

        const primes = Object.keys(this.prime_factors).map(Number).sort((a, b) => a - b); // Sorted primes
        for (let i = 0; i < primes.length; i++) {
            const prime = primes[i];
            const count = this.prime_factors[prime];

            // Use golden ratio for hue, but offset by index to further differentiate
            h += (0.618033988749895 * prime) + (i * 0.1);
            h %= 1; // Keep hue within 0-1

            s += count / (1 + Math.log(Math.max(...primes))); // Saturation based on count
            v += 1 / (1 + Math.log(prime));  // Value based on prime
        }
         s = Math.min(1, s); //Ensure s and v stays between 0 and 1
         v = Math.min(1,v);

        // Convert HSV to RGB
        let r, g, b;
        const i = Math.floor(h * 6);
        const f = h * 6 - i;
        const p = v * (1 - s);
        const q = v * (1 - s * f);
        const t = v * (1 - s * (1 - f));

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
    setReactionTime(time) {
        this.lastReactionTime = time;
    }

    getHaloColor() {
        //  alpha based on time since last reaction
        const timeSinceReaction = performance.now() - this.lastReactionTime;
        const alpha = Math.max(0, 1 - timeSinceReaction / 1000); // Fade out over 1 second
        return [...this.color, alpha]; // Return color with alpha
    }

    getHaloSize() {
        const timeSinceReaction = performance.now() - this.lastReactionTime;
        const size = 0.4 + 0.1 * Math.log2(this.number)
        const factor = Math.max(0, 1 - timeSinceReaction/500); // Shrink
        return size * (1 + factor);

    }
}