// molecule.js
export class PrimeMolecule {
    constructor(number, position) {
        this.number = number;
        this.position = position;
        this.velocity = [0, 0, 0];
        this.acceleration = [0, 0, 0];
        this.prime_factors = this.factorize(number);
        this.mass = Math.log2(number) * 2;
        this.charge = this.calculate_charge();
        this.color = this.generate_color();
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
            return [0.5, 0.5, 0.5];
        }

        const color = [0, 0, 0];
        const maxPrime = Math.max(...Object.keys(this.prime_factors).map(Number)); // Convert to numbers

        for (const prime in this.prime_factors) {
            const count = this.prime_factors[prime];
            const hue = (prime * 0.618033988749895) % 1.0;
            let h = hue * 6.0;
            const c = count / (1 + Math.log(maxPrime));
            const x = c * (1 - Math.abs((h % 2) - 1));

            let rgb;
            if (0 <= h && h < 1) {
                rgb = [c, x, 0];
            } else if (1 <= h && h < 2) {
                rgb = [x, c, 0];
            } else if (2 <= h && h < 3) {
                rgb = [0, c, x];
            } else if (3 <= h && h < 4) {
                rgb = [0, x, c];
            } else if (4 <= h && h < 5) {
                rgb = [x, 0, c];
            } else {
                rgb = [c, 0, x];
            }

            color[0] += rgb[0];
            color[1] += rgb[1];
            color[2] += rgb[2];
        }

        // Normalize and make more vivid
        for (let i = 0; i < 3; i++) {
            color[i] = Math.min(Math.max(color[i], 0), 1); // Clip
            color[i] = color[i] * 0.7 + 0.3; // Vividness
        }

        return color;
    }
}