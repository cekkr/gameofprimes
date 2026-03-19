const GOLDEN_RATIO = 0.618033988749895;

const PRIME_VALUES = [2];
const PRIME_INDEX_CACHE = new Map([[2, 1]]);
const FACTOR_CACHE = new Map();
const HSV_CACHE = new Map();

function randomNormal() {
    let u = 0;
    let v = 0;
    while (u === 0) u = Math.random();
    while (v === 0) v = Math.random();
    return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
}

function clamp(value, minValue, maxValue) {
    return Math.max(minValue, Math.min(maxValue, value));
}

function isPrimeCandidate(value) {
    if (value < 2) {
        return false;
    }
    if (value % 2 === 0) {
        return value === 2;
    }

    const limit = Math.floor(Math.sqrt(value));
    for (let divisor = 3; divisor <= limit; divisor += 2) {
        if (value % divisor === 0) {
            return false;
        }
    }
    return true;
}

export function primeIndex(primeValue) {
    const prime = Math.floor(Math.abs(primeValue));
    const cached = PRIME_INDEX_CACHE.get(prime);
    if (cached !== undefined) {
        return cached;
    }

    let candidate = PRIME_VALUES[PRIME_VALUES.length - 1] + 1;
    while (PRIME_VALUES[PRIME_VALUES.length - 1] < prime) {
        if (isPrimeCandidate(candidate)) {
            PRIME_VALUES.push(candidate);
            PRIME_INDEX_CACHE.set(candidate, PRIME_VALUES.length);
        }
        candidate += 1;
    }

    const resolved = PRIME_INDEX_CACHE.get(prime);
    if (resolved === undefined) {
        throw new Error(`${prime} is not a prime number`);
    }
    return resolved;
}

export function factorizeInteger(value) {
    const n = Math.max(2, Math.floor(Math.abs(value)));
    const cached = FACTOR_CACHE.get(n);
    if (cached) {
        return { ...cached };
    }

    const factors = {};
    let remaining = n;
    let divisor = 2;

    while (remaining > 1) {
        while (remaining % divisor === 0) {
            factors[divisor] = (factors[divisor] || 0) + 1;
            remaining /= divisor;
        }
        divisor += 1;
        if (divisor * divisor > remaining) {
            if (remaining > 1) {
                factors[remaining] = (factors[remaining] || 0) + 1;
            }
            break;
        }
    }

    if (n < 200000) {
        FACTOR_CACHE.set(n, { ...factors });
        if (FACTOR_CACHE.size > 5000) {
            FACTOR_CACHE.clear();
        }
    }

    return factors;
}

function recursiveIndexScalar(indexValue, depth) {
    if (depth <= 0 || indexValue < 2) {
        return indexValue;
    }

    const factors = factorizeInteger(indexValue);
    const entries = Object.entries(factors);
    if (entries.length === 0) {
        return indexValue;
    }

    let nextIndex = 0.0;
    let totalWeight = 0.0;
    for (const [factorPrime, exponent] of entries) {
        const prime = Number(factorPrime);
        nextIndex += primeIndex(prime) * exponent;
        totalWeight += exponent;
    }

    if (totalWeight <= 0.0) {
        return indexValue;
    }

    const averaged = Math.max(2, Math.round(nextIndex / totalWeight));
    return indexValue + 0.5 * recursiveIndexScalar(averaged, depth - 1);
}

function normalizeVector3(vector) {
    const norm = Math.sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
    if (norm <= 1e-8) {
        return [0, 0, 0];
    }
    return [vector[0] / norm, vector[1] / norm, vector[2] / norm];
}

export class PrimeMolecule {
    constructor(number, position) {
        this.id = null;
        this.number = Math.max(2, Math.floor(Math.abs(number)));
        this.position = [
            Number(position?.[0] ?? 0),
            Number(position?.[1] ?? 0),
            Number(position?.[2] ?? 0)
        ];

        // Start with non-zero motion to avoid static dead states.
        this.velocity = [randomNormal() * 0.9, randomNormal() * 0.45, randomNormal() * 0.9];
        this.acceleration = [0, 0, 0];
        this.inertial_force = [randomNormal(), randomNormal() * 0.35, randomNormal()];

        this.angularVelocity = [
            (Math.random() - 0.5) * 0.2,
            (Math.random() - 0.5) * 0.2,
            (Math.random() - 0.5) * 0.2
        ];

        this._primeFactors = null;
        this.prime_indices = {};
        this._mass = 1;
        this.even_index_weight = 0.0;
        this.odd_index_weight = 0.0;
        this.fractal_signature = [0.0, 0.0, 0.0];
        this._charge = 0.0;
        this._color = [0.5, 0.5, 0.5];

        this.lastReactionTime = -Infinity;
        this.parentIds = [];
        this.reactionType = null;

        this.refreshDerivedProperties();
    }

    get prime_factors() {
        if (!this._primeFactors) {
            this._primeFactors = factorizeInteger(this.number);
        }
        return this._primeFactors;
    }

    get mass() {
        return this._mass;
    }

    get charge() {
        return this._charge;
    }

    get color() {
        return this._color;
    }

    set color(value) {
        if (Array.isArray(value) && value.length === 3) {
            this._color = [
                clamp(Number(value[0]), 0, 1),
                clamp(Number(value[1]), 0, 1),
                clamp(Number(value[2]), 0, 1)
            ];
        }
    }

    refreshDerivedProperties() {
        this._primeFactors = factorizeInteger(this.number);
        this.prime_indices = {};

        for (const key of Object.keys(this._primeFactors)) {
            const prime = Number(key);
            this.prime_indices[prime] = primeIndex(prime);
        }

        this._mass = this.calculateMassByPrimeIndex();
        const [evenWeight, oddWeight] = this.calculateParityWeights();
        this.even_index_weight = evenWeight;
        this.odd_index_weight = oddWeight;
        this.fractal_signature = this.buildFractalSignature(3);
        this._charge = this.calculateCharge();
        this._color = this.generateColor();
    }

    calculateMassByPrimeIndex() {
        let totalMass = 0;
        for (const [primeKey, exponent] of Object.entries(this.prime_factors)) {
            const prime = Number(primeKey);
            const indexValue = this.prime_indices[prime] || primeIndex(prime);
            totalMass += indexValue * exponent;
        }
        return Math.max(1, totalMass);
    }

    calculateParityWeights() {
        let evenWeight = 0.0;
        let oddWeight = 0.0;

        for (const [primeKey, exponent] of Object.entries(this.prime_factors)) {
            const prime = Number(primeKey);
            const indexValue = this.prime_indices[prime] || primeIndex(prime);
            const contribution = indexValue * exponent;
            if (indexValue % 2 === 0) {
                evenWeight += contribution;
            } else {
                oddWeight += contribution;
            }
        }

        return [evenWeight, oddWeight];
    }

    buildFractalSignature(depth = 3) {
        const signature = [0.0, 0.0, 0.0];

        for (const [primeKey, exponent] of Object.entries(this.prime_factors)) {
            const prime = Number(primeKey);
            const indexValue = this.prime_indices[prime] || primeIndex(prime);
            const recursiveScalar = recursiveIndexScalar(indexValue, depth);
            const weighted = recursiveScalar * exponent;

            signature[0] += Math.sin(weighted * 0.31);
            signature[1] += Math.cos(weighted * 0.53);
            signature[2] += Math.sin(weighted * 0.17 + indexValue);
        }

        return normalizeVector3(signature);
    }

    calculateCharge() {
        let charge = 0.0;
        for (const [primeKey, exponent] of Object.entries(this.prime_factors)) {
            const prime = Number(primeKey);
            const indexValue = this.prime_indices[prime] || primeIndex(prime);
            const sign = indexValue % 2 === 0 ? 1.0 : -1.0;
            charge += sign * indexValue * exponent;
        }

        const fractalBias = this.fractal_signature[0] - this.fractal_signature[1];
        return charge / (1.0 + this._mass) + 0.35 * fractalBias;
    }

    generateColor() {
        if (Object.keys(this.prime_factors).length === 0) {
            return [0.5, 0.5, 0.5];
        }

        let h = 0;
        let s = 0;
        let v = 0;

        const sortedPrimes = Object.keys(this.prime_factors).map(Number).sort((a, b) => a - b);
        for (let i = 0; i < sortedPrimes.length; i++) {
            const prime = sortedPrimes[i];
            const exponent = this.prime_factors[prime] || 0;
            const indexValue = this.prime_indices[prime] || primeIndex(prime);

            h = (h + GOLDEN_RATIO * indexValue + i * 0.083) % 1;
            s += exponent / (1 + Math.sqrt(this._mass));
            v += 1 / (1 + Math.log(indexValue + 1));
        }

        s = clamp(s, 0.30, 1.0);
        v = clamp(v, 0.45, 1.0);

        const cacheKey = `${h.toFixed(4)}_${s.toFixed(4)}_${v.toFixed(4)}`;
        const cached = HSV_CACHE.get(cacheKey);
        if (cached) {
            return [...cached];
        }

        const i = Math.floor(h * 6);
        const f = h * 6 - i;
        const p = v * (1 - s);
        const q = v * (1 - s * f);
        const t = v * (1 - s * (1 - f));

        let r;
        let g;
        let b;
        switch (i % 6) {
            case 0:
                r = v;
                g = t;
                b = p;
                break;
            case 1:
                r = q;
                g = v;
                b = p;
                break;
            case 2:
                r = p;
                g = v;
                b = t;
                break;
            case 3:
                r = p;
                g = q;
                b = v;
                break;
            case 4:
                r = t;
                g = p;
                b = v;
                break;
            default:
                r = v;
                g = p;
                b = q;
                break;
        }

        const color = [clamp(r, 0, 1), clamp(g, 0, 1), clamp(b, 0, 1)];
        if (HSV_CACHE.size > 4096) {
            HSV_CACHE.clear();
        }
        HSV_CACHE.set(cacheKey, [...color]);
        return color;
    }

    setReactionTime(time = performance.now()) {
        this.lastReactionTime = time;
    }

    getHaloColor() {
        const now = typeof performance !== 'undefined' ? performance.now() : Date.now();
        const timeSinceReaction = now - this.lastReactionTime;
        const alpha = Math.max(0, 1 - timeSinceReaction / 1200);
        return [...this.color, alpha];
    }

    getHaloSize() {
        const now = typeof performance !== 'undefined' ? performance.now() : Date.now();
        const timeSinceReaction = now - this.lastReactionTime;
        const baseSize = 0.22 + 0.035 * Math.log1p(this.mass);
        const pulse = Math.max(0, 1 - timeSinceReaction / 700);
        return baseSize * (1 + pulse);
    }

    dispose() {
        this._primeFactors = null;
        this.prime_indices = {};
        this.fractal_signature = [0, 0, 0];
        this.parentIds = [];
    }

    static clearCaches() {
        FACTOR_CACHE.clear();
        HSV_CACHE.clear();
    }
}
