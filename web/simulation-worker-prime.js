import { PrimeMolecule, primeIndex } from './molecule.js';

const DEFAULT_CONSTANTS = {
    interactionDistance: 24.0,
    reactionDistance: 2.4,
    minDistance: 0.35,
    maxForce: 12.0,
    damping: 0.988,
    fluidDrag: 0.08,
    thermalStrength: 0.12,
    chargeCoupling: 0.35,
    collisionRepulsionDistance: 1.15,
    collisionRepulsionStrength: 0.80,
    sharedPrimeRepulsion: 4.2,
    evenIndexAttraction: 2.2,
    oddIndexRepulsion: 1.8,
    mixedParityOrbital: 1.4,
    fractalPolarization: 2.2,
    fractalTorsion: 1.3,
    inertialForceStrength: 2.0,
    inertialForcePersistence: 0.982,
    inertialForceJitter: 0.42,
    maxSpeed: 12.0,
    maxMolecules: 5000
};

function clamp(value, minValue, maxValue) {
    return Math.max(minValue, Math.min(maxValue, value));
}

function randomNormal() {
    let u = 0;
    let v = 0;
    while (u === 0) u = Math.random();
    while (v === 0) v = Math.random();
    return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
}

function vectorMagnitude(vector) {
    return Math.sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
}

function normalizeVector(vector) {
    const norm = vectorMagnitude(vector);
    if (norm <= 1e-8) {
        return [0, 0, 0];
    }
    return [vector[0] / norm, vector[1] / norm, vector[2] / norm];
}

function dot(a, b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function cross(a, b) {
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    ];
}

function smallestPrimeFactor(value) {
    const n = Math.max(2, Math.floor(Math.abs(value)));
    if (n % 2 === 0) {
        return 2;
    }
    for (let divisor = 3; divisor * divisor <= n; divisor += 2) {
        if (n % divisor === 0) {
            return divisor;
        }
    }
    return n;
}

function countSharedPrimeFactors(factorsA, factorsB) {
    let count = 0;
    for (const key of Object.keys(factorsA)) {
        if (Object.prototype.hasOwnProperty.call(factorsB, key)) {
            count += Math.min(factorsA[key], factorsB[key]);
        }
    }
    return count;
}

class PrimeChemistrySimulation {
    constructor(size, moleculeCount, maxNumber, timeScale = 1.0) {
        this.size = Number(size) || 100;
        this.maxNumber = Math.max(32, Number(maxNumber) || 2000);
        this.timeScale = clamp(Number(timeScale) || 1.0, 0.1, 10.0);
        this.temperature = 1.0;
        this.isPaused = false;

        this.constants = { ...DEFAULT_CONSTANTS };

        this.molecules = [];
        this.nextId = 1;
        this.reactionCount = 0;
        this.stepCount = 0;

        this.initializeMolecules(Math.max(2, Number(moleculeCount) || 200));
    }

    initializeMolecules(count) {
        const primes = this.generatePrimes(Math.max(32, Math.floor(this.maxNumber * 0.3)));
        const compounds = this.generateCompounds(Math.max(16, Math.floor(this.maxNumber * 0.2)), primes);
        const numberPool = [...primes, ...compounds];

        for (let i = 0; i < count; i++) {
            const position = [
                (Math.random() - 0.5) * this.size,
                (Math.random() - 0.5) * this.size,
                (Math.random() - 0.5) * this.size
            ];

            const sampledNumber = numberPool[Math.floor(Math.random() * numberPool.length)] ||
                (2 + Math.floor(Math.random() * Math.min(this.maxNumber, 1000)));

            const molecule = new PrimeMolecule(sampledNumber, position);
            molecule.id = `mol-${this.nextId++}`;

            this.molecules.push(molecule);
        }
    }

    generatePrimes(limit) {
        const max = Math.max(2, Math.floor(limit));
        const sieve = Array(max + 1).fill(true);
        sieve[0] = false;
        sieve[1] = false;

        for (let i = 2; i * i <= max; i++) {
            if (!sieve[i]) {
                continue;
            }
            for (let j = i * i; j <= max; j += i) {
                sieve[j] = false;
            }
        }

        const primes = [];
        for (let n = 2; n <= max; n++) {
            if (sieve[n]) {
                primes.push(n);
            }
        }
        return primes;
    }

    generateCompounds(limit, primes) {
        const max = Math.max(4, Math.floor(limit));
        const primeSet = new Set(primes);
        const compounds = [];

        for (let n = 4; n <= max; n++) {
            if (!primeSet.has(n)) {
                compounds.push(n);
            }
        }

        return compounds.length > 0 ? compounds : [4, 6, 8, 9, 10, 12];
    }

    setPauseState(paused) {
        this.isPaused = Boolean(paused);
    }

    setTimeScale(value) {
        this.timeScale = clamp(Number(value) || 1.0, 0.1, 10.0);
    }

    setTemperature(value) {
        this.temperature = clamp(Number(value) || 1.0, 0.05, 4.0);
    }

    addRandomMolecules(count = 20) {
        const toAdd = Math.max(1, Math.min(400, Math.floor(count)));
        for (let i = 0; i < toAdd; i++) {
            const number = 2 + Math.floor(Math.random() * this.maxNumber);
            const molecule = new PrimeMolecule(number, [
                (Math.random() - 0.5) * this.size,
                (Math.random() - 0.5) * this.size,
                (Math.random() - 0.5) * this.size
            ]);
            molecule.id = `mol-${this.nextId++}`;
            this.molecules.push(molecule);
        }
    }

    computeDt() {
        return clamp(this.timeScale * 0.02, 0.002, 0.22);
    }

    updateInertialForce(molecule, dt) {
        const persistence = clamp(this.constants.inertialForcePersistence, 0.0, 0.9995);
        const jitter = Math.max(0.0, this.constants.inertialForceJitter);

        for (let axis = 0; axis < 3; axis++) {
            molecule.inertial_force[axis] =
                molecule.inertial_force[axis] * persistence + randomNormal() * jitter * Math.sqrt(dt);
        }
        molecule.inertial_force[1] *= 0.35;

        const targetMagnitude = this.constants.inertialForceStrength *
            (0.45 + 0.12 * Math.log1p(Math.max(molecule.mass, 1.0)));

        let currentMagnitude = vectorMagnitude(molecule.inertial_force);
        if (currentMagnitude < 1e-8) {
            molecule.inertial_force = [randomNormal(), randomNormal() * 0.35, randomNormal()];
            currentMagnitude = vectorMagnitude(molecule.inertial_force);
        }

        const scale = targetMagnitude / Math.max(currentMagnitude, 1e-8);
        molecule.inertial_force[0] *= scale;
        molecule.inertial_force[1] *= scale;
        molecule.inertial_force[2] *= scale;
    }

    buildSpatialHash() {
        const cellSize = this.constants.interactionDistance;
        const halfSize = this.size * 0.5;
        const grid = new Map();
        const coordinates = new Array(this.molecules.length);

        for (let i = 0; i < this.molecules.length; i++) {
            const molecule = this.molecules[i];
            const cx = Math.floor((molecule.position[0] + halfSize) / cellSize);
            const cy = Math.floor((molecule.position[1] + halfSize) / cellSize);
            const cz = Math.floor((molecule.position[2] + halfSize) / cellSize);
            const key = `${cx}|${cy}|${cz}`;

            coordinates[i] = [cx, cy, cz];
            if (!grid.has(key)) {
                grid.set(key, []);
            }
            grid.get(key).push(i);
        }

        return { grid, coordinates };
    }

    calculatePairForce(mol1, mol2, direction, distance) {
        const softened = distance + 0.4;
        const invSq = 1.0 / (softened * softened);
        const invSoft = 1.0 / (softened * Math.sqrt(softened));
        const totalMass = Math.max(1e-6, mol1.mass + mol2.mass);

        let sharedOverlap = 0.0;
        for (const [primeKey, countA] of Object.entries(mol1.prime_factors)) {
            const countB = mol2.prime_factors[primeKey] || 0;
            if (countB <= 0) {
                continue;
            }
            const prime = Number(primeKey);
            const indexValue = mol1.prime_indices[prime] || primeIndex(prime);
            sharedOverlap += indexValue * Math.min(countA, countB);
        }

        const denominator = (1.0 + mol1.mass) * (1.0 + mol2.mass);
        const evenCoupling = (mol1.even_index_weight * mol2.even_index_weight) / denominator;
        const oddCoupling = (mol1.odd_index_weight * mol2.odd_index_weight) / denominator;
        const mixedCoupling = (
            mol1.even_index_weight * mol2.odd_index_weight +
            mol1.odd_index_weight * mol2.even_index_weight
        ) / denominator;

        const force = [0, 0, 0];

        const sharedScale = (this.constants.sharedPrimeRepulsion * sharedOverlap / (1.0 + totalMass)) * invSq;
        force[0] -= direction[0] * sharedScale;
        force[1] -= direction[1] * sharedScale;
        force[2] -= direction[2] * sharedScale;

        const evenScale = this.constants.evenIndexAttraction * evenCoupling * invSq;
        force[0] += direction[0] * evenScale;
        force[1] += direction[1] * evenScale;
        force[2] += direction[2] * evenScale;

        const oddScale = this.constants.oddIndexRepulsion * oddCoupling * invSoft;
        force[0] -= direction[0] * oddScale;
        force[1] -= direction[1] * oddScale;
        force[2] -= direction[2] * oddScale;

        const orbitalAxisRaw = cross(direction, [0, 1, 0]);
        const orbitalAxis = normalizeVector(orbitalAxisRaw);
        const orbitalScale = this.constants.mixedParityOrbital * mixedCoupling * invSoft;

        force[0] += orbitalAxis[0] * orbitalScale;
        force[1] += orbitalAxis[1] * orbitalScale;
        force[2] += orbitalAxis[2] * orbitalScale;

        const fractalAlignment = dot(mol1.fractal_signature, mol2.fractal_signature);
        const polarizationScale = (-fractalAlignment) * this.constants.fractalPolarization * invSq;
        force[0] += direction[0] * polarizationScale;
        force[1] += direction[1] * polarizationScale;
        force[2] += direction[2] * polarizationScale;

        const torsionCross = cross(mol1.fractal_signature, mol2.fractal_signature);
        const torsionScalar = dot(torsionCross, [0, 1, 0]);
        const torsionScale = this.constants.fractalTorsion * torsionScalar * invSq;
        force[0] += orbitalAxis[0] * torsionScale;
        force[1] += orbitalAxis[1] * torsionScale;
        force[2] += orbitalAxis[2] * torsionScale;

        const chargeScale = -this.constants.chargeCoupling * (mol1.charge * mol2.charge) * invSq / (1.0 + totalMass);
        force[0] += direction[0] * chargeScale;
        force[1] += direction[1] * chargeScale;
        force[2] += direction[2] * chargeScale;

        const repulsionDistance = this.constants.collisionRepulsionDistance;
        if (distance < repulsionDistance) {
            const repulsion = this.constants.collisionRepulsionStrength * (1.0 - distance / repulsionDistance);
            force[0] -= direction[0] * repulsion;
            force[1] -= direction[1] * repulsion;
            force[2] -= direction[2] * repulsion;
        }

        const maxForce = this.constants.maxForce;
        force[0] = clamp(force[0], -maxForce, maxForce);
        force[1] = clamp(force[1], -maxForce, maxForce);
        force[2] = clamp(force[2], -maxForce, maxForce);

        return force;
    }

    maybeReact(mol1, mol2, distance, dt, now) {
        if (distance > this.constants.reactionDistance) {
            return null;
        }

        const sharedCount = countSharedPrimeFactors(mol1.prime_factors, mol2.prime_factors);
        const relativeVelocity = [
            mol1.velocity[0] - mol2.velocity[0],
            mol1.velocity[1] - mol2.velocity[1],
            mol1.velocity[2] - mol2.velocity[2]
        ];
        const relativeSpeed = vectorMagnitude(relativeVelocity);

        const baseProbability = 0.028 + sharedCount * 0.022 + Math.min(0.09, 0.018 * relativeSpeed);
        const scaledProbability = clamp(
            baseProbability * Math.sqrt(dt / 0.02) * (0.65 + 0.35 * this.temperature),
            0.0,
            0.55
        );

        if (Math.random() >= scaledProbability) {
            return null;
        }

        const totalMass = Math.max(1e-6, mol1.mass + mol2.mass);
        const center = [
            (mol1.position[0] * mol1.mass + mol2.position[0] * mol2.mass) / totalMass,
            (mol1.position[1] * mol1.mass + mol2.position[1] * mol2.mass) / totalMass,
            (mol1.position[2] * mol1.mass + mol2.position[2] * mol2.mass) / totalMass
        ];
        const combinedVelocity = [
            (mol1.velocity[0] * mol1.mass + mol2.velocity[0] * mol2.mass) / totalMass,
            (mol1.velocity[1] * mol1.mass + mol2.velocity[1] * mol2.mass) / totalMass,
            (mol1.velocity[2] * mol1.mass + mol2.velocity[2] * mol2.mass) / totalMass
        ];

        const sharedPrimes = Object.keys(mol1.prime_factors)
            .filter(key => Object.prototype.hasOwnProperty.call(mol2.prime_factors, key))
            .map(Number)
            .sort((a, b) => a - b);

        const products = [];

        if (sharedPrimes.length > 0 && Math.random() < 0.72) {
            const emittedPrime = sharedPrimes[0];
            const rawFusion = Math.floor((mol1.number * mol2.number) / Math.max(2, emittedPrime));
            const fusedNumber = this.clampNumber(rawFusion);

            const fused = new PrimeMolecule(fusedNumber, center);
            fused.velocity = [
                combinedVelocity[0] + randomNormal() * 0.08,
                combinedVelocity[1] + randomNormal() * 0.04,
                combinedVelocity[2] + randomNormal() * 0.08
            ];
            fused.parentIds = [mol1.id, mol2.id];
            fused.reactionType = 'fusion';
            fused.setReactionTime(now);
            products.push(fused);

            if (Math.random() < 0.55) {
                const offsetDirection = normalizeVector([randomNormal(), randomNormal() * 0.4, randomNormal()]);
                const fragmentPos = [
                    center[0] + offsetDirection[0] * 0.6,
                    center[1] + offsetDirection[1] * 0.6,
                    center[2] + offsetDirection[2] * 0.6
                ];
                const fragment = new PrimeMolecule(this.clampNumber(emittedPrime), fragmentPos);
                fragment.velocity = [
                    combinedVelocity[0] - offsetDirection[0] * 0.8,
                    combinedVelocity[1] - offsetDirection[1] * 0.4,
                    combinedVelocity[2] - offsetDirection[2] * 0.8
                ];
                fragment.parentIds = [mol1.id, mol2.id];
                fragment.reactionType = 'emission';
                fragment.setReactionTime(now);
                products.push(fragment);
            }
        } else {
            const sum = mol1.number + mol2.number;
            const diff = Math.abs(mol1.number - mol2.number);
            const fragmentPrime = smallestPrimeFactor(Math.max(2, diff || sum));

            const productA = this.clampNumber(Math.floor(sum * 0.5));
            const productB = this.clampNumber(Math.max(2, diff + fragmentPrime));

            const tangent = normalizeVector(cross(combinedVelocity, [0, 1, 0]));
            const tangentScale = 0.45;

            const first = new PrimeMolecule(productA, [
                center[0] + tangent[0] * 0.4,
                center[1] + tangent[1] * 0.4,
                center[2] + tangent[2] * 0.4
            ]);
            first.velocity = [
                combinedVelocity[0] + tangent[0] * tangentScale + randomNormal() * 0.08,
                combinedVelocity[1] + tangent[1] * tangentScale * 0.6 + randomNormal() * 0.05,
                combinedVelocity[2] + tangent[2] * tangentScale + randomNormal() * 0.08
            ];
            first.parentIds = [mol1.id, mol2.id];
            first.reactionType = 'exchange';
            first.setReactionTime(now);

            const second = new PrimeMolecule(productB, [
                center[0] - tangent[0] * 0.4,
                center[1] - tangent[1] * 0.4,
                center[2] - tangent[2] * 0.4
            ]);
            second.velocity = [
                combinedVelocity[0] - tangent[0] * tangentScale + randomNormal() * 0.08,
                combinedVelocity[1] - tangent[1] * tangentScale * 0.6 + randomNormal() * 0.05,
                combinedVelocity[2] - tangent[2] * tangentScale + randomNormal() * 0.08
            ];
            second.parentIds = [mol1.id, mol2.id];
            second.reactionType = 'exchange';
            second.setReactionTime(now);

            products.push(first, second);

            if (Math.random() < 0.22) {
                const catalyst = new PrimeMolecule(this.clampNumber(fragmentPrime), center);
                catalyst.velocity = [
                    combinedVelocity[0] + randomNormal() * 0.15,
                    combinedVelocity[1] + randomNormal() * 0.08,
                    combinedVelocity[2] + randomNormal() * 0.15
                ];
                catalyst.parentIds = [mol1.id, mol2.id];
                catalyst.reactionType = 'catalyst';
                catalyst.setReactionTime(now);
                products.push(catalyst);
            }
        }

        return products;
    }

    clampNumber(value) {
        return clamp(Math.floor(Math.max(2, value)), 2, this.maxNumber);
    }

    integrateMolecule(molecule, force, dt) {
        const mass = Math.max(1.0, molecule.mass);
        const thermalSigma = this.constants.thermalStrength * this.temperature * Math.sqrt(dt) / Math.sqrt(mass + 1.0);
        const inertialScale = dt / Math.max(1.0, Math.sqrt(mass));

        for (let axis = 0; axis < 3; axis++) {
            molecule.velocity[axis] += randomNormal() * thermalSigma;
            molecule.velocity[axis] += (force[axis] / mass) * dt;
            molecule.velocity[axis] += molecule.inertial_force[axis] * inertialScale;
        }

        const dragMultiplier = Math.max(0.0, 1.0 - this.constants.fluidDrag * dt);
        molecule.velocity[0] *= dragMultiplier;
        molecule.velocity[1] *= dragMultiplier;
        molecule.velocity[2] *= dragMultiplier;

        molecule.velocity[0] *= this.constants.damping;
        molecule.velocity[1] *= this.constants.damping;
        molecule.velocity[2] *= this.constants.damping;

        const speed = vectorMagnitude(molecule.velocity);
        if (speed > this.constants.maxSpeed) {
            const scale = this.constants.maxSpeed / speed;
            molecule.velocity[0] *= scale;
            molecule.velocity[1] *= scale;
            molecule.velocity[2] *= scale;
        }

        molecule.position[0] += molecule.velocity[0] * dt;
        molecule.position[1] += molecule.velocity[1] * dt;
        molecule.position[2] += molecule.velocity[2] * dt;

        this.enforceBoundaries(molecule);
    }

    enforceBoundaries(molecule) {
        const halfSize = this.size * 0.5;
        const rebound = -0.82;

        for (let axis = 0; axis < 3; axis++) {
            if (molecule.position[axis] < -halfSize) {
                molecule.position[axis] = -halfSize;
                molecule.velocity[axis] *= rebound;
            } else if (molecule.position[axis] > halfSize) {
                molecule.position[axis] = halfSize;
                molecule.velocity[axis] *= rebound;
            }
        }
    }

    step() {
        if (this.isPaused || this.molecules.length === 0) {
            return [];
        }

        this.stepCount += 1;
        const dt = this.computeDt();
        this.temperature = clamp(this.temperature + Math.sin(this.stepCount * 0.01) * 0.002, 0.05, 4.0);

        for (const molecule of this.molecules) {
            this.updateInertialForce(molecule, dt);
        }

        const forces = this.molecules.map(() => [0, 0, 0]);
        const removedIndices = new Set();
        const producedMolecules = [];

        const interactionDistanceSq = this.constants.interactionDistance * this.constants.interactionDistance;
        const { grid, coordinates } = this.buildSpatialHash();

        for (let i = 0; i < this.molecules.length; i++) {
            if (removedIndices.has(i)) {
                continue;
            }

            const mol1 = this.molecules[i];
            const [cx, cy, cz] = coordinates[i];

            for (let dx = -1; dx <= 1; dx++) {
                for (let dy = -1; dy <= 1; dy++) {
                    for (let dz = -1; dz <= 1; dz++) {
                        const neighborKey = `${cx + dx}|${cy + dy}|${cz + dz}`;
                        const neighbors = grid.get(neighborKey);
                        if (!neighbors) {
                            continue;
                        }

                        for (const j of neighbors) {
                            if (j <= i || removedIndices.has(j) || removedIndices.has(i)) {
                                continue;
                            }

                            const mol2 = this.molecules[j];
                            const delta = [
                                mol2.position[0] - mol1.position[0],
                                mol2.position[1] - mol1.position[1],
                                mol2.position[2] - mol1.position[2]
                            ];
                            const distanceSq = dot(delta, delta);
                            if (distanceSq > interactionDistanceSq) {
                                continue;
                            }

                            const safeDistance = Math.max(Math.sqrt(distanceSq), this.constants.minDistance);
                            const direction = [
                                delta[0] / safeDistance,
                                delta[1] / safeDistance,
                                delta[2] / safeDistance
                            ];

                            const pairForce = this.calculatePairForce(mol1, mol2, direction, safeDistance);
                            forces[i][0] += pairForce[0];
                            forces[i][1] += pairForce[1];
                            forces[i][2] += pairForce[2];
                            forces[j][0] -= pairForce[0];
                            forces[j][1] -= pairForce[1];
                            forces[j][2] -= pairForce[2];

                            const products = this.maybeReact(mol1, mol2, safeDistance, dt, performance.now());
                            if (products && products.length > 0) {
                                removedIndices.add(i);
                                removedIndices.add(j);
                                producedMolecules.push(...products);
                                this.reactionCount += 1;
                            }
                        }
                    }
                }
            }
        }

        for (let i = 0; i < this.molecules.length; i++) {
            if (removedIndices.has(i)) {
                continue;
            }
            this.integrateMolecule(this.molecules[i], forces[i], dt);
        }

        const removedIds = [];
        if (removedIndices.size > 0) {
            const survivors = [];
            for (let i = 0; i < this.molecules.length; i++) {
                const molecule = this.molecules[i];
                if (removedIndices.has(i)) {
                    removedIds.push(molecule.id);
                } else {
                    survivors.push(molecule);
                }
            }
            this.molecules = survivors;
        }

        for (const molecule of producedMolecules) {
            molecule.id = `mol-${this.nextId++}`;
            this.enforceBoundaries(molecule);
            this.molecules.push(molecule);
        }

        if (this.molecules.length > this.constants.maxMolecules) {
            const overflow = this.molecules.length - this.constants.maxMolecules;
            const trimmed = this.molecules.splice(this.constants.maxMolecules, overflow);
            for (const removed of trimmed) {
                removedIds.push(removed.id);
            }
        }

        return removedIds;
    }

    serializeMolecule(molecule) {
        return {
            id: molecule.id,
            number: molecule.number,
            position: [...molecule.position],
            velocity: [...molecule.velocity],
            prime_factors: { ...molecule.prime_factors },
            mass: molecule.mass,
            charge: molecule.charge,
            color: [...molecule.color],
            angularVelocity: [...molecule.angularVelocity],
            lastReactionTime: molecule.lastReactionTime,
            parentIds: molecule.parentIds ? [...molecule.parentIds] : [],
            reactionType: molecule.reactionType || null,
            inertial_force: [...molecule.inertial_force],
            even_index_weight: molecule.even_index_weight,
            odd_index_weight: molecule.odd_index_weight,
            fractal_signature: [...molecule.fractal_signature]
        };
    }

    createUpdatePayload(removedIds = []) {
        return {
            type: 'update',
            temperature: this.temperature,
            reactionCount: this.reactionCount,
            molecules: {
                molecules: this.molecules.map(molecule => this.serializeMolecule(molecule)),
                removedIds
            }
        };
    }
}

let simulation = null;
let stepInProgress = false;
let queuedStep = false;

function sendUpdate(removedIds = []) {
    if (!simulation) {
        return;
    }
    postMessage(simulation.createUpdatePayload(removedIds));
}

function runStep() {
    if (!simulation) {
        return;
    }
    if (stepInProgress) {
        queuedStep = true;
        return;
    }

    stepInProgress = true;
    try {
        const removedIds = simulation.step();
        sendUpdate(removedIds);
    } catch (error) {
        postMessage({
            type: 'error',
            message: `Simulation step failed: ${error.message}`,
            stack: error.stack
        });
    } finally {
        stepInProgress = false;
        if (queuedStep) {
            queuedStep = false;
            queueMicrotask(runStep);
        }
    }
}

globalThis.onmessage = function onWorkerMessage(event) {
    const message = event.data || {};

    try {
        switch (message.type) {
            case 'init': {
                const size = Number(message.size) || 100;
                const moleculeCount = Number(message.moleculeCount) || 200;
                const maxNumber = Number(message.maxNumber) || 2000;
                const timeScale = Number(message.timeScale) || 1.0;

                simulation = new PrimeChemistrySimulation(size, moleculeCount, maxNumber, timeScale);
                stepInProgress = false;
                queuedStep = false;
                sendUpdate([]);
                break;
            }

            case 'step':
                runStep();
                break;

            case 'set_temperature':
                if (simulation) {
                    simulation.setTemperature(message.value);
                }
                break;

            case 'set_timescale':
                if (simulation) {
                    simulation.setTimeScale(message.value);
                }
                break;

            case 'set_pause':
                if (simulation) {
                    simulation.setPauseState(message.isPaused);
                    if (simulation.isPaused) {
                        sendUpdate([]);
                    }
                }
                break;

            case 'add_molecules':
                if (simulation) {
                    simulation.addRandomMolecules(message.count || 20);
                    sendUpdate([]);
                }
                break;

            case 'set_visualization':
                // UI-only option; accepted for compatibility.
                break;

            case 'cleanup':
                simulation = null;
                stepInProgress = false;
                queuedStep = false;
                postMessage({ type: 'cleanup_complete' });
                break;

            default:
                postMessage({ type: 'error', message: `Unknown message type: ${message.type}` });
        }
    } catch (error) {
        postMessage({
            type: 'error',
            message: `Worker error while handling '${message.type}': ${error.message}`,
            stack: error.stack
        });
    }
};
