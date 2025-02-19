import * as THREE from 'three';
import {
    OrbitControls
} from 'three/addons/controls/OrbitControls.js';

let scene, camera, renderer, controls;

let simulationData = {
    molecules: [],
    temperature: 0
}; // Holds the latest data

let paused = false;
let autoStep = true;
let showVectors = false;
const infoDiv = document.getElementById('info');
const moleculeInfoDiv = document.getElementById('moleculeInfo');
let selectedMolecule = null;

// Store molecule meshes with their IDs
const moleculeMeshes = new Map();
let simulationWorker;

function initializeSimulation(size = 10, moleculeCount = 10, maxNumber = 100, timeScale = 1) {
    // Termina il worker esistente se presente
    if (simulationWorker) {
        simulationWorker.terminate();
    }

    try {
        // Crea il nuovo worker con gestione errori
        console.log("Inizializzazione worker...");
        if (false) {
            simulationWorker = new Worker(new URL('./worker.js', import.meta.url), {
                type: 'module'
            });
        } else {
            simulationWorker = new Worker(new URL('./minimal-worker.js', import.meta.url), {
                type: 'module'
            });
        }

        // Aggiungi gestore errori
        simulationWorker.onerror = function(error) {
            console.error("Worker error:", error);
            alert("Errore nel worker di simulazione: " + error.message);
        };

        // Imposta la gestione messaggi
        simulationWorker.onmessage = function(event) {
            console.log("Ricevuto messaggio dal worker:", event.data);
            handleWorkerMessage(event);
        };

        // Invia il messaggio di inizializzazione
        console.log("Invio messaggio init al worker");
        simulationWorker.postMessage({
            type: 'init',
            size: size,
            moleculeCount: moleculeCount * size,
            maxNumber: maxNumber,
            timeScale: timeScale
        });

        console.log("Messaggio init inviato");
    } catch (error) {
        console.error("Errore nell'inizializzazione del worker:", error);
        alert("Impossibile inizializzare la simulazione: " + error.message);

        // Crea molecole di test per debug
        createTestMolecules();
    }
}

// Aggiungi un timer per verificare che il worker stia funzionando
function checkWorkerStatus() {
    console.log("Verifica stato worker...");
    if (simulationData.molecules.length === 0) {
        console.warn("Nessuna molecola ricevuta dal worker dopo 5 secondi");
        console.log("Creo molecole di test...");
        createTestMolecules();
    }
}

// Funzione per creare molecole di test
function createTestMolecules() {
    console.log("Creazione molecole di test");
    const testMolecules = [];

    for (let i = 0; i < 10; i++) {
        testMolecules.push({
            id: `test-${i}`,
            number: 10 + i,
            position: [
                (Math.random() - 0.5) * 8,
                (Math.random() - 0.5) * 8,
                (Math.random() - 0.5) * 8
            ],
            velocity: [0, 0, 0],
            prime_factors: {
                2: 1,
                5: 1
            },
            mass: 5,
            charge: 0,
            color: [Math.random(), Math.random(), Math.random()],
            angularVelocity: [0, 0, 0],
            lastReactionTime: -1
        });
    }

    simulationData.molecules = testMolecules;
    updateMoleculeMeshes();
}

function handleWorkerMessage(event) {
    console.log("Received worker message:", event.data);

    if (event.data.type === 'update') {
        // Update temperature
        simulationData.temperature = event.data.temperature;

        // Debug message structure
        console.log("Molecules data structure:", event.data.molecules);

        // Process molecule updates - check structure carefully
        if (event.data.molecules) {
            if (Array.isArray(event.data.molecules)) {
                // Old format - direct array
                console.log("Using older API format (direct array)");
                simulationData.molecules = event.data.molecules;
            } else if (event.data.molecules.molecules) {
                // New format with molecules and removedIds
                console.log("Using new API format (molecules + removedIds)");
                simulationData.molecules = event.data.molecules.molecules;

                // Handle removed molecules
                if (event.data.molecules.removedIds) {
                    console.log("Processing removedIds:", event.data.molecules.removedIds);
                    event.data.molecules.removedIds.forEach(id => {
                        const mesh = moleculeMeshes.get(id);
                        if (mesh) {
                            scene.remove(mesh);
                            if (mesh.userData.halo) {
                                scene.remove(mesh.userData.halo);
                            }
                            if (mesh.userData.velocityArrow) {
                                scene.remove(mesh.userData.velocityArrow);
                            }
                            moleculeMeshes.delete(id);
                        }
                    });
                }
            } else {
                console.error("Unknown molecule data format:", event.data.molecules);
            }

            console.log(`Molecules to render: ${simulationData.molecules.length}`);

            // Update molecule meshes
            updateMoleculeMeshes();
        } else {
            console.warn("No molecules data in update message");
        }
    } else if (event.data.type === 'cleanup_complete') {
        console.log('Worker cleanup completed');
    }
}

// Modifica la funzione updateMoleculeMeshes nel main.js
// Modifica la funzione updateMoleculeMeshes nel main.js
function updateMoleculeMeshes() {
    // Process each molecule
    simulationData.molecules.forEach(molData => {
        let mesh = moleculeMeshes.get(molData.id);

        // Create new mesh if it doesn't exist
        if (!mesh) {
            // Usa la nuova funzione per creare mesh avanzati
            mesh = createEnhancedMoleculeMesh(molData);
            scene.add(mesh);
            moleculeMeshes.set(molData.id, mesh);
        } else {
            // Aggiorna i dati
            mesh.userData.moleculeData = molData;

            // Aggiorna posizione
            mesh.position.set(...molData.position);

            // Aggiorna rotazione se presente
            if (molData.angularVelocity) {
                const axis = new THREE.Vector3(...molData.angularVelocity).normalize();
                const angle = new THREE.Vector3(...molData.angularVelocity).length();
                const quaternion = new THREE.Quaternion().setFromAxisAngle(axis, angle);
                mesh.quaternion.multiplyQuaternions(quaternion, mesh.quaternion);
            }

            // Aggiungi o aggiorna vettori di velocità se richiesto
            updateVelocityVector(mesh, molData, showVectors);
        }

        // Aggiorna effetti di reazione (alone, particelle, glow)
        updateHalo(mesh, molData);
    });
}

// Helper per i vettori di velocità
function updateVelocityVector(mesh, molData, show) {
    // Rimuovi il vettore esistente se presente
    if (mesh.userData.velocityArrow) {
        scene.remove(mesh.userData.velocityArrow);
        mesh.userData.velocityArrow = null;
    }

    // Aggiungi nuovo vettore se richiesto
    if (show && molData.velocity) {
        const dir = new THREE.Vector3(...molData.velocity);
        const origin = new THREE.Vector3(...molData.position);
        const length = dir.length();
        const hex = 0xffff00; // Giallo

        const arrowHelper = new THREE.ArrowHelper(
            dir.normalize(),
            origin,
            length,
            hex
        );

        scene.add(arrowHelper);
        mesh.userData.velocityArrow = arrowHelper;
    }
}

// Modifica la funzione animate nel main.js
function animate(timestamp) {
    // Calcola delta time per animazioni
    if (!lastFrameTime) lastFrameTime = timestamp;
    const deltaTime = timestamp - lastFrameTime;
    lastFrameTime = timestamp;

    requestAnimationFrame(animate);
    controls.update();

    // Request a simulation step from the worker (if not paused)
    if (autoStep && !paused) {
        simulationWorker.postMessage({
            type: 'step'
        });
    }

    // Aggiorna gli effetti delle particelle
    updateParticleEffects(deltaTime);

    // Aggiorna le visualizzazioni avanzate
    updateEnhancedVisuals(deltaTime, timestamp);

    // Update info display
    infoDiv.innerText = `Molecules: ${simulationData.molecules.length}\n` +
        `Temperature: ${simulationData.temperature.toFixed(2)}\n` +
        `${paused ? 'PAUSED' : 'RUNNING'}\n` +
        `${autoStep ? 'AUTO' : 'MANUAL'}`;
    render();
}

// Aggiungi variabile per il timestamp
let lastFrameTime = 0;

// Pulisci effetti quando necessario (cleanup)
function cleanupSimulation() {
    if (simulationWorker) {
        // Request cleanup from worker
        simulationWorker.postMessage({
            type: 'cleanup'
        });

        // Set a timeout to force termination if worker doesn't respond
        setTimeout(() => {
            simulationWorker.terminate();
            simulationWorker = null;
        }, 500);
    }

    // Pulisci effetti e mesh
    moleculeMeshes.forEach(mesh => {
        scene.remove(mesh);

        // Pulisci anche tutti gli effetti associati
        const effects = [
            'halo', 'glowEffect', 'particleSystem',
            'velocityArrow', 'elementalGlow', 'orbitalRing',
            'atomicStructure'
        ];

        effects.forEach(effectType => {
            if (mesh.userData[effectType]) {
                scene.remove(mesh.userData[effectType]);
            }
        });
    });

    moleculeMeshes.clear();
}

// Aggiungi informazioni miglioratei sulle molecole nel pannello info
function onClick(event) {
    event.preventDefault();

    const mouse = new THREE.Vector2();
    mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
    mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;

    const raycaster = new THREE.Raycaster();
    raycaster.setFromCamera(mouse, camera);

    // Convert Map to array for raycasting
    const meshArray = Array.from(moleculeMeshes.values());
    const intersects = raycaster.intersectObjects(meshArray);

    if (intersects.length > 0) {
        // Get the closest intersection
        const closestIntersection = intersects[0];
        let molIntersected = closestIntersection.object;

        // If we clicked on a part of a complex molecule (like ring or atom)
        // find the parent molecule
        if (!molIntersected.userData.moleculeData) {
            for (const mesh of meshArray) {
                if (mesh.userData.elementalGlow === molIntersected ||
                    mesh.userData.orbitalRing === molIntersected ||
                    (mesh.userData.atomicStructure &&
                        mesh.userData.atomicStructure.children.includes(molIntersected))) {
                    molIntersected = mesh;
                    break;
                }
            }
        }

        // Make sure we have molecule data
        if (molIntersected.userData.moleculeData) {
            selectedMolecule = molIntersected.userData.moleculeData;
            const visualProps = molIntersected.userData.visualProperties;

            // Display enhanced molecule info
            moleculeInfoDiv.style.display = 'block';
            moleculeInfoDiv.innerHTML = `
          <h3>Molecule ${selectedMolecule.id}</h3>
          <div class="molecule-info-container">
            <div class="basic-properties">
              <p><strong>Number:</strong> ${selectedMolecule.number}</p>
              <p><strong>Prime Factors:</strong> ${formatPrimeFactors(selectedMolecule.prime_factors)}</p>
              <p><strong>Mass:</strong> ${selectedMolecule.mass.toFixed(2)}</p>
              <p><strong>Charge:</strong> ${selectedMolecule.charge.toFixed(2)}</p>
              <p><strong>Velocity:</strong> ${vectorMagnitude(selectedMolecule.velocity).toFixed(2)}</p>
            </div>
            
            <div class="advanced-properties">
              <p><strong>Category:</strong> ${visualProps.category}</p>
              <p><strong>Stability:</strong> ${visualProps.isStable ? 'Stable' : 'Unstable'}</p>
              <p><strong>Energy:</strong> ${visualProps.energy.toFixed(2)}</p>
              <p><strong>Complexity:</strong> ${visualProps.complexity}</p>
              <p><strong>Largest Prime:</strong> ${visualProps.largestPrime}</p>
            </div>
          </div>
          
          <div class="position-data">
            <p><strong>Position:</strong> 
              (${selectedMolecule.position[0].toFixed(2)}, 
               ${selectedMolecule.position[1].toFixed(2)}, 
               ${selectedMolecule.position[2].toFixed(2)})</p>
          </div>
        `;

            // Evidenzia la molecola selezionata
            highlightSelectedMolecule(molIntersected);
        }
    } else {
        selectedMolecule = null;
        moleculeInfoDiv.style.display = 'none';
        clearMoleculeHighlight();
    }
}

// Funzione per evidenziare la molecola selezionata
function highlightSelectedMolecule(mesh) {
    // Rimuovi evidenziazione precedente
    clearMoleculeHighlight();

    // Crea un effetto di evidenziazione
    const highlightGeometry = new THREE.SphereGeometry(
        mesh.geometry.parameters.radius * 1.1,
        32,
        32
    );

    const highlightMaterial = new THREE.MeshBasicMaterial({
        color: 0xffffff,
        wireframe: true,
        transparent: true,
        opacity: 0.3
    });

    const highlightMesh = new THREE.Mesh(highlightGeometry, highlightMaterial);
    highlightMesh.position.copy(mesh.position);
    scene.add(highlightMesh);

    // Salva l'evidenziazione
    currentHighlight = highlightMesh;
}

// Variabile per tenere traccia dell'evidenziazione corrente
let currentHighlight = null;

// Funzione per rimuovere l'evidenziazione
function clearMoleculeHighlight() {
    if (currentHighlight) {
        scene.remove(currentHighlight);
        currentHighlight = null;
    }
}

// Aggiungi stili CSS per il div di informazioni molecola
function addMoleculeInfoStyles() {
    const style = document.createElement('style');
    style.textContent = `
      #moleculeInfo {
        position: absolute;
        bottom: 20px;
        right: 20px;
        background-color: rgba(0, 0, 0, 0.7);
        color: white;
        padding: 15px;
        border-radius: 10px;
        font-family: Arial, sans-serif;
        max-width: 400px;
        display: none;
      }
      
      #moleculeInfo h3 {
        margin-top: 0;
        color: #88ccff;
        border-bottom: 1px solid #555;
        padding-bottom: 5px;
      }
      
      .molecule-info-container {
        display: flex;
        justify-content: space-between;
      }
      
      .basic-properties, .advanced-properties {
        flex: 1;
      }
      
      .position-data {
        margin-top: 10px;
        font-size: 0.9em;
        color: #aaa;
      }
      
      #moleculeInfo p {
        margin: 5px 0;
      }
      
      #moleculeInfo strong {
        color: #aaffaa;
      }
    `;
    document.head.appendChild(style);
}

function init() {
    const spaceDimension = 10;

    // Scene setup
    scene = new THREE.Scene();
    camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
    camera.position.z = spaceDimension;

    renderer = new THREE.WebGLRenderer({
        antialias: true
    });
    renderer.setSize(window.innerWidth, window.innerHeight);
    document.body.appendChild(renderer.domElement);

    // Orbit Controls
    controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.1;

    // Lighting
    const ambientLight = new THREE.AmbientLight(0x404040);
    scene.add(ambientLight);
    const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
    directionalLight.position.set(1, 1, 1).normalize();
    scene.add(directionalLight);

    // Axes helper (for debugging)
    const axesHelper = new THREE.AxesHelper(5);
    scene.add(axesHelper);

    // Boundary box
    const boxSize = spaceDimension; // Same as Python simulation size
    const boxGeometry = new THREE.BoxGeometry(boxSize, boxSize, boxSize);
    const wireframe = new THREE.WireframeGeometry(boxGeometry);
    const line = new THREE.LineSegments(wireframe);
    line.material.color.set(0x808080);
    scene.add(line);

    initializeSimulation(spaceDimension);

    addMoleculeInfoStyles();

    // Event Listeners
    window.addEventListener('resize', onWindowResize, false);
    document.addEventListener('keydown', onKeyDown, false);
    renderer.domElement.addEventListener('click', onClick, false);

    animate();
}

function render() {
    renderer.render(scene, camera);
}

function onWindowResize() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
}

function onKeyDown(event) {
    switch (event.code) {
        case 'Space':
            paused = !paused;
            break;
        case 'KeyA':
            autoStep = !autoStep;
            break;
        case 'KeyV':
            showVectors = !showVectors;
            updateMoleculeMeshes();
            break;
        case 'KeyR': // Reset Camera
            controls.reset();
            break;
        case 'ArrowRight':
            if (paused) {
                simulationWorker.postMessage({
                    type: 'step'
                });
            }
            break;
    }
}

function formatPrimeFactors(factors) {
    let str = "";
    for (const prime in factors) {
        str += `${prime}^${factors[prime]} `;
    }
    return str.trim();
}

function vectorMagnitude(vec) {
    return Math.sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

// Attach cleanup to page events
window.addEventListener('beforeunload', cleanupSimulation);
window.addEventListener('visibilitychange', () => {
    if (document.visibilityState === 'hidden') {
        // Pause simulation when tab is not visible
        paused = true;
    } else {
        // Keep it paused if it was manually paused
        if (autoStep) {
            paused = false;
        }
    }
});

///
///
///

// Migliore gestione degli effetti di reazione
function updateHalo(mesh, molData) {
    // Remove old halo and effects if exist
    if (mesh.userData.halo) {
        scene.remove(mesh.userData.halo);
        mesh.userData.halo = null;
    }

    if (mesh.userData.glowEffect) {
        scene.remove(mesh.userData.glowEffect);
        mesh.userData.glowEffect = null;
    }

    if (mesh.userData.particleSystem) {
        scene.remove(mesh.userData.particleSystem);
        mesh.userData.particleSystem = null;
    }

    // Add reaction effects if recently reacted
    if (molData.lastReactionTime > 0) {
        const timeSinceReaction = performance.now() - molData.lastReactionTime;
        const maxEffectDuration = 2000; // Effects last 2 seconds

        if (timeSinceReaction < maxEffectDuration) {
            // Calculate effect intensity (1.0 to 0.0 over duration)
            const intensity = 1 - (timeSinceReaction / maxEffectDuration);
            addReactionEffects(mesh, molData, intensity);
        }
    }
}

function addReactionEffects(mesh, molData, intensity) {
    // Extract base properties
    const baseColor = new THREE.Color(...molData.color);
    const position = new THREE.Vector3(...molData.position);
    const size = mesh.geometry.parameters.radius;

    // 1. Outer halo (pulsating sphere)
    const haloSize = size * (2.0 + Math.sin(performance.now() * 0.01) * 0.3);
    const haloGeometry = new THREE.SphereGeometry(haloSize, 32, 16);
    const haloMaterial = new THREE.MeshBasicMaterial({
        color: baseColor.clone().offsetHSL(0, 0, 0.2), // Slightly lighter
        transparent: true,
        opacity: 0.3 * intensity,
        side: THREE.BackSide,
        blending: THREE.AdditiveBlending
    });

    const halo = new THREE.Mesh(haloGeometry, haloMaterial);
    halo.position.copy(position);
    scene.add(halo);
    mesh.userData.halo = halo;

    // 2. Emission particles for more complex reactions
    if (Object.keys(molData.prime_factors).length > 2) {
        createParticleSystem(mesh, molData, intensity);
    }

    // 3. Inner glow effect for all reactions
    addGlowEffect(mesh, molData, intensity);
}

function createParticleSystem(mesh, molData, intensity) {
    const particleCount = 15 + Math.floor(intensity * 30);
    const particleGeometry = new THREE.BufferGeometry();
    const particleMaterial = new THREE.PointsMaterial({
        color: new THREE.Color(...molData.color).offsetHSL(0, 0.5, 0.2),
        size: 0.1 + 0.1 * intensity,
        transparent: true,
        opacity: 0.6 * intensity,
        blending: THREE.AdditiveBlending,
        sizeAttenuation: true
    });

    // Create particle positions
    const positions = new Float32Array(particleCount * 3);
    const velocities = [];
    const size = mesh.geometry.parameters.radius;

    for (let i = 0; i < particleCount; i++) {
        // Random positions within a sphere
        const theta = Math.random() * Math.PI * 2;
        const phi = Math.acos(2 * Math.random() - 1);
        const radius = size * (1 + Math.random() * intensity);

        const x = radius * Math.sin(phi) * Math.cos(theta);
        const y = radius * Math.sin(phi) * Math.sin(theta);
        const z = radius * Math.cos(phi);

        positions[i * 3] = mesh.position.x + x;
        positions[i * 3 + 1] = mesh.position.y + y;
        positions[i * 3 + 2] = mesh.position.z + z;

        // Store velocity for animation
        velocities.push(new THREE.Vector3(x, y, z).normalize().multiplyScalar(0.01 + 0.02 * Math.random()));
    }

    particleGeometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));

    const particleSystem = new THREE.Points(particleGeometry, particleMaterial);
    particleSystem.userData = {
        velocities: velocities,
        creationTime: performance.now(),
        update: function(deltaTime) {
            const positions = particleSystem.geometry.attributes.position.array;
            const elapsed = performance.now() - this.creationTime;
            const lifePercentage = Math.min(1, elapsed / 2000);

            // Update particle positions based on velocity and time
            for (let i = 0; i < particleCount; i++) {
                positions[i * 3] += this.velocities[i].x * deltaTime;
                positions[i * 3 + 1] += this.velocities[i].y * deltaTime;
                positions[i * 3 + 2] += this.velocities[i].z * deltaTime;

                // Fade out based on life percentage
                particleSystem.material.opacity = (1 - lifePercentage) * 0.6 * intensity;
            }

            particleSystem.geometry.attributes.position.needsUpdate = true;
        }
    };

    scene.add(particleSystem);
    mesh.userData.particleSystem = particleSystem;
}

function addGlowEffect(mesh, molData, intensity) {
    // Create a slightly larger, emissive sphere inside the molecule
    const glowGeometry = new THREE.SphereGeometry(
        mesh.geometry.parameters.radius * 0.9,
        32, 16
    );

    // Determine glow color based on prime factors
    const primeCount = Object.keys(molData.prime_factors).length;
    let glowColor;

    if (primeCount === 1) {
        // Pure prime - bright blue glow
        glowColor = new THREE.Color(0.2, 0.4, 1.0);
    } else if (primeCount === 2) {
        // Two primes - purple glow
        glowColor = new THREE.Color(0.8, 0.2, 1.0);
    } else {
        // Multiple primes - orange-red glow
        glowColor = new THREE.Color(1.0, 0.4, 0.1);
    }

    const glowMaterial = new THREE.MeshBasicMaterial({
        color: glowColor,
        transparent: true,
        opacity: 0.7 * intensity,
        blending: THREE.AdditiveBlending
    });

    const glowMesh = new THREE.Mesh(glowGeometry, glowMaterial);
    glowMesh.position.copy(mesh.position);
    scene.add(glowMesh);
    mesh.userData.glowEffect = glowMesh;
}

// Aggiungi questa funzione all'animate loop per aggiornare gli effetti particellari
function updateParticleEffects(deltaTime) {
    moleculeMeshes.forEach(mesh => {
        if (mesh.userData.particleSystem) {
            mesh.userData.particleSystem.userData.update(deltaTime);
        }
    });
}

// Sistema avanzato di visualizzazione molecole
function createEnhancedMoleculeMesh(molData) {
    // Analizza le proprietà della molecola per determinarne l'aspetto
    const visualProperties = analyzeMoleculeProperties(molData);

    // Crea il materiale base con le caratteristiche calcolate
    const material = createMoleculeMaterial(visualProperties);

    // Crea la geometria appropriata
    const geometry = createMoleculeGeometry(visualProperties);

    // Crea il mesh
    const mesh = new THREE.Mesh(geometry, material);
    mesh.position.set(...molData.position);
    mesh.userData.moleculeData = molData;
    mesh.userData.visualProperties = visualProperties;

    // Aggiungi effetti extra in base alla composizione della molecola
    addExtraDetails(mesh, visualProperties);

    return mesh;
}

function analyzeMoleculeProperties(molData) {
    const properties = {
        size: 0.2 + 0.1 * Math.log2(molData.number),
        baseColor: new THREE.Color(...molData.color),
        primeFactors: Object.keys(molData.prime_factors).map(Number),
        complexity: Object.keys(molData.prime_factors).length,
        largestPrime: Math.max(...Object.keys(molData.prime_factors).map(Number)),
        mass: molData.mass,
        charge: molData.charge,
        isStable: isStableMolecule(molData),
        energy: calculateMoleculeEnergy(molData),
        category: categorizeMolecule(molData)
    };

    // Calcola luminosità e saturazione basate sulle proprietà
    properties.luminosity = calculateLuminosity(properties);
    properties.saturation = calculateSaturation(properties);
    properties.roughness = calculateRoughness(properties);
    properties.metalness = calculateMetalness(properties);

    return properties;
}

function isStableMolecule(molData) {
    // Molecole più stabili: numeri con pochi fattori primi, o potenze di un primo
    const factorCount = Object.keys(molData.prime_factors).length;
    const isPrimePower = factorCount === 1;

    if (isPrimePower) return true;

    // Verifica se i numeri hanno fattori primi "compatibili"
    const primes = Object.keys(molData.prime_factors).map(Number);
    if (primes.includes(2) && primes.includes(3)) return true;
    if (primes.includes(5) && primes.includes(2)) return true;

    return factorCount <= 2; // Pochi fattori = più stabile
}

function calculateMoleculeEnergy(molData) {
    // L'energia è correlata alla somma dei fattori primi
    let energy = 0;
    for (const [prime, exponent] of Object.entries(molData.prime_factors)) {
        energy += Number(prime) * exponent;
    }
    return energy / (1 + Math.log(molData.number));
}

function categorizeMolecule(molData) {
    // Categorizzazione basata sulla composizione dei fattori primi
    const primes = Object.keys(molData.prime_factors).map(Number);

    if (primes.length === 1) {
        return 'elemental'; // Elementare (numero primo)
    }

    if (primes.every(p => p < 10)) {
        return 'stable'; // Stabile (fattori primi piccoli)
    }

    if (primes.some(p => p > 50)) {
        return 'exotic'; // Esotico (fattori primi grandi)
    }

    if (primes.length > 3) {
        return 'complex'; // Complesso (molti fattori)
    }

    return 'standard'; // Standard
}

function calculateLuminosity(props) {
    // Molecole più leggere sono più luminose
    const baseLuminosity = 0.5;

    if (props.category === 'elemental') {
        // Elementi puri sono più luminosi
        return Math.min(0.9, baseLuminosity + 0.3);
    }

    if (props.category === 'exotic') {
        // Elementi esotici hanno luminosità variabile
        return baseLuminosity + 0.2 * Math.sin(props.largestPrime * 0.1);
    }

    if (props.isStable) {
        // Molecole stabili hanno luminosità media-alta
        return baseLuminosity + 0.15;
    }

    // Modifica in base alla carica
    return baseLuminosity + props.charge * 0.1;
}

function calculateSaturation(props) {
    // Molecole con più fattori primi hanno minore saturazione
    const baseSaturation = 0.7;
    const complexityFactor = Math.max(0, 0.2 - props.complexity * 0.05);

    if (props.category === 'elemental') {
        // Elementi puri sono molto saturi
        return Math.min(1.0, baseSaturation + 0.3);
    }

    if (props.category === 'exotic') {
        // Elementi esotici hanno alta saturazione
        return Math.min(1.0, baseSaturation + 0.2);
    }

    return Math.max(0.3, baseSaturation + complexityFactor);
}

function calculateRoughness(props) {
    // Molecole complesse sono più "ruvide"
    const baseRoughness = 0.5;

    if (props.category === 'elemental') {
        // Elementi puri sono lisci
        return baseRoughness - 0.3;
    }

    if (props.category === 'exotic') {
        // Elementi esotici sono molto ruvidi
        return Math.min(1.0, baseRoughness + 0.3);
    }

    return baseRoughness + 0.1 * props.complexity;
}

function calculateMetalness(props) {
    // Molecole con certi numeri primi sono più metalliche
    const hasMetallicPrime = props.primeFactors.some(p => [2, 3, 5, 7, 13, 19].includes(p));

    if (props.category === 'elemental' && hasMetallicPrime) {
        return 0.7; // Alto metallismo per elementi puri metallici
    }

    if (props.category === 'stable' && hasMetallicPrime) {
        return 0.5; // Metallismo medio per composti stabili metallici  
    }

    return 0.2; // Basso metallismo per altri
}

function createMoleculeMaterial(props) {
    // Prepara il colore con luminosità e saturazione calcolate
    const enhancedColor = props.baseColor.clone();

    // Aggiusta HSL
    let hsl = {};
    enhancedColor.getHSL(hsl);
    enhancedColor.setHSL(
        hsl.h,
        props.saturation,
        props.luminosity
    );

    // Definisci emissività in base alla categoria
    const emissiveIntensity = {
        'elemental': 0.3,
        'exotic': 0.5,
        'complex': 0.1,
        'stable': 0.05,
        'standard': 0
    } [props.category];

    const emissiveColor = props.baseColor.clone().offsetHSL(0.1, 0, 0);
    emissiveColor.multiplyScalar(emissiveIntensity);

    let material;

    // Usa materiali diversi in base alla categoria
    if (props.category === 'elemental' || props.category === 'exotic') {
        // MeshStandardMaterial per elementi speciali
        material = new THREE.MeshStandardMaterial({
            color: enhancedColor,
            emissive: emissiveColor,
            roughness: props.roughness,
            metalness: props.metalness,
            envMapIntensity: 1.0,
            emissiveIntensity: emissiveIntensity,
        });
    } else {
        // MeshPhongMaterial per molecole normali (più leggero)
        material = new THREE.MeshPhongMaterial({
            color: enhancedColor,
            emissive: emissiveColor,
            specular: new THREE.Color(0xffffff),
            shininess: 30 * (1 - props.roughness),
            reflectivity: props.metalness
        });
    }

    return material;
}

function createMoleculeGeometry(props) {
    const baseSize = props.size;

    // Scegli la geometria in base alla categoria
    switch (props.category) {
        case 'elemental':
            // Sfere perfette per elementi puri
            return new THREE.SphereGeometry(baseSize, 32, 32);

        case 'exotic':
            // Forme complesse per elementi esotici
            if (props.largestPrime > 50) {
                return new THREE.OctahedronGeometry(baseSize, 1);
            } else {
                return new THREE.IcosahedronGeometry(baseSize, 1);
            }

        case 'complex':
            // Forme irregolari per molecole complesse
            return new THREE.DodecahedronGeometry(baseSize, 0);

        case 'stable':
            // Forme regolari per molecole stabili
            return new THREE.SphereGeometry(baseSize, 24, 24);

        default:
            // Sfere semplici per il resto
            return new THREE.SphereGeometry(baseSize, 16, 16);
    }
}

function addExtraDetails(mesh, props) {
    // Aggiungi effetti visivi extra in base alle proprietà

    if (props.category === 'elemental') {
        // Aggiungi un alone permanente per elementi puri
        addElementalGlow(mesh, props);
    }

    if (props.category === 'exotic') {
        // Aggiungi particelle orbitanti per elementi esotici
        addOrbitalRing(mesh, props);
    }

    if (props.isStable && props.complexity > 1) {
        // Aggiungi struttura atomica per molecole stabili
        addAtomicStructure(mesh, props);
    }
}

function addElementalGlow(mesh, props) {
    // Crea un alone permanente per elementi puri
    const haloSize = props.size * 1.2;
    const haloGeometry = new THREE.SphereGeometry(haloSize, 32, 16);
    const haloMaterial = new THREE.MeshBasicMaterial({
      color: props.baseColor.clone().offsetHSL(0, -0.5, 0.3),
      transparent: true,
      opacity: 0.15,
      side: THREE.BackSide,
      blending: THREE.AdditiveBlending
    });
    
    const halo = new THREE.Mesh(haloGeometry, haloMaterial);
    halo.position.copy(mesh.position);
    
    // Salviamo i parametri di pulsazione nel halo stesso
    halo.userData.pulseParams = {
      baseSize: haloSize,
      speed: 0.001,
      amplitude: 0.1,
      phase: Math.random() * Math.PI * 2
    };
    
    // Funzione di aggiornamento rivista per accedere correttamente a scale
    halo.userData.updateFunction = function(time) {
      const p = this.userData.pulseParams;
      const scale = 1 + Math.sin(time * p.speed + p.phase) * p.amplitude;
      this.scale.set(scale, scale, scale);
    };
    
    scene.add(halo);
    mesh.userData.elementalGlow = halo;
}
  
// Aggiorna la funzione per gestire correttamente l'aggiornamento
function updateEnhancedVisuals(deltaTime, time) {
    moleculeMeshes.forEach(mesh => {
      // Aggiorna effetti elementali - CORRETTO
      if (mesh.userData.elementalGlow) {
        const halo = mesh.userData.elementalGlow;
        halo.position.copy(mesh.position);
        halo.userData.updateFunction.call(halo, time);
      }
      
      // Aggiorna anelli orbitali
      if (mesh.userData.orbitalRing) {
        const ring = mesh.userData.orbitalRing;
        ring.userData.updateFunction.call(ring, deltaTime);
      }
      
      // Aggiorna strutture atomiche 
      if (mesh.userData.atomicStructure) {
        const structure = mesh.userData.atomicStructure;
        structure.userData.updateFunction.call(structure, time);
      }
    });
}
  
  // Aggiorna anche la funzione orbitalRing per coerenza
function addOrbitalRing(mesh, props) {
    // Crea un anello orbitale per elementi esotici
    const ringRadius = props.size * 1.8;
    const ringGeometry = new THREE.TorusGeometry(ringRadius, 0.05, 8, 50);
    
    // Colore complementare
    const ringColor = props.baseColor.clone();
    let hsl = {};
    ringColor.getHSL(hsl);
    ringColor.setHSL((hsl.h + 0.5) % 1, hsl.s, hsl.l);
    
    const ringMaterial = new THREE.MeshBasicMaterial({
      color: ringColor,
      transparent: true,
      opacity: 0.7,
      blending: THREE.AdditiveBlending
    });
    
    const ring = new THREE.Mesh(ringGeometry, ringMaterial);
    
    // Randomizza l'orientamento
    ring.rotation.x = Math.random() * Math.PI;
    ring.rotation.y = Math.random() * Math.PI;
    
    // Aggiungi parametri di rotazione
    ring.userData.rotationSpeed = {
      x: (Math.random() - 0.5) * 0.01,
      y: (Math.random() - 0.5) * 0.01,
      z: (Math.random() - 0.5) * 0.01
    };
    
    // Coerente con le altre funzioni di aggiornamento
    ring.userData.updateFunction = function(deltaTime) {
      this.rotation.x += this.userData.rotationSpeed.x * deltaTime;
      this.rotation.y += this.userData.rotationSpeed.y * deltaTime;
      this.rotation.z += this.userData.rotationSpeed.z * deltaTime;
    };
    
    mesh.add(ring);
    mesh.userData.orbitalRing = ring;
}

function addAtomicStructure(mesh, props) {
    // Crea punti orbitanti che rappresentano la struttura atomica
    const atomRadius = props.size * 0.25;
    const atomGeometry = new THREE.SphereGeometry(atomRadius, 8, 8);
    
    // Crea un gruppo per contenere la struttura
    const atomicStructure = new THREE.Group();
    
    // Aggiungi un atomo per ogni fattore primo
    props.primeFactors.forEach((prime, index) => {
      // Calcola posizione orbitale
      const distance = props.size * 1.2;
      const angleOffset = (index / props.primeFactors.length) * Math.PI * 2;
      
      const x = Math.cos(angleOffset) * distance;
      const y = Math.sin(angleOffset) * distance * 0.5;
      const z = Math.sin(angleOffset + Math.PI/2) * distance * 0.7;
      
      // Colore basato sul numero primo
      const atomColor = new THREE.Color().setHSL(
        (prime % 12) / 12,
        0.8,
        0.6
      );
      
      const atomMaterial = new THREE.MeshBasicMaterial({
        color: atomColor,
        transparent: true,
        opacity: 0.8
      });
      
      const atom = new THREE.Mesh(atomGeometry, atomMaterial);
      atom.position.set(x, y, z);
      
      // Aggiungi parametri di orbita
      atom.userData.orbit = {
        distance: distance,
        speed: 0.0005 * prime % 5,
        offset: angleOffset,
        vertical: (index % 3 - 1) * 0.3
      };
      
      atomicStructure.add(atom);
    });
    
    // Salviamo la funzione di aggiornamento nei userData della struttura atomica
    atomicStructure.userData.updateFunction = function(time) {
      const atoms = this.children;
      for (let i = 0; i < atoms.length; i++) {
        const atom = atoms[i];
        const orbit = atom.userData.orbit;
        const angle = time * orbit.speed + orbit.offset;
        
        atom.position.x = Math.cos(angle) * orbit.distance;
        atom.position.y = Math.sin(angle) * orbit.distance * 0.5 + orbit.vertical * orbit.distance;
        atom.position.z = Math.sin(angle + Math.PI/2) * orbit.distance * 0.7;
      }
    };
    
    mesh.add(atomicStructure);
    mesh.userData.atomicStructure = atomicStructure;
}

///
///
///

init();