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
        if(false){
            simulationWorker = new Worker(new URL('./worker.js', import.meta.url), {
                type: 'module'
            });
        }
        else {
            simulationWorker = new Worker(new URL('./minimal-worker.js', import.meta.url), { type: 'module' });
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

    // Clear scene
    moleculeMeshes.forEach(mesh => {
        scene.remove(mesh);
        if (mesh.userData.halo) scene.remove(mesh.userData.halo);
        if (mesh.userData.velocityArrow) scene.remove(mesh.userData.velocityArrow);
    });

    moleculeMeshes.clear();
}

function updateMoleculeMeshes() {
    // Process each molecule
    simulationData.molecules.forEach(molData => {
        let mesh = moleculeMeshes.get(molData.id);

        // Create new mesh if it doesn't exist
        if (!mesh) {
            const geometry = new THREE.SphereGeometry(0.2 + 0.1 * Math.log2(molData.number), 32, 16);
            const material = new THREE.MeshPhongMaterial({
                color: new THREE.Color(...molData.color)
            });
            mesh = new THREE.Mesh(geometry, material);
            mesh.userData.moleculeData = molData;
            scene.add(mesh);
            moleculeMeshes.set(molData.id, mesh);
        } else {
            // Update existing mesh
            mesh.userData.moleculeData = molData;

            // Update material color if needed
            mesh.material.color.setRGB(...molData.color);

            // Update geometry if number has changed
            if (mesh.geometry.parameters.radius !== 0.2 + 0.1 * Math.log2(molData.number)) {
                mesh.geometry.dispose();
                mesh.geometry = new THREE.SphereGeometry(0.2 + 0.1 * Math.log2(molData.number), 32, 16);
            }
        }

        // Update position
        mesh.position.set(...molData.position);

        // Apply Rotation
        if (molData.angularVelocity) {
            const axis = new THREE.Vector3(...molData.angularVelocity).normalize();
            const angle = new THREE.Vector3(...molData.angularVelocity).length();
            const quaternion = new THREE.Quaternion().setFromAxisAngle(axis, angle);
            mesh.quaternion.multiplyQuaternions(quaternion, mesh.quaternion);
        }

        // Update or create velocity vector
        if (showVectors && molData.velocity) {
            // Remove old arrow if it exists
            if (mesh.userData.velocityArrow) {
                scene.remove(mesh.userData.velocityArrow);
            }

            const dir = new THREE.Vector3(...molData.velocity);
            const origin = new THREE.Vector3(...molData.position);
            const length = dir.length();
            const hex = 0xffff00; // Yellow

            const arrowHelper = new THREE.ArrowHelper(
                dir.normalize(),
                origin,
                length,
                hex
            );
            scene.add(arrowHelper);
            mesh.userData.velocityArrow = arrowHelper;
        } else if (!showVectors && mesh.userData.velocityArrow) {
            // Remove arrow if vectors are disabled
            scene.remove(mesh.userData.velocityArrow);
            mesh.userData.velocityArrow = null;
        }

        // Update or create halo
        updateHalo(mesh, molData);
    });
}

function updateHalo(mesh, molData) {
    // Remove old halo if it exists
    if (mesh.userData.halo) {
        scene.remove(mesh.userData.halo);
        mesh.userData.halo = null;
    }

    // Add Halo (if recently reacted)
    if (molData.lastReactionTime > 0) {
        const haloColor = new THREE.Color(...molData.color);
        const timeSinceReaction = performance.now() - molData.lastReactionTime;
        const alpha = Math.max(0, 1 - timeSinceReaction / 1000); // Fade out over 1 second

        if (alpha > 0) { // Only add if visible
            const radius = mesh.geometry.parameters.radius;
            const haloGeometry = new THREE.SphereGeometry(
                radius * (1.5 + 0.5 * alpha),
                32,
                16
            );
            const haloMaterial = new THREE.MeshBasicMaterial({
                color: haloColor,
                transparent: true,
                opacity: alpha,
            });
            const halo = new THREE.Mesh(haloGeometry, haloMaterial);
            halo.position.copy(mesh.position);
            scene.add(halo);
            mesh.userData.halo = halo;
        }
    }
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

    // Event Listeners
    window.addEventListener('resize', onWindowResize, false);
    document.addEventListener('keydown', onKeyDown, false);
    renderer.domElement.addEventListener('click', onClick, false);

    animate();
}

function animate() {
    requestAnimationFrame(animate);
    controls.update();

    // Request a simulation step from the worker (if not paused)
    if (autoStep && !paused) {
        simulationWorker.postMessage({
            type: 'step'
        });
    }

    // Update info display
    infoDiv.innerText = `Molecules: ${simulationData.molecules.length}\n` +
        `Temperature: ${simulationData.temperature.toFixed(2)}\n` +
        `${paused ? 'PAUSED' : 'RUNNING'}\n` +
        `${autoStep ? 'AUTO' : 'MANUAL'}`;
    render();
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

        // If we clicked on a halo, find its parent molecule
        if (!molIntersected.userData.moleculeData) {
            for (const mesh of meshArray) {
                if (mesh.userData.halo === molIntersected) {
                    molIntersected = mesh;
                    break;
                }
            }
        }

        // Make sure we have molecule data
        if (molIntersected.userData.moleculeData) {
            selectedMolecule = molIntersected.userData.moleculeData;

            // Display molecule info
            moleculeInfoDiv.style.display = 'block';
            moleculeInfoDiv.innerHTML = `
        Number: ${selectedMolecule.number}<br>
        Prime Factors: ${formatPrimeFactors(selectedMolecule.prime_factors)}<br>
        Mass: ${selectedMolecule.mass.toFixed(2)}<br>
        Charge: ${selectedMolecule.charge.toFixed(2)}<br>
        Position: (${selectedMolecule.position[0].toFixed(2)}, 
                  ${selectedMolecule.position[1].toFixed(2)}, 
                  ${selectedMolecule.position[2].toFixed(2)})<br>
        Velocity: ${vectorMagnitude(selectedMolecule.velocity).toFixed(2)}
      `;
        }
    } else {
        selectedMolecule = null;
        moleculeInfoDiv.style.display = 'none';
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

init();