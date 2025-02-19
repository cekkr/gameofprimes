import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

let scene, camera, renderer, controls;
let molecules = [];  // Array to hold Three.js objects representing molecules
let simulationWorker;
let simulationData = { molecules: [], temperature: 0 }; // Holds the latest data, initialize temperature
let paused = false;
let autoStep = true;
let showVectors = false;
const infoDiv = document.getElementById('info');
const moleculeInfoDiv = document.getElementById('moleculeInfo');
let selectedMolecule = null;


function init() {
    const simulationSize = 25

    // Scene setup
    scene = new THREE.Scene();
    camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
    camera.position.z = simulationSize * 2;

    renderer = new THREE.WebGLRenderer({ antialias: true });
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
    const boxSize = simulationSize; // Same as Python simulation size
    const boxGeometry = new THREE.BoxGeometry(boxSize, boxSize, boxSize);
    const wireframe = new THREE.WireframeGeometry(boxGeometry);
    const line = new THREE.LineSegments(wireframe);
    line.material.color.set(0x808080);
    scene.add(line);


    // Initialize Simulation Worker
    simulationWorker = new Worker('simulation.js', { type: "module" }); // Add type: "module"
    simulationWorker.onmessage = handleWorkerMessage;

    // Initial simulation setup (send parameters to worker)
     simulationWorker.postMessage({
        type: 'init',
        size: simulationSize,
        moleculeCount: 300,
        maxNumber: 100,
        timeScale: 0.1, // Send timeScale
    });


    // Event Listeners
    window.addEventListener('resize', onWindowResize, false);
    document.addEventListener('keydown', onKeyDown, false);
    renderer.domElement.addEventListener('click', onClick, false);

    animate();
}


function handleWorkerMessage(event) {
    if (event.data.type === 'update') {
        //console.log("Received update from worker:", event.data); // Debugging
        simulationData = event.data;
        updateMoleculeMeshes(); // Update Three.js objects based on worker data
    }
}


function updateMoleculeMeshes() {
    //console.log("Updating molecule meshes. simulationData:", simulationData); // Debugging

    // Remove old molecules
    molecules.forEach(molecule => {
        scene.remove(molecule);
    });
    molecules = []; // Clear the array

    // Create/update molecules based on simulationData
    simulationData.molecules.forEach(molData => {
       // console.log("Processing molecule data:", molData); // Debugging
        const geometry = new THREE.SphereGeometry(0.2 + 0.1 * Math.log2(molData.number), 32, 16);
        const material = new THREE.MeshPhongMaterial({ color: new THREE.Color(...molData.color) });
        const sphere = new THREE.Mesh(geometry, material);
        sphere.position.set(...molData.position);
        sphere.userData.moleculeData = molData; // Store simulation data with the mesh
        scene.add(sphere);
        molecules.push(sphere);

        if (showVectors && molData.velocity) {
            const dir = new THREE.Vector3(...molData.velocity);
            const origin = new THREE.Vector3(...molData.position);
            const length = dir.length();
            const hex = 0xffff00; // Yellow

            const arrowHelper = new THREE.ArrowHelper(dir.normalize(), origin, length, hex);
            scene.add(arrowHelper);
            molecules.push(arrowHelper);
        }
    });
}



function animate() {
    requestAnimationFrame(animate);
    controls.update();

    // Request a simulation step from the worker (if not paused)
     if (autoStep && !paused) {
        simulationWorker.postMessage({ type: 'step' });
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
            // Force an update to show/hide vectors immediately.
            updateMoleculeMeshes();
            break;
        case 'KeyR': // Reset Camera
            controls.reset();
            break;
        case 'ArrowRight':
            if (paused) {
                simulationWorker.postMessage({ type: 'step' });
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
    const intersects = raycaster.intersectObjects(molecules.filter(obj => obj instanceof THREE.Mesh));

    if (intersects.length> 0) {
        // Get the closest intersection (if multiple)
        const closestIntersection = intersects[0];
        selectedMolecule = closestIntersection.object.userData.moleculeData;

        // Display molecule info
        moleculeInfoDiv.style.display = 'block';
        moleculeInfoDiv.innerHTML = `
            Number: ${selectedMolecule.number}<br>
            Prime Factors: ${formatPrimeFactors(selectedMolecule.prime_factors)}<br>
            Mass: ${selectedMolecule.mass.toFixed(2)}<br>
            Charge: ${selectedMolecule.charge.toFixed(2)}<br>
            Position: (${selectedMolecule.position[0].toFixed(2)}, ${selectedMolecule.position[1].toFixed(2)}, ${selectedMolecule.position[2].toFixed(2)})<br>
            Velocity: ${vectorMagnitude(selectedMolecule.velocity).toFixed(2)}
        `;
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

init();