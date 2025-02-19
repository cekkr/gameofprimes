import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

// Core simulation state
let scene, camera, renderer, controls;
let simulationWorker;
let init_done = false;
let lastFrameTime = 0;
let currentHighlight = null;

// Simulation parameters
const simulationData = {
  molecules: [],
  temperature: 0,
  reactionCount: 0
};

// UI state
let paused = false;
let autoStep = true;
let showVectors = false;
let visualizationMode = 'simple';
let selectedMolecule = null;

// DOM elements
const moleculeMeshes = new Map();
let infoDiv, moleculeInfoDiv;

/**
 * Main initialization function
 */
function init() {
  const spaceDimension = 150;
  
  // Scene setup
  scene = new THREE.Scene();
  camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
  camera.position.z = spaceDimension;

  renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setSize(window.innerWidth, window.innerHeight);
  document.body.appendChild(renderer.domElement);

  // Setup UI elements
  setupUserInterface();

  // Setup controls
  controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.dampingFactor = 0.25;

  // Add lighting
  setupLighting();
  
  // Add coordinate axes and boundary box
  setupSceneHelpers(spaceDimension);

  // Initialize simulation
  initializeSimulation(spaceDimension, 2, 100, 1);

  // Make the simulation state globally accessible
  makeGloballyAccessible();
  
  // Bind events
  bindEvents();
  
  // Start animation loop
  requestAnimationFrame(animate);
}

/**
 * Setup scene lighting
 */
function setupLighting() {
  const ambientLight = new THREE.AmbientLight(0x404040);
  scene.add(ambientLight);
  
  const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
  directionalLight.position.set(1, 1, 1).normalize();
  scene.add(directionalLight);
}

/**
 * Setup scene helpers (axes, boundary box)
 */
function setupSceneHelpers(size) {
  // Axes helper
  const axesHelper = new THREE.AxesHelper(5);
  scene.add(axesHelper);

  // Boundary box
  const boxGeometry = new THREE.BoxGeometry(size, size, size);
  const wireframe = new THREE.WireframeGeometry(boxGeometry);
  const line = new THREE.LineSegments(wireframe);
  line.material.color.set(0x808080);
  scene.add(line);
}

/**
 * Setup UI elements
 */
function setupUserInterface() {
  // Create DOM containers
  createInfoContainers();
  
  // Add CSS styles
  addSimulationStyles();
  
  // Create control panel
  createControlPanel();
  
  // Expose control functions globally
  window.pauseSimulation = togglePause;
  window.toggleVectors = toggleVectorDisplay;
  window.updateControlsState = updateControlsState;
}

/**
 * Create info containers for UI
 */
function createInfoContainers() {
  // Info panel
  infoDiv = document.createElement('div');
  infoDiv.id = 'info';
  document.body.appendChild(infoDiv);
  
  // Molecule info panel
  moleculeInfoDiv = document.createElement('div');
  moleculeInfoDiv.id = 'moleculeInfo';
  document.body.appendChild(moleculeInfoDiv);
}

/**
 * Initialize the simulation
 */
function initializeSimulation(size, moleculeCount = 5, maxNumber = 200, timeScale = 1.0) {
  console.log("Initializing simulation...");
  
  // Terminate existing worker if present
  if (simulationWorker) {
    simulationWorker.terminate();
  }

  try {
    // Create worker
    simulationWorker = new Worker(new URL('./simulation-worker-gpu.js', import.meta.url), {
      type: 'module'
    });
    
    // Setup message handler
    simulationWorker.onmessage = handleWorkerMessage;
    
    // Setup error handler
    simulationWorker.onerror = function(error) {
      console.error("Worker error:", error);
      displayErrorMessage(`Simulation error: ${error.message}`);
    };
    
    // Calculate total molecules based on size
    const totalMolecules = Math.floor(moleculeCount * size);
    
    // Initialize UI elements
    updateTimeScaleUI(timeScale);
    
    // Reset pause state
    paused = false;
    updateControlsState();
    
    // Send initialization message to worker
    simulationWorker.postMessage({
      type: 'init',
      size: size,
      moleculeCount: totalMolecules,
      maxNumber: maxNumber,
      timeScale: timeScale
    });
    
    // Set timer to mark initialization as complete
    setTimeout(() => {
      init_done = true;
    }, 500);
    
    // Add worker status check
    setTimeout(checkWorkerStatus, 5000);
    
  } catch (error) {
    console.error("Critical error initializing worker:", error);
    displayErrorMessage(`Initialization error: ${error.message}`);
    createTestMolecules(); // Create test molecules for UI to function
  }
}

/**
 * Handle messages from simulation worker
 */
function handleWorkerMessage(event) {
  const data = event.data;
  
  switch (data.type) {
    case 'update':
        console.log("main.js mols update:", data);
        handleUpdateMessage(data);
        break;
      
    case 'cleanup_complete':
      console.log('Worker cleanup completed');
      break;
      
    case 'error':
      console.error('Worker error:', data.message);
      break;
      
    default:
      console.warn('Unknown message type from worker:', data.type);
  }
}

/**
 * Handle update messages from worker
 */
function handleUpdateMessage(data) {
  // Update temperature
  simulationData.temperature = data.temperature;
  
  // Update reaction count if provided
  if (data.reactionCount !== undefined) {
    simulationData.reactionCount = data.reactionCount;
  }
  
  // Handle molecule updates
  if (data.molecules) {
    const updatedMolecules = data.molecules.molecules || [];
    const removedIds = data.molecules.removedIds || [];
    
    // Process molecule sizes
    for (let mol of updatedMolecules) {
      mol.size = calculateMoleculeSize(mol);      
    }
    
    // Remove deleted molecules
    if (removedIds.length > 0) {
      removedIds.forEach(id => {
        const mesh = moleculeMeshes.get(id);
        if (mesh) {
          removeMoleculeAndEffects(mesh);
          moleculeMeshes.delete(id);
        }
      });
    }
    
    // Update simulation data
    simulationData.molecules = updatedMolecules;
    
    // Update molecule meshes
    updateMoleculeMeshes();
  }
}

/**
 * Calculate molecule size based on its properties
 */
function calculateMoleculeSize(molData) {
  // Base size
  const baseSize = 1.0;
  
  // Size factors
  const numberFactor = Math.log2(molData.number) * 0.05;
  const factorCount = Object.keys(molData.prime_factors).length;
  const complexityFactor = factorCount * 0.06;
  
  // Total exponents factor
  const totalExponents = Object.values(molData.prime_factors).reduce((sum, exp) => sum + exp, 0);
  const exponentFactor = totalExponents * 0.04;
  
  // Final size
  return baseSize + numberFactor + complexityFactor + exponentFactor;
}

/**
 * Update or create molecule meshes
 */
function updateMoleculeMeshes() {
  simulationData.molecules.forEach(molData => {
    let mesh = moleculeMeshes.get(molData.id);
    
    // Create new mesh if needed
    if (!mesh) {
      mesh = createMoleculeMesh(molData);
      scene.add(mesh);
      moleculeMeshes.set(molData.id, mesh);
    } else {
      // Update existing mesh
      mesh.userData.moleculeData = molData;
      mesh.position.set(...molData.position);
      
      // Update velocity vectors if enabled
      updateVelocityVector(mesh, molData, showVectors);
    }
    
    // Update visual effects
    updateHalo(mesh, molData);
  });
}

/**
 * Create a molecule mesh
 */
function createMoleculeMesh(molData) {
  // Analyze properties
  const properties = analyzeEnhancedMoleculeProperties(molData);
  
  // Create geometry
  const geometry = new THREE.SphereGeometry(molData.size, 16, 16);
  
  // Create material
  const material = new THREE.MeshPhongMaterial({
    color: new THREE.Color(...molData.color),
    shininess: 30
  });
  
  // Create mesh
  const mesh = new THREE.Mesh(geometry, material);
  mesh.position.set(...molData.position);
  
  // Store molecule data
  mesh.userData.moleculeData = molData;
  mesh.userData.visualProperties = properties;
  
  return mesh;
}

/**
 * Analyze molecule properties
 */
function analyzeEnhancedMoleculeProperties(molData) {
  // Extract factor information
  const factorCount = Object.keys(molData.prime_factors).length;
  const factors = Object.entries(molData.prime_factors);
  const largestPrime = Math.max(...Object.keys(molData.prime_factors).map(Number));
  const smallestPrime = Math.min(...Object.keys(molData.prime_factors).map(Number));
  
  // Calculate total exponent
  const totalExponent = factors.reduce((sum, [_prime, exp]) => sum + exp, 0);
  
  // Property analysis
  const isPrimePure = factorCount === 1 && factors[0][1] === 1;
  const isPrimePower = factorCount === 1 && factors[0][1] > 1;
  const hasOnlySmallPrimes = Object.keys(molData.prime_factors).every(prime => parseInt(prime) < 10);
  const hasLargePrimes = Object.keys(molData.prime_factors).some(prime => parseInt(prime) > 50);
  
  // Categorization
  let category, stability;
  
  // Determine category
  if (isPrimePure) {
    if (largestPrime < 10) category = 'elementary';
    else if (largestPrime > 50) category = 'exotic';
    else category = 'prime';
  } else if (isPrimePower) {
    category = 'power';
  } else if (hasOnlySmallPrimes && factorCount <= 2) {
    category = 'stable';
  } else if (hasLargePrimes || factorCount > 3) {
    category = 'complex';
  } else {
    category = 'standard';
  }
  
  // Determine stability
  if (isPrimePure || isPrimePower) {
    stability = 'high';
  } else if (hasOnlySmallPrimes && totalExponent < 4) {
    stability = 'medium';
  } else if (hasLargePrimes || totalExponent > 6) {
    stability = 'low';
  } else {
    stability = 'normal';
  }
  
  // Use provided color or generate one
  let color = molData.color || generateColorFromFactors(largestPrime, factorCount);
  
  return {
    factorCount,
    largestPrime,
    smallestPrime,
    totalExponent,
    isPrimePure,
    isPrimePower,
    hasLargePrimes,
    hasOnlySmallPrimes,
    category,
    stability,
    energy: molData.number / (1 + Math.log(molData.number)),
    complexity: factorCount * totalExponent,
    primeFactors: factors,
    color
  };
}

/**
 * Generate color based on prime factors
 */
function generateColorFromFactors(largestPrime, factorCount) {
  const hue = (largestPrime % 360) / 360;
  const saturation = 0.5 + (factorCount * 0.1);
  const lightness = 0.5;
  return hslToRgb(hue, saturation, lightness);
}

/**
 * Convert HSL to RGB
 */
function hslToRgb(h, s, l) {
  let r, g, b;
  
  if (s === 0) {
    r = g = b = l;
  } else {
    const hue2rgb = (p, q, t) => {
      if (t < 0) t += 1;
      if (t > 1) t -= 1;
      if (t < 1/6) return p + (q - p) * 6 * t;
      if (t < 1/2) return q;
      if (t < 2/3) return p + (q - p) * (2/3 - t) * 6;
      return p;
    };
    
    const q = l < 0.5 ? l * (1 + s) : l + s - l * s;
    const p = 2 * l - q;
    
    r = hue2rgb(p, q, h + 1/3);
    g = hue2rgb(p, q, h);
    b = hue2rgb(p, q, h - 1/3);
  }
  
  return [r, g, b];
}

/**
 * Display error message
 */
function displayErrorMessage(message) {
  const errorDiv = document.getElementById('errorMessage');
  if (errorDiv) {
    const errorTextEl = document.getElementById('errorText');
    if (errorTextEl) {
      errorTextEl.textContent = message;
    } else {
      errorDiv.innerHTML = `<p>${message}</p>`;
    }
    errorDiv.style.display = 'block';
    
    // Auto-hide after 10 seconds
    setTimeout(() => {
      errorDiv.style.display = 'none';
    }, 10000);
  } else {
    // Fallback
    alert(`Simulation error: ${message}`);
  }
}

/**
 * Remove molecule mesh and all associated effects
 */
function removeMoleculeAndEffects(mesh) {
  // Remove main mesh
  scene.remove(mesh);
  
  // Remove all associated effects
  const effectTypes = [
    'halo', 'glowEffect', 'particleSystem',
    'velocityArrow', 'elementalGlow', 'orbitalRing',
    'atomicStructure'
  ];
  
  effectTypes.forEach(effectType => {
    if (mesh.userData[effectType]) {
      scene.remove(mesh.userData[effectType]);
      mesh.userData[effectType] = null;
    }
  });
  
  // Clean up geometry and materials
  if (mesh.material) {
    if (Array.isArray(mesh.material)) {
      mesh.material.forEach(mat => mat.dispose());
    } else {
      mesh.material.dispose();
    }
  }
  
  if (mesh.geometry) {
    mesh.geometry.dispose();
  }
}

/**
 * Update info display
 */
function updateInfoDisplay() {
  if (!infoDiv) return;
  
  infoDiv.innerHTML = `
    <div class="info-section">
      <h3>Simulation</h3>
      <p>Molecules: <strong>${simulationData.molecules.length}</strong></p>
      <p>Temperature: <strong>${simulationData.temperature.toFixed(2)}</strong></p>
      <p>Reactions: <strong>${simulationData.reactionCount || 0}</strong></p>
      <p>Status: <strong>${paused ? 'PAUSED' : 'RUNNING'}</strong></p>
    </div>

    <div class="controls-info">
      <p>Space: Pause | A: Auto | V: Vectors</p>
      <p>Click: Select molecule</p>
    </div>
  `;
}

/**
 * Add CSS styles
 */
function addSimulationStyles() {
  const styleElement = document.createElement('style');
  styleElement.textContent = `
    /* General styles */
    body {
      margin: 0;
      overflow: hidden;
      font-family: 'Arial', sans-serif;
    }

    /* Info panel */
    #info {
      position: absolute;
      top: 10px;
      left: 10px;
      background-color: rgba(0, 0, 0, 0.6);
      color: #fff;
      padding: 15px;
      border-radius: 8px;
      font-size: 14px;
      max-width: 250px;
      z-index: 100;
      pointer-events: none;
    }

    .info-section {
      margin-bottom: 10px;
    }

    .info-section h3 {
      margin: 0 0 8px 0;
      font-size: 16px;
      color: #88ccff;
      border-bottom: 1px solid rgba(255, 255, 255, 0.3);
      padding-bottom: 4px;
    }

    .info-section p {
      margin: 4px 0;
      font-size: 13px;
    }

    .info-section strong {
      color: #aaffaa;
    }

    .controls-info {
      margin-top: 15px;
      font-size: 12px;
      color: #aaaaaa;
    }

    /* Molecule info panel */
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
      z-index: 100;
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
      flex-wrap: wrap;
    }

    .basic-properties, .advanced-properties {
      flex: 1;
      min-width: 160px;
      padding-right: 10px;
    }

    .position-data {
      margin-top: 10px;
      font-size: 0.9em;
      color: #aaa;
    }

    #moleculeInfo p {
      margin: 5px 0;
      font-size: 13px;
    }

    #moleculeInfo strong {
      color: #aaffaa;
    }

    /* Control panel */
    #controlPanel {
      position: absolute;
      top: 10px;
      right: 10px;
      background-color: rgba(0, 0, 0, 0.7);
      padding: 15px;
      border-radius: 8px;
      z-index: 100;
      width: 180px;
    }

    #controlPanel h3 {
      margin-top: 0;
      color: #88ccff;
      border-bottom: 1px solid rgba(255, 255, 255, 0.3);
      padding-bottom: 4px;
      font-size: 16px;
      margin-bottom: 12px;
    }

    .control-group {
      margin-bottom: 12px;
    }

    .control-label {
      display: block;
      color: #fff;
      margin-bottom: 5px;
      font-size: 12px;
    }

    .slider-container {
      display: flex;
      align-items: center;
    }

    .control-slider {
      flex-grow: 1;
      margin-right: 8px;
      height: 5px;
      -webkit-appearance: none;
      appearance: none;
      background: #555;
      outline: none;
      border-radius: 3px;
    }

    .control-slider::-webkit-slider-thumb {
      -webkit-appearance: none;
      appearance: none;
      width: 15px;
      height: 15px;
      border-radius: 50%;
      background: #88ccff;
      cursor: pointer;
    }

    .control-slider::-moz-range-thumb {
      width: 15px;
      height: 15px;
      border-radius: 50%;
      background: #88ccff;
      cursor: pointer;
    }

    .slider-value {
      color: #fff;
      min-width: 30px;
      font-size: 12px;
      text-align: right;
    }

    .control-button {
      background-color: #2a2a2a;
      color: white;
      border: none;
      border-radius: 4px;
      padding: 6px 12px;
      margin-right: 8px;
      margin-bottom: 8px;
      cursor: pointer;
      font-size: 12px;
      transition: background-color 0.2s;
    }

    .control-button:hover {
      background-color: #3a3a3a;
    }

    .control-button.active {
      background-color: #88ccff;
      color: #000;
    }

    /* Tooltip */
    .tooltip {
      position: relative;
      display: inline-block;
    }

    .tooltip .tooltiptext {
      visibility: hidden;
      background-color: rgba(0, 0, 0, 0.8);
      color: #fff;
      text-align: center;
      border-radius: 4px;
      padding: 5px 10px;
      position: absolute;
      z-index: 101;
      bottom: 125%;
      left: 50%;
      transform: translateX(-50%);
      opacity: 0;
      transition: opacity 0.3s;
      font-size: 11px;
      white-space: nowrap;
    }

    .tooltip:hover .tooltiptext {
      visibility: visible;
      opacity: 1;
    }
    
    /* Error message */
    #errorMessage {
      position: absolute;
      top: 50%;
      left: 50%;
      transform: translate(-50%, -50%);
      background-color: rgba(220, 0, 0, 0.8);
      color: white;
      padding: 20px;
      border-radius: 8px;
      z-index: 200;
      display: none;
      text-align: center;
      max-width: 80%;
    }
    
    #errorMessage button {
      margin-top: 15px;
      padding: 8px 16px;
      background: #fff;
      color: #a00;
      border: none;
      border-radius: 4px;
      cursor: pointer;
    }
  `;

  document.head.appendChild(styleElement);
}

/**
 * Create control panel
 */
function createControlPanel() {
  const panel = document.createElement('div');
  panel.id = 'controlPanel';
  panel.innerHTML = `
    <h3>Controls</h3>
    
    <div class="control-group">
      <label class="control-label" for="temperatureSlider">Temperature</label>
      <div class="slider-container">
        <input type="range" min="0.5" max="2.0" step="0.1" value="1.0" class="control-slider" id="temperatureSlider">
        <span class="slider-value" id="temperatureValue">1.0</span>
      </div>
    </div>
    
    <div class="control-group">
      <label class="control-label" for="timeScaleSlider">Simulation Speed</label>
      <div class="slider-container">
        <input type="range" min="0.2" max="10.0" step="0.1" value="0.1" class="control-slider" id="timeScaleSlider">
        <span class="slider-value" id="timeScaleValue">0.1</span>
      </div>
    </div>
    
    <div class="control-group">
      <button id="pauseButton" class="control-button">Pause</button>
      <button id="resetButton" class="control-button">Reset</button>
      <button id="addMoleculesButton" class="control-button tooltip">
        Add Molecules
        <span class="tooltiptext">Adds 20 new molecules</span>
      </button>
      <button id="toggleVectorsButton" class="control-button">Vectors</button>
    </div>
    
    <div class="control-group">
      <label class="control-label">Visualization</label>
      <button id="simpleModeButton" class="control-button active">Basic</button>
      <button id="detailedModeButton" class="control-button">Detailed</button>
    </div>
  `;

  // Add to DOM
  document.body.appendChild(panel);

  // Add error message container
  const errorDiv = document.createElement('div');
  errorDiv.id = 'errorMessage';
  errorDiv.innerHTML = `
    <p id="errorText">Simulation error</p>
    <button id="dismissError">Close</button>
  `;
  document.body.appendChild(errorDiv);
}

/**
 * Update UI controls state
 */
function updateControlsState() {
  // Update pause button
  const pauseButton = document.getElementById('pauseButton');
  if (pauseButton) {
    pauseButton.classList.toggle('active', paused);
    pauseButton.textContent = paused ? 'Resume' : 'Pause';
  }
  
  // Update vectors button
  const vectorsButton = document.getElementById('toggleVectorsButton');
  if (vectorsButton) {
    vectorsButton.classList.toggle('active', showVectors);
  }
  
  // Update visualization mode buttons
  const simpleButton = document.getElementById('simpleModeButton');
  const detailedButton = document.getElementById('detailedModeButton');
  if (simpleButton && detailedButton) {
    if (visualizationMode === 'simple') {
      simpleButton.classList.add('active');
      detailedButton.classList.remove('active');
    } else {
      simpleButton.classList.remove('active');
      detailedButton.classList.add('active');
    }
  }
}

/**
 * Bind event handlers
 */
function bindEvents() {
  // Temperature slider
  const tempSlider = document.getElementById('temperatureSlider');
  const tempValue = document.getElementById('temperatureValue');
  if (tempSlider && tempValue) {
    tempSlider.addEventListener('input', function() {
      const value = parseFloat(this.value).toFixed(1);
      tempValue.textContent = value;
      if (simulationWorker) {
        simulationWorker.postMessage({
          type: 'set_temperature',
          value: parseFloat(value)
        });
      }
    });
  }

  // Time scale slider
  const timeSlider = document.getElementById('timeScaleSlider');
  const timeValue = document.getElementById('timeScaleValue');
  if (timeSlider && timeValue) {
    timeSlider.addEventListener('input', function() {
      const value = parseFloat(this.value).toFixed(1);
      timeValue.textContent = value;
      if (simulationWorker) {
        simulationWorker.postMessage({
          type: 'set_timescale',
          value: parseFloat(value)
        });
        window.currentTimeScale = parseFloat(value);
      }
    });
  }

  // Pause button
  const pauseButton = document.getElementById('pauseButton');
  if (pauseButton) {
    pauseButton.addEventListener('click', togglePause);
  }

  // Reset button
  const resetButton = document.getElementById('resetButton');
  if (resetButton) {
    resetButton.addEventListener('click', resetSimulation);
  }

  // Add molecules button
  const addButton = document.getElementById('addMoleculesButton');
  if (addButton) {
    addButton.addEventListener('click', addMolecules);
  }

  // Vectors button
  const vectorsButton = document.getElementById('toggleVectorsButton');
  if (vectorsButton) {
    vectorsButton.addEventListener('click', toggleVectorDisplay);
  }

  // Visualization mode buttons
  const simpleButton = document.getElementById('simpleModeButton');
  const detailedButton = document.getElementById('detailedModeButton');
  if (simpleButton && detailedButton) {
    simpleButton.addEventListener('click', () => setVisualizationMode('simple'));
    detailedButton.addEventListener('click', () => setVisualizationMode('detailed'));
  }

  // Dismiss error button
  const dismissButton = document.getElementById('dismissError');
  if (dismissButton) {
    dismissButton.addEventListener('click', () => {
      document.getElementById('errorMessage').style.display = 'none';
    });
  }

  // Keyboard controls
  document.addEventListener('keydown', onKeyDown);
  
  // Mouse controls
  renderer.domElement.addEventListener('click', onClick);
  
  // Window resize
  window.addEventListener('resize', onWindowResize);
  
  // Page visibility
  window.addEventListener('visibilitychange', onVisibilityChange);
  
  // Page unload
  window.addEventListener('beforeunload', cleanupSimulation);
}

/**
 * Animation loop
 */
function animate(timestamp) {
  // Calculate delta time
  const deltaTime = lastFrameTime ? timestamp - lastFrameTime : 16.6;
  lastFrameTime = timestamp;

  requestAnimationFrame(animate);
  
  // Update controls
  if (controls) controls.update();

  // Step simulation if not paused
  if (autoStep && !paused && init_done) {
    if (simulationWorker) {
      simulationWorker.postMessage({ type: 'step' });
    }
  }

  // Update visual effects
  updateEnhancedVisuals(deltaTime, timestamp);
  
  // Update info display
  updateInfoDisplay();

  // Render scene
  renderer.render(scene, camera);
}

/**
 * Update special visual effects
 */
function updateEnhancedVisuals(deltaTime, time) {
  moleculeMeshes.forEach(mesh => {
    // Update elemental glow effects
    if (mesh.userData.elementalGlow) {
      const halo = mesh.userData.elementalGlow;
      halo.position.copy(mesh.position);
      if (halo.userData.updateFunction) {
        halo.userData.updateFunction.call(halo, time);
      }
    }

    // Update orbital rings
    if (mesh.userData.orbitalRing) {
      const ring = mesh.userData.orbitalRing;
      if (ring.userData.updateFunction) {
        ring.userData.updateFunction.call(ring, deltaTime);
      }
    }

    // Update atomic structures
    if (mesh.userData.atomicStructure) {
      const structure = mesh.userData.atomicStructure;
      if (structure.userData.updateFunction) {
        structure.userData.updateFunction.call(structure, time);
      }
    }
    
    // Update particle systems
    if (mesh.userData.particleSystem) {
      const particles = mesh.userData.particleSystem;
      if (particles.userData.updateFunction) {
        particles.userData.updateFunction.call(particles, deltaTime);
      }
    }
  });
}

/**
 * Update velocity vector visualization
 */
function updateVelocityVector(mesh, molData, show) {
  // Remove existing arrow if present
  if (mesh.userData.velocityArrow) {
    scene.remove(mesh.userData.velocityArrow);
    mesh.userData.velocityArrow = null;
  }

  // Add new vector if enabled
  if (show && molData.velocity) {
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
  }
}

/**
 * Update molecule halo/reaction effects
 */
function updateHalo(mesh, molData) {
  // Remove old effects
  removeHaloEffects(mesh);

  // Add reaction effects if molecule recently reacted
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

/**
 * Remove halo effects
 */
function removeHaloEffects(mesh) {
  const effectTypes = ['halo', 'glowEffect', 'particleSystem'];
  
  effectTypes.forEach(effectType => {
    if (mesh.userData[effectType]) {
      scene.remove(mesh.userData[effectType]);
      mesh.userData[effectType] = null;
    }
  });
}

/**
 * Add reaction visual effects
 */
function addReactionEffects(mesh, molData, intensity) {
  // Extract properties
  const baseColor = new THREE.Color(...molData.color);
  const position = new THREE.Vector3(...molData.position);
  const size = mesh.geometry.parameters.radius;

  // Outer halo (pulsating sphere)
  addOuterHalo(mesh, baseColor, position, size, intensity);
  
  // Particle emission for complex reactions
  if (Object.keys(molData.prime_factors).length > 2) {
    createParticleSystem(mesh, molData, intensity);
  }
  
  // Inner glow
  addGlowEffect(mesh, molData, intensity);
}

/**
 * Add outer halo effect
 */
function addOuterHalo(mesh, baseColor, position, size, intensity) {
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
}

/**
 * Create particle system for reactions
 */
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

/**
 * Add inner glow effect
 */
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

/**
 * Add elemental glow effect
 */
function addElementalGlowEffect(mesh, properties) {
  // Create a slightly larger, emissive sphere to create glow effect
  const glowSize = mesh.geometry.parameters.radius * 1.2;
  const glowGeometry = new THREE.SphereGeometry(glowSize, 32, 16);
  
  // Get base color
  const baseColor = properties.color ? 
    new THREE.Color(...properties.color) : 
    new THREE.Color(0.5, 0.5, 1.0);
  
  // Choose glow color based on properties
  let glowColor;
  if (properties.category === 'elemental') {
    // Pure elements glow with bright blue/cyan
    glowColor = new THREE.Color(0.3, 0.7, 1.0);
  } else if (properties.category === 'prime') {
    // Prime molecules glow with purple hue
    glowColor = new THREE.Color(0.7, 0.3, 1.0);
  } else {
    // Default glow color - shift from base color
    glowColor = new THREE.Color().copy(baseColor);
    let glowHsl = {};
    glowColor.getHSL(glowHsl);
    glowColor.setHSL(
      glowHsl.h,
      Math.max(0, glowHsl.s - 0.5),
      Math.min(1, glowHsl.l + 0.3)
    );
  }
  
  const glowMaterial = new THREE.MeshBasicMaterial({
    color: glowColor,
    transparent: true,
    opacity: 0.15,
    side: THREE.BackSide,
    blending: THREE.AdditiveBlending
  });
  
  const glowMesh = new THREE.Mesh(glowGeometry, glowMaterial);
  
  // Setup pulsation parameters
  glowMesh.userData.pulseParams = {
    baseSize: glowSize,
    speed: 0.002 + Math.random() * 0.001,
    amplitude: 0.1 + (properties.energy / 1000),
    phase: Math.random() * Math.PI * 2
  };
  
  // Update function for the glow effect
  glowMesh.userData.updateFunction = function(time) {
    const p = this.userData.pulseParams;
    const scale = 1 + Math.sin(time * p.speed + p.phase) * p.amplitude;
    this.scale.set(scale, scale, scale);
  };
  
  mesh.add(glowMesh);
  mesh.userData.elementalGlow = glowMesh;
}

/**
 * Toggle pause state
 */
function togglePause() {
  paused = !paused;
  updateControlsState();
}

/**
 * Toggle vector display
 */
function toggleVectorDisplay() {
  showVectors = !showVectors;
  updateControlsState();
  updateMoleculeMeshes();
}

/**
 * Reset simulation
 */
function resetSimulation() {
  const spaceDimension = 10;
  const timeSlider = document.getElementById('timeScaleSlider');
  const moleculeCount = timeSlider ? parseInt(timeSlider.value) * 5 : 5;
  const currentTimeScale = window.currentTimeScale || 1.0;
  
  initializeSimulation(spaceDimension, moleculeCount, 150, currentTimeScale);
}

/**
 * Add new molecules
 */
function addMolecules() {
  if (simulationWorker) {
    simulationWorker.postMessage({
      type: 'add_molecules',
      count: 20
    });
  }
}

/**
 * Set visualization mode
 */
function setVisualizationMode(mode) {
  visualizationMode = mode;
  
  // Update UI
  const simpleButton = document.getElementById('simpleModeButton');
  const detailedButton = document.getElementById('detailedModeButton');
  
  if (simpleButton && detailedButton) {
    if (mode === 'simple') {
      simpleButton.classList.add('active');
      detailedButton.classList.remove('active');
    } else {
      simpleButton.classList.remove('active');
      detailedButton.classList.add('active');
    }
  }
  
  // Notify worker
  if (simulationWorker) {
    simulationWorker.postMessage({
      type: 'set_visualization',
      mode: mode
    });
  }
  
  // Update visuals
  updateMoleculeMeshes();
}

/**
 * Update timeScale UI
 */
function updateTimeScaleUI(value) {
  const timeSlider = document.getElementById('timeScaleSlider');
  const timeValue = document.getElementById('timeScaleValue');
  
  if (timeSlider && timeValue) {
    timeSlider.value = value;
    timeValue.textContent = value.toFixed(1);
  }
}

/**
 * Create test molecules when worker fails
 */
function createTestMolecules() {
  console.log("Creating test molecules");
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
      velocity: [
        (Math.random() - 0.5) * 0.05,
        (Math.random() - 0.5) * 0.05,
        (Math.random() - 0.5) * 0.05
      ],
      prime_factors: {
        2: i % 2 + 1,
        3: Math.floor(i / 3) + 1,
        5: i % 3
      },
      mass: 4 + i,
      charge: (Math.random() - 0.5) * 2,
      color: [
        0.2 + Math.random() * 0.8,
        0.2 + Math.random() * 0.8,
        0.2 + Math.random() * 0.8
      ],
      angularVelocity: [
        (Math.random() - 0.5) * 0.01,
        (Math.random() - 0.5) * 0.01,
        (Math.random() - 0.5) * 0.01
      ],
      lastReactionTime: -1
    });
  }

  simulationData.molecules = testMolecules;
  updateMoleculeMeshes();
}

/**
 * Handle window resize
 */
function onWindowResize() {
  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(window.innerWidth, window.innerHeight);
}

/**
 * Handle key down events
 */
function onKeyDown(event) {
  switch (event.code) {
    case 'Space':
      togglePause();
      break;
      
    case 'KeyA':
      autoStep = !autoStep;
      break;
      
    case 'KeyV':
      toggleVectorDisplay();
      break;
      
    case 'KeyR': // Reset Camera
      controls.reset();
      break;

    case 'KeyS':
        simulationWorker.postMessage({ type: 'step' });
        console.log("Requested forced step");
        break;
      
    case 'ArrowRight':
      if (paused && simulationWorker) {
        simulationWorker.postMessage({
          type: 'step'
        });
      }
      break;
  }
}

/**
 * Handle mouse click events
 */
function onClick(event) {
  event.preventDefault();

  const mouse = new THREE.Vector2();
  mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
  mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;

  selectMoleculeAtPosition(mouse.x, mouse.y);
}

/**
 * Handle visibility change (tab switching)
 */
function onVisibilityChange() {
  if (document.visibilityState === 'hidden') {
    // Pause simulation when tab is not visible
    paused = true;
  } else {
    // Keep it paused if it was manually paused
    if (autoStep) {
      paused = false;
    }
  }
}

/**
 * Select molecule at specified position
 */
function selectMoleculeAtPosition(x, y) {
  const raycaster = new THREE.Raycaster();
  const mouse = new THREE.Vector2(x, y);

  raycaster.setFromCamera(mouse, camera);
  const meshArray = Array.from(moleculeMeshes.values());
  const intersects = raycaster.intersectObjects(meshArray);

  if (intersects.length > 0) {
    // Get the closest intersection
    const mol = intersects[0].object;
    
    if (mol.userData.moleculeData) {
      selectedMolecule = mol.userData.moleculeData;
      updateMoleculeInfoPanel(selectedMolecule, mol.userData.visualProperties);
      moleculeInfoDiv.style.display = 'block';
      highlightSelectedMolecule(mol);
    }
  } else {
    selectedMolecule = null;
    moleculeInfoDiv.style.display = 'none';
    clearMoleculeHighlight();
  }
}

/**
 * Update molecule info panel with selected molecule data
 */
function updateMoleculeInfoPanel(molecule, properties) {
  if (!molecule || !moleculeInfoDiv) return;

  // Format prime factors for display
  const formattedFactors = formatPrimeFactors(molecule.prime_factors);

  // Descriptive text for categories and stability
  const categoryDescriptions = {
    'elementary': 'Elementary (small pure prime)',
    'exotic': 'Exotic (large pure prime)',
    'prime': 'Prime number',
    'power': 'Prime power',
    'stable': 'Stable (simple compound)',
    'complex': 'Complex structure',
    'standard': 'Standard'
  };

  const stabilityDescriptions = {
    'high': 'High',
    'medium': 'Medium',
    'normal': 'Normal',
    'low': 'Low'
  };

  const categoryText = categoryDescriptions[properties.category] || properties.category;
  const stabilityText = stabilityDescriptions[properties.stability] || properties.stability;

  // Update UI
  moleculeInfoDiv.innerHTML = `
    <h3>Molecule ${molecule.number}</h3>
    <div class="molecule-info-container">
      <div class="basic-properties">
        <p><strong>Number:</strong> ${molecule.number}</p>
        <p><strong>Prime Factors:</strong> ${formattedFactors}</p>
        <p><strong>Mass:</strong> ${molecule.mass.toFixed(2)}</p>
        <p><strong>Charge:</strong> ${molecule.charge.toFixed(2)}</p>
        <p><strong>Velocity:</strong> ${vectorMagnitude(molecule.velocity).toFixed(2)}</p>
      </div>
      
      <div class="advanced-properties">
        <p><strong>Category:</strong> ${categoryText}</p>
        <p><strong>Stability:</strong> ${stabilityText}</p>
        <p><strong>Energy:</strong> ${properties.energy.toFixed(2)}</p>
        <p><strong>Complexity:</strong> ${properties.complexity.toFixed(1)}</p>
        <p><strong>Factor count:</strong> ${properties.factorCount}</p>
      </div>
    </div>
    
    <div class="position-data">
      <p><strong>Position:</strong> 
        (${molecule.position[0].toFixed(2)}, 
         ${molecule.position[1].toFixed(2)}, 
         ${molecule.position[2].toFixed(2)})</p>
    </div>
  `;
}

/**
 * Highlight selected molecule
 */
function highlightSelectedMolecule(mesh) {
  // Remove previous highlight
  clearMoleculeHighlight();

  // Create highlight effect
  const highlightGeometry = new THREE.SphereGeometry(
    mesh.geometry.parameters.radius * 1.1,
    32, 32
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

  // Store reference to current highlight
  currentHighlight = highlightMesh;
}

/**
 * Clear highlight from previously selected molecule
 */
function clearMoleculeHighlight() {
  if (currentHighlight) {
    scene.remove(currentHighlight);
    currentHighlight = null;
  }
}

/**
 * Format prime factors for display
 */
function formatPrimeFactors(factors) {
  let str = "";
  for (const prime in factors) {
    str += `${prime}^${factors[prime]} `;
  }
  return str.trim();
}

/**
 * Calculate vector magnitude
 */
function vectorMagnitude(vec) {
  return Math.sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

/**
 * Check if worker is functioning
 */
function checkWorkerStatus() {
  if (simulationData.molecules.length === 0) {
    console.warn("No molecules received after 5 seconds");
    createTestMolecules();
  }
}

/**
 * Cleanup simulation
 */
function cleanupSimulation() {
  if (simulationWorker) {
    // Request cleanup from worker
    simulationWorker.postMessage({
      type: 'cleanup'
    });

    // Force termination after timeout
    setTimeout(() => {
      simulationWorker.terminate();
      simulationWorker = null;
    }, 500);
  }

  // Clean up meshes and effects
  moleculeMeshes.forEach(mesh => {
    scene.remove(mesh);
    removeAllEffects(mesh);
  });

  moleculeMeshes.clear();
}

/**
 * Remove all effects from a mesh
 */
function removeAllEffects(mesh) {
  const effectTypes = [
    'halo', 'glowEffect', 'particleSystem',
    'velocityArrow', 'elementalGlow', 'orbitalRing',
    'atomicStructure'
  ];

  effectTypes.forEach(effectType => {
    if (mesh.userData[effectType]) {
      scene.remove(mesh.userData[effectType]);
      mesh.userData[effectType] = null;
    }
  });
}

/**
 * Make important variables globally accessible
 */
function makeGloballyAccessible() {
  window.scene = scene;
  window.camera = camera;
  window.renderer = renderer;
  window.controls = controls;
  window.paused = paused;
  window.showVectors = showVectors;
  window.autoStep = autoStep;
  window.visualizationMode = visualizationMode;
  window.simulationData = simulationData;
  window.moleculeMeshes = moleculeMeshes;
  window.simulationWorker = simulationWorker;
  window.currentTimeScale = 0.1;
}

// Initialize the simulation when script loads
init();

//import './advanced-simulation-worker.js';