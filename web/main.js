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

// Variabile globale per il timeScale
window.currentTimeScale = 0.1;

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
        
        // Nascondi automaticamente dopo 10 secondi
        setTimeout(() => {
            errorDiv.style.display = 'none';
        }, 10000);
    } else {
        // Fallback se l'elemento non esiste
        alert(`Errore nella simulazione: ${message}`);
    }
}

// 2. Gestione messaggi migliorata
function handleWorkerMessage(event) {
    // Debug
    console.log(`Ricevuto messaggio dal worker di tipo: ${event.data.type}`);

    if (event.data.type === 'update') {
        // Aggiorna temperatura
        simulationData.temperature = event.data.temperature;

        // Gestione molecole
        if (event.data.molecules) {
            const updatedMolecules = event.data.molecules.molecules || [];
            const removedIds = event.data.molecules.removedIds || [];

            // Debug
            console.log(`Main: ricevute ${updatedMolecules.length} molecole, ${removedIds.length} rimosse`);

            // Aggiorna dati simulazione
            simulationData.molecules = updatedMolecules;
            simulationData.reactionCount = event.data.reactionCount || 0;

            // Rimuovi mesh di molecole non più presenti
            if (removedIds.length > 0) {
                removedIds.forEach(id => {
                    const mesh = moleculeMeshes.get(id);
                    if (mesh) {
                        // Rimuovi mesh e tutti gli effetti associati
                        removeMoleculeAndEffects(mesh);
                        moleculeMeshes.delete(id);
                    }
                });
            }

            // Aggiorna o crea mesh per molecole
            updateMoleculeMeshes();
        } else {
            console.warn("Messaggio update ricevuto senza dati molecole");
        }
    } else if (event.data.type === 'cleanup_complete') {
        console.log('Pulizia worker completata');
    } else if (event.data.type === 'error') {
        console.error('Errore nel worker:', event.data.message);
        displayErrorMessage(event.data.message);
    }
}

// 3. Funzione migliorata per rimuovere mesh e effetti
function removeMoleculeAndEffects(mesh) {
    // Rimuovi mesh principale
    scene.remove(mesh);

    // Rimuovi tutti gli effetti associati
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

    // Assicurati che il materiale venga correttamente deallocato
    if (mesh.material) {
        if (Array.isArray(mesh.material)) {
            mesh.material.forEach(mat => mat.dispose());
        } else {
            mesh.material.dispose();
        }
    }

    // Libera la geometria
    if (mesh.geometry) {
        mesh.geometry.dispose();
    }
}

// 5. Calcolo dimensione proporzionale ai fattori primi
function calculateMoleculeSize(molData) {
    // Base size
    const baseSize = 0.2;

    // Fattore di scala basato sul numero
    const numberFactor = Math.log2(molData.number) * 0.05;

    // Fattore basato sulla complessità (numero di fattori primi diversi)
    const complexityFactor = Object.keys(molData.prime_factors).length * 0.06;

    // Fattore basato sull'esponente totale (somma di tutti gli esponenti)
    const totalExponents = Object.values(molData.prime_factors).reduce((sum, exp) => sum + exp, 0);
    const exponentFactor = totalExponents * 0.04;

    // Calcola dimensione finale con maggiore variabilità
    return baseSize + numberFactor + complexityFactor + exponentFactor;
}

// 6. Analisi avanzata delle proprietà
function analyzeEnhancedMoleculeProperties(molData) {
    // Estrai informazioni dai fattori primi
    const factorCount = Object.keys(molData.prime_factors).length;
    const factors = Object.entries(molData.prime_factors);
    const largestPrime = Math.max(...Object.keys(molData.prime_factors).map(Number));
    const smallestPrime = Math.min(...Object.keys(molData.prime_factors).map(Number));

    // Calcola esponente totale
    const totalExponent = factors.reduce((sum, [_prime, exp]) => sum + exp, 0);

    // Verifica se è un numero primo puro (un solo fattore con esponente 1)
    const isPrimePure = factorCount === 1 && factors[0][1] === 1;

    // Verifica se è una potenza di un primo
    const isPrimePower = factorCount === 1 && factors[0][1] > 1;

    // Verifica se contiene solo primi piccoli (<10)
    const hasOnlySmallPrimes = Object.keys(molData.prime_factors)
        .every(prime => parseInt(prime) < 10);

    // Verifica se contiene grandi primi
    const hasLargePrimes = Object.keys(molData.prime_factors)
        .some(prime => parseInt(prime) > 50);

    // Categorizzazione
    let category;
    if (isPrimePure) {
        if (largestPrime < 10) {
            category = 'elementary'; // Primo puro piccolo (elementare)
        } else if (largestPrime > 50) {
            category = 'exotic'; // Primo puro grande (esotico)
        } else {
            category = 'prime'; // Primo puro medio
        }
    } else if (isPrimePower) {
        category = 'power'; // Potenza di un primo
    } else if (hasOnlySmallPrimes && factorCount <= 2) {
        category = 'stable'; // Molecola semplice stabile
    } else if (hasLargePrimes || factorCount > 3) {
        category = 'complex'; // Molecola complessa
    } else {
        category = 'standard'; // Molecola standard
    }

    // Stabilità basata su vari fattori
    let stability;
    if (isPrimePure || isPrimePower) {
        stability = 'high';
    } else if (hasOnlySmallPrimes && totalExponent < 4) {
        stability = 'medium';
    } else if (hasLargePrimes || totalExponent > 6) {
        stability = 'low';
    } else {
        stability = 'normal';
    }

    // CORREZIONE: Assegna il colore dai dati della molecola
    // Se il colore non è presente, generiamo un colore basato sui fattori primi
    let color;
    if (molData.color) {
        color = molData.color;
    } else {
        // Genera un colore basato sui fattori primi
        const hue = (largestPrime % 360) / 360;
        const saturation = 0.5 + (factorCount * 0.1);
        const lightness = 0.5;
        color = hslToRgb(hue, saturation, lightness);
    }

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
        color: color  // CORREZIONE: Aggiungi il colore all'output
    };
}

// Funzione di supporto per convertire HSL a RGB
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

// 8. UI informativa migliorata
function updateInfoDisplay() {
    // Statistiche simulazione
    infoDiv.innerHTML = `
      <div class="info-section">
        <h3>Simulazione</h3>
        <p>Molecole: <strong>${simulationData.molecules.length}</strong></p>
        <p>Temperatura: <strong>${simulationData.temperature.toFixed(2)}</strong></p>
        <p>Reazioni: <strong>${simulationData.reactionCount || 0}</strong></p>
        <p>Stato: <strong>${paused ? 'PAUSA' : 'ATTIVO'}</strong></p>
      </div>
      
      <div class="controls-info">
        <p>Spazio: Pausa | A: Auto | V: Vettori</p>
        <p>Click: Seleziona molecola</p>
      </div>
    `;
}

// Aggiungi queste funzioni e poi chiamale in init()
// ===== CORREZIONI AL CODICE DI SIMULAZIONE MOLECOLARE =====
// Le modifiche principali riguardano:
// 1. Correzione controlli dell'interfaccia utente
// 2. Rimozione animazioni superflue
// 3. Ottimizzazione dei messaggi tra worker e thread principale

// ===== 1. CORREZIONI ALL'INTERFACCIA UTENTE =====

// Problema principale: le funzioni di controllo non sono globali
// Soluzione: rendere le funzioni di controllo globali ed esporle a window

// MODIFICARE: Riga 1120 circa - Aggiungi queste righe all'inizio della funzione setupUserInterface
function setupUserInterface() {
    // Rendi globali le funzioni principali di controllo
    window.handleControlMessage = handleControlMessage;
    window.updateControlsState = updateControlsState;
    window.pauseSimulation = () => {
        paused = !paused;
        updateControlsState();
    };
    window.toggleVectors = () => {
        showVectors = !showVectors;
        updateControlsState();
        updateMoleculeMeshes();
    };
    
    // Resto del codice rimane invariato
    createInfoContainers();
    addSimulationStyles();
    createControlPanel();
}

// MODIFICARE: Riga 1390 circa - Sostituire bindControlEvents con questa versione
function bindControlEvents() {
    // Slider temperatura
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

    // Slider velocità
    const timeSlider = document.getElementById('timeScaleSlider');
    const timeValue = document.getElementById('timeScaleValue');
    
    if (timeSlider && timeValue) {
        // Imposta valore iniziale
        timeValue.textContent = timeSlider.value;
        
        timeSlider.addEventListener('input', function() {
            const value = parseFloat(this.value).toFixed(1);
            timeValue.textContent = value;
            if (simulationWorker) {
                const timeScaleValue = parseFloat(value);
                console.log(`Impostazione timeScale: ${timeScaleValue}`);
                simulationWorker.postMessage({
                    type: 'set_timescale',
                    value: timeScaleValue
                });
                
                // Aggiorniamo anche una variabile globale per tenere traccia del timeScale attuale
                window.currentTimeScale = timeScaleValue;
            }
        });
    }

    // Pulsante pausa
    const pauseButton = document.getElementById('pauseButton');
    if (pauseButton) {
        pauseButton.addEventListener('click', function() {
            window.pauseSimulation();
        });
    }

    // Pulsante reset
    const resetButton = document.getElementById('resetButton');
    if (resetButton) {
        resetButton.addEventListener('click', function() {
            const spaceDimension = 10;
            const timeSlider = document.getElementById('timeScaleSlider');
            const molecolePerUnit = timeSlider ? parseInt(timeSlider.value) * 5 : 5;
            initializeEnhancedSimulation(spaceDimension, molecolePerUnit);
        });
    }

    // Pulsante aggiungi molecole
    const addButton = document.getElementById('addMoleculesButton');
    if (addButton) {
        addButton.addEventListener('click', function() {
            if (simulationWorker) {
                simulationWorker.postMessage({
                    type: 'add_molecules',
                    count: 20
                });
            }
        });
    }

    // Pulsante vettori
    const vectorsButton = document.getElementById('toggleVectorsButton');
    if (vectorsButton) {
        vectorsButton.addEventListener('click', function() {
            window.toggleVectors();
        });
    }

    // Pulsanti modalità di visualizzazione
    const simpleButton = document.getElementById('simpleModeButton');
    const detailedButton = document.getElementById('detailedModeButton');
    
    if (simpleButton && detailedButton) {
        simpleButton.addEventListener('click', function() {
            setVisualizationMode('simple');
            simpleButton.classList.add('active');
            detailedButton.classList.remove('active');
        });

        detailedButton.addEventListener('click', function() {
            setVisualizationMode('detailed');
            detailedButton.classList.add('active');
            simpleButton.classList.remove('active');
        });
    }

    // Pulsante di chiusura errore
    const dismissButton = document.getElementById('dismissError');
    if (dismissButton) {
        dismissButton.addEventListener('click', function() {
            document.getElementById('errorMessage').style.display = 'none';
        });
    }
    
    // Aggiungi event listener per click e tasti
    document.addEventListener('keydown', onKeyDown);
    renderer.domElement.addEventListener('click', onClick);
}

// MODIFICARE: Sovrascrivere la funzione init per chiamare setupUserInterface dopo la creazione del renderer
function init() {
    const spaceDimension = 10;

    // Inizializza simulazione migliorata
    initializeEnhancedSimulation(spaceDimension);


    // Scene setup
    scene = new THREE.Scene();
    camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
    camera.position.z = spaceDimension;

    renderer = new THREE.WebGLRenderer({
        antialias: true
    });
    renderer.setSize(window.innerWidth, window.innerHeight);
    document.body.appendChild(renderer.domElement);

    // Inizializza UI prima dei controlli per garantire che gli elementi siano creati
    setupUserInterface();

    // Orbit Controls
    controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.25;

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

    // Usa setTimeout per garantire che gli elementi DOM siano creati
    setTimeout(() => {
        bindControlEvents();
    }, 100);

    // Usa setTimeout per verificare lo stato del worker dopo 5 secondi
    setTimeout(checkWorkerStatus, 5000);
    
    // Imposta lo stato iniziale dei controlli
    updateControlsState();

    // Applica le correzioni di emergenza
    ensureGlobalVariables();
    
    // Ricollegamento controlli
    setTimeout(connectAllControlEvents, 200);
    
    // Sostituisci animate
    window.requestAnimationFrame = window.requestAnimationFrame || 
                                 window.mozRequestAnimationFrame ||
                                 window.webkitRequestAnimationFrame || 
                                 window.msRequestAnimationFrame;
    
   requestAnimationFrame(fixedAnimate);
}

// ===== 2. RIMOZIONE EFFETTI VISIVI SUPERFLUI =====

// MODIFICARE: Sostituire createEnhancedMoleculeMesh con questa versione semplificata
function createEnhancedMoleculeMesh(molData) {
    // Analisi essenziale delle proprietà per visualizzazione
    const properties = {
        number: molData.number,
        prime_factors: molData.prime_factors,
        factorCount: Object.keys(molData.prime_factors).length,
        color: molData.color
    };

    // Scegli una dimensione semplice basata sul numero
    const size = 0.2 + 0.1 * Math.log2(molData.number);

    // Usa una geometria semplice: sfera per tutti
    const geometry = new THREE.SphereGeometry(size, 16, 16);

    // Crea un materiale semplice con il colore della molecola
    const material = new THREE.MeshPhongMaterial({
        color: new THREE.Color(...molData.color),
        shininess: 30
    });

    // Crea mesh
    const mesh = new THREE.Mesh(geometry, material);
    mesh.position.set(...molData.position);
    mesh.userData.moleculeData = molData;
    mesh.userData.visualProperties = properties;

    return mesh;
}

// MODIFICARE: Sostituire la funzione updateMoleculeMeshes
function updateMoleculeMeshes() {
    // Process each molecule
    simulationData.molecules.forEach(molData => {
        let mesh = moleculeMeshes.get(molData.id);

        // Create new mesh if it doesn't exist
        if (!mesh) {
            mesh = createEnhancedMoleculeMesh(molData);
            scene.add(mesh);
            moleculeMeshes.set(molData.id, mesh);
        } else {
            // Aggiorna i dati
            mesh.userData.moleculeData = molData;

            // Aggiorna posizione
            mesh.position.set(...molData.position);

            // Aggiungi o aggiorna vettori di velocità se richiesto
            updateVelocityVector(mesh, molData, showVectors);
        }
    });
}

// MODIFICARE: Sostituire la funzione animate per rimuovere aggiornamenti superflui
function animate(timestamp) {
    // Delta time per animazioni fluide
    const deltaTime = lastFrameTime ? timestamp - lastFrameTime : 16.6;
    lastFrameTime = timestamp;

    requestAnimationFrame(animate);
    controls.update();

    // Simulation step
    if (autoStep && !paused && init_done) {
        if (simulationWorker) {
            simulationWorker.postMessage({
                type: 'step'
            });
        }
    }

    // Aggiorna la UI
    updateInfoDisplay();

    // Rendering
    renderer.render(scene, camera);
}

// ===== 3. MIGLIORIE ALLA COMUNICAZIONE CON IL WORKER =====

let init_done = false;

// MODIFICARE: Sovrascrivere initializeEnhancedSimulation
function initializeEnhancedSimulation(size = 10, molecolePerUnit = 5, maxNumber = 150, timeScale = 0.1) {
    console.log("Inizializzazione simulazione avanzata...");

    // Termina worker esistente se presente
    if (simulationWorker) {
        simulationWorker.terminate();
    }

    try {
        // Crea il worker con gestione errori
        simulationWorker = new Worker(new URL('./advanced-simulation-worker.js', import.meta.url), { // advanced-simulation-worker.js
            type: 'module'
        });

        // Gestione errori
        simulationWorker.onerror = function(error) {
            console.error("Errore nel worker:", error);
            displayErrorMessage(`Errore nella simulazione: ${error.message}`);
        };       

        // Gestione messaggi
        simulationWorker.onmessage = handleWorkerMessage;

        // Inizializza parametri simulazione
        const totalMolecules = Math.floor(molecolePerUnit * size);
        console.log(`Inizializzazione con ${totalMolecules} molecole...`);

        // Imposta il valore iniziale della UI
        const timeSlider = document.getElementById('timeScaleSlider');
        if (timeSlider) {
            timeSlider.value = timeScale;
            const timeValue = document.getElementById('timeScaleValue');
            if (timeValue) {
                timeValue.textContent = timeScale.toFixed(1);
            }
        }

        // Imposta pausa a false all'inizializzazione
        paused = false;
        updateControlsState();

        simulationWorker.postMessage({
            type: 'init',
            size: size,
            moleculeCount: totalMolecules,
            maxNumber: maxNumber,
            timeScale: timeScale
        });

        setTimeout(()=>{
            init_done = true;
        }, 500)

        // Debug
        console.log("Messaggio di inizializzazione inviato al worker");
    } catch (error) {
        console.error("Errore critico nell'inizializzazione del worker:", error);
        displayErrorMessage(`Errore nell'inizializzazione: ${error.message}`);
        createTestMolecules(); // Crea molecole di test per far funzionare la UI
    }
}

// Funzione per verificare che le altre funzioni importanti siano accessibili
function initGlobalHandlers() {
    // Assicura che le funzioni di gestione eventi siano accessibili a livello globale
    window.onKeyDown = function(event) {
        switch (event.code) {
            case 'Space':
                // Usa pauseSimulation per gestire la pausa
                window.pauseSimulation();
                break;
            case 'KeyA':
                autoStep = !autoStep;
                break;
            case 'KeyV':
                // Usa toggleVectors
                window.toggleVectors();
                break;
            case 'KeyR': // Reset Camera
                controls.reset();
                break;
            case 'ArrowRight':
                if (paused && simulationWorker) {
                    simulationWorker.postMessage({
                        type: 'step'
                    });
                }
                break;
        }
    };
    
    window.onClick = function(event) {
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

            // Make sure we have molecule data
            if (molIntersected.userData.moleculeData) {
                selectedMolecule = molIntersected.userData.moleculeData;
                const visualProps = molIntersected.userData.visualProperties;

                // Display simplified molecule info
                moleculeInfoDiv.style.display = 'block';
                moleculeInfoDiv.innerHTML = `
                  <h3>Molecola ${selectedMolecule.number}</h3>
                  <div class="molecule-info-container">
                    <div class="basic-properties">
                      <p><strong>Numero:</strong> ${selectedMolecule.number}</p>
                      <p><strong>Fattori Primi:</strong> ${formatPrimeFactors(selectedMolecule.prime_factors)}</p>
                      <p><strong>Massa:</strong> ${selectedMolecule.mass.toFixed(2)}</p>
                      <p><strong>Velocità:</strong> ${vectorMagnitude(selectedMolecule.velocity).toFixed(2)}</p>
                    </div>
                  </div>
                  
                  <div class="position-data">
                    <p><strong>Posizione:</strong> 
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
    };
}

// Aggiungi questa riga alla fine dell'init
document.addEventListener('DOMContentLoaded', function() {
    initGlobalHandlers();
});

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
            // CORREZIONE: Aggiungi colori variabili
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

// Nuova funzione per aggiornare lo stato di tutti i controlli
function updateControlsState() {
    // Aggiorna stato pulsante pausa
    const pauseButton = document.getElementById('pauseButton');
    if (pauseButton) {
        pauseButton.classList.toggle('active', paused);
        pauseButton.textContent = paused ? 'Riprendi' : 'Pausa';
    }
    
    // Aggiorna stato pulsante vettori
    const vectorsButton = document.getElementById('toggleVectorsButton');
    if (vectorsButton) {
        vectorsButton.classList.toggle('active', showVectors);
    }
    
    // Aggiorna modalità visualizzazione
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
            // Correzione pausa con spazio
            console.log("SPAZIO - Toggle pausa");
            window.paused = !window.paused;
            const pauseButton = document.getElementById('pauseButton');
            if (pauseButton) {
                pauseButton.classList.toggle('active', window.paused);
                pauseButton.textContent = window.paused ? 'Riprendi' : 'Pausa';
            }
            break;
        case 'KeyA':
            window.autoStep = !window.autoStep;
            break;
        case 'KeyV':
            // Correzione vettori con V
            console.log("V - Toggle vettori");
            window.showVectors = !window.showVectors;
            const vectorsButton = document.getElementById('toggleVectorsButton');
            if (vectorsButton) {
                vectorsButton.classList.toggle('active', window.showVectors);
            }
            if (typeof updateMoleculeMeshes === 'function') {
                updateMoleculeMeshes();
            }
            break;
        case 'KeyR': // Reset Camera
            if (window.controls) window.controls.reset();
            break;
        case 'ArrowRight':
            if (window.paused && window.simulationWorker) {
                window.simulationWorker.postMessage({
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
        const z = Math.sin(angleOffset + Math.PI / 2) * distance * 0.7;

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
            atom.position.z = Math.sin(angle + Math.PI / 2) * orbit.distance * 0.7;
        }
    };

    mesh.add(atomicStructure);
    mesh.userData.atomicStructure = atomicStructure;
}

///
///
///

function createAdvancedMoleculeVisualization(properties) {
    // Determine the base material type based on category
    const isMaterialComplex = ['exotic', 'elemental', 'prime'].includes(properties.category);
    
    // CORREZIONE: Assicurati che il colore sia sempre disponibile
    // Create base color from properties.color
    const baseColor = properties.color ? 
        new THREE.Color(...properties.color) : 
        new THREE.Color(0.5, 0.7, 0.9); // Colore di fallback blu
    
    // Calculate HSL adjustments
    let hsl = {};
    baseColor.getHSL(hsl);
    
    // Adjust saturation based on stability
    const saturationMap = {
        'high': Math.min(1.0, hsl.s + 0.2),
        'medium': hsl.s,
        'normal': Math.max(0.4, hsl.s - 0.1),
        'low': Math.max(0.3, hsl.s - 0.2)
    };
    const saturation = saturationMap[properties.stability] || hsl.s;
    
    // Adjust luminosity based on category
    const luminosityMap = {
        'elemental': Math.min(0.8, hsl.l + 0.15),
        'exotic': Math.min(0.7, hsl.l + 0.1),
        'prime': Math.min(0.75, hsl.l + 0.12),
        'power': Math.min(0.7, hsl.l + 0.05),
        'complex': Math.max(0.3, hsl.l - 0.1),
        'stable': hsl.l,
        'standard': hsl.l
    };
    const luminosity = luminosityMap[properties.category] || hsl.l;
    
    // Create enhanced color
    const enhancedColor = new THREE.Color().copy(baseColor).setHSL(hsl.h, saturation, luminosity);
    
    // Determine emission properties
    const emissiveIntensity = properties.isPrimePure ? 0.5 : 
                             properties.hasLargePrimes ? 0.3 :
                             properties.category === 'exotic' ? 0.4 :
                             properties.isPrimePower ? 0.25 :
                             properties.stability === 'high' ? 0.2 : 0.1;
    
    // Calculate emissive color (slight hue shift from base)
    const emissiveColor = new THREE.Color().copy(baseColor);
    let emissiveHsl = {};
    emissiveColor.getHSL(emissiveHsl);
    emissiveColor.setHSL(
        (emissiveHsl.h + 0.1) % 1,
        emissiveHsl.s,
        emissiveHsl.l + 0.1
    );
    emissiveColor.multiplyScalar(emissiveIntensity);
    
    let material;
    
    // Create appropriate material based on molecule type
    if (isMaterialComplex) {
        // Use physically-based rendering for complex materials
        const roughness = properties.category === 'elemental' ? 0.2 : 
                         properties.category === 'exotic' ? 0.5 :
                         properties.hasLargePrimes ? 0.6 : 0.4;
        
        const metalness = properties.category === 'elemental' && properties.largestPrime < 10 ? 0.8 :
                         properties.category === 'prime' ? 0.5 :
                         properties.hasLargePrimes ? 0.2 : 0.3;
        
        material = new THREE.MeshStandardMaterial({
            color: enhancedColor,
            emissive: emissiveColor,
            roughness: roughness,
            metalness: metalness,
            envMapIntensity: 1.0,
            flatShading: properties.complexity > 3
        });
        
        // Add subtle normal mapping for texture to complex molecules
        if (properties.complexity > 2) {
            // Simulating bump/normal through vertex displacement would go here
            // For actual implementation, would need to generate or load normal maps
        }
    } else {
        // Use simpler material for standard molecules
        const shininess = properties.stability === 'high' ? 80 :
                         properties.stability === 'medium' ? 60 : 30;
        
        material = new THREE.MeshPhongMaterial({
            color: enhancedColor,
            emissive: emissiveColor,
            specular: new THREE.Color(0xffffff),
            shininess: shininess,
            reflectivity: properties.isPrimePower ? 0.8 : 0.5,
            flatShading: properties.complexity > 4
        });
    }
    
    // Apply transparency for unstable molecules
    if (properties.stability === 'low') {
        material.transparent = true;
        material.opacity = 0.8;
    }
    
    // Special case for very large prime numbers
    if (properties.largestPrime > 100) {
        material.wireframe = true;
        material.wireframeLinewidth = 2;
    }
    
    return material;
}

// Helper function to select appropriate geometry based on molecule properties
function selectMoleculeGeometry(properties, size) {
    // Base geometry selection on molecule properties
    switch (properties.category) {
        case 'elemental':
            // Pure elements are perfect spheres with high detail
            return new THREE.SphereGeometry(size, 32, 32);
            
        case 'exotic':
            // Exotic molecules have crystalline structures
            if (properties.largestPrime > 50) {
                return new THREE.IcosahedronGeometry(size, 1);
            } else {
                return new THREE.OctahedronGeometry(size, 1);
            }
            
        case 'prime':
            // Prime molecules have regular but distinct shapes
            return new THREE.DodecahedronGeometry(size, 0);
            
        case 'power':
            // Power molecules (powers of primes) have geometric perfection
            return new THREE.TorusKnotGeometry(size * 0.7, size * 0.3, 64, 8, 2, 3);
            
        case 'complex':
            // Complex molecules have irregular surfaces
            if (properties.factorCount > 3) {
                // Very complex molecules have more chaotic geometry
                return new THREE.TetrahedronGeometry(size, 1);
            }
            return new THREE.BoxGeometry(size, size, size, 2, 2, 2);
            
        case 'stable':
            // Stable molecules are smooth but not perfect spheres
            return new THREE.SphereGeometry(size, 16, 16);
            
        default:
            // Standard molecules use basic spheres
            return new THREE.SphereGeometry(size, 12, 12);
    }
}

// Function to add special visual effects based on molecule properties
function addSpecialEffects(mesh, properties) {
    // Add glowing effect for elemental molecules
    if (properties.category === 'elemental' || properties.category === 'prime') {
        addElementalGlowEffect(mesh, properties);
    }
    
    // Add orbital rings for exotic molecules
    if (properties.category === 'exotic' || 
       (properties.hasLargePrimes && properties.factorCount > 2)) {
        addOrbitalRingEffect(mesh, properties);
    }
    
    // Add particle system for complex or unstable molecules
    if (properties.category === 'complex' || properties.stability === 'low') {
        addParticleFieldEffect(mesh, properties);
    }
    
    // Add atomic structure visualization for stable compound molecules
    if (properties.factorCount > 1 && 
       (properties.stability === 'high' || properties.stability === 'medium')) {
        addAtomicStructureEffect(mesh, properties);
    }
}

// Implementation of elemental glow effect
function addElementalGlowEffect(mesh, properties) {
    // Create a slightly larger, emissive sphere to create glow effect
    const glowSize = mesh.geometry.parameters.radius * 1.2;
    const glowGeometry = new THREE.SphereGeometry(glowSize, 32, 16);
    
    // Create base color from molData.color or use a default color
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

// Implementation of orbital ring effect
function addOrbitalRingEffect(mesh, properties) {
    // Create an orbital ring
    const ringRadius = mesh.geometry.parameters.radius * 1.8;
    const tubeRadius = 0.05 + (properties.factorCount * 0.01);
    const ringGeometry = new THREE.TorusGeometry(ringRadius, tubeRadius, 8, 50);
    
    // Create base color from molData.color or use a default color
    const baseColor = properties.color ? 
        new THREE.Color(...properties.color) : 
        new THREE.Color(0.5, 0.5, 1.0);
        
    // Choose complementary color for the ring
    const ringColor = new THREE.Color().copy(baseColor);
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
    
    // Randomize orientation
    ring.rotation.x = Math.random() * Math.PI;
    ring.rotation.y = Math.random() * Math.PI;
    
    // Setup rotation parameters
    ring.userData.rotationSpeed = {
        x: (Math.random() - 0.5) * 0.01 * properties.energy / 100,
        y: (Math.random() - 0.5) * 0.01 * properties.energy / 120,
        z: (Math.random() - 0.5) * 0.01 * properties.energy / 140
    };
    
    // Update function for the ring
    ring.userData.updateFunction = function(deltaTime) {
        this.rotation.x += this.userData.rotationSpeed.x * deltaTime;
        this.rotation.y += this.userData.rotationSpeed.y * deltaTime;
        this.rotation.z += this.userData.rotationSpeed.z * deltaTime;
    };
    
    mesh.add(ring);
    mesh.userData.orbitalRing = ring;
}

// Implementation of particle field effect for complex/unstable molecules
function addParticleFieldEffect(mesh, properties) {
    // Number of particles based on complexity and energy
    const particleCount = Math.floor(10 + (properties.complexity * 5) + (properties.energy / 20));
    
    // Create particle system
    const particleGeometry = new THREE.BufferGeometry();
    
    // Create base color from molData.color or use a default color
    const baseColor = properties.color ? 
        new THREE.Color(...properties.color) : 
        new THREE.Color(0.5, 0.5, 1.0);
        
    // Get a color for particles based on base color with HSL shift
    const particleColor = new THREE.Color().copy(baseColor);
    let particleHsl = {};
    particleColor.getHSL(particleHsl);
    particleColor.setHSL(
        (particleHsl.h + 0.1) % 1,
        Math.min(1, particleHsl.s + 0.3),
        Math.min(1, particleHsl.l + 0.2)
    );
    
    const particleMaterial = new THREE.PointsMaterial({
        color: particleColor,
        size: 0.05 + (properties.energy / 500),
        transparent: true,
        opacity: 0.6,
        blending: THREE.AdditiveBlending,
        sizeAttenuation: true
    });
    
    // Create particle positions
    const positions = new Float32Array(particleCount * 3);
    const velocities = [];
    const radius = mesh.geometry.parameters.radius * 1.5;
    
    for (let i = 0; i < particleCount; i++) {
        // Distribute particles in a shell around the molecule
        const theta = Math.random() * Math.PI * 2;
        const phi = Math.acos(2 * Math.random() - 1);
        const r = radius * (0.8 + Math.random() * 0.2);
        
        const x = r * Math.sin(phi) * Math.cos(theta);
        const y = r * Math.sin(phi) * Math.sin(theta);
        const z = r * Math.cos(phi);
        
        positions[i * 3] = x;
        positions[i * 3 + 1] = y;
        positions[i * 3 + 2] = z;
        
        // Create velocity vector for animation
        const speed = 0.002 + (properties.energy / 10000);
        velocities.push(new THREE.Vector3(
            (Math.random() - 0.5) * speed,
            (Math.random() - 0.5) * speed,
            (Math.random() - 0.5) * speed
        ));
    }
    
    particleGeometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    
    const particles = new THREE.Points(particleGeometry, particleMaterial);
    
    // Setup update function
    particles.userData = {
        velocities: velocities,
        radius: radius,
        updateFunction: function(deltaTime) {
            const positions = this.geometry.attributes.position.array;
            
            for (let i = 0; i < particleCount; i++) {
                // Update position
                positions[i * 3] += this.velocities[i].x * deltaTime;
                positions[i * 3 + 1] += this.velocities[i].y * deltaTime;
                positions[i * 3 + 2] += this.velocities[i].z * deltaTime;
                
                // Keep particles within a shell - bounce when hitting boundary
                const x = positions[i * 3];
                const y = positions[i * 3 + 1];
                const z = positions[i * 3 + 2];
                const distance = Math.sqrt(x*x + y*y + z*z);
                
                if (distance > this.radius * 1.2 || distance < this.radius * 0.8) {
                    // Reverse direction if hitting boundary
                    this.velocities[i].x *= -1;
                    this.velocities[i].y *= -1;
                    this.velocities[i].z *= -1;
                }
            }
            
            this.geometry.attributes.position.needsUpdate = true;
        }
    };
    
    mesh.add(particles);
    mesh.userData.particleSystem = particles;
}

// Implementation of atomic structure visualization
function addAtomicStructureEffect(mesh, properties) {
    // Create a group to hold all atoms
    const atomicStructure = new THREE.Group();
    
    // Calculate atom attributes based on prime factors
    const primeFactors = properties.primeFactors;
    const atomRadius = mesh.geometry.parameters.radius * 0.2;
    
    // Create an atom for each prime factor
    primeFactors.forEach((prime, index) => {
        // Calculate orbital parameters
        const orbitRadius = mesh.geometry.parameters.radius * (1.2 + (index * 0.1));
        const orbitSpeed = 0.0003 * (50 / Math.max(prime, 10)); // Larger primes orbit slower
        const orbitOffset = (index / primeFactors.length) * Math.PI * 2;
        const orbitTilt = (index % 3 - 1) * 0.3; // Tilt different orbits
        
        // Calculate atom color based on prime value
        const atomColor = new THREE.Color().setHSL(
            (prime % 12) / 12, // Hue cycles through color wheel
            0.8,
            0.6
        );
        
        // Create atom geometry and material
        const atomGeometry = new THREE.SphereGeometry(atomRadius, 12, 12);
        const atomMaterial = new THREE.MeshBasicMaterial({
            color: atomColor,
            transparent: true,
            opacity: 0.9
        });
        
        const atom = new THREE.Mesh(atomGeometry, atomMaterial);
        
        // Position atom initially
        const angle = orbitOffset;
        atom.position.x = Math.cos(angle) * orbitRadius;
        atom.position.y = Math.sin(angle) * orbitRadius * 0.5 + orbitTilt * orbitRadius;
        atom.position.z = Math.sin(angle + Math.PI/2) * orbitRadius * 0.7;
        
        // Store orbital parameters
        atom.userData.orbit = {
            radius: orbitRadius,
            speed: orbitSpeed,
            offset: orbitOffset,
            tilt: orbitTilt
        };
        
        atomicStructure.add(atom);
    });
    
    // Setup update function
    atomicStructure.userData.updateFunction = function(time) {
        this.children.forEach(atom => {
            const orbit = atom.userData.orbit;
            const angle = time * orbit.speed + orbit.offset;
            
            // Update position based on orbital parameters
            atom.position.x = Math.cos(angle) * orbit.radius;
            atom.position.y = Math.sin(angle) * orbit.radius * 0.5 + orbit.tilt * orbit.radius;
            atom.position.z = Math.sin(angle + Math.PI/2) * orbit.radius * 0.7;
        });
    };
    
    mesh.add(atomicStructure);
    mesh.userData.atomicStructure = atomicStructure;
}

///
///
///

// Funzioni per gestire stili e interfaccia utente della simulazione

function addSimulationStyles() {
    const styleElement = document.createElement('style');
    styleElement.textContent = `
      /* Stili generali */
      body {
        margin: 0;
        overflow: hidden;
        font-family: 'Arial', sans-serif;
      }
  
      /* Pannello informazioni */
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
  
      /* Pannello dettagli molecola */
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
  
      /* Pannello controlli */
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
  
      /* Effetto tooltip */
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
      
      /* Messaggio di errore */
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

function createControlPanel() {
    // Aggiungi elemento panel
    const panel = document.createElement('div');
    panel.id = 'controlPanel';
    panel.innerHTML = `
      <h3>Controlli</h3>
      
      <div class="control-group">
        <label class="control-label" for="temperatureSlider">Temperatura</label>
        <div class="slider-container">
          <input type="range" min="0.5" max="2.0" step="0.1" value="1.0" class="control-slider" id="temperatureSlider">
          <span class="slider-value" id="temperatureValue">1.0</span>
        </div>
      </div>
      
      <div class="control-group">
        <label class="control-label" for="timeScaleSlider">Velocità simulazione</label>
        <div class="slider-container">
          <input type="range" min="0.2" max="2.0" step="0.1" value="0.1" class="control-slider" id="timeScaleSlider">
          <span class="slider-value" id="timeScaleValue">0.1</span>
        </div>
      </div>
      
      <div class="control-group">
        <button id="pauseButton" class="control-button">Pausa</button>
        <button id="resetButton" class="control-button">Reset</button>
        <button id="addMoleculesButton" class="control-button tooltip">
          Aggiungi molecole
          <span class="tooltiptext">Aggiunge 20 nuove molecole</span>
        </button>
        <button id="toggleVectorsButton" class="control-button">Vettori</button>
      </div>
      
      <div class="control-group">
        <label class="control-label">Visualizzazione</label>
        <button id="simpleModeButton" class="control-button active">Base</button>
        <button id="detailedModeButton" class="control-button">Dettagliata</button>
      </div>
    `;

    // Aggiungi al DOM
    document.body.appendChild(panel);

    // Aggiungi elementi per messaggi di errore
    const errorDiv = document.createElement('div');
    errorDiv.id = 'errorMessage';
    errorDiv.innerHTML = `
      <p id="errorText">Errore nella simulazione</p>
      <button id="dismissError">Chiudi</button>
    `;
    document.body.appendChild(errorDiv);

    // PARTE CRITICA: collega immediatamente gli eventi
    connectAllControlEvents();
}

// 2. NUOVA FUNZIONE: collega TUTTI gli eventi direttamente 
function connectAllControlEvents() {
    console.log("EMERGENZA: Collegamento diretto eventi controlli");
    
    // PAUSA BUTTON
    const pauseButton = document.getElementById('pauseButton');
    if (pauseButton) {
        pauseButton.onclick = function() {
            console.log("CLIC PAUSA");
            window.paused = !window.paused;
            this.classList.toggle('active', window.paused);
            this.textContent = window.paused ? 'Riprendi' : 'Pausa';
        };
    }
    
    // VETTORI BUTTON  
    const vectorsButton = document.getElementById('toggleVectorsButton');
    if (vectorsButton) {
        vectorsButton.onclick = function() {
            console.log("CLIC VETTORI");
            window.showVectors = !window.showVectors;
            this.classList.toggle('active', window.showVectors);
            if (typeof updateMoleculeMeshes === 'function') {
                updateMoleculeMeshes();
            }
        };
    }
    
    // RESET BUTTON
    const resetButton = document.getElementById('resetButton');
    if (resetButton) {
        resetButton.onclick = function() {
            console.log("CLIC RESET");
            const spaceDimension = 10;
            if (typeof initializeEnhancedSimulation === 'function') {
                initializeEnhancedSimulation(spaceDimension, 5);
            }
        };
    }
    
    // AGGIUNGI MOLECOLE
    const addButton = document.getElementById('addMoleculesButton');
    if (addButton) {
        addButton.onclick = function() {
            console.log("CLIC AGGIUNGI MOLECOLE");
            if (window.simulationWorker) {
                window.simulationWorker.postMessage({
                    type: 'add_molecules',
                    count: 20
                });
            }
        };
    }
    
    // MODALITÀ VISUALIZZAZIONE
    const simpleButton = document.getElementById('simpleModeButton');
    const detailedButton = document.getElementById('detailedModeButton');
    
    if (simpleButton) {
        simpleButton.onclick = function() {
            console.log("CLIC VISUALIZZAZIONE SEMPLICE");
            window.visualizationMode = 'simple';
            simpleButton.classList.add('active');
            if (detailedButton) detailedButton.classList.remove('active');
            if (window.simulationWorker) {
                window.simulationWorker.postMessage({
                    type: 'set_visualization',
                    mode: 'simple'
                });
            }
            if (typeof updateMoleculeMeshes === 'function') {
                updateMoleculeMeshes();
            }
        };
    }
    
    if (detailedButton) {
        detailedButton.onclick = function() {
            console.log("CLIC VISUALIZZAZIONE DETTAGLIATA");
            window.visualizationMode = 'detailed';
            detailedButton.classList.add('active');
            if (simpleButton) simpleButton.classList.remove('active');
            if (window.simulationWorker) {
                window.simulationWorker.postMessage({
                    type: 'set_visualization',
                    mode: 'detailed'
                });
            }
            if (typeof updateMoleculeMeshes === 'function') {
                updateMoleculeMeshes();
            }
        };
    }
    
    // SLIDER TEMPERATURA
    const tempSlider = document.getElementById('temperatureSlider');
    const tempValue = document.getElementById('temperatureValue');
    
    if (tempSlider && tempValue) {
        tempSlider.oninput = function() {
            console.log("MODIFICA TEMPERATURA");
            const value = parseFloat(this.value).toFixed(1);
            tempValue.textContent = value;
            if (window.simulationWorker) {
                window.simulationWorker.postMessage({
                    type: 'set_temperature',
                    value: parseFloat(value)
                });
            }
        };
    }
    
    // SLIDER VELOCITÀ
    const timeSlider = document.getElementById('timeScaleSlider');
    const timeValue = document.getElementById('timeScaleValue');
    
    if (timeSlider && timeValue) {
        timeSlider.oninput = function() {
            console.log("MODIFICA TIMESCALE");
            const value = parseFloat(this.value).toFixed(1);
            timeValue.textContent = value;
            if (window.simulationWorker) {
                window.simulationWorker.postMessage({
                    type: 'set_timescale',
                    value: parseFloat(value)
                });
                window.currentTimeScale = parseFloat(value);
            }
        };
    }
    
    // CHIUDI ERRORE
    const dismissButton = document.getElementById('dismissError');
    if (dismissButton) {
        dismissButton.onclick = function() {
            document.getElementById('errorMessage').style.display = 'none';
        };
    }
}

// 4. ASSICURIAMO CHE TUTTE LE VARIABILI FONDAMENTALI SIANO GLOBALI
function ensureGlobalVariables() {
    if (typeof window.paused === 'undefined') window.paused = false;
    if (typeof window.showVectors === 'undefined') window.showVectors = false; 
    if (typeof window.autoStep === 'undefined') window.autoStep = true;
    if (typeof window.visualizationMode === 'undefined') window.visualizationMode = 'simple';
    if (typeof window.scene === 'undefined' && typeof scene !== 'undefined') window.scene = scene;
    if (typeof window.camera === 'undefined' && typeof camera !== 'undefined') window.camera = camera;
    if (typeof window.renderer === 'undefined' && typeof renderer !== 'undefined') window.renderer = renderer;
    if (typeof window.controls === 'undefined' && typeof controls !== 'undefined') window.controls = controls;
    if (typeof window.simulationData === 'undefined' && typeof simulationData !== 'undefined') window.simulationData = simulationData;
    if (typeof window.moleculeMeshes === 'undefined' && typeof moleculeMeshes !== 'undefined') window.moleculeMeshes = moleculeMeshes;
    if (typeof window.simulationWorker === 'undefined' && typeof simulationWorker !== 'undefined') window.simulationWorker = simulationWorker;
}

// 5. SOSTITUIAMO LA FUNZIONE animate CON UNA VERSIONE PIÙ ROBUSTA
function fixedAnimate(timestamp) {
    // Garantisci che le variabili siano accessibili
    ensureGlobalVariables();
    
    // Delta time per animazioni fluide
    const deltaTime = window.lastFrameTime ? timestamp - window.lastFrameTime : 16.6;
    window.lastFrameTime = timestamp;

    requestAnimationFrame(fixedAnimate);
    
    if (window.controls) window.controls.update();

    // Simulation step
    if (window.autoStep && !window.paused && init_done) {
        if (window.simulationWorker) {
            window.simulationWorker.postMessage({
                type: 'step'
            });
        }
    }

    // Aggiorna la UI
    if (typeof updateInfoDisplay === 'function') {
        updateInfoDisplay();
    }

    // Rendering
    if (window.renderer && window.scene && window.camera) {
        window.renderer.render(window.scene, window.camera);
    }
}

// Assicura che questo script venga eseguito all'avvio della pagina
document.addEventListener('DOMContentLoaded', function() {
    console.log("EMERGENZA: Applicazione correzioni UI");
    ensureGlobalVariables();
    // Ricollegamento controlli se la pagina è già caricata
    setTimeout(connectAllControlEvents, 200);
});

function setVisualizationMode(mode) {
    visualizationMode = mode;
    if (simulationWorker) {
        simulationWorker.postMessage({
            type: 'set_visualization',
            mode: mode
        });
    }
    updateMoleculeMeshes();
}

// Implementazione di funzioni mancanti per il pannello di controllo
let visualizationMode = 'simple'; // Default mode

function setupSimulationStructure() {
    // Aggiungi i contenitori DOM necessari
    createInfoContainers();

    // Aggiungi gli stili CSS
    addSimulationStyles();

    // Crea il pannello di controllo
    createControlPanel();

    // Inizializza gesture handler per mobile/touch
    initTouchControls();

    // Imposta la gestione del ridimensionamento della finestra
    setupWindowResize();
}

function createInfoContainers() {
    // Crea il container info se non esiste
    if (!document.getElementById('info')) {
        const infoDiv = document.createElement('div');
        infoDiv.id = 'info';
        document.body.appendChild(infoDiv);
    }

    // Crea il container per i dettagli molecola se non esiste
    if (!document.getElementById('moleculeInfo')) {
        const moleculeInfoDiv = document.createElement('div');
        moleculeInfoDiv.id = 'moleculeInfo';
        document.body.appendChild(moleculeInfoDiv);
    }
}

function handleControlMessage(messageType, value) {
    if (!simulationWorker) return;

    switch (messageType) {
        case 'temperature':
            simulationWorker.postMessage({
                type: 'set_temperature',
                value: parseFloat(value)
            });
            break;

        case 'timescale':
            simulationWorker.postMessage({
                type: 'set_timescale',
                value: parseFloat(value)
            });
            break;

        case 'reset':
            const spaceDimension = 10;
            let molecolePerUnit = 5;
            let currentTimeScale = 1.0;
            const timeSlider = document.getElementById('timeScaleSlider');
            if (timeSlider) {
                currentTimeScale = parseFloat(timeSlider.value);
                molecolePerUnit = Math.round(currentTimeScale * 5);
            }
            // Passa esplicitamente il timeScale all'inizializzazione
            initializeEnhancedSimulation(spaceDimension, molecolePerUnit, 150, currentTimeScale);
            break;

        case 'add_molecules':
            simulationWorker.postMessage({
                type: 'add_molecules',
                count: value || 20
            });
            break;

        case 'toggle_pause':
            paused = !paused;
            const pauseButton = document.getElementById('pauseButton');
            if (pauseButton) {
                pauseButton.classList.toggle('active', paused);
                pauseButton.textContent = paused ? 'Riprendi' : 'Pausa';
            }
            break;

        case 'toggle_vectors':
            showVectors = !showVectors;
            const vectorsButton = document.getElementById('toggleVectorsButton');
            if (vectorsButton) {
                vectorsButton.classList.toggle('active', showVectors);
            }
            updateMoleculeMeshes();
            break;

        case 'set_visualization':
            setVisualizationMode(value);
            break;
    }
}

function initTouchControls() {
    // Implementa controlli touch per dispositivi mobile
    let touchStartX = 0;
    let touchStartY = 0;
    let touchStartDistance = 0;

    // Gestisci pinch-to-zoom
    renderer.domElement.addEventListener('touchstart', function(event) {
        if (event.touches.length === 1) {
            touchStartX = event.touches[0].clientX;
            touchStartY = event.touches[0].clientY;
        } else if (event.touches.length === 2) {
            // Calcola distanza iniziale per pinch
            const dx = event.touches[0].clientX - event.touches[1].clientX;
            const dy = event.touches[0].clientY - event.touches[1].clientY;
            touchStartDistance = Math.sqrt(dx * dx + dy * dy);
        }
    });

    renderer.domElement.addEventListener('touchmove', function(event) {
        if (event.touches.length === 2) {
            event.preventDefault(); // Previeni scroll

            // Calcola nuova distanza per pinch
            const dx = event.touches[0].clientX - event.touches[1].clientX;
            const dy = event.touches[0].clientY - event.touches[1].clientY;
            const newDistance = Math.sqrt(dx * dx + dy * dy);

            // Applica zoom
            if (touchStartDistance > 0) {
                const delta = touchStartDistance - newDistance;
                camera.position.z += delta * 0.01;
                camera.position.z = Math.max(1, Math.min(30, camera.position.z));
                touchStartDistance = newDistance;
            }
        } else if (event.touches.length === 1) {
            // Rotazione camera
            const touchX = event.touches[0].clientX;
            const touchY = event.touches[0].clientY;
            const deltaX = touchX - touchStartX;
            const deltaY = touchY - touchStartY;

            controls.rotateLeft(deltaX * 0.005);
            controls.rotateUp(deltaY * 0.005);

            touchStartX = touchX;
            touchStartY = touchY;
        }
    });

    // Doppio tap per selezionare molecola
    let lastTap = 0;
    renderer.domElement.addEventListener('touchend', function(event) {
        const currentTime = new Date().getTime();
        const tapLength = currentTime - lastTap;
        if (tapLength < 300 && tapLength > 0) {
            // Doppio tap
            const touch = event.changedTouches[0];
            const rect = renderer.domElement.getBoundingClientRect();
            const x = (touch.clientX - rect.left) / rect.width * 2 - 1;
            const y = -((touch.clientY - rect.top) / rect.height) * 2 + 1;

            selectMoleculeAtPosition(x, y);
        }
        lastTap = currentTime;
    });
}

function selectMoleculeAtPosition(x, y) {
    const raycaster = new THREE.Raycaster();
    const mouse = new THREE.Vector2(x, y);

    raycaster.setFromCamera(mouse, camera);
    const meshArray = Array.from(moleculeMeshes.values());
    const intersects = raycaster.intersectObjects(meshArray);

    if (intersects.length > 0) {
        const mesh = intersects[0].object;
        selectedMolecule = mesh.userData.moleculeData;
        updateMoleculeInfoPanel(selectedMolecule, mesh.userData.visualProperties);
        moleculeInfoDiv.style.display = 'block';
        highlightSelectedMolecule(mesh);
    } else {
        selectedMolecule = null;
        moleculeInfoDiv.style.display = 'none';
        clearMoleculeHighlight();
    }
}

function setupWindowResize() {
    window.addEventListener('resize', function() {
        camera.aspect = window.innerWidth / window.innerHeight;
        camera.updateProjectionMatrix();
        renderer.setSize(window.innerWidth, window.innerHeight);

        // Riposiziona UI se necessario
        positionUIElements();
    });
}

function positionUIElements() {
    // Adatta posizione UI in base a dimensione schermo
    const isMobile = window.innerWidth < 768;

    if (isMobile) {
        // Adatta per mobile
        const infoDiv = document.getElementById('info');
        if (infoDiv) {
            infoDiv.style.maxWidth = '120px';
            infoDiv.style.fontSize = '12px';
        }

        const controlPanel = document.getElementById('controlPanel');
        if (controlPanel) {
            controlPanel.style.width = '140px';
        }

        const moleculeInfo = document.getElementById('moleculeInfo');
        if (moleculeInfo) {
            moleculeInfo.style.maxWidth = '280px';
            moleculeInfo.style.bottom = '10px';
            moleculeInfo.style.right = '10px';
        }
    } else {
        // Ripristina per desktop
        const infoDiv = document.getElementById('info');
        if (infoDiv) {
            infoDiv.style.maxWidth = '250px';
            infoDiv.style.fontSize = '14px';
        }

        const controlPanel = document.getElementById('controlPanel');
        if (controlPanel) {
            controlPanel.style.width = '180px';
        }

        const moleculeInfo = document.getElementById('moleculeInfo');
        if (moleculeInfo) {
            moleculeInfo.style.maxWidth = '400px';
            moleculeInfo.style.bottom = '20px';
            moleculeInfo.style.right = '20px';
        }
    }
}

function updateMoleculeInfoPanel(molecule, properties) {
    if (!molecule || !moleculeInfoDiv) return;

    // Crea rappresentazione dei fattori primi formattata
    const formattedFactors = formatPrimeFactors(molecule.prime_factors);

    // Determina categoria e stabilità in modo descrittivo
    const categoryDescriptions = {
        'elementary': 'Elementare (primo puro piccolo)',
        'exotic': 'Esotico (primo puro grande)',
        'prime': 'Primo (numero primo)',
        'power': 'Potenza (potenza di un primo)',
        'stable': 'Stabile (composto semplice)',
        'complex': 'Complesso (struttura articolata)',
        'standard': 'Standard'
    };

    const stabilityDescriptions = {
        'high': 'Alta',
        'medium': 'Media',
        'normal': 'Normale',
        'low': 'Bassa'
    };

    const categoryText = categoryDescriptions[properties.category] || properties.category;
    const stabilityText = stabilityDescriptions[properties.stability] || properties.stability;

    // Aggiorna UI con informazioni dettagliate
    moleculeInfoDiv.innerHTML = `
    <h3>Molecola ${molecule.number}</h3>
    <div class="molecule-info-container">
      <div class="basic-properties">
        <p><strong>Numero:</strong> ${molecule.number}</p>
        <p><strong>Fattori Primi:</strong> ${formattedFactors}</p>
        <p><strong>Massa:</strong> ${molecule.mass.toFixed(2)}</p>
        <p><strong>Carica:</strong> ${molecule.charge.toFixed(2)}</p>
        <p><strong>Velocità:</strong> ${vectorMagnitude(molecule.velocity).toFixed(2)}</p>
      </div>
      
      <div class="advanced-properties">
        <p><strong>Categoria:</strong> ${categoryText}</p>
        <p><strong>Stabilità:</strong> ${stabilityText}</p>
        <p><strong>Energia:</strong> ${properties.energy.toFixed(2)}</p>
        <p><strong>Complessità:</strong> ${properties.complexity.toFixed(1)}</p>
        <p><strong>N° fattori:</strong> ${properties.factorCount}</p>
      </div>
    </div>
    
    <div class="position-data">
      <p><strong>Posizione:</strong> 
        (${molecule.position[0].toFixed(2)}, 
         ${molecule.position[1].toFixed(2)}, 
         ${molecule.position[2].toFixed(2)})</p>
    </div>
  `;
}

///
///
///

init();
