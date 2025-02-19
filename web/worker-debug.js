// Crea un nuovo file chiamato worker-debug.js

// Questo è un wrapper per il worker che consente di catturare errori di importazione
console.log("Inizializzazione worker-debug.js");

// Cattura errori globali
self.onerror = function(message, source, lineno, colno, error) {
  console.error("Worker error:", {message, source, lineno, colno, error});
  self.postMessage({
    type: 'error',
    message: message,
    source: source,
    lineno: lineno,
    colno: colno,
    stack: error ? error.stack : null
  });
  return true;
};

// Funzione per verificare che i moduli siano raggiungibili
async function testImports() {
  try {
    console.log("Worker: Verifico import di rules.js");
    const rulesModule = await import('./rules.js');
    console.log("Worker: rules.js importato correttamente", rulesModule);
    
    console.log("Worker: Verifico import di molecule.js");
    const moleculeModule = await import('./molecule.js');
    console.log("Worker: molecule.js importato correttamente", moleculeModule);
    
    return true;
  } catch (error) {
    console.error("Worker: Errore importazione moduli:", error);
    self.postMessage({
      type: 'error',
      message: 'Errore importazione moduli',
      details: error.message,
      stack: error.stack
    });
    return false;
  }
}

// Intercetta i messaggi prima di importare il worker reale
self.addEventListener('message', async function(event) {
  if (event.data.type === 'init') {
    console.log("Worker-debug: Ricevuto messaggio init");
    
    // Verifica gli import
    const importsOk = await testImports();
    if (!importsOk) {
      return; // Termina se gli import falliscono
    }
    
    // Prova a importare il worker reale
    try {
      const { default: initRealWorker } = await import('./worker.js');
      console.log("Worker: worker.js importato correttamente");
      
      // Crea alcune molecole di test
      self.postMessage({
        type: 'update',
        molecules: {
          molecules: Array.from({length: 5}, (_, i) => ({
            id: `test-${i}`,
            number: 10 + i,
            position: [
              (Math.random() - 0.5) * 8,
              (Math.random() - 0.5) * 8,
              (Math.random() - 0.5) * 8
            ],
            velocity: [0, 0, 0],
            prime_factors: {2: 1, 5: 1},
            mass: 5,
            charge: 0,
            color: [Math.random(), Math.random(), Math.random()],
            angularVelocity: [0, 0, 0],
            lastReactionTime: -1
          })),
          removedIds: []
        },
        temperature: 1.0
      });
    } catch (error) {
      console.error("Worker: Errore caricamento worker principale:", error);
      self.postMessage({
        type: 'error',
        message: 'Errore caricamento worker principale',
        details: error.message,
        stack: error.stack
      });
    }
  }
}, false);

console.log("Worker-debug.js inizializzato");