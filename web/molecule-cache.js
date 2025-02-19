/**
 * MoleculeCentricCache - Sistema di caching avanzato per molecole
 * Mantiene le relazioni tra molecole e include un meccanismo di selezione casuale
 * per garantire che, nel lungo termine, vengano considerate tutte le possibili interazioni
 */
class MoleculeCentricCache {
    constructor() {
      // Mappa di ID molecola alle sue relazioni
      this.moleculeRelationships = new Map();
      
      // Soglia di obsolescenza relazioni (ms)
      this.stalenessThreshold = 1000*30;
      
      // Massimo numero di relazioni per molecola
      this.maxRelationshipsPerMolecule = 20;
      
      // Probabilità di considerare molecole casuali
      this.randomInteractionProbability = 0.05;
      
      // Cache delle interazioni recenti per evitare duplicati
      this.recentInteractions = new Set();
      
      // Contatore per tracciare l'efficienza del caching
      this.cacheHits = 0;
      this.cacheMisses = 0;
    }
    
    // Restituisce tutte le relazioni per una molecola
    getRelationshipsForMolecule(molecule) {
      if (!molecule || !molecule.id) return [];
      return this.moleculeRelationships.get(molecule.id) || [];
    }
    
    // Restituisce tutte le relazioni per una lista di molecole
    getRelationshipsForMolecules(molecules) {
      const result = {};
      
      for (const mol of molecules) {
        if (!mol.id) continue;
        
        const relationships = this.getRelationshipsForMolecule(mol);
        if (relationships.length > 0) {
          result[mol.id] = relationships.map(rel => ({
            otherId: rel.otherId,
            distance: rel.distance,
            lastUpdated: rel.lastUpdated
          }));
          this.cacheHits += relationships.length;
        } else {
          this.cacheMisses++;
        }
      }
      
      return result;
    }
    
    // Ottieni molecole senza cache da una lista, considerando anche molecole casuali
    getUncachedMolecules(molecules) {
      // Molecole che devono essere considerate perché non sono nella cache
      const uncachedMolecules = molecules.filter(mol => 
        !mol.id || !this.moleculeRelationships.has(mol.id) || 
        this.moleculeRelationships.get(mol.id).length === 0
      );
      
      // Considera anche un sottoinsieme casuale di molecole già nella cache
      // per garantire che, col tempo, tutte le interazioni vengano considerate
      const cachedMolecules = molecules.filter(mol => 
        mol.id && this.moleculeRelationships.has(mol.id) && 
        this.moleculeRelationships.get(mol.id).length > 0
      );
      
      // Seleziona un sottoinsieme casuale di molecole già nella cache
      const randomSelection = cachedMolecules.filter(() => 
        Math.random() < this.randomInteractionProbability
      );
      
      // Unisci le molecole non nella cache con la selezione casuale
      return [...uncachedMolecules, ...randomSelection];
    }
    
    // Crea una nuova relazione tra molecole
    createRelationship(mol1, mol2, distance, timestamp) {
      if (!mol1.id || !mol2.id || mol1.id === mol2.id) return;
      
      // Crea chiave di interazione unica
      const interactionKey = this._getInteractionKey(mol1.id, mol2.id);
      
      // Evita di ricreare interazioni troppo recenti
      if (this.recentInteractions.has(interactionKey)) return;
      
      // Verifica se la relazione esiste già
      if (this.hasRelationship(mol1.id, mol2.id)) {
        // Aggiorna relazione esistente
        this.updateRelationship(mol1.id, mol2.id, distance, timestamp);
        return;
      }
      
      // Crea nuove relazioni in entrambe le direzioni
      this._ensureRelationshipsList(mol1.id);
      this._ensureRelationshipsList(mol2.id);
      
      // Verifica se è necessario eliminare relazioni
      this._pruneRelationshipsIfNeeded(mol1.id);
      this._pruneRelationshipsIfNeeded(mol2.id);
      
      // Crea relazione bidirezionale
      const rel1 = {
        otherId: mol2.id,
        distance: distance,
        lastUpdated: timestamp
      };
      
      const rel2 = {
        otherId: mol1.id,
        distance: distance,
        lastUpdated: timestamp
      };
      
      this.moleculeRelationships.get(mol1.id).push(rel1);
      this.moleculeRelationships.get(mol2.id).push(rel2);
      
      // Aggiungi alle interazioni recenti
      this.recentInteractions.add(interactionKey);
      
      // Limita la dimensione delle interazioni recenti
      if (this.recentInteractions.size > 1000) {
        // Rimuovi un elemento casuale (approssimazione efficiente)
        const keys = Array.from(this.recentInteractions);
        this.recentInteractions.delete(keys[0]);
      }
    }
    
    // Genera una chiave di interazione unica
    _getInteractionKey(mol1Id, mol2Id) {
      return mol1Id < mol2Id ? `${mol1Id}:${mol2Id}` : `${mol2Id}:${mol1Id}`;
    }
    
    // Verifica se esiste una relazione
    hasRelationship(mol1Id, mol2Id) {
      if (!this.moleculeRelationships.has(mol1Id)) return false;
      
      return this.moleculeRelationships.get(mol1Id).some(rel => 
        rel.otherId === mol2Id
      );
    }
    
    // Aggiorna una relazione esistente
    updateRelationship(mol1Id, mol2Id, distance, timestamp) {
      if (!this.hasRelationship(mol1Id, mol2Id)) return;
      
      // Aggiorna da mol1 a mol2
      const rels1 = this.moleculeRelationships.get(mol1Id);
      const rel1Index = rels1.findIndex(rel => rel.otherId === mol2Id);
      
      if (rel1Index >= 0) {
        rels1[rel1Index].distance = distance;
        rels1[rel1Index].lastUpdated = timestamp;
      }
      
      // Aggiorna da mol2 a mol1
      if (this.moleculeRelationships.has(mol2Id)) {
        const rels2 = this.moleculeRelationships.get(mol2Id);
        const rel2Index = rels2.findIndex(rel => rel.otherId === mol1Id);
        
        if (rel2Index >= 0) {
          rels2[rel2Index].distance = distance;
          rels2[rel2Index].lastUpdated = timestamp;
        }
      }
    }
    
    // Rimuove una relazione
    removeRelationship(mol1Id, mol2Id) {
      // Rimuove da mol1 a mol2
      if (this.moleculeRelationships.has(mol1Id)) {
        const rels1 = this.moleculeRelationships.get(mol1Id);
        const filteredRels1 = rels1.filter(rel => rel.otherId !== mol2Id);
        this.moleculeRelationships.set(mol1Id, filteredRels1);
      }
      
      // Rimuove da mol2 a mol1
      if (this.moleculeRelationships.has(mol2Id)) {
        const rels2 = this.moleculeRelationships.get(mol2Id);
        const filteredRels2 = rels2.filter(rel => rel.otherId !== mol1Id);
        this.moleculeRelationships.set(mol2Id, filteredRels2);
      }
      
      // Rimuovi dall'elenco delle interazioni recenti
      const interactionKey = this._getInteractionKey(mol1Id, mol2Id);
      this.recentInteractions.delete(interactionKey);
    }
    
    // Rimuove tutte le relazioni per una molecola
    removeAllRelationshipsForMolecule(molId) {
      if (!this.moleculeRelationships.has(molId)) return;
      
      // Ottieni lista delle altre molecole
      const relationships = this.moleculeRelationships.get(molId);
      const otherIds = relationships.map(rel => rel.otherId);
      
      // Rimuovi questa molecola dalle relazioni di tutte le altre molecole
      for (const otherId of otherIds) {
        if (this.moleculeRelationships.has(otherId)) {
          const otherRels = this.moleculeRelationships.get(otherId);
          const filteredOtherRels = otherRels.filter(rel => rel.otherId !== molId);
          this.moleculeRelationships.set(otherId, filteredOtherRels);
          
          // Rimuovi dalle interazioni recenti
          const interactionKey = this._getInteractionKey(molId, otherId);
          this.recentInteractions.delete(interactionKey);
        }
      }
      
      // Rimuovi le relazioni di questa molecola
      this.moleculeRelationships.delete(molId);
    }
    
    // Pulisci relazioni obsolete
    cleanupStaleRelationships() {
      const now = performance.now();
      
      for (const [molId, relationships] of this.moleculeRelationships.entries()) {
        // Filtra le relazioni obsolete
        const freshRelationships = relationships.filter(rel => 
          now - rel.lastUpdated < this.stalenessThreshold
        );
        
        if (freshRelationships.length !== relationships.length) {
          // Aggiorna con le sole relazioni fresche
          this.moleculeRelationships.set(molId, freshRelationships);
        }
      }
      
      // Rimuovi liste di relazioni vuote
      for (const [molId, relationships] of this.moleculeRelationships.entries()) {
        if (relationships.length === 0) {
          this.moleculeRelationships.delete(molId);
        }
      }
      
      // Pulizia periodica dell'insieme delle interazioni recenti
      if (this.recentInteractions.size > 500) {
        this.recentInteractions.clear();
      }
    }
    
    // Cancella tutte le relazioni
    clear() {
      this.moleculeRelationships.clear();
      this.recentInteractions.clear();
      this.cacheHits = 0;
      this.cacheMisses = 0;
    }
    
    // Ottiene statistiche del cache
    getStats() {
      const totalRequests = this.cacheHits + this.cacheMisses;
      const hitRate = totalRequests > 0 ? this.cacheHits / totalRequests : 0;
      
      return {
        totalRelationships: this._countTotalRelationships(),
        uniqueMolecules: this.moleculeRelationships.size,
        hitRate: hitRate,
        recentInteractionsSize: this.recentInteractions.size
      };
    }
    
    // Conta il numero totale di relazioni
    _countTotalRelationships() {
      let count = 0;
      for (const relationships of this.moleculeRelationships.values()) {
        count += relationships.length;
      }
      return count;
    }
    
    // Helper per assicurarsi che una molecola abbia una lista di relazioni
    _ensureRelationshipsList(molId) {
      if (!this.moleculeRelationships.has(molId)) {
        this.moleculeRelationships.set(molId, []);
      }
    }
    
    // Helper per potare relazioni in eccesso, mantenendo un mix di nuove e vecchie
    _pruneRelationshipsIfNeeded(molId) {
      if (!this.moleculeRelationships.has(molId)) return;
      
      const relationships = this.moleculeRelationships.get(molId);
      if (relationships.length >= this.maxRelationshipsPerMolecule) {
        // Ordina per recenza (più recenti prima)
        relationships.sort((a, b) => b.lastUpdated - a.lastUpdated);
        
        // Strategia ibrida: mantieni alcune relazioni recenti e alcune casuali
        const recentCount = Math.floor(this.maxRelationshipsPerMolecule * 0.8);
        const recentRelationships = relationships.slice(0, recentCount);
        
        // Seleziona alcune relazioni casuali dalle rimanenti
        const remainingRelationships = relationships.slice(recentCount);
        const randomRelationships = [];
        
        const randomCount = this.maxRelationshipsPerMolecule - recentCount;
        for (let i = 0; i < Math.min(randomCount, remainingRelationships.length); i++) {
          // Scegli un indice casuale dalle relazioni rimanenti
          const randomIndex = Math.floor(Math.random() * remainingRelationships.length);
          randomRelationships.push(remainingRelationships[randomIndex]);
          // Rimuovi per evitare duplicati
          remainingRelationships.splice(randomIndex, 1);
        }
        
        // Combina relazioni recenti e casuali
        this.moleculeRelationships.set(molId, [...recentRelationships, ...randomRelationships]);
      }
    }
  }

  /**
 * InteractionCache - Versione migliorata del sistema di caching per interazioni fisiche
 * Include sistema di selezione casuale per garantire copertura completa nel lungo termine
 */
class InteractionCache {
    constructor(lifetime) {
        this.cache = new Map();
        this.lifetime = lifetime || 10;
        this.invalidationSets = new Map(); // Maps molecule IDs to sets of interaction keys
        this.lastRandomSelections = new Map(); // Tiene traccia delle ultime selezioni casuali
        this.randomSelectionRate = 0.05; // Probabilità di considerare interazioni non in cache
        this.hitCount = 0;
        this.missCount = 0;
    }
  
    getInteractionKey(mol1Id, mol2Id) {
        // Create consistent key regardless of molecule order
        return mol1Id < mol2Id 
            ? `${mol1Id}:${mol2Id}` 
            : `${mol2Id}:${mol1Id}`;
    }
  
    hasInteraction(key) {
        const hasKey = this.cache.has(key);
        if (hasKey) this.hitCount++;
        else this.missCount++;
        return hasKey;
    }
  
    getInteraction(key) {
        const entry = this.cache.get(key);
        if (entry) {
            entry.age = 0; // Reset age when accessed
            entry.accessCount = (entry.accessCount || 0) + 1; // Incrementa contatore accessi
            return entry.data;
        }
        return null;
    }
  
    storeInteraction(key, data) {
        // Store interaction data with age tracking
        this.cache.set(key, {
            data: data,
            age: 0,
            accessCount: 1, // Inizializza contatore accessi
            createdAt: performance.now()
        });
        
        // Track which molecules are involved in this interaction
        const [mol1Id, mol2Id] = key.split(':');
        
        // Add to invalidation map for mol1
        if (!this.invalidationSets.has(mol1Id)) {
            this.invalidationSets.set(mol1Id, new Set());
        }
        this.invalidationSets.get(mol1Id).add(key);
        
        // Add to invalidation map for mol2
        if (!this.invalidationSets.has(mol2Id)) {
            this.invalidationSets.set(mol2Id, new Set());
        }
        this.invalidationSets.get(mol2Id).add(key);
    }
  
    invalidateForMolecule(molId) {
        // Invalidate all cached interactions involving this molecule
        const keys = this.invalidationSets.get(molId);
        if (keys) {
            for (const key of keys) {
                this.cache.delete(key);
            }
            this.invalidationSets.delete(molId);
        }
    }
  
    invalidateForMolecules(molecules) {
        // Invalidate cache for multiple molecules
        for (const mol of molecules) {
            this.invalidateForMolecule(mol.id);
        }
    }
  
    ageEntries() {
        // Increase age of all cache entries
        const now = performance.now();
        for (const [key, entry] of this.cache.entries()) {
            entry.age++;
            
            // Remove entries that exceed lifetime
            if (entry.age > this.lifetime) {
                this.cache.delete(key);
                
                // Clean up invalidation sets
                const [mol1Id, mol2Id] = key.split(':');
                const set1 = this.invalidationSets.get(mol1Id);
                if (set1) set1.delete(key);
                
                const set2 = this.invalidationSets.get(mol2Id);
                if (set2) set2.delete(key);
            }
        }
    }
  
    clear() {
        this.cache.clear();
        this.invalidationSets.clear();
        this.lastRandomSelections.clear();
        this.hitCount = 0;
        this.missCount = 0;
    }
  
    getSerializableCache() {
        // Create a serializable version of the cache for sending to sub-workers
        const result = {};
        for (const [key, entry] of this.cache.entries()) {
            if (entry.age <= this.lifetime / 2) { // Only share fresh entries
                result[key] = entry.data;
            }
        }
        return result;
    }
  
    importCache(serializedCache) {
        // Import cache from serialized format (from main worker)
        if (!serializedCache) return;
        
        for (const [key, data] of Object.entries(serializedCache)) {
            this.storeInteraction(key, data);
        }
    }
    
    // Nuovo metodo: determina se dovremmo considerare un'interazione casuale
    shouldConsiderRandomInteraction(mol1Id, mol2Id) {
        const key = this.getInteractionKey(mol1Id, mol2Id);
        
        // Se l'interazione è già nella cache, non selezionarla casualmente
        if (this.cache.has(key)) return false;
        
        // Verifica quando è stata l'ultima selezione casuale per queste molecole
        const now = performance.now();
        const lastTime1 = this.lastRandomSelections.get(mol1Id) || 0;
        const lastTime2 = this.lastRandomSelections.get(mol2Id) || 0;
        
        // Evita di selezionare troppo frequentemente le stesse molecole
        if (now - lastTime1 < 1000 || now - lastTime2 < 1000) return false;
        
        // Selezione casuale con probabilità configurata
        if (Math.random() < this.randomSelectionRate) {
            // Registra il timestamp della selezione
            this.lastRandomSelections.set(mol1Id, now);
            this.lastRandomSelections.set(mol2Id, now);
            return true;
        }
        
        return false;
    }
    
    // Ottieni statistiche del cache
    getStats() {
        const totalRequests = this.hitCount + this.missCount;
        const hitRate = totalRequests > 0 ? this.hitCount / totalRequests : 0;
        
        return {
            cacheSize: this.cache.size,
            invalidationSetsSize: this.invalidationSets.size,
            hitRate: hitRate,
            randomSelectionsCount: this.lastRandomSelections.size
        };
    }
}
  
export { MoleculeCentricCache, InteractionCache };