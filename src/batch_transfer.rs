use h3o::CellIndex;
use std::collections::HashMap;
use std::sync::{Mutex, OnceLock};
use crate::constants::LAYER_COUNT;

/// Unique identifier for a cell-layer combination using string pattern
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct CellLayerKey {
    pub cell_id: CellIndex,
    pub layer_index: usize,
}

impl CellLayerKey {
    pub fn new(cell_id: CellIndex, layer_index: usize) -> Self {
        Self { cell_id, layer_index }
    }
    
    /// Create string representation for debugging
    pub fn to_string(&self) -> String {
        format!("{}_{}", self.cell_id, self.layer_index)
    }
}

/// Two-level hash structure for tracking transfers: from_key -> to_key -> volume
type TransferMap = HashMap<CellLayerKey, HashMap<CellLayerKey, f64>>;

/// Singleton tracker for aggregating all volume and energy transfers
pub struct BatchTransferTracker {
    transfers: Mutex<TransferMap>,
}

impl BatchTransferTracker {
    /// Get the global singleton instance
    pub fn instance() -> &'static BatchTransferTracker {
        static INSTANCE: OnceLock<BatchTransferTracker> = OnceLock::new();
        INSTANCE.get_or_init(|| BatchTransferTracker {
            transfers: Mutex::new(HashMap::new()),
        })
    }
    
    /// Register a transfer from one cell-layer to another
    pub fn register_transfer(
        &self,
        from_key: CellLayerKey,
        to_key: CellLayerKey,
        volume: f64,
    ) {
        let mut transfers = self.transfers.lock().unwrap();
        
        // Get or create the from_key entry
        let from_map = transfers.entry(from_key).or_insert_with(HashMap::new);
        
        // Add to existing transfer or create new one
        *from_map.entry(to_key).or_insert(0.0) += volume;
    }
    
    /// Get all aggregated transfers and clear the tracker
    pub fn consume_transfers(&self) -> TransferMap {
        let mut transfers = self.transfers.lock().unwrap();
        std::mem::take(&mut *transfers)
    }
    
    /// Clear all pending transfers without returning them
    pub fn clear(&self) {
        let mut transfers = self.transfers.lock().unwrap();
        transfers.clear();
    }
    
    /// Get count of pending transfer operations for debugging
    pub fn pending_count(&self) -> usize {
        let transfers = self.transfers.lock().unwrap();
        transfers.values().map(|to_map| to_map.len()).sum()
    }
}

/// Performance timer for measuring step timing
#[derive(Debug, Default)]
pub struct StepTimer {
    pub leveling_time_ms: f64,
    pub mixing_time_ms: f64,
    pub ejection_time_ms: f64,
    pub convection_time_ms: f64,
    pub total_step_time_ms: f64,
    pub transfer_commit_time_ms: f64,
}

impl StepTimer {
    pub fn new() -> Self {
        Self::default()
    }
    
    /// Print timing statistics
    pub fn print_stats(&self, step: u32) {
        println!(
            "⏱️ Step {} timing: total={:.2}ms (leveling={:.2}ms, mixing={:.2}ms, convection={:.2}ms, transfers={:.2}ms)",
            step,
            self.total_step_time_ms,
            self.leveling_time_ms,
            self.mixing_time_ms,
            self.convection_time_ms,
            self.transfer_commit_time_ms
        );
    }
    
    /// Get total step time in seconds
    pub fn total_step_time_s(&self) -> f64 {
        self.total_step_time_ms / 1000.0
    }
}

/// Transfer record for accounting system
#[derive(Debug, Clone)]
pub struct TransferRecord {
    pub from_key: CellLayerKey,
    pub to_key: CellLayerKey,
    pub volume: f64,
    pub step: u32,
}

impl TransferRecord {
    pub fn new(from_key: CellLayerKey, to_key: CellLayerKey, volume: f64, step: u32) -> Self {
        Self { from_key, to_key, volume, step }
    }
}

/// Singleton accounting system to track all volume transfers
pub struct TransferAccountingSystem {
    records: Mutex<Vec<TransferRecord>>,
    current_step: Mutex<u32>,
}

impl TransferAccountingSystem {
    /// Get the global singleton instance
    pub fn instance() -> &'static TransferAccountingSystem {
        static INSTANCE: OnceLock<TransferAccountingSystem> = OnceLock::new();
        INSTANCE.get_or_init(|| TransferAccountingSystem {
            records: Mutex::new(Vec::new()),
            current_step: Mutex::new(0),
        })
    }

    /// Set the current simulation step
    pub fn set_step(&self, step: u32) {
        let mut current_step = self.current_step.lock().unwrap();
        *current_step = step;
    }

    /// Record a volume transfer
    pub fn record_transfer(&self, from_key: CellLayerKey, to_key: CellLayerKey, volume: f64) {
        if volume <= 0.0 {
            return; // Don't record zero or negative transfers
        }

        let current_step = *self.current_step.lock().unwrap();
        let record = TransferRecord::new(from_key, to_key, volume, current_step);
        
        let mut records = self.records.lock().unwrap();
        records.push(record);
    }

    /// Get all transfer records
    pub fn get_all_records(&self) -> Vec<TransferRecord> {
        let records = self.records.lock().unwrap();
        records.clone()
    }

    /// Get transfer records for a specific step
    pub fn get_records_for_step(&self, step: u32) -> Vec<TransferRecord> {
        let records = self.records.lock().unwrap();
        records.iter()
            .filter(|r| r.step == step)
            .cloned()
            .collect()
    }

    /// Get total volume transferred in a specific step
    pub fn get_total_volume_for_step(&self, step: u32) -> f64 {
        let records = self.records.lock().unwrap();
        records.iter()
            .filter(|r| r.step == step)
            .map(|r| r.volume)
            .sum()
    }

    /// Clear all records
    pub fn clear(&self) {
        let mut records = self.records.lock().unwrap();
        records.clear();
    }

    /// Get record count
    pub fn record_count(&self) -> usize {
        let records = self.records.lock().unwrap();
        records.len()
    }

    /// Flatten all current step transfers into a consolidated from->to volume mapping
    pub fn flatten_transfers_for_step(&self, step: u32) -> HashMap<CellLayerKey, HashMap<CellLayerKey, f64>> {
        let records = self.records.lock().unwrap();
        let mut flattened: HashMap<CellLayerKey, HashMap<CellLayerKey, f64>> = HashMap::new();
        
        for record in records.iter().filter(|r| r.step == step) {
            let from_map = flattened.entry(record.from_key.clone()).or_insert_with(HashMap::new);
            let existing_volume = from_map.entry(record.to_key.clone()).or_insert(0.0);
            *existing_volume += record.volume;
        }
        
        flattened
    }

    /// Clear records for a specific step (after processing)
    pub fn clear_step(&self, step: u32) {
        let mut records = self.records.lock().unwrap();
        records.retain(|r| r.step != step);
    }
}

/// Utility function to measure execution time of a closure
pub fn time_operation<F, R>(operation: F) -> (R, f64)
where
    F: FnOnce() -> R,
{
    let start = std::time::Instant::now();
    let result = operation();
    let duration_ms = start.elapsed().as_secs_f64() * 1000.0;
    (result, duration_ms)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_cell_layer_key_creation() {
        let cell_id = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let key = CellLayerKey::new(cell_id, 1);
        
        assert_eq!(key.cell_id, cell_id);
        assert_eq!(key.layer_index, 1);
        assert_eq!(key.to_string(), format!("{}_{}", cell_id, 1));
    }
    
    #[test]
    fn test_volume_aggregation() {
        let tracker = BatchTransferTracker::instance();
        tracker.clear(); // Start clean
        
        let cell_a = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let cell_b = CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap();
        let key_a = CellLayerKey::new(cell_a, 0);
        let key_b = CellLayerKey::new(cell_b, 1);
        
        // Register multiple volumes to same destination
        tracker.register_transfer(key_a.clone(), key_b.clone(), 100.0);
        tracker.register_transfer(key_a.clone(), key_b.clone(), 50.0);
        
        let transfers = tracker.consume_transfers();
        let from_map = transfers.get(&key_a).unwrap();
        let aggregated_volume = from_map.get(&key_b).unwrap();
        
        assert_eq!(*aggregated_volume, 150.0);
    }
    
    #[test]
    fn test_batch_transfer_tracker() {
        let tracker = BatchTransferTracker::instance();
        tracker.clear(); // Start clean
        
        let cell_a = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let cell_b = CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap();
        let key_a = CellLayerKey::new(cell_a, 0);
        let key_b = CellLayerKey::new(cell_b, 1);
        
        // Register two transfers
        tracker.register_transfer(key_a.clone(), key_b.clone(), 100.0);
        tracker.register_transfer(key_a.clone(), key_b.clone(), 50.0);
        
        assert_eq!(tracker.pending_count(), 1); // Should aggregate into single transfer
        
        let transfers = tracker.consume_transfers();
        assert_eq!(transfers.len(), 1);
        
        let from_map = transfers.get(&key_a).unwrap();
        let aggregated_volume = from_map.get(&key_b).unwrap();
        
        assert_eq!(*aggregated_volume, 150.0);
        
        assert_eq!(tracker.pending_count(), 0); // Should be empty after consume
    }
}