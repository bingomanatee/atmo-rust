use crate::asth_sim::VolumeEnergyTransfer;
use crate::asthenosphere::AsthenosphereCell;
use crate::planet::Planet;
use crate::plate::Plate;
use crate::sim::Sim;
use bincode;
use h3o::CellIndex;
use rayon::prelude::*;
use rocksdb::Error as RocksDbError;
use rocksdb::{DB, Options, WriteBatch};
use std::collections::{HashMap, HashSet};
use uuid::Uuid;

pub struct RockStore {
    db: DB,
}

impl RockStore {
    pub fn open(path: &str) -> Result<Self, rocksdb::Error> {
        let mut opts = Options::default();
        opts.create_if_missing(true);
        let db = DB::open(&opts, path)?;
        Ok(RockStore { db })
    }

    // Helper to generate keys with prefixes
    fn key_asth(id: &CellIndex, step: u32) -> Vec<u8> {
        let mut key = b"asth:".to_vec();
        key.extend_from_slice(step.to_string().as_bytes());
        key.extend_from_slice(b":");
        key.extend_from_slice(id.to_string().as_bytes());
        key
    }
    // Helper to generate keys with prefixes
    fn key_sim(id: &Uuid) -> Vec<u8> {
        let mut key = b"sim:".to_vec();
        key.extend(id.as_bytes());
        key
    }

    fn key_planet(id: &Uuid) -> Vec<u8> {
        let mut key = b"planet:".to_vec();
        key.extend(id.as_bytes());
        key
    }
    fn key_transfer(id: &Uuid, step: u32) -> Vec<u8> {
        let mut key = b"transfer:".to_vec();
        key.extend_from_slice(step.to_string().as_bytes());
        key.extend_from_slice(b":");
        key.extend(id.as_bytes());
        key
    }

    fn key_plate(id: &Uuid) -> Vec<u8> {
        let mut key = b"plate:".to_vec();
        key.extend(id.as_bytes());
        key
    }

    // Save Sim
    pub fn put_sim(&self, sim: &Sim) -> Result<(), rocksdb::Error> {
        let key = Self::key_sim(&sim.id);
        let value = bincode::serialize(sim).unwrap();
        self.db.put(key, value)
    }

    // Get AsthenosphereCell
    pub fn get_asth(
        &self,
        cell: CellIndex,
        step: u32,
    ) -> Result<Option<AsthenosphereCell>, rocksdb::Error> {
        let key = Self::key_asth(&cell, step);
        match self.db.get(key)? {
            Some(value) => {
                let a: AsthenosphereCell = bincode::deserialize(&value).unwrap();
                Ok(Some(a))
            }
            None => Ok(None),
        }
    }
    // Save AsthenosphereCell
    pub fn put_asth(&self, a: &AsthenosphereCell) -> Result<(), rocksdb::Error> {
        let key = Self::key_asth(&a.cell, a.step);
        let value = bincode::serialize(a).unwrap();
        self.db.put(key, value)
    }

    // Get Sim
    pub fn get_sim(&self, id: &Uuid) -> Result<Option<Sim>, rocksdb::Error> {
        let key = Self::key_sim(id);
        match self.db.get(key)? {
            Some(value) => {
                let sim: Sim = bincode::deserialize(&value).unwrap();
                Ok(Some(sim))
            }
            None => Ok(None),
        }
    }

    // Save Planet
    pub fn put_planet(&self, planet: &Planet) -> Result<(), rocksdb::Error> {
        let key = Self::key_planet(&planet.id);
        let value = bincode::serialize(planet).unwrap();
        self.db.put(key, value)
    }

    // Get Planet
    pub fn get_planet(&self, id: &Uuid) -> Result<Option<Planet>, rocksdb::Error> {
        let key = Self::key_planet(id);
        match self.db.get(key)? {
            Some(value) => {
                let planet: Planet = bincode::deserialize(&value).unwrap();
                Ok(Some(planet))
            }
            None => Ok(None),
        }
    }

    // Save Plate
    pub fn put_plate(&self, plate: &Plate) -> Result<(), rocksdb::Error> {
        let key = Self::key_plate(&plate.id);
        let value = bincode::serialize(plate).unwrap();
        self.db.put(key, value)
    }

    // Get Plate
    pub fn get_plate(&self, id: &Uuid) -> Result<Option<Plate>, rocksdb::Error> {
        let key = Self::key_plate(id);
        match self.db.get(key)? {
            Some(value) => {
                let plate: Plate = bincode::deserialize(&value).unwrap();
                Ok(Some(plate))
            }
            None => Ok(None),
        }
    }

    pub fn each_plate<F, E>(&self, mut callback: F) -> Result<(), E>
    where
        F: FnMut(Plate, Box<[u8]>) -> Result<(), E>,
        E: From<rocksdb::Error> + From<bincode::Error>,
    {
        let prefix = b"plate:";
        let iter = self.db.prefix_iterator(prefix);

        for item in iter {
            let (key, value) = item.map_err(E::from)?;
            if !key.starts_with(prefix) {
                break;
            }
            let plate: Plate = bincode::deserialize(&value).map_err(E::from)?;
            callback(plate, key)?;
        }
        Ok(())
    }

    pub fn put_asth_batch(&self, cells: &[AsthenosphereCell]) {
        let mut batch = WriteBatch::default();
        for cell in cells {
            let key = RockStore::key_asth(&cell.cell, cell.step);
            let value = bincode::serialize(cell).expect("serialize failed");
            batch.put(key, value);
        }
        self.db.write(batch);
    }

    pub fn batch_write_transfers(&self, transfers: Vec<VolumeEnergyTransfer>) {
        let mut batch = WriteBatch::default();
        for transfer in transfers {
            let key = RockStore::key_transfer(&transfer.id, transfer.step);
            let value = bincode::serialize(&transfer).expect("serialize failure");
            batch.put(key, value);
        }
        self.db.write(batch);
    }

    /// Retrieve all transfers for a specific step using prefix iteration
    pub fn each_transfer_for_step<F>(&self, step: u32, mut callback: F) -> Result<(), String>
    where
        F: FnMut(VolumeEnergyTransfer, Vec<u8>) -> Result<(), String>,
    {
        // Create prefix that matches the key_transfer format: "transfer:{step}:"
        let mut prefix = b"transfer:".to_vec();
        prefix.extend_from_slice(step.to_string().as_bytes());
        prefix.extend_from_slice(b":");

        let iter = self.db.prefix_iterator(&prefix);

        for item in iter {
            let (key, value) = item.map_err(|e| format!("RocksDB error: {}", e))?;
            if !key.starts_with(&prefix) {
                break;
            }
            let transfer: VolumeEnergyTransfer = bincode::deserialize(&value)
                .map_err(|e| format!("Deserialization error: {}", e))?;
            callback(transfer, key.to_vec())?;
        }
        Ok(())
    }

    /// Process all volume transfers for a given step using parallel processing
    pub fn transfer_volume(&self, step: u32) {
        // Collect all transfers for this step
        let mut transfers = Vec::new();
        let mut transfer_keys = Vec::new();

        let _ = self.each_transfer_for_step(step, |transfer, key| {
            transfers.push(transfer);
            transfer_keys.push(key);
            Ok(())
        });

        if transfers.is_empty() {
            return; // No transfers to process
        }

        // Group transfers by both source and target cells to aggregate changes
        let mut changes: HashMap<CellIndex, (f64, f64)> = HashMap::new(); // (volume_delta, energy_delta)

        for transfer in &transfers {
            let source_entry = changes.entry(transfer.from_cell).or_insert((0.0, 0.0));
            source_entry.0 -= transfer.volume;
            source_entry.1 -= transfer.energy;

            let target_entry = changes.entry(transfer.to_cell).or_insert((0.0, 0.0));
            target_entry.0 += transfer.volume;
            target_entry.1 += transfer.energy;
        }

        // Process all affected cells in parallel
        let updated_cells: Vec<AsthenosphereCell> = changes
            .par_iter()
            .filter_map(|(&cell_index, change)| {
                if let Ok(Some(mut cell)) = self.get_asth(cell_index, step) {
                    let (vol, energy) = change;
                    cell.volume += vol; 
                    cell.energy_k += energy; 
                    Some(cell)
                } else {
                    None
                }
            })
            .collect();

        // Batch write all updated cells
        if !updated_cells.is_empty() {
            self.put_asth_batch(&updated_cells);
        }

        // Clean up processed transfers
        let mut batch = WriteBatch::default();
        for key in transfer_keys {
            batch.delete(key);
        }
        self.db.write(batch);
    }

    /// Store a single VolumeEnergyTransfer record in RocksDB.
    pub fn put_transfer(&self, transfer: &VolumeEnergyTransfer) -> Result<(), RocksDbError> {
        let serialized = bincode::serialize(transfer).expect("serialize failed");
        let key = RockStore::key_transfer(&transfer.id, transfer.step);
        self.db.put(key, serialized)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::asthenosphere::AsthenosphereCell;
    use h3o::{CellIndex, Resolution};
    use tempfile::tempdir;
    use uuid::Uuid;

    #[test]
    fn test_transfer_volume() {
        // Create temporary RocksDB directory
        let temp_dir = tempdir().expect("failed to create temp dir");
        let db_path = temp_dir.path().to_str().unwrap().to_string();

        // Open RockStore
        let store = RockStore::open(&db_path).expect("failed to open RockStore");

        let step = 1u32;

        // Create test cells
        let cell1 = CellIndex::try_from(0x85283473fffffff).unwrap(); // Example H3 cell
        let cell2 = CellIndex::try_from(0x85283447fffffff).unwrap(); // Another H3 cell

        let mut asth_cell1 = AsthenosphereCell {
            cell: cell1,
            volume: 1000.0,
            energy_k: 2000.0,
            step,
            neighbors: vec![cell2],
        };

        let mut asth_cell2 = AsthenosphereCell {
            cell: cell2,
            volume: 500.0,
            energy_k: 1000.0,
            step,
            neighbors: vec![cell1],
        };

        // Store initial cells
        store.put_asth(&asth_cell1).expect("failed to store cell1");
        store.put_asth(&asth_cell2).expect("failed to store cell2");

        // Create a transfer from cell1 to cell2
        let transfer = VolumeEnergyTransfer {
            id: Uuid::new_v4(),
            step,
            from_cell: cell1,
            to_cell: cell2,
            volume: 100.0,
            energy: 200.0,
        };

        // Store the transfer
        store
            .put_transfer(&transfer)
            .expect("failed to store transfer");

        // Execute transfer_volume
        store.transfer_volume(step);

        // Verify the results
        let updated_cell1 = store
            .get_asth(cell1, step)
            .expect("failed to get cell1")
            .expect("cell1 not found");
        let updated_cell2 = store
            .get_asth(cell2, step)
            .expect("failed to get cell2")
            .expect("cell2 not found");

        // Cell1 should have lost volume and energy
        assert_eq!(updated_cell1.volume, 900.0); // 1000.0 - 100.0
        assert_eq!(updated_cell1.energy_k, 1800.0); // 2000.0 - 200.0

        // Cell2 should have gained volume and energy
        assert_eq!(updated_cell2.volume, 600.0); // 500.0 + 100.0
        assert_eq!(updated_cell2.energy_k, 1200.0); // 1000.0 + 200.0

        // Verify that the transfer was deleted
        let mut transfer_count = 0;
        let _ = store.each_transfer_for_step(step, |_, _| {
            transfer_count += 1;
            Ok(())
        });
        assert_eq!(transfer_count, 0, "Transfer should have been deleted");
    }
}
