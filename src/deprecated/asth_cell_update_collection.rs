use crate::asthenosphere::AsthenosphereCell;
use crate::h3_utils::H3Utils;
use crate::rock_store::RockStore;
use h3o::{CellIndex, Resolution};
use rayon::iter::IntoParallelRefIterator;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use crate::deprecated::asth_cell_update::AsthCellUpdate;

pub struct AsthCellUpdateCollection<'a> {
    store: &'a RockStore,
    asth_cell_updates: HashMap<CellIndex, AsthCellUpdate>,
    resolution: Resolution,
    step: u32,
}

impl<'a> AsthCellUpdateCollection<'a> {
    pub fn new(store: &'a RockStore, resolution: Resolution, step: u32) -> Self {
        AsthCellUpdateCollection {
            store,
            asth_cell_updates: HashMap::new(),
            resolution,
            step,
        }
    }

    pub fn insert_update(&mut self, cell_update: AsthCellUpdate) -> bool {
        let cell_index = cell_update.cell.id;

        if self.asth_cell_updates.contains_key(&cell_index) {
            return false;
        }

        if cell_update.cell.step != self.step {
            eprintln!(
                "Warning: Cell step {} doesn't match collection step {} for cell {}",
                cell_update.cell.step, self.step, cell_index
            );
            return false;
        }

        self.asth_cell_updates.insert(cell_index, cell_update);
        true
    }

    pub fn insert_asth_cell(&mut self, cell: AsthenosphereCell) -> bool {
        self.insert_update(AsthCellUpdate::new(cell))
    }

    pub fn insert_id(&mut self, cell_index: CellIndex) -> bool {
        if self.asth_cell_updates.contains_key(&cell_index) {
            return false;
        }

        if let Ok(Some(cell)) = self.store.get_asth(cell_index, self.step) {
            self.insert_update(AsthCellUpdate::new(cell))
        } else {
            false
        }
    }

    fn revise_update(&mut self, cell_index: CellIndex, cell_update: AsthCellUpdate) {
        if cell_update.cell.step != self.step {
            eprintln!(
                "Warning: Cell step {} doesn't match collection step {} for cell {}",
                cell_update.cell.step, self.step, cell_index
            );
            return;
        }

        self.asth_cell_updates.insert(cell_index, cell_update);
    }

    pub fn load_from_store(&mut self, cell_index: CellIndex) -> bool {
        if self.asth_cell_updates.contains_key(&cell_index) {
            return true;
        }

        self.insert_id(cell_index)
    }

    pub fn reload_from_store(&mut self, cell_index: CellIndex) -> bool {
        if let Ok(Some(cell)) = self.store.get_asth(cell_index, self.step) {
            self.revise_update(cell_index, AsthCellUpdate::new(cell));
            return true;
        }

        false
    }

    pub fn load_all(&mut self) {
        let res = self.resolution;
        let mut store = self.store;
        let ids: Vec<CellIndex> = H3Utils::iter_cells_with_base(res)
            .map(|(cell, _)| cell)
            .collect();

        let loaded_map: HashMap<CellIndex, AsthCellUpdate> = ids
            .par_iter()
            .filter_map(|&id| match store.get_asth(id, self.step) {
                Ok(None) => None,
                Err(_) => None,
                Ok(Some(cell)) => Some((id, AsthCellUpdate::new(cell))),
            })
            .collect();
        self.asth_cell_updates.extend(loaded_map);
    }

    pub fn save(&mut self) -> Result<(), String> {
        let cells_to_save: Vec<AsthenosphereCell> = self
            .asth_cell_updates
            .values()
            .filter_map(|cell_update| {
                if !cell_update.is_dirty() {
                    return None;
                }
                println!("cell has changed: {:?}", cell_update);
                let new_cell = cell_update.update();
                Some(new_cell)
            })
            .collect();

        
        if !cells_to_save.is_empty() {
            self.store.put_asth_batch(&cells_to_save);
            for cell in cells_to_save {
                self.reload_from_store(cell.id);
            }
        }
        Ok(())
    }

    pub fn get(&self, cell_index: &CellIndex) -> Option<&AsthCellUpdate> {
        self.asth_cell_updates.get(cell_index).clone()
    }

    pub fn contains(&self, cell_index: &CellIndex) -> bool {
        self.asth_cell_updates.contains_key(cell_index)
    }

    pub fn len(&self) -> usize {
        self.asth_cell_updates.len()
    }

    pub fn is_empty(&self) -> bool {
        self.asth_cell_updates.is_empty()
    }

    pub fn cell_indices(&self) -> impl Iterator<Item = &CellIndex> {
        self.asth_cell_updates.keys()
    }

    pub fn cells(&self) -> impl Iterator<Item = &AsthenosphereCell> {
        self.asth_cell_updates.values().map(|update| &update.cell)
    }

    pub fn clear(&mut self) {
        self.asth_cell_updates.clear();
    }
}

