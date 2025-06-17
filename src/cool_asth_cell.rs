use crate::asth_constants::{
    AVG_VOLUME_TO_ADD, CELL_ENERGY_EQUILIBRIUM, CELL_ENERGY_START, COOLING_RATE, K_PER_VOLUME,
    MAX_SUNK_TEMP, STANDARD_STEPS,
};
use crate::asthenosphere::AsthenosphereCell;
use crate::rock_store::RockStore;
use h3o::CellIndex;
use rand::Rng;
use std::sync::Once;

pub struct ProcessResult {
    pub volume_added: f64,
    pub volume_removed: f64,
    pub new_volume: f64,
    pub energy_k: f64,
    pub cell: AsthenosphereCell,
}

static COOL_RATE_PRINTED: Once = Once::new();

pub fn cool_asth_cell(
    store: &RockStore,
    l2_cell: CellIndex,
    step: u32,
    total_mio_years: u32,
) -> Result<ProcessResult, String> {
    let cool_rate = (CELL_ENERGY_EQUILIBRIUM / CELL_ENERGY_START).powf(1.0 / (total_mio_years as f64));

    // Print cool rate only once, even when called in parallel
    COOL_RATE_PRINTED.call_once(|| {
        println!("cool rate: {}", cool_rate);
    });
    match store.get_asth(l2_cell, step - 1) {
        // to tune the system we are ONLY considering cooling rate;
        Ok(Some(old_cell)) => {
            let new_energy = old_cell.energy_k as f64 * cool_rate;

            let new_cell = AsthenosphereCell {
                step,
                volume: old_cell.volume,
                energy_k: new_energy,
                ..old_cell.clone()
            };
            
            Ok(ProcessResult {
                volume_added: 0.0,
                volume_removed: 0.0,
                new_volume: new_cell.volume,
                energy_k: new_energy,
                cell: new_cell,
            })
        }
        Ok(None) => Err(format!("Cell not found: {}, {}", l2_cell, step - 1)),
        Err(e) => Err(format!("Error fetching cell: {:?}", e)),
    }
}