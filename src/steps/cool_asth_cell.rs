use crate::constants::{
    AVG_VOLUME_TO_ADD, CELL_JOULES_EQUILIBRIUM, CELL_JOULES_START, COOLING_RATE, JOULES_PER_KM3,
    LAYER_COUNT, MAX_SUNK_TEMP, STANDARD_STEPS,
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
    let cool_rate = (CELL_JOULES_EQUILIBRIUM / CELL_JOULES_START).powf(1.1 / (total_mio_years as f64));

    // Print cool rate only once, even when called in parallel
    COOL_RATE_PRINTED.call_once(|| {
        println!("cool rate: {}", cool_rate);
    });
    match store.get_asth(l2_cell, step - 1) {
        // to tune the system we are ONLY considering cooling rate;
        Ok(Some(old_cell)) => {
            let mut new_cell = old_cell.clone();
            new_cell.step = step;
            
            // Apply cooling to all layers
            for layer_idx in 0..LAYER_COUNT {
                new_cell.energy_layers[layer_idx] *= cool_rate;
            }
            
            // Use surface layer for result reporting
            let surface_layer = LAYER_COUNT - 1;
            let new_energy = new_cell.energy_layers[surface_layer];
            let new_volume = new_cell.volume_layers[surface_layer];
            
            Ok(ProcessResult {
                volume_added: 0.0,
                volume_removed: 0.0,
                new_volume,
                energy_k: new_energy,
                cell: new_cell,
            })
        }
        Ok(None) => Err(format!("Cell not found: {}, {}", l2_cell, step - 1)),
        Err(e) => Err(format!("Error fetching cell: {:?}", e)),
    }
}