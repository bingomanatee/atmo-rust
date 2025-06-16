use crate::asthenosphere::AsthenosphereCell;
use crate::rock_store::RockStore;
use h3o::CellIndex;
use rand::Rng;
use crate::asth_constants::{AVG_VOLUME_TO_ADD, COOLING_RATE, K_PER_VOLUME};

pub struct ProcessResult {
    pub volume_added: f64,
    pub volume_removed: f64,
    pub new_volume: f64,
    pub energy_k: f64,
    pub cell: AsthenosphereCell
}

pub fn process_cell(
    store: &RockStore,
    l2_cell: CellIndex,
    to_mio_years: u32,
    total_mio_years: u32,
) -> Result<ProcessResult, String> {
    match store.get_asth(l2_cell, to_mio_years - 1) {
        // to tune the system we are ONLY considering cooling rate;
        Ok(Some(old_cell)) => {
            let new_energy = old_cell.energy_k as f64 * *COOLING_RATE;

            let new_cell = AsthenosphereCell {
                step: to_mio_years,
                volume: old_cell.volume,
                energy_k: new_energy,
                ..old_cell.clone()
            };

           /* store
                .put_asth(&new_cell)
                .map_err(|e| format!("Failed to save cell: {:?}", e))?;
*/
            Ok(ProcessResult {
                volume_added: 0.0,
                volume_removed: 0.0,
                new_volume: new_cell.volume,
                energy_k: new_energy,
                cell: new_cell
            })
        }
        Ok(None) => Err(format!("Cell not found: {}, {}", l2_cell, to_mio_years - 1)),
        Err(e) => Err(format!("Error fetching cell: {:?}", e)),
    }
}
pub fn process_cell_old(
    store: &RockStore,
    l2_cell: CellIndex,
    to_mio_years: u32,
    total_mio_years: u32,
) -> Result<ProcessResult, String> {
    match store.get_asth(l2_cell, to_mio_years - 1) {
        Ok(Some(old_cell)) => {
            //    let progress = to_mio_years as f64 / total_mio_years as f64;
            let average = AVG_VOLUME_TO_ADD; // we were doing a gradient input - now it is constant
            // - (AE_SPAN * progress);
            let mut rng = rand::rng();
            let rand_multiplier: f64 = rng.random_range(0.5..=1.5);
            let volume_to_add = average * rand_multiplier;
            let sunk_volume = old_cell.sunk_volume();
            let new_volume = old_cell.volume + volume_to_add - sunk_volume;

            let old_energy = *COOLING_RATE * old_cell.energy_k as f64;
            let added_energy = K_PER_VOLUME * volume_to_add;
            let mut total_energy = old_energy + added_energy;
            if sunk_volume > 0.0 {
                total_energy -= old_energy * sunk_volume / old_cell.volume;
            }

            let new_cell = AsthenosphereCell {
                step: to_mio_years,
                volume: new_volume,
                energy_k: total_energy,
                ..old_cell.clone()
            };
           /* store
                .put_asth(&new_cell)
                .map_err(|e| format!("Failed to save cell: {:?}", e))?;
*/
            Ok(ProcessResult {
                volume_added: volume_to_add,
                volume_removed: sunk_volume,
                new_volume,
                energy_k: total_energy,
                cell: new_cell
            })
        }
        Ok(None) => Err(format!("Cell not found: {}, {}", l2_cell, to_mio_years - 1)),
        Err(e) => Err(format!("Error fetching cell: {:?}", e)),
    }
}

