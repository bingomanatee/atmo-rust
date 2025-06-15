use crate::asthenosphere::AsthenosphereCell;
use crate::rock_store::RockStore;
use crate::asth_sim::{AVERAGE_ENERGY_AT_START, AE_SPAN, K_PER_VOLUME, COOLING_RATE};
use h3o::CellIndex;
use rand::Rng;

pub struct ProcessResult {
    pub volume_added: f32,
    pub volume_removed: f32,
    pub new_volume: f32,
    pub total_energy: f32,
}

pub fn process_cell(
    store: &RockStore,
    l2_cell: CellIndex,
    to_mio_years: u32,
    total_mio_years: u32,
) -> Result<ProcessResult, String> {
    match store.get_asth(l2_cell, to_mio_years - 1) {
        Ok(Some(old_cell)) => {
            let progress = to_mio_years as f32 / total_mio_years as f32;
            let average = AVERAGE_ENERGY_AT_START - (AE_SPAN * progress);
            let mut rng = rand::thread_rng();
            let rand_multiplier: f32 = rng.gen_range(0.5..=1.5);
            let volume_to_add = average * rand_multiplier;
            let sunk_volume = old_cell.sunk_volume();
            let new_volume = old_cell.volume + volume_to_add - sunk_volume;

            let old_energy = COOLING_RATE * old_cell.energy_k as f32;
            let added_energy = K_PER_VOLUME * volume_to_add;
            let mut total_energy = old_energy + added_energy;
            if sunk_volume > 0.0 {
                total_energy *= (old_cell.volume - sunk_volume) / old_cell.volume;
            }

            let new_cell = AsthenosphereCell {
                step: to_mio_years,
                volume: new_volume,
                energy_k: total_energy as u32,
                ..old_cell.clone()
            };
            store.put_asth(&new_cell).map_err(|e| format!("Failed to save cell: {:?}", e))?;

            Ok(ProcessResult {
                volume_added: volume_to_add,
                volume_removed: sunk_volume,
                new_volume,
                total_energy,
            })
        }
        Ok(None) => Err(format!("Cell not found: {}, {}", l2_cell, to_mio_years - 1)),
        Err(e) => Err(format!("Error fetching cell: {:?}", e)),
    }
}
