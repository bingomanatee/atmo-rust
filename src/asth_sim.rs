use crate::asth_process_cell::{ProcessResult, process_cell};
use crate::asthenosphere::{ASTH_RES, AsthenosphereCell, CellsForPlanetArgs};
use crate::h30_utils::H3Utils;
use crate::planet::Planet;
use crate::rock_store::RockStore;
use h3o::{CellIndex, Resolution};
use rand::Rng;
use rayon::prelude::*;
use std::collections::HashMap;

pub struct AsthSim<'a> {
    store: RockStore,
    planet: &'a Planet,
    mio_years: u32,
}

pub struct AsthSimArgs<'a> {
    pub db_path: String,
    pub planet: &'a Planet,
    pub mio_years: u32,
    pub resolution: Resolution,
}

pub struct VolumeStats {
    pub volume_added: f64,
    pub volume_removed: f64,
    pub total_volume: f64,
    pub mean_energy: f64,
    pub mean_volume: f64,
}

impl<'a> AsthSim<'a> {
    pub fn new(args: AsthSimArgs<'a>) -> Self {
        let AsthSimArgs {
            db_path,
            planet,
            mio_years,
            resolution,
        } = args;

        let store = RockStore::open(&db_path).expect("cannot make store");

        AsthenosphereCell::cells_for_planet(CellsForPlanetArgs {
            planet: planet.clone(),
            res: Some(ASTH_RES),
            energy_per_volume: 200.0,
            on_cell: |asth_cell: AsthenosphereCell| {
                match store.put_asth(&asth_cell) {
                    Ok(_) => {
                        // Successfully stored
                    }
                    Err(e) => {
                        eprintln!("Failed to store cell {:?}: {:?}", asth_cell.cell, e);
                        // Optionally handle error, e.g., panic!("DB write failed");
                    }
                }
                asth_cell
            },
        });

        AsthSim {
            store,
            planet: &planet,
            mio_years,
        }
    }

    pub fn run(&self) {
        // Print table header
        println!(
            "{:<6} | {:<13} | {:<14} | {:<13} | {:<11} | {:<11}",
            "Step", "Volume Added", "Volume Removed", "Total Volume", "Mean Energy", "Mean_volume"
        );
        println!("{}", "-".repeat(65));

        for current in 1..=self.mio_years {
            let stats = self.advance(current);
            if current < 10 || (current % 10) == 0 {
                println!(
                    "{:<6} | {:<13.2} | {:<14.2} | {:<13.2} | {:<11.2}, {:<11.2}",
                    current,
                    stats.volume_added,
                    stats.volume_removed,
                    stats.total_volume,
                    stats.mean_energy,
                    stats.mean_volume
                );
            }
        }
    }

    fn advance(&self, to_mio_years: u32) -> VolumeStats {
        // Collect all cells at ASTH_RES into a Vec<CellIndex>
        let cells: Vec<CellIndex> = {
            let mut v = Vec::new();
            H3Utils::iter_at(ASTH_RES, |cell| v.push(cell));
            v
        };

        // Process cells in parallel using Rayon
        let results: Result<Vec<ProcessResult>, String> = cells
            .par_iter()
            .map(|l2_cell| process_cell(&self.store, *l2_cell, to_mio_years, self.mio_years))
            .collect();

        let results = match results {
            Ok(res) => res,
            Err(e) => panic!("Error processing cell: {}", e),
        };

        let mut volume_added = 0.0f64;
        let mut volume_removed = 0.0f64;
        let mut total_volume = 0.0f64;
        let mut total_energy = 0.0f64;

        let mut updated_cells = Vec::with_capacity(results.len());
    
        for r in &results {
            volume_added += r.volume_added;
            volume_removed += r.volume_removed;
            total_volume += r.new_volume;
            total_energy += r.energy_k;
            
            updated_cells.push(r.cell.clone());
        }

        // Batch write all updated cells at once
        self.store.put_asth_batch(&updated_cells);

        let cell_count = results.len() as f64;

        VolumeStats {
            volume_added,
            volume_removed: volume_removed,
            total_volume: total_volume,
            mean_energy: if cell_count > 0.0 {
                total_energy / cell_count
            } else {
                0.0
            },
            mean_volume: total_volume as f64 / cell_count,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::asth_constants::STANDARD_STEPS;
    use crate::planet::{EARTH_RADIUS_KM, PlanetParams, RHO_EARTH};
    use tempfile::tempdir;
    use uuid::Uuid;

    // You need to provide a minimal Planet for testing.
    // Assuming Planet implements Clone and Default or you can create a dummy.
    // Replace this with your actual Planet creation logic.
    fn create_test_planet() -> Planet {
        // TODO: Replace with your actual Planet initialization
        Planet::new(PlanetParams {
            radius: EARTH_RADIUS_KM,
            mantle_density_gcm3: Some(RHO_EARTH),
            sim_id: Uuid::new_v4(),
        })
    }

    #[test]
    fn test_run_and_report_mean_energy() {
        // Create a temporary directory for the RocksDB store
        let temp_dir = tempdir().expect("failed to create temp dir");
        let db_path = temp_dir.path().to_str().unwrap().to_string();

        let planet = create_test_planet();
        let mio_years = STANDARD_STEPS;
        let resolution = ASTH_RES;

        // Create the simulation
        let sim = AsthSim::new(AsthSimArgs {
            db_path,
            planet: &planet,
            mio_years,
            resolution,
        });

        // Run the simulation for all steps
        sim.run();

        let mut energy_sums: HashMap<u32, f32> = HashMap::new();
        let mut counts: HashMap<u32, usize> = HashMap::new();

        H3Utils::iter_at(ASTH_RES, |l2_cell| {
            for step in 1..mio_years {
                if let Ok(Some(cell)) = sim.store.get_asth(l2_cell, step) {
                    *energy_sums.entry(step).or_insert(0.0) += cell.energy_k as f32;
                    *counts.entry(step).or_insert(0) += 1;
                }
            }
        });

        // Now compute mean energy per step and print
        for step in 1..mio_years {
            let total_energy = energy_sums.get(&step).copied().unwrap_or(0.0);
            let cell_count = counts.get(&step).copied().unwrap_or(0);
            let mean_energy = if cell_count == 0 {
                0.0
            } else {
                total_energy / cell_count as f32
            };
        }
    }
}
