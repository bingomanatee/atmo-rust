use rayon::prelude::*;
use crate::asth_process_cell::{process_cell, ProcessResult};
use crate::asthenosphere::{ASTH_RES, AsthenosphereCell, CellsForPlanetArgs};
use crate::h30_utils::H3Utils;
use crate::planet::Planet;
use crate::rock_store::RockStore;
use h3o::CellIndex;
use rand::Rng;

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
    pub volume_added: f32,
    pub volume_removed: f32,
    pub total_volume: f32,
    pub mean_energy: f32,
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
            on_cell: |cell: AsthenosphereCell| {
                store.put_asth(&cell);
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
            "{:<6} | {:<13} | {:<14} | {:<13} | {:<11}",
            "Step", "Volume Added", "Volume Removed", "Total Volume", "Mean Energy/cell"
        );
        println!("{}", "-".repeat(65));

        for current in 1..=self.mio_years {
            let stats = self.advance(current);
            if (current % 5) == 0 {
                println!(
                    "{:<6} | {:<13.2} | {:<14.2} | {:<13.2} | {:<11.2}",
                    current, stats.volume_added, stats.volume_removed, stats.total_volume, stats.mean_energy
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
        let results: Vec<ProcessResult> = cells
            .par_iter()
            .filter_map(|l2_cell| process_cell(&self.store, *l2_cell, to_mio_years, self.mio_years).ok())
            .collect();
        
        let volume_added = results.iter().map(|r| r.volume_added).sum();
        let volume_removed = results.iter().map(|r| r.volume_removed).sum();
        let total_volume = results.iter().map(|r| r.new_volume).sum();
        let total_energy: f32 = results.iter().map(|r| r.total_energy).sum();
        let cell_count = results.len() as f32;

        let mean_energy = if cell_count > 0.0 {
            total_energy / cell_count
        } else {
            0.0
                    };

        VolumeStats {
            volume_added,
            volume_removed,
            total_volume,
            mean_energy,
        }
    }
}

pub const AVERAGE_ENERGY_AT_START: f32 = 2.0;
pub const AVERAGE_ENERGY_AT_END: f32 = 0.5;

pub const AE_SPAN: f32 = AVERAGE_ENERGY_AT_START - AVERAGE_ENERGY_AT_END;
pub const K_PER_VOLUME: f32 = 1000.0;
pub const COOLING_RATE: f32 = 0.95;
#[cfg(test)]
mod tests {
    use super::*;
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
        let mio_years = 500;
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

        // For each step, collect all cells and compute mean energy
        for step in 1..mio_years {
            let mut total_energy = 0.0f32;
            let mut cell_count = 0usize;

            H3Utils::iter_at(ASTH_RES, |l2_cell| {
                if let Ok(Some(cell)) = sim.store.get_asth(l2_cell, step) {
                    total_energy += cell.energy_k as f32;
                    cell_count += 1;
                }
            });

            let mean_energy = if cell_count == 0 {
                0.0
            } else {
                total_energy / cell_count as f32
            };

          // println!("Step {}: mean energy = {:.2}", step, mean_energy);
        }
    }
}
