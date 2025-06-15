use crate::asthenosphere::{ASTH_RES, AsthenosphereCell, CellsForPlanetArgs};
use crate::h30_utils::H3Utils;
use crate::planet::Planet;
use crate::rock_store::RockStore;
use h3o::{CellIndex, Resolution};
use noise::{NoiseFn, Perlin};

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
        for current in 1..self.mio_years {
            self.advance(current);
        }
    }

    fn advance(&self, to_mio_years: u32) {
        let per = Perlin::new(100 + to_mio_years);
        H3Utils::iter_at(ASTH_RES, |l2_cell| {
            match self.store.get_asth(l2_cell, to_mio_years - 1) {
                Ok(Some(old_cell)) => {
                    let [x, y, z] = old_cell.location.to_array();
                    let rand = per.get([x as f64, y as f64, z as f64]);
                    let progress = to_mio_years as f32 / self.mio_years as f32;
                    let average = AVERAGE_ENEGY_AT_START as f32 - (AE_SPAN as f32 * progress);
                    let volume_to_add = average * (2.0 + rand as f32) / 3.0;
                    let sunk_volume = old_cell.sunk_volume();
                    let mut new_volume = old_cell.volume + volume_to_add - sunk_volume;

                    let mut new_energy = ((0.98 * old_cell.energy_k as f32 * old_cell.volume)
                        + (K_PER_VOLUME * volume_to_add))
                        / (volume_to_add + old_cell.volume);

                    if sunk_volume > 0.0 {
                        new_energy *= (old_cell.volume - sunk_volume) / old_cell.volume;
                    }

                    // @TODO: remove some mass to the lithosphere.

                    let new_cell = AsthenosphereCell {
                        step: to_mio_years,
                        volume: new_volume,
                        energy_k: new_energy as u32,
                        ..old_cell.clone()
                    };
                    self.store
                        .put_asth(&new_cell)
                        .expect("cannot save new cell")
                }
                Ok(None) => {
                    panic!("cannot find cell {}, {}", l2_cell, to_mio_years - 1)
                }
                Err(e) => {
                    panic!("Error fetching asth cell: {:?}", e);
                }
            }
        });
    }
}

const AVERAGE_ENEGY_AT_START: u32 = 20;
const AVERAGE_ENERGY_AT_END: u32 = 5;

const AE_SPAN: u32 = AVERAGE_ENEGY_AT_START - AVERAGE_ENERGY_AT_END;
const K_PER_VOLUME: f32 = 1800.0;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::planet::{EARTH_RADIUS_KM, PlanetParams, RHO_EARTH};
    use h3o::Resolution;
    use std::path::PathBuf;
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
        let mio_years = 100;
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
            let mut sum_energy = 0.0f32;
            let mut count = 0usize;

            H3Utils::iter_at(ASTH_RES, |l2_cell| {
                if let Ok(Some(cell)) = sim.store.get_asth(l2_cell, step) {
                    sum_energy += cell.energy_k as f32;
                    count += 1;
                }
            });

            let mean_energy = if count == 0 {
                0.0
            } else {
                sum_energy / count as f32
            };

            println!("Step {}: mean energy = {:.2}", step, mean_energy);
        }
    }
}
