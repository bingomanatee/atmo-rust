use crate::asth_constants::{AVG_STARTING_VOLUME, CELL_ENERGY_START, K_PER_VOLUME, STANDARD_STEPS};
use crate::asthenosphere::{ASTH_RES, AsthenosphereCell, CellsForPlanetArgs};
use crate::cool_asth_cell::{ProcessResult, cool_asth_cell};
use crate::h3o_utils::H3Utils;
use crate::planet::Planet;
use crate::plate::PLATE_RESOLUTION;
use crate::rock_store::RockStore;
use h3o::{CellIndex, Resolution};
use rand::Rng;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use uuid::Uuid;
use crate::level_cells::level_cells;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct VolumeEnergyTransfer {
    pub id: Uuid,
    pub step: u32,
    pub from_cell: CellIndex,
    pub to_cell: CellIndex,
    pub volume: f64,
    pub energy: f64,
}

pub struct AsthSim<'a> {
    store: RockStore,
    planet: &'a Planet,
    mio_years: u32,
    resolution: Resolution
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

        AsthenosphereCell::initial_cells_for_planet(CellsForPlanetArgs {
            planet: planet.clone(),
            res: Some(ASTH_RES),
            energy_per_volume: K_PER_VOLUME,
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
            resolution
        }
    }

    pub fn run(&self, debug: bool) {
        // Print table header
        if debug {
            println!(
                "{:<6} | {:<13} | {:<14} | {:<13} | {:<11} | {:<11}",
                "Step",
                "Volume Added",
                "Volume Removed",
                "Total Volume",
                "Mean Energy",
                "Mean_volume"
            );
            println!("{}", "-".repeat(65));
        }

        for current in 1..=self.mio_years {
            let stats = self.cool_cells(current);
            if debug && (current < 10 || (current % 25) == 0) {
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
            self.level(current);
            self.add_random_material(current);
        }
    }
    fn add_energy(old_volume: f64, old_energy: f64, added_volume: f64) -> f64 {
        (old_volume * old_energy + added_volume * CELL_ENERGY_START) / (old_volume + added_volume)
    }
    fn cool_cells(&self, step: u32) -> VolumeStats {
        // Collect all cells at ASTH_RES into a Vec<CellIndex>
        let cells: Vec<CellIndex> = {
            let mut v = Vec::new();
            H3Utils::iter_at(ASTH_RES, |cell| v.push(cell));
            v
        };

        // Process cells in parallel using Rayon
        let results: Result<Vec<ProcessResult>, String> = cells
            .par_iter()
            .map(|l2_cell| cool_asth_cell(&self.store, *l2_cell, step, self.mio_years))
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
            volume_removed,
            total_volume,
            mean_energy: total_energy / cell_count,
            mean_volume: total_volume / cell_count,
        }
    }

    /// Main leveling function that processes all cells at PLATE_RESOLUTION for the given step.
    fn level(&self, current_step: u32) {
       let base_cells: Vec<CellIndex> = CellIndex::base_cells()
           .collect();
        
        base_cells
            .par_iter()
            .map(|base_cell| level_cells(*base_cell, &self.store, current_step, self.resolution))
            .collect();
    }

    /// Returns a vector of neighbors with volume less than the root cell's volume.
    fn lower_volume_neighbors(
        &self,
        root_cell: &AsthenosphereCell,
        current_step: u32,
    ) -> Vec<AsthenosphereCell> {
        root_cell
            .neighbors
            .iter()
            .filter_map(
                |neighbor_idx| match self.store.get_asth(*neighbor_idx, current_step) {
                    Ok(Some(neighbor_cell)) if neighbor_cell.volume < root_cell.volume => {
                        Some(neighbor_cell)
                    }
                    _ => None,
                },
            )
            .collect()
    }

    const LEVEL_SCALE: f64 = 0.001; // Much more subtle: 0.1% per step for realistic geological flow

    /**
    writes transfers for energy/volume based on volume difference between cells
    */
    fn create_volume_transfer_vectors(
        &self,
        root_cell: &AsthenosphereCell,
        lower_volume_neighbors: &[AsthenosphereCell],
        current_step: u32,
    ) {
        let mut gaps: HashMap<CellIndex, f64> = HashMap::new();
        let mut total_gap: f64 = 0.0;

        for neighbor in lower_volume_neighbors {
            let gap = (root_cell.volume - neighbor.volume) * AsthSim::LEVEL_SCALE;
            gaps.insert(neighbor.cell, gap);
            total_gap += gap;
        }

        let amt_to_move = total_gap.min(root_cell.volume / 4.0); // never move more than 25% of the amount

        let mut transfers = Vec::with_capacity(lower_volume_neighbors.len());

        for (idx, gap) in gaps {
            let volume = gap * amt_to_move / total_gap;
            let energy = root_cell.energy_k * volume / root_cell.volume;

            let transfer = VolumeEnergyTransfer {
                volume,
                energy,
                id: Uuid::new_v4(),
                from_cell: root_cell.cell,
                to_cell: idx,
                step: current_step,
            };
            transfers.push(transfer);
        }

        self.store.batch_write_transfers(transfers);
    }

    /// Add or subtract random material to each cell using balanced random distribution
    /// Material added is at CELL_ENERGY_START temperature
    /// Material removed reduces energy proportionally
    fn add_random_material(&self, current_step: u32) {
        use rand::SeedableRng;
        use rand::seq::SliceRandom;

        // Collect all cells at ASTH_RES into a Vec<CellIndex>
        let mut cells: Vec<CellIndex> = {
            let mut v = Vec::new();
            H3Utils::iter_at(ASTH_RES, |cell| v.push(cell));
            v
        };

        // Create deterministic RNG with step-based seed for reproducible results
        let mut rng = rand::rngs::StdRng::seed_from_u64(current_step as u64 * 12345);

        // Shuffle the cells for random distribution
        cells.shuffle(&mut rng);

        // Process cells sequentially to maintain deterministic RNG state
        let mut updated_cells = Vec::new();

        for &cell_index in &cells {
            // Get cell from previous step (current_step - 1) to evolve it to current_step
            let previous_step = if current_step == 0 {
                0
            } else {
                current_step - 1
            };

            if let Ok(Some(mut cell)) = self.store.get_asth(cell_index, previous_step) {
                // Generate random number from 1 to 20
                let random_action = rng.random_range(1..=20);

                match random_action {
                    1 => {
                        // Add small amount of material (1/400 of AVG_STARTING_VOLUME)
                        let volume_to_add = AVG_STARTING_VOLUME / 400.0;
                        let old_volume = cell.volume;
                        let new_volume = old_volume + volume_to_add;

                        // Calculate weighted average energy
                        let old_energy = cell.energy_k;
                        let added_energy = volume_to_add * (CELL_ENERGY_START / AVG_STARTING_VOLUME);
                        let new_energy = (old_energy * old_volume + added_energy * volume_to_add) / new_volume;

                        cell.volume = new_volume;
                        cell.energy_k = new_energy;
                        cell.step = current_step;
                        updated_cells.push(cell);
                    },
                    2 => {
                        // Add large amount of material (1/100 of AVG_STARTING_VOLUME)
                        let volume_to_add = AVG_STARTING_VOLUME / 100.0;
                        let old_volume = cell.volume;
                        let new_volume = old_volume + volume_to_add;

                        // Calculate weighted average energy
                        let old_energy = cell.energy_k;
                        let added_energy = volume_to_add * (CELL_ENERGY_START / AVG_STARTING_VOLUME);
                        let new_energy = (old_energy * old_volume + added_energy * volume_to_add) / new_volume;

                        cell.volume = new_volume;
                        cell.energy_k = new_energy;
                        cell.step = current_step;
                        updated_cells.push(cell);
                    },
                    19 => {
                        // Remove small amount of material (1/400 of AVG_STARTING_VOLUME)
                        let volume_to_remove = AVG_STARTING_VOLUME / 400.0;
                        if cell.volume > volume_to_remove {
                            let old_volume = cell.volume;
                            let new_volume = old_volume - volume_to_remove;

                            // Scale energy proportionally
                            cell.energy_k *= new_volume / old_volume;
                            cell.volume = new_volume;
                            cell.step = current_step;
                            updated_cells.push(cell);
                        }
                    },
                    20 => {
                        // Remove large amount of material (1/100 of AVG_STARTING_VOLUME)
                        let volume_to_remove = AVG_STARTING_VOLUME / 100.0;
                        if cell.volume > volume_to_remove {
                            let old_volume = cell.volume;
                            let new_volume = old_volume - volume_to_remove;

                            // Scale energy proportionally
                            cell.energy_k *= new_volume / old_volume;
                            cell.volume = new_volume;
                            cell.step = current_step;
                            updated_cells.push(cell);
                        }
                    },
                    _ => {
                        // No change for numbers 3-18 (16 out of 20 cases = 80% no change)
                        cell.step = current_step;
                        updated_cells.push(cell);
                    }
                }
            }
        }

        // Batch write all updated cells
        if !updated_cells.is_empty() {
            self.store.put_asth_batch(&updated_cells);
        }
    }

    /// Returns Vec of (root_cell, neighbor_cell) pairs where neighbor volume < root volume at given step.
    /// This is a test helper to validate leveling flow detection.
    pub fn level_flows(&self, current_step: u32) -> Vec<(AsthenosphereCell, AsthenosphereCell)> {
        let mut flows = Vec::new();

        H3Utils::iter_at(PLATE_RESOLUTION, |idx| {
            if let Ok(Some(root_cell)) = self.store.get_asth(idx, current_step) {
                let neighbors: Vec<AsthenosphereCell> = root_cell
                    .neighbors
                    .iter()
                    .filter_map(|idx2| match self.store.get_asth(*idx2, current_step) {
                        Ok(Some(neighbor_cell)) if neighbor_cell.volume < root_cell.volume => {
                            Some(neighbor_cell)
                        }
                        _ => None,
                    })
                    .collect();

                for neighbor in neighbors {
                    flows.push((root_cell.clone(), neighbor));
                }
            }
        });

        flows
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::planet::EARTH;
    use tempfile::tempdir;

    #[test]
    fn test_png_voronoi_export_500() {
        use crate::png_exporter::PngExporter;

        // Create temporary RocksDB directory
        let temp_dir = tempdir().expect("failed to create temp dir");
        let db_path = temp_dir.path().to_str().unwrap().to_string();

        // Create a Planet
        let planet = EARTH.clone();

        // Create AsthSim using the proper constructor
        let sim = AsthSim::new(AsthSimArgs {
            db_path,
            planet: &planet,
            mio_years: 500,
            resolution: ASTH_RES,
        });

        // Run simulation for 500 steps with enhanced material changes
        println!("Running 500-step simulation for voronoi...");
        sim.run(true);
        let target = "vis/asthenosphere_pngs_500_final";

        // Create PNG exporter and generate Voronoi images
        println!("Generating PNG images with Voronoi patterns...");
        let mut png_exporter = PngExporter::new(800, 400, planet.clone());

        // Export every 25th step (21 images total: 0, 25, 50, 75, ..., 500)
        let result = png_exporter.export_asthenosphere_pngs(
            &sim.store, 0,   // start_step
            500, // end_step
            25,  // step_interval - render every 25th step
            target,
        );

        match result {
            Ok(_) => println!("Successfully created PNG images in {} directory", target),
            Err(e) => eprintln!("Failed to create PNGs: {}", e),
        }
    }
}
