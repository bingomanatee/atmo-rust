use crate::asth_constants::{STANDARD_STEPS, AVG_STARTING_VOLUME, CELL_ENERGY_START};
use crate::asth_process_cell::{ProcessResult, process_cell};
use crate::asthenosphere::{ASTH_RES, AsthenosphereCell, CellsForPlanetArgs};
use crate::geoconverter::GeoCellConverter;
use crate::h30_utils::H3Utils;
use crate::planet::Planet;
use crate::plate::PLATE_RESOLUTION;
use crate::rock_store::RockStore;
use h3o::{CellIndex, Resolution};
use noise::{NoiseFn, Perlin};
use rand::Rng;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use uuid::Uuid;

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
            if current < 10 || (current % 25) == 0 {
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

    /// Main leveling function that processes all cells at PLATE_RESOLUTION for the given step.
    fn level(&self, current_step: u32) {
        H3Utils::iter_at(PLATE_RESOLUTION, |idx| {
            if let Ok(Some(root_cell)) = self.store.get_asth(idx, current_step) {
                let lower_volume_neighbors = self.lower_volume_neighbors(&root_cell, current_step);
                self.create_volume_transfer_vectors(&root_cell, &lower_volume_neighbors, current_step);
            }
        });
        
        self.store.transfer_volume(current_step);
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
            .filter_map(|neighbor_idx| match self.store.get_asth(*neighbor_idx, current_step) {
                Ok(Some(neighbor_cell)) if neighbor_cell.volume < root_cell.volume => Some(neighbor_cell),
                _ => None,
            })
            .collect()
    }

    const LEVEL_SCALE: f64 = 0.33; // Much more subtle: 0.1% per step for realistic geological flow

    /**
    writes transfers for energy/volume based on volume difference between cells
    */
    fn create_volume_transfer_vectors(
        &self,
        root_cell: &AsthenosphereCell,
        lower_volume_neighbors: &[AsthenosphereCell],
        current_step: u32
    ) {
        let mut gaps: HashMap<CellIndex, f64> = HashMap::new();
        let mut total_gap : f64 = 0.0;
        
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

    /// Add or subtract random material to each cell based on Perlin noise
    /// Material added is at CELL_ENERGY_START temperature
    /// Material removed reduces energy proportionally
    fn add_random_material(&self, current_step: u32) {
        // Collect all cells at ASTH_RES into a Vec<CellIndex>
        let cells: Vec<CellIndex> = {
            let mut v = Vec::new();
            H3Utils::iter_at(ASTH_RES, |cell| v.push(cell));
            v
        };

        // Create Perlin noise generator with step-based seed - one per frame
        let perlin = Perlin::new(current_step * 5);
        let gc = GeoCellConverter::new(self.planet.radius_km as f64, ASTH_RES);

        // Much more subtle material change: ±0.5% of AVG_STARTING_VOLUME for gradual evolution
        let max_material_change = AVG_STARTING_VOLUME * 0.005;

        // Process cells in parallel using Rayon
        let updated_cells: Vec<AsthenosphereCell> = cells
            .par_iter()
            .filter_map(|&cell_index| {
                // Get cell from previous step (current_step - 1) to evolve it to current_step
                let previous_step = if current_step == 0 { 0 } else { current_step - 1 };
                if let Ok(Some(mut cell)) = self.store.get_asth(cell_index, previous_step) {
                    // Use the same Perlin generator for all cells in this step
                    // This creates coherent patterns that evolve frame by frame

                    // Get cell location for noise sampling
                    let location = gc.cell_to_vec3(cell_index);
                    let [x, y, z] = location.normalize().to_array();

                    // Scale up the coordinates to create more detailed noise patterns
                    // Higher scale = more detailed/frequent variations
                    let noise_scale = 8.0; // Moderate scale for balanced patterns
                    let scaled_x = x as f64 * noise_scale;
                    let scaled_y = y as f64 * noise_scale;
                    let scaled_z = z as f64 * noise_scale;

                    // Sample Perlin noise (-1.0 to 1.0) with scaled coordinates using frame-wide seed
                    let raw_noise = perlin.get([scaled_x, scaled_y, scaled_z]);

                    // Apply exponential transformation for moderate spikes
                    // This creates localized changes but not too extreme
                    let exponential_noise = if raw_noise > 0.0 {
                        raw_noise.powf(2.0) // Square for moderate peaks
                    } else {
                        -((-raw_noise).powf(2.0)) // Square for moderate valleys
                    };

                    // Calculate material change with exponential spikes
                    let material_change = exponential_noise * max_material_change;

                    if material_change > 0.0 {
                        // Adding material at CELL_ENERGY_START temperature
                        let added_energy = material_change * (CELL_ENERGY_START / AVG_STARTING_VOLUME);
                        cell.volume += material_change;
                        cell.energy_k += added_energy;
                    } else if material_change < 0.0 {
                        // Removing material - reduce energy proportionally
                        let volume_to_remove = -material_change;
                        if cell.volume > volume_to_remove {
                            let energy_ratio = cell.energy_k / cell.volume;
                            let energy_to_remove = volume_to_remove * energy_ratio;

                            cell.volume -= volume_to_remove;
                            cell.energy_k -= energy_to_remove;
                        }
                        // If we would remove more volume than exists, do nothing
                    }

                    // Update the step to current_step so it gets stored correctly
                    cell.step = current_step;
                    Some(cell)
                } else {
                    None
                }
            })
            .collect();

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
        use h3o::Resolution;
        use tempfile::tempdir;
        use uuid::Uuid;

        #[test]
        #[ignore] // Skip this test - not a good test yet
        fn test_level_flows_with_real_rockstore() {
            // Step: 1
            let step = 1;

            // Create temporary RocksDB directory
            let temp_dir = tempdir().expect("failed to create temp dir");
            let db_path = temp_dir.path().to_str().unwrap().to_string();

            // Create a Planet (use your existing helper or dummy)
            let planet = EARTH.clone();

            // Open RockStore
            let store = RockStore::open(&db_path).expect("failed to open RockStore");

            // Find the first L2 cell at PLATE_RESOLUTION
            let mut root_idx_opt = None;
            H3Utils::iter_at(PLATE_RESOLUTION, |cell| {
                if root_idx_opt.is_none() {
                    root_idx_opt = Some(cell);
                }
            });
            let root_idx = root_idx_opt.expect("No cells found at PLATE_RESOLUTION");

            // For the root cell neighbors, get the first two neighbors (if available)
            let neighbors = {
                let mut n = Vec::new();
                if let Ok(Some(root_cell)) = store.get_asth(root_idx, step) {
                    n = root_cell.neighbors.clone();
                }
                n
            };

            // If neighbors are empty, create dummy neighbors by getting neighbors of root_idx from H3Utils
            let neighbors = if neighbors.is_empty() {
                let mut n = Vec::new();
                H3Utils::iter_at(PLATE_RESOLUTION, |cell| {
                    if cell != root_idx && n.len() < 2 {
                        n.push(cell);
                    }
                });
                n
            } else {
                neighbors
            };

            // Create root cell with volume 10.0 and neighbors
            let root_cell = AsthenosphereCell {
                cell: root_idx,
                volume: 10.0,
                energy_k: 100.0,
                neighbors: neighbors.clone(),
                ..AsthenosphereCell::default()
            };

            // Insert root cell into store at step
            store.put_asth(&root_cell).expect("failed to put root cell");

            // Insert neighbors with volumes: one lower, one higher
            if neighbors.len() >= 2 {
                let neighbor_low = AsthenosphereCell {
                    cell: neighbors[0],
                    volume: 5.0,
                    energy_k: 50.0,
                    neighbors: vec![],
                    ..Default::default()
                };
                let neighbor_high = AsthenosphereCell {
                    cell: neighbors[1],
                    volume: 15.0,
                    energy_k: 150.0,
                    neighbors: vec![],
                    ..Default::default()
                };
                store
                    .put_asth(&neighbor_low)
                    .expect("failed to put low neighbor");
                store
                    .put_asth(&neighbor_high)
                    .expect("failed to put high neighbor");

                // Create sim with this store
                let sim = AsthSim {
                    store,
                    planet: &planet,
                    mio_years: 1,
                };

                // Call level_flows at step
                let flows = sim.level_flows(step);

                // We expect exactly one flow: root -> neighbor_low
                assert_eq!(flows.len(), 1);
                assert_eq!(flows[0].0.cell, root_idx);
                assert_eq!(flows[0].1.cell, neighbors[0]);
            } else {
                panic!("Not enough neighbors found for root cell");
            }
        }

        #[test]
        fn test_add_random_material() {
            // Create temporary RocksDB directory
            let temp_dir = tempdir().expect("failed to create temp dir");
            let db_path = temp_dir.path().to_str().unwrap().to_string();

            // Create a Planet
            let planet = EARTH.clone();

            // Open RockStore
            let store = RockStore::open(&db_path).expect("failed to open RockStore");

            let step = 1u32;

            // Create a test cell at ASTH_RES (Resolution::Two)
            let cell_index = CellIndex::try_from(0x85283473fffffff).unwrap();
            // Ensure the cell is at the correct resolution
            let cell_index = cell_index.parent(ASTH_RES).unwrap_or(cell_index);

            let initial_volume = 1000.0;
            let initial_energy = 2000.0;

            let test_cell = AsthenosphereCell {
                cell: cell_index,
                volume: initial_volume,
                energy_k: initial_energy,
                step,
                neighbors: vec![],
            };

            // Store the cell
            store.put_asth(&test_cell).expect("failed to store test cell");

            // Create AsthSim
            let sim = AsthSim {
                store,
                planet: &planet,
                mio_years: 1,
            };

            // Call add_random_material
            sim.add_random_material(step);

            // Retrieve the updated cell
            let updated_cell = sim.store.get_asth(cell_index, step)
                .expect("failed to get updated cell")
                .expect("updated cell not found");

            // Calculate the expected noise value for verification
            let perlin = Perlin::new(step);
            let gc = GeoCellConverter::new(sim.planet.radius_km as f64, ASTH_RES);
            let location = gc.cell_to_vec3(cell_index);
            let [x, y, z] = location.normalize().to_array();

            // Use the same scaling as in add_random_material
            let noise_scale = 8.0;
            let scaled_x = x as f64 * noise_scale;
            let scaled_y = y as f64 * noise_scale;
            let scaled_z = z as f64 * noise_scale;
            let noise_val = perlin.get([scaled_x, scaled_y, scaled_z]);
            let max_material_change = AVG_STARTING_VOLUME * 0.005;
            let material_change = noise_val * max_material_change;

            // Verify that the cell was modified (volume and energy should have changed)
            // Since we're using Perlin noise with a fixed seed (step), the result should be deterministic
            assert_ne!(updated_cell.volume, initial_volume, "Volume should have changed");
            assert_ne!(updated_cell.energy_k, initial_energy, "Energy should have changed");

            // Verify the volume change matches the expected material change
            let actual_volume_change = updated_cell.volume - initial_volume;
            assert!((actual_volume_change - material_change).abs() < 0.001,
                "Volume change {} should match material change {}", actual_volume_change, material_change);

            // Verify that the changes are within expected bounds (±1% of AVG_STARTING_VOLUME)
            let max_change = AVG_STARTING_VOLUME * 0.01;
            let volume_change = (updated_cell.volume - initial_volume).abs();
            assert!(volume_change <= max_change, "Volume change {} exceeds maximum {}", volume_change, max_change);

            // Verify that energy changes are reasonable
            if updated_cell.volume > initial_volume {
                // Material was added, energy should increase
                assert!(updated_cell.energy_k > initial_energy, "Energy should increase when material is added");
            } else {
                // Material was removed, energy should decrease
                assert!(updated_cell.energy_k < initial_energy, "Energy should decrease when material is removed");
            }
        }

        #[test]
        fn test_500_step_gif_animation() {
            use crate::gif_exporter::GifExporter;

            // Create temporary RocksDB directory
            let temp_dir = tempdir().expect("failed to create temp dir");
            let db_path = temp_dir.path().to_str().unwrap().to_string();

            // Create a Planet
            let planet = EARTH.clone();

            // Open RockStore and create AsthSim
            let store = RockStore::open(&db_path).expect("failed to open RockStore");
            let sim = AsthSim {
                store,
                planet: &planet,
                mio_years: 500, // 500 million years = 500 steps
            };

            // Initialize asthenosphere cells (step 0)
            AsthenosphereCell::cells_for_planet(CellsForPlanetArgs {
                planet: planet.clone(),
                res: Some(ASTH_RES),
                energy_per_volume: 200.0,
                on_cell: |asth_cell: AsthenosphereCell| {
                    match sim.store.put_asth(&asth_cell) {
                        Ok(_) => {},
                        Err(e) => eprintln!("Failed to store cell {:?}: {:?}", asth_cell.cell, e),
                    }
                    asth_cell
                },
            });

            // Run simulation for 500 steps
            println!("Running 500-step simulation...");
            for step in 1..=500 {
                if step % 50 == 0 {
                    println!("Step {}/500", step);
                }
                sim.level(step);
                sim.add_random_material(step);
            }

            // Create GIF exporter and generate animation
            println!("Generating animated GIF...");
            let gif_exporter = GifExporter::new(800, 400, planet.clone());

            // Export every 5th step (100 frames total: 0, 5, 10, 15, ..., 500)
            let result = gif_exporter.export_asthenosphere_animation(
                &sim.store,
                0,    // start_step
                500,  // end_step
                5,    // step_interval - render every 5th step
                "asthenosphere_500_steps.gif"
            );

            match result {
                Ok(_) => println!("Successfully created asthenosphere_500_steps.gif with 100 frames"),
                Err(e) => eprintln!("Failed to create GIF: {}", e),
            }
        }

        #[test]
        fn test_png_voronoi_export_500() {
            use crate::png_exporter::PngExporter;

            // Create temporary RocksDB directory
            let temp_dir = tempdir().expect("failed to create temp dir");
            let db_path = temp_dir.path().to_str().unwrap().to_string();

            // Create a Planet
            let planet = EARTH.clone();

            // Open RockStore and create AsthSim
            let store = RockStore::open(&db_path).expect("failed to open RockStore");
            let sim = AsthSim {
                store,
                planet: &planet,
                mio_years: 500, // 500 million years = 500 steps
            };

            // Initialize asthenosphere cells (step 0)
            AsthenosphereCell::cells_for_planet(CellsForPlanetArgs {
                planet: planet.clone(),
                res: Some(ASTH_RES),
                energy_per_volume: 200.0,
                on_cell: |asth_cell: AsthenosphereCell| {
                    match sim.store.put_asth(&asth_cell) {
                        Ok(_) => {},
                        Err(e) => eprintln!("Failed to store cell {:?}: {:?}", asth_cell.cell, e),
                    }
                    asth_cell
                },
            });

            // Run simulation for 500 steps with enhanced material changes
            println!("Running 500-step simulation with enhanced exponential spikes and per-cell seeding...");
            for step in 1..=500 {
                if step % 50 == 0 {
                    println!("Step {}/500", step);
                }
                // Only do random material changes - this will create evolving patterns
                sim.add_random_material(step);
            }

            // Create PNG exporter and generate Voronoi images
            println!("Generating PNG images with Voronoi patterns...");
            let mut png_exporter = PngExporter::new(800, 400, planet.clone());

            // Export every 25th step (21 images total: 0, 25, 50, 75, ..., 500)
            let result = png_exporter.export_asthenosphere_pngs(
                &sim.store,
                0,    // start_step
                500,  // end_step
                25,   // step_interval - render every 25th step
                "asthenosphere_pngs_500_final"
            );

            match result {
                Ok(_) => println!("Successfully created PNG images in asthenosphere_pngs_500_final/ directory"),
                Err(e) => eprintln!("Failed to create PNGs: {}", e),
            }
        }

        #[test]
        fn test_png_voronoi_export_simple() {
            use crate::png_exporter::PngExporter;

            // Create temporary RocksDB directory
            let temp_dir = tempdir().expect("failed to create temp dir");
            let db_path = temp_dir.path().to_str().unwrap().to_string();

            // Create a Planet
            let planet = EARTH.clone();

            // Open RockStore and create AsthSim
            let store = RockStore::open(&db_path).expect("failed to open RockStore");
            let sim = AsthSim {
                store,
                planet: &planet,
                mio_years: 500, // 500 million years = 500 steps
            };

            // Initialize asthenosphere cells (step 0)
            AsthenosphereCell::cells_for_planet(CellsForPlanetArgs {
                planet: planet.clone(),
                res: Some(ASTH_RES),
                energy_per_volume: 200.0,
                on_cell: |asth_cell: AsthenosphereCell| {
                    match sim.store.put_asth(&asth_cell) {
                        Ok(_) => {},
                        Err(e) => eprintln!("Failed to store cell {:?}: {:?}", asth_cell.cell, e),
                    }
                    asth_cell
                },
            });

            // Run the full 500-step simulation using the proper run() method
            println!("Running 500-step simulation with leveling and random material changes...");
            sim.run();

            // Create PNG exporter and generate Voronoi images
            println!("Generating PNG images with Voronoi patterns...");
            let mut png_exporter =      PngExporter::new(800, 400, planet.clone());

            // Sub-sample key steps from the 500-step simulation for detailed analysis
            // Export steps: 0, 1, 2, 3, 4, 5, 10, 15, 20, 25, 50, 75, 100, 150, 200, 250, 300, 400, 500
            let result = png_exporter.export_asthenosphere_pngs(
                &sim.store,
                0,    // start_step
                500,  // end_step
                25,   // step_interval - render every 25th step (plus we'll get 0, 25, 50, 75, etc.)
                "asthenosphere_pngs_simple"
            );

            match result {
                Ok(_) => println!("Successfully created PNG images (sub-sampled from 500 steps) in asthenosphere_pngs_simple/ directory"),
                Err(e) => eprintln!("Failed to create PNGs: {}", e),
            }
        }
    }
