use crate::asthenosphere::{ASTH_RES, AsthenosphereCell, CellsForPlanetArgs};
use crate::asthenosphere_linked::AsthenosphereCellLinked;
use crate::constants::{ANOMALY_DECAY_RATE, ANOMALY_ENERGY_AMOUNT, ANOMALY_SPAWN_CHANCE, ANOMALY_VOLUME_AMOUNT, AVG_STARTING_VOLUME_KM_3, CELL_JOULES_EQUILIBRIUM, CELL_JOULES_START, JOULES_PER_KM3, LEVEL_AMT};
use crate::planet::Planet;
use crate::png_exporter::PngExporter;
use crate::rock_store::RockStore;
use h3o::CellIndex;
use rand::Rng;
use std::cell::RefCell;
use std::collections::HashMap;
use std::fs;
use std::rc::Rc;

pub struct ASLParams {
    pub planet: Planet,
    pub steps: u64,
    pub store_path: String,
    pub visualize: bool,
    pub vis_freq: u64,
    pub debug: bool,
}

pub struct AsthSimLinked {
    cells: HashMap<CellIndex, Rc<RefCell<AsthenosphereCellLinked>>>,
    current_step: u64,
    steps: u64,
    planet: Planet,
    store: RockStore,
    vis_freq: u64,
    visualize: bool,
}

impl AsthSimLinked {
    pub fn new(config: ASLParams) -> AsthSimLinked {
        let mut sim = match RockStore::open(&config.store_path) {
            Ok(store) => AsthSimLinked {
                cells: HashMap::new(),
                store,
                planet: config.planet.clone(),
                current_step: 0,
                steps: config.steps,
                visualize: config.visualize,
                vis_freq: config.vis_freq,
            },
            Err(e) => {
                panic!("failed to load store: {:?}", e)
            }
        };

        let init_config = CellsForPlanetArgs {
            planet: config.planet.clone(),
            on_cell: |cell| {
                let linked_cell = AsthenosphereCellLinked::from_cell(&cell);
                sim.cells
                    .insert(cell.id, Rc::new(RefCell::new(linked_cell)));
                cell
            },
            res: ASTH_RES,
            joules_per_km3: CELL_JOULES_START / AVG_STARTING_VOLUME_KM_3,
            seed: 42,
        };

        AsthenosphereCell::initial_cells_for_planet(init_config);
        sim
    }

    pub fn run_step(&mut self) {
        println!("  🔥 Cooling cells...");
        self.cool_cells();

        println!("  📏 Leveling cells...");
        self.level_cells();

        println!("  🌀 Processing anomalies...");
        self.process_anomalies();

        println!("  ⚡ Spawning new anomalies...");
        self.spawn_anomalies();

        println!("  ⏭️  Advancing to next step...");
        self.advance_to_next_step();

        // Export visualization every 5 frames if enabled
        if self.visualize && self.current_step % self.vis_freq == 0 {
            println!("  🖼️  Exporting visualization...");
            self.export_visualization();
        }
    }

    fn export_visualization(&self) {
        println!("    📁 Creating output directory...");
        // Create output directory if it doesn't exist
        let output_dir = "vis/asth_sim_linked";
        if let Err(_) = fs::create_dir_all(output_dir) {
            eprintln!(
                "Warning: Could not create visualization directory {}",
                output_dir
            );
            return;
        }

        println!("    📊 Collecting {} cells for export...", self.cells.len());
        // Collect current cell data for visualization with simplified energy values for speed
        let cells_for_export: Vec<(CellIndex, AsthenosphereCell)> = self
            .cells
            .iter()
            .map(|(&cell_id, cell)| {
                let borrowed_cell = cell.borrow();
                let mut simplified_cell = borrowed_cell.cell.clone();
                // Round energy to reduce floating point precision for faster processing
                simplified_cell.energy_j = (simplified_cell.energy_j / 100.0).round() * 100.0;
                (cell_id, simplified_cell)
            })
            .collect();

        println!("    🎨 Creating PNG exporter (900x450)...");
        // Create PNG exporter with moderate resolution for faster rendering
        let mut exporter = PngExporter::new(900, 450, self.planet.clone());

        println!("    🖌️  Rendering Voronoi image...");
        let image = exporter.render_voronoi_image_from_cells(&cells_for_export);

        println!("    💾 Saving image file...");
        // Save the image
        let filename = format!("{}/step_{:04}.png", output_dir, self.current_step);
        match image.save(&filename) {
            Ok(_) => println!("    ✅ Exported visualization: {}", filename),
            Err(e) => eprintln!(
                "    ❌ Failed to export visualization for step {}: {}",
                self.current_step, e
            ),
        }
    }

    fn advance_to_next_step(&mut self) {
        // 1. Save each cell to disk using batch save
        let cells_to_save: Vec<AsthenosphereCell> = self
            .cells
            .values()
            .map(|cell| {
                let mut cell_data = cell.borrow().cell.clone();
                cell_data.step = self.current_step as u32;
                cell_data
            })
            .collect();

        if !cells_to_save.is_empty() {
            self.store.put_asth_batch(&cells_to_save);
        }

        // 2. Increment current_step
        self.current_step += 1;

        // 3. Collect cell IDs that have next cells
        let cell_ids_with_next: Vec<CellIndex> = self
            .cells
            .iter()
            .filter_map(|(&cell_id, current_cell)| {
                if current_cell.borrow().next.is_some() {
                    Some(cell_id)
                } else {
                    None
                }
            })
            .collect();

        // 4. For each cell with a next: promote next to current, trim tails, add new next
        for cell_id in cell_ids_with_next {
            if let Some(current_cell) = self.cells.get(&cell_id) {
                let next_cell = {
                    let current_borrowed = current_cell.borrow();
                    current_borrowed.next.as_ref().map(|nc| Rc::clone(nc))
                };

                if let Some(next_cell) = next_cell {
                    // Update the step on the next cell
                    next_cell.borrow_mut().cell.step = self.current_step as u32;

                    // Trim the tail (remove previous links)
                    next_cell.borrow_mut().unlink_all_prev();

                    // Update hash to point to the next cell (now current)
                    self.cells.insert(cell_id, next_cell.clone());

                    // Add a new next cell
                    AsthenosphereCellLinked::add(&next_cell);
                }
            }
        }
    }

    pub fn cool_cells(&self) {
        let cool_rate =
            (CELL_JOULES_EQUILIBRIUM / CELL_JOULES_START).powf(0.7 / (self.steps as f64));
        for (_id, cell) in &self.cells {
            AsthenosphereCellLinked::add(cell);

            // Access the next cell through borrowing
            if let Some(next_cell) = &cell.borrow().next {
                let mut next = next_cell.borrow_mut();
                next.cell.energy_j *= cool_rate;
            }
        }
    }

    pub fn level_cells(&self) {
        // Level between spatial neighbors, reading from current states and writing to next states
        for (_id, cell_a) in &self.cells {
            let neighbors = &cell_a.borrow().cell.neighbors;
            for neighbor_id in neighbors {
                if let Some(cell_b) = self.cells.get(neighbor_id) {
                    self.level_between_neighbors(cell_a, cell_b);
                }
            }
        }
    }

    fn level_between_neighbors(
        &self,
        cell_a: &Rc<RefCell<AsthenosphereCellLinked>>,
        cell_b: &Rc<RefCell<AsthenosphereCellLinked>>,
    ) {
        // PARADIGM: Read current cells, write to next cells
        // This allows us to read consistent state while writing changes for next step
        let current_a = cell_a.borrow();
        let current_b = cell_b.borrow();

        let volume_a = current_a.cell.volume;
        let volume_b = current_b.cell.volume;
        let volume_diff = volume_a - volume_b;

        if volume_diff.abs() < f64::EPSILON {
            return; // No significant difference
        }

        // Only flow downhill - from higher volume to lower volume
        if volume_diff > 0.0 {
            // Moderate exponential leveling: smooth gradual scaling instead of sharp threshold
            // Use volume difference relative to average starting volume as the baseline
            let relative_diff = volume_diff / AVG_STARTING_VOLUME_KM_3;
            

            // Calculate transfer amount with moderate exponential scaling
            let base_transfer = volume_diff * LEVEL_AMT;
            let transfer_amount = base_transfer;

            // Calculate energy per unit volume for source cell
            let energy_per_volume_a = if volume_a > 0.0 {
                current_a.cell.energy_j / volume_a
            } else {
                0.0
            };

            let energy_to_transfer = transfer_amount * energy_per_volume_a;

            // Write changes to next states only - both must exist
            if let (Some(next_a), Some(next_b)) = (&current_a.next, &current_b.next) {
                let mut next_a_b = next_a.borrow_mut();
                let mut next_b_b = next_b.borrow_mut();
                next_a_b.cell.volume -= transfer_amount;
                next_b_b.cell.volume += transfer_amount;
                next_a_b.cell.energy_j -= energy_to_transfer;
                next_b_b.cell.energy_j += energy_to_transfer;
            }
        }
    }

    pub fn process_anomalies(&self) {
        // PARADIGM: Read current cells, write to next cells
        // Decay existing anomalies and apply their effects
        for (_id, cell) in &self.cells {
            let current = cell.borrow();

            // Only process if there are anomalies and next cell exists
            if (current.cell.anomaly_energy.abs() > 1e-6
                || current.cell.anomaly_volume.abs() > 1e-6)
                && current.next.is_some()
            {
                if let Some(next_cell) = &current.next {
                    let mut next = next_cell.borrow_mut();

                    // Apply anomaly effects to actual energy and volume
                    next.cell.energy_j += current.cell.anomaly_energy;
                    next.cell.volume += current.cell.anomaly_volume;

                    // Decay anomalies by 3% per step - the anomaly fields themselves decay
                    next.cell.anomaly_energy =
                        current.cell.anomaly_energy * (1.0 - ANOMALY_DECAY_RATE);
                    next.cell.anomaly_volume =
                        current.cell.anomaly_volume * (1.0 - ANOMALY_DECAY_RATE);

                    // Set to zero if negligible (less than +/-1.0)
                    if next.cell.anomaly_energy.abs() < 1.0 {
                        next.cell.anomaly_energy = 0.0;
                    }
                    if next.cell.anomaly_volume.abs() < 1.0 {
                        next.cell.anomaly_volume = 0.0;
                    }
                }
            }
        }
    }

    pub fn spawn_anomalies(&self) {
        let mut rng = rand::thread_rng();

        // 15% chance to spawn an anomaly each cycle
        if rng.random::<f64>() < ANOMALY_SPAWN_CHANCE {
            // Select a random cell
            let cell_indices: Vec<_> = self.cells.keys().collect();
            if let Some(&random_cell_id) =
                cell_indices.get((rng.random::<f64>() * cell_indices.len() as f64) as usize)
            {
                if let Some(target_cell) = self.cells.get(random_cell_id) {
                    self.apply_anomaly_to_cell_and_neighbors(target_cell, &mut rng);
                }
            }
        }
    }

    fn apply_anomaly_to_cell_and_neighbors(
        &self,
        target_cell: &Rc<RefCell<AsthenosphereCellLinked>>,
        rng: &mut impl Rng,
    ) {
        // PARADIGM: Read current cells, write to next cells
        let current = target_cell.borrow();
        let mut rand = rand::rng();
        let sign = if rng.random::<bool>() { 1.0 } else { -1.0 };
        let size = rand.random_range(0.01..0.1);
        let volume = ANOMALY_VOLUME_AMOUNT * sign * size;
        // Apply to target cell
        if let Some(next_cell) = &current.next {
            let mut next = next_cell.borrow_mut();

            // Randomly add or remove energy/volume (positive or negative anomaly)
         
            next.cell.anomaly_volume += volume;
            next.cell.anomaly_energy += volume * JOULES_PER_KM3;
        }

        // Apply smaller effects to neighbors (50% of main effect)
        let neighbors = current.cell.neighbors.clone();
        drop(current); // Release borrow before iterating neighbors
        for neighbor_id in neighbors {
            if let Some(neighbor_cell) = self.cells.get(&neighbor_id) {
                let neighbor_current = neighbor_cell.borrow();
                if let Some(neighbor_next) = &neighbor_current.next {
                    let mut neighbor_next_mut = neighbor_next.borrow_mut();
                    neighbor_next_mut.cell.anomaly_volume += volume / 6.0;
                    neighbor_next_mut.cell.anomaly_energy += volume * JOULES_PER_KM3 / 6.0;
                }
            }
        }
    }

    /// Get the number of cells in the simulation
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }

    /// Run the simulation for a specified number of steps
    pub fn run_simulation(&mut self, steps: u64) {
        println!(
            "🌍 Starting AsthSimLinked simulation for {} steps...",
            steps
        );
        println!("📊 Initial state: {} cells loaded", self.cells.len());

        for step in 1..=steps {
            println!("🔄 Running step {} of {}:", step, steps);
            self.run_step();

            if step % 10 == 0 {
                println!("✅ Completed step {} of {}", step, steps);
            }

            if self.visualize && step % 5 == 0 {
                println!("🖼️  Exported visualization for step {}", self.current_step);
            }
        }

        if self.visualize {
            println!(
                "🎉 Simulation completed! {} PNG files exported to vis/asth_sim_linked/",
                (steps / 5) + 1
            );
        } else {
            println!("🎉 Simulation completed!");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::EARTH;

    #[test]
    fn test_asth_sim_linked_initial_cells_count() {
        let config = ASLParams {
            planet: EARTH.clone(),
            steps: 10,
            store_path: String::from("/tmp/test_rocksdb_store"),
            visualize: false,
            vis_freq: 0,
            debug: false,
        };

        let sim = AsthSimLinked::new(config);

        assert!(sim.cells.len() > 0, "Expected some cells to be initialized");
    }

    #[test]
    fn test_simulation_run_through() {
        println!("🚀 Testing AsthSimLinked full simulation run-through...");

        let config = ASLParams {
            planet: EARTH.clone(),
            steps: 25, // Run for 25 steps to see multiple visualization exports
            store_path: String::from("/tmp/test_asth_sim_linked"),
            visualize: true,
            vis_freq: 0,
            debug: false,
        };

        let mut sim = AsthSimLinked::new(config);

        // Run the simulation for 25 steps (will export PNGs at steps 5, 10, 15, 20, 25)
        sim.run_simulation(25);

        // Verify final state
        assert!(
            sim.cells.len() > 0,
            "Expected cells to remain after simulation"
        );
        assert_eq!(sim.current_step, 25, "Expected to reach step 25");

        println!("✅ Simulation run-through test completed successfully!");
    }

    #[test]
    fn test_level_cells() {
        let config = ASLParams {
            planet: EARTH.clone(),
            steps: 10,
            store_path: String::from("/tmp/test_rocksdb_store_level"),
            visualize: false,
            vis_freq: 0,
            debug: false,
        };

        let mut sim = AsthSimLinked::new(config);

        // Find two cells that are neighbors of each other
        let mut cell_a = None;
        let mut cell_b = None;

        for (_id, cell) in &sim.cells {
            let neighbors = &cell.borrow().cell.neighbors;
            if !neighbors.is_empty() {
                if let Some(neighbor_id) = neighbors.first() {
                    if let Some(neighbor_cell) = sim.cells.get(neighbor_id) {
                        cell_a = Some(cell);
                        cell_b = Some(neighbor_cell);
                        break;
                    }
                }
            }
        }

        if let (Some(cell_a), Some(cell_b)) = (cell_a, cell_b) {
            // Set different volumes to test levelling first
            cell_a.borrow_mut().cell.volume = 1000.0;
            cell_a.borrow_mut().cell.energy_j = 2000.0;
            cell_b.borrow_mut().cell.volume = 500.0;
            cell_b.borrow_mut().cell.energy_j = 1000.0;

            // Create next cells for both - these will copy the current modified values
            let _next_a = AsthenosphereCellLinked::add(cell_a);
            let _next_b = AsthenosphereCellLinked::add(cell_b);

            let volume_before_a = cell_a.borrow().cell.volume;
            let volume_before_b = cell_b.borrow().cell.volume;
            let energy_before_a = cell_a.borrow().cell.energy_j;
            let energy_before_b = cell_b.borrow().cell.energy_j;

            // Run levelling - this reads from current cells and writes to next cells
            sim.level_cells();

            // PARADIGM: Check the NEXT cells where changes are written
            // The level_cells method reads current state but writes changes to next state
            let volume_after_a = if let Some(next) = &cell_a.borrow().next {
                next.borrow().cell.volume
            } else {
                cell_a.borrow().cell.volume // fallback to current if no next
            };
            let volume_after_b = if let Some(next) = &cell_b.borrow().next {
                next.borrow().cell.volume
            } else {
                cell_b.borrow().cell.volume // fallback to current if no next
            };
            let energy_after_a = if let Some(next) = &cell_a.borrow().next {
                next.borrow().cell.energy_j
            } else {
                cell_a.borrow().cell.energy_j
            };
            let energy_after_b = if let Some(next) = &cell_b.borrow().next {
                next.borrow().cell.energy_j
            } else {
                cell_b.borrow().cell.energy_j
            };

            // Check that volume flowed downhill (from higher to lower) by a significant amount
            if volume_before_a > volume_before_b {
                let volume_decrease = volume_before_a - volume_after_a;
                let volume_increase = volume_after_b - volume_before_b;

                assert!(
                    volume_decrease > 10.0,
                    "Source volume should decrease by at least 10.0 units, decreased by: {}",
                    volume_decrease
                );
                assert!(
                    volume_increase > 10.0,
                    "Dest volume should increase by at least 10.0 units, increased by: {}",
                    volume_increase
                );
                assert!(
                    (volume_decrease - volume_increase).abs() < 1e-6,
                    "Volume transfer should be equal: decrease={}, increase={}",
                    volume_decrease,
                    volume_increase
                );
            }

            // Check that total volume and energy are conserved
            let total_volume_before = volume_before_a + volume_before_b;
            let total_volume_after = volume_after_a + volume_after_b;
            let total_energy_before = energy_before_a + energy_before_b;
            let total_energy_after = energy_after_a + energy_after_b;

            assert!(
                (total_volume_before - total_volume_after).abs() < 1e-6,
                "Volume should be conserved"
            );
            assert!(
                (total_energy_before - total_energy_after).abs() < 1e-6,
                "Energy should be conserved"
            );
        } else {
            panic!("Could not find neighboring cells for test");
        }
    }
}
