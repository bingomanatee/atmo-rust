use crate::asthenosphere::{ASTH_RES, AsthenosphereCell, CellsForPlanetArgs};
use crate::asthenosphere_linked::AsthenosphereCellLinked;
use crate::constants::{
    ANOMALY_DECAY_RATE, ANOMALY_SPAWN_CHANCE, ANOMALY_VOLUME_AMOUNT, AVG_STARTING_VOLUME_KM_3,
    CELL_JOULES_EQUILIBRIUM, CELL_JOULES_START, JOULES_PER_KM3, LEVEL_AMT,
};
use crate::planet::Planet;
use crate::png_exporter::PngExporter;
use crate::rock_store::RockStore;
use h3o::CellIndex;
use rand::Rng;
use std::cell::RefCell;
use std::collections::HashMap;
use std::fs;
use std::rc::Rc;
use uuid::Uuid;
use crate::asl_leveller::AslMapper;

struct AsthChange {
    energy_j: f64,
    volume: f64,
    energy_j_change: f64,
    volume_change: f64,
    neighbors: Vec<CellIndex>,
}

impl AsthChange {
    pub fn new(source: &AsthenosphereCell) -> AsthChange {
        AsthChange {
            energy_j: source.energy_j,
            volume: source.volume,
            energy_j_change: 0.0,
            volume_change: 0.0,
            neighbors: source.neighbors.clone(),
        }
    }
}
pub struct ASLParams {
    pub planet: Planet,
    pub steps: u64,
    pub store_path: String,
    pub visualize: bool,
    pub vis_freq: u64,
    pub debug: bool,
}

pub struct AsthSimLinked {
    pub(crate) cells: HashMap<CellIndex, Rc<RefCell<AsthenosphereCellLinked>>>,
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

    fn debug_print(&self, message: &str) {
        if self.debug {
            println!("{}", message);
        }
    }

    pub fn run_step(&mut self) {
        self.debug_print("  📏 Leveling cells...");
        self.level_cells();

        self.debug_print("  🔥 Cooling cells...");
        self.cool_cells();

        self.debug_print("  🌀 Processing anomalies...");
        self.process_anomalies();

        self.debug_print("  ⚡ Spawning new anomalies...");
        self.spawn_anomalies();

        self.debug_print("  ⏭️  Advancing to next step...");
        self.advance();

        if self.visualize && self.current_step % self.vis_freq == 0 {
            self.debug_print("  🖼️  Exporting visualization...");
            self.export_visualization();
        }
    }

    fn export_visualization(&self) {
        self.debug_print("    📁 Creating output directory...");
        let output_dir = "vis/asth_sim_linked";
        if let Err(_) = fs::create_dir_all(output_dir) {
            eprintln!(
                "Warning: Could not create visualization directory {}",
                output_dir
            );
            return;
        }

        self.debug_print(&format!("    📊 Collecting {} cells for export...", self.cells.len()));
        let cells_for_export: Vec<(CellIndex, AsthenosphereCell)> = self
            .cells
            .iter()
            .map(|(&cell_id, cell)| {
                let borrowed_cell = cell.borrow();
                let mut simplified_cell = borrowed_cell.cell.clone();
                simplified_cell.energy_j = (simplified_cell.energy_j / 100.0).round() * 100.0;
                (cell_id, simplified_cell)
            })
            .collect();

        self.debug_print("    🎨 Creating PNG exporter (900x450)...");
        let mut exporter = PngExporter::new(900, 450, self.planet.clone());

        self.debug_print("    🖌️  Rendering Voronoi image...");
        let image = exporter.render_voronoi_image_from_cells(&cells_for_export);

        self.debug_print("    💾 Saving image file...");
        let filename = format!("{}/step_{:04}.png", output_dir, self.current_step);
        match image.save(&filename) {
            Ok(_) => self.debug_print(&format!("    ✅ Exported visualization: {}", filename)),
            Err(e) => eprintln!(
                "    ❌ Failed to export visualization for step {}: {}",
                self.current_step, e
            ),
        }
    }

    fn advance(&mut self) {
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

        self.current_step += 1;

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

        for cell_id in cell_ids_with_next {
            if let Some(current_cell) = self.cells.get(&cell_id) {
                let next_cell = {
                    let current_borrowed = current_cell.borrow();
                    current_borrowed.next.as_ref().map(|nc| Rc::clone(nc))
                };

                if let Some(next_cell) = next_cell {
                    next_cell.borrow_mut().cell.step = self.current_step as u32;
                    next_cell.borrow_mut().unlink_all_prev();
                    self.cells.insert(cell_id, next_cell.clone());
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

            if let Some(next_cell) = &cell.borrow().next {
                let mut next = next_cell.borrow_mut();
                next.cell.energy_j *= cool_rate;
            }
        }
    }

    pub fn level_cells(&self) {
        for (id,_) in &self.cells {
           let mapper = AslMapper::new(self, id);
            mapper.level();
        }
    }

    pub fn process_anomalies(&self) {
        for (_id, cell) in &self.cells {
            let current = cell.borrow();

            if (current.cell.anomaly_energy.abs() > 1e-6
                || current.cell.anomaly_volume.abs() > 1e-6)
                && current.next.is_some()
            {
                if let Some(next_cell) = &current.next {
                    let mut next = next_cell.borrow_mut();

                    next.cell.energy_j += current.cell.anomaly_energy;
                    next.cell.volume += current.cell.anomaly_volume;

                    next.cell.anomaly_energy =
                        current.cell.anomaly_energy * (1.0 - ANOMALY_DECAY_RATE);
                    next.cell.anomaly_volume =
                        current.cell.anomaly_volume * (1.0 - ANOMALY_DECAY_RATE);

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

        if rng.random::<f64>() < ANOMALY_SPAWN_CHANCE {
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
        let current = target_cell.borrow();
        let mut rand = rand::rng();
        let sign = if rng.random::<bool>() { 1.0 } else { -1.0 };
        let size = rand.random_range(0.01..0.1);
        let volume = ANOMALY_VOLUME_AMOUNT * sign * size;

        if let Some(next_cell) = &current.next {
            let mut next = next_cell.borrow_mut();

            next.cell.anomaly_volume += volume;
            next.cell.anomaly_energy += volume * JOULES_PER_KM3;
        }

        let neighbors = current.cell.neighbors.clone();
        drop(current);
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

    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }

    pub fn run_simulation(&mut self, steps: u64) {
        self.debug_print(&format!(
            "🌍 Starting AsthSimLinked simulation for {} steps...",
            steps
        ));
        self.debug_print(&format!("📊 Initial state: {} cells loaded", self.cells.len()));

        for step in 1..=steps {
            self.debug_print(&format!("🔄 Running step {} of {}:", step, steps));
            self.run_step();

            if step % 10 == 0 {
                println!("✅ Completed step {} of {}", step, steps);
            }

            if self.visualize && step % 5 == 0 {
                self.debug_print(&format!("🖼️  Exported visualization for step {}", self.current_step));
            }
        }

        if self.visualize {
            self.debug_print(&format!(
                "🎉 Simulation completed! {} PNG files exported to vis/asth_sim_linked/",
                (steps / 5) + 1
            ));
        } else {
            self.debug_print("🎉 Simulation completed!");
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
            steps: 25,
            store_path: String::from("/tmp/test_asth_sim_linked"),
            visualize: true,
            vis_freq: 0,
            debug: false,
        };

        let mut sim = AsthSimLinked::new(config);

        sim.run_simulation(25);

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
            cell_a.borrow_mut().cell.volume = 1000.0;
            cell_a.borrow_mut().cell.energy_j = 2000.0;
            cell_b.borrow_mut().cell.volume = 500.0;
            cell_b.borrow_mut().cell.energy_j = 1000.0;

            let _next_a = AsthenosphereCellLinked::add(cell_a);
            let _next_b = AsthenosphereCellLinked::add(cell_b);

            let volume_before_a = cell_a.borrow().cell.volume;
            let volume_before_b = cell_b.borrow().cell.volume;
            let energy_before_a = cell_a.borrow().cell.energy_j;
            let energy_before_b = cell_b.borrow().cell.energy_j;

            sim.level_cells();

            let volume_after_a = if let Some(next) = &cell_a.borrow().next {
                next.borrow().cell.volume
            } else {
                cell_a.borrow().cell.volume
            };
            let volume_after_b = if let Some(next) = &cell_b.borrow().next {
                next.borrow().cell.volume
            } else {
                cell_b.borrow().cell.volume
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
