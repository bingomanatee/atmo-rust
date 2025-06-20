use crate::asthenosphere::{AsthenosphereCell, CellsForPlanetArgs};
use crate::asth_cell_next::AsthenosphereCellNext;
use crate::binary_pair::BinaryPair;
use crate::constants::{
    ANOMALY_SPAWN_CHANCE, JOULES_PER_KM3,
};
use crate::planet::Planet;
use crate::png_exporter::PngExporter;
use crate::rock_store::RockStore;
use h3o::CellIndex;
use rand::Rng;
use std::collections::HashMap;
use std::fs;
use uuid::Uuid;

/// Configuration properties for creating a new SimNext simulation
pub struct SimNextProps {
    pub planet: Planet,
    pub store: RockStore,
    pub visualize: bool,
    pub vis_freq: u32,
    pub debug: bool,
    pub seed: u64,
    pub anomaly_freq: f64,
    pub resolution: h3o::Resolution,
    pub joules_per_km3: f64,
}

impl SimNextProps {
    /// Create default props with required parameters
    pub fn new(planet: Planet, store: RockStore) -> Self {
        Self {
            planet,
            store,
            visualize: true,
            vis_freq: 10,
            debug: true,
            seed: 42,
            anomaly_freq: ANOMALY_SPAWN_CHANCE * 5.0 / 4000.0,
            resolution: crate::asthenosphere::ASTH_RES,
            joules_per_km3: JOULES_PER_KM3,
        }
    }

    /// Builder pattern methods for customization
    pub fn with_visualization(mut self, visualize: bool, vis_freq: u32) -> Self {
        self.visualize = visualize;
        self.vis_freq = vis_freq;
        self
    }

    pub fn with_debug(mut self, debug: bool) -> Self {
        self.debug = debug;
        self
    }

    pub fn with_seed(mut self, seed: u64) -> Self {
        self.seed = seed;
        self
    }

    pub fn with_anomaly_freq(mut self, anomaly_freq: f64) -> Self {
        self.anomaly_freq = anomaly_freq;
        self
    }

    pub fn with_resolution(mut self, resolution: h3o::Resolution) -> Self {
        self.resolution = resolution;
        self
    }

    pub fn with_joules_per_km3(mut self, joules_per_km3: f64) -> Self {
        self.joules_per_km3 = joules_per_km3;
        self
    }
}

/// Fast asthenosphere simulation without RefCell complexity
pub struct SimNext {
    /// All simulation cells - direct HashMap access, no RefCell overhead
    pub cells: HashMap<CellIndex, AsthenosphereCellNext>,

    /// Binary pairs for levelling (precomputed for efficiency)
    pub binary_pairs: Vec<BinaryPair>,

    /// Simulation metadata
    pub step: u32,
    pub planet: Planet,
    pub store: RockStore,
    pub png_exporter: Option<PngExporter>,
    pub id: String,

    /// Visualization settings
    pub visualize: bool,
    pub vis_freq: u32,
    pub debug: bool,
}

impl SimNext {
    /// Create a new simulation with configuration props
    pub fn new(props: SimNextProps) -> Self {
        let id = Uuid::new_v4().to_string();

        let mut sim = Self {
            cells: HashMap::new(),
            binary_pairs: Vec::new(),
            step: 0,
            planet: props.planet,
            store: props.store,
            png_exporter: None,
            id,
            visualize: props.visualize,
            vis_freq: props.vis_freq,
            debug: props.debug,
        };

        // Initialize cells with props configuration
        sim.initialize_cells_with_props(props.seed, props.anomaly_freq, props.resolution, props.joules_per_km3);
        
        // Setup visualization if enabled
        if sim.visualize {
            sim.setup_visualization();
        }
        
        sim
    }


    /// Helper method for debug printing
    fn debug_print(&self, message: &str) {
        if self.debug {
            println!("{}", message);
        }
    }

    /// Initialize cells for the planet using provided configuration
    pub fn initialize_cells_with_props(&mut self, seed: u64, anomaly_freq: f64, resolution: h3o::Resolution, joules_per_km3: f64) {
        println!("🌍 Initializing cells for simulation...");

        let mut asth_cells = Vec::new();

        let args = CellsForPlanetArgs {
            planet: self.planet.clone(),
            on_cell: |cell| {
                asth_cells.push(cell.clone());
                cell
            },
            res: resolution,
            joules_per_km3,
            seed: seed as u32,
            anomaly_freq,
        };

        AsthenosphereCell::initial_cells_for_planet(args);
        println!("🌍 Generated {} asthenosphere cells", asth_cells.len());

        // Convert to simulation cells
        self.cells = asth_cells
            .into_iter()
            .map(|asth_cell| {
                let sim_cell = AsthenosphereCellNext::new(asth_cell);
                (sim_cell.cell_index(), sim_cell)
            })
            .collect();

        // Generate binary pairs
        self.generate_binary_pairs();

        println!("🌍 Initialized {} cells with {} binary pairs",
                 self.cells.len(), self.binary_pairs.len());
    }

    /// Initialize cells for the planet (backward compatibility)
    pub fn initialize_cells(&mut self) {
        self.initialize_cells_with_props(42, ANOMALY_SPAWN_CHANCE * 5.0 / 4000.0, crate::asthenosphere::ASTH_RES, JOULES_PER_KM3);
    }

    /// Setup visualization components
    fn setup_visualization(&mut self) {
        self.png_exporter = Some(PngExporter::new(720, 480, self.planet.clone()));

        // Create visualization directory
        let output_dir = "vis/sim_next";
        if let Err(e) = fs::create_dir_all(output_dir) {
            eprintln!("Warning: Could not create visualization directory {}: {}", output_dir, e);
        } else {
            println!("📁 Created visualization directory: {}", output_dir);
        }
    }

    /// Generate binary pairs for efficient levelling
    fn generate_binary_pairs(&mut self) {
        println!("🔗 Generating binary pairs...");

        // Use a HashMap to store unique pairs
        let mut pairs_map = HashMap::new();

        // Process each cell using its precomputed neighbors
        for (&cell_index, cell) in &self.cells {
            // Use the precomputed neighbors list
            for &neighbor_index in &cell.cell.neighbors {
                let pair = BinaryPair::new(cell_index, neighbor_index);
                pairs_map.insert(pair.to_string_id(), pair);
            }
        }

        // Convert HashMap values to Vector
        self.binary_pairs = pairs_map.into_values().collect();
    }

    /// Run a single simulation step
    pub fn run_step(&mut self) {
        self.debug_print(&format!("🔄 Running step {}", self.step + 1));

        // 1. Level cells (equilibrate volumes/energies)
        self.level_cells();

        // 2. Cool cells
        self.cool_cells();

        // 3. Process existing anomalies and potentially add new ones
        self.process_anomalies();
        self.try_add_new_anomaly();

        // 4. Commit next state to current state
        self.commit_step();

        self.step += 1;

        // 5. Export visualization if needed
        if self.visualize && self.step % self.vis_freq == 0 {
            self.debug_print(&format!("🖼️ Exporting visualization for step {}", self.step));
            self.export_visualization();
        }
    }

    /// Export visualization of current simulation state
    fn export_visualization(&mut self) {
        // Only export if visualization is enabled and exporter exists
        if !self.visualize || self.png_exporter.is_none() {
            return;
        }

        self.debug_print("  📁 Preparing visualization export...");
        let output_dir = "vis/sim_next";

        // Convert simulation cells to AsthenosphereCell for rendering
        let cells_for_export: Vec<(CellIndex, AsthenosphereCell)> = self
            .cells
            .iter()
            .map(|(&cell_id, cell)| {
                // Convert to AsthenosphereCell for export
                let mut export_cell = cell.cell.clone();
                // Round values for better visualization
                export_cell.energy_j = (export_cell.energy_j / 100.0).round() * 100.0;
                (cell_id, export_cell)
            })
            .collect();

        self.debug_print(&format!("  🖌️ Rendering Voronoi image with {} cells...", cells_for_export.len()));

        // Use the PNG exporter to render the image
        let image = self.png_exporter.as_mut().unwrap()
            .render_voronoi_image_from_cells(&cells_for_export);

        // Save the image to file
        self.debug_print("  💾 Saving image file...");
        let filename = format!("{}/step_{:04}.png", output_dir, self.step);
        match image.save(&filename) {
            Ok(_) => self.debug_print(&format!("  ✅ Exported visualization: {}", filename)),
            Err(e) => eprintln!("  ❌ Failed to export visualization for step {}: {}", self.step, e),
        }
    }

    /// Level cells using binary pairs (simple, fast version)
    fn level_cells(&mut self) {
        self.debug_print(&format!("⚖️ Levelling {} pairs", self.binary_pairs.len()));

        // Clone pairs to avoid borrowing issues
        let pairs = self.binary_pairs.clone();
        for pair in &pairs {
            self.level_binary_pair(pair);
        }

        self.debug_print("⚖️ Completed levelling");
    }

    /// Level a single binary pair using conservative volume transfer
    fn level_binary_pair(&mut self, pair: &BinaryPair) {
        // Determine which cell has higher volume first
        let (higher_id, lower_id) = {
            let vol_a = match self.cells.get(&pair.cell_a) {
                Some(cell) => cell.next_cell.volume,
                None => return,
            };
            let vol_b = match self.cells.get(&pair.cell_b) {
                Some(cell) => cell.next_cell.volume,
                None => return,
            };

            // Early termination if volume difference is too small
            if (vol_a - vol_b).abs() < 1.0 {
                return;
            }

            if vol_a > vol_b {
                (pair.cell_a, pair.cell_b)
            } else {
                (pair.cell_b, pair.cell_a)
            }
        };

        // Clone the pair and perform transfer operations
        let (mut higher_cell, mut lower_cell) = (self.cells[&higher_id].clone(), self.cells[&lower_id].clone());

        // Calculate equilibrium and transfer amount
        let total_volume = higher_cell.next_cell.volume + lower_cell.next_cell.volume;
        let equilibrium_volume = total_volume / 2.0;
        let excess = higher_cell.next_cell.volume - equilibrium_volume;

        // Scale transfer for multiple neighbors (assume ~6 neighbors per cell)
        let transfer_amount = excess * 0.16; // 1.0 / 6.0 ≈ 0.16

        if transfer_amount < 0.1 {
            return;
        }

        // Perform downstream transfer from higher to lower on clones
        higher_cell.next_cell.transfer_volume(transfer_amount, &mut lower_cell.next_cell);

        // Update the HashMap with the modified cells
        self.cells.insert(higher_id, higher_cell);
        self.cells.insert(lower_id, lower_cell);
    }



    /// Cool all cells
    fn cool_cells(&mut self) {
        for cell in self.cells.values_mut() {
            cell.cool();
        }
    }

    /// Process existing anomalies in all cells
    fn process_anomalies(&mut self) {
        // Clone cell IDs to avoid borrowing conflicts
        let cell_ids: Vec<CellIndex> = self.cells.keys().cloned().collect();

        for cell_id in cell_ids {
            // Check if this cell has a strong anomaly that should propagate
            let (should_propagate, neighbor_ids, volume_per_neighbor) = {
                let cell = &self.cells[&cell_id];
                if cell.next_cell.anomaly_volume.abs() > 10.0 && !cell.cell.neighbors.is_empty() {
                    let volume_per_neighbor = cell.next_cell.anomaly_volume / cell.cell.neighbors.len() as f64;
                    (true, cell.cell.neighbors.clone(), volume_per_neighbor)
                } else {
                    (false, vec![], 0.0)
                }
            };

            // Process the anomaly for this cell
            if let Some(cell) = self.cells.get_mut(&cell_id) {
                cell.process_anomaly();
            }

            // Propagate to neighbors immediately
            if should_propagate {
                for neighbor_id in neighbor_ids {
                    if let Some(neighbor_cell) = self.cells.get_mut(&neighbor_id) {
                        neighbor_cell.volume_from_anomaly(volume_per_neighbor);
                    }
                }
            }
        }
    }

    /// Try to add a new anomaly based on spawn chance
    fn try_add_new_anomaly(&mut self) {
        let mut rng = rand::rng();

        // Check spawn chance
        if rng.random::<f64>() > ANOMALY_SPAWN_CHANCE {
            return; // No anomaly this step
        }

        // Select a random cell
        let cell_ids: Vec<CellIndex> = self.cells.keys().cloned().collect();
        if cell_ids.is_empty() {
            return;
        }

        let random_index = rng.random_range(0..cell_ids.len());
        let cell_id = cell_ids[random_index];

        // Try to add anomaly (will fail if cell already has one)
        if let Some(cell) = self.cells.get_mut(&cell_id) {
            if cell.add_anomaly() {
                self.debug_print(&format!("🌋 New anomaly spawned at cell {}", cell_id));
            }
        }
    }

    /// Commit next state to current state
    fn commit_step(&mut self) {
        for cell in self.cells.values_mut() {
            cell.commit_step();
        }
    }

    /// Get cell count
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }

    /// Get current step
    pub fn current_step(&self) -> u32 {
        self.step
    }

    /// Calculate simulation statistics
    fn calculate_statistics(&self) -> SimulationStats {
        let volumes: Vec<f64> = self.cells.values().map(|cell| cell.next_cell.volume).collect();
        let energies: Vec<f64> = self.cells.values().map(|cell| cell.next_cell.energy_j).collect();

        let total_volume: f64 = volumes.iter().sum();
        let total_energy: f64 = energies.iter().sum();
        let cell_count = volumes.len() as f64;

        let avg_volume = total_volume / cell_count;
        let avg_energy = total_energy / cell_count;

        let volume_variance: f64 = volumes.iter()
            .map(|v| (v - avg_volume).powi(2))
            .sum::<f64>() / cell_count;
        let energy_variance: f64 = energies.iter()
            .map(|e| (e - avg_energy).powi(2))
            .sum::<f64>() / cell_count;

        SimulationStats {
            total_volume,
            total_energy,
            avg_volume,
            avg_energy,
            volume_std_dev: volume_variance.sqrt(),
            energy_std_dev: energy_variance.sqrt(),
            cell_count: cell_count as usize,
        }
    }

    /// Print statistics in a formatted way
    fn print_statistics(&self, stats: &SimulationStats) {
        println!("  Cells: {}", stats.cell_count);
        println!("  Total Volume: {:.2} km³", stats.total_volume);
        println!("  Total Energy: {:.2e} J", stats.total_energy);
        println!("  Avg Volume/Cell: {:.2} km³", stats.avg_volume);
        println!("  Avg Energy/Cell: {:.2e} J", stats.avg_energy);
        println!("  Volume Std Dev: {:.2} km³", stats.volume_std_dev);
        println!("  Energy Std Dev: {:.2e} J", stats.energy_std_dev);
    }

    /// Run complete simulation with statistics
    pub fn run_simulation(&mut self, steps: u32) {
        println!("🚀 Starting simulation for {} steps", steps);

        // Calculate initial statistics
        let initial_stats = self.calculate_statistics();
        println!("📊 Initial Statistics:");
        self.print_statistics(&initial_stats);

        // Export initial state if visualization is enabled
        if self.visualize {
            self.export_visualization();
        }

        // Run simulation steps
        for _ in 1..=steps {
            self.run_step();

            if self.step % 50 == 0 {
                println!("⏳ Completed step {}/{}", self.step, steps);
            }
        }

        // Calculate final statistics
        let final_stats = self.calculate_statistics();
        println!("📊 Final Statistics:");
        self.print_statistics(&final_stats);

        // Print comparison
        println!("📈 Changes:");
        println!("  Volume: {:.2} → {:.2} ({:+.2}%)",
                 initial_stats.total_volume, final_stats.total_volume,
                 ((final_stats.total_volume - initial_stats.total_volume) / initial_stats.total_volume) * 100.0);
        println!("  Energy: {:.2e} → {:.2e} ({:+.2}%)",
                 initial_stats.total_energy, final_stats.total_energy,
                 ((final_stats.total_energy - initial_stats.total_energy) / initial_stats.total_energy) * 100.0);
        println!("  Volume Std Dev: {:.2} → {:.2}",
                 initial_stats.volume_std_dev, final_stats.volume_std_dev);
        println!("  Energy Std Dev: {:.2e} → {:.2e}",
                 initial_stats.energy_std_dev, final_stats.energy_std_dev);

        if self.visualize {
            println!("🖼️ Visualization images saved to vis/sim_next/");
        }
    }
}

#[derive(Debug, Clone)]
struct SimulationStats {
    total_volume: f64,
    total_energy: f64,
    avg_volume: f64,
    avg_energy: f64,
    volume_std_dev: f64,
    energy_std_dev: f64,
    cell_count: usize,
}
