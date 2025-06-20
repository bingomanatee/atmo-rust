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
use uuid::Uuid;

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
}

impl SimNext {
    /// Create a new simulation
    pub fn new(planet: Planet, store: RockStore) -> Self {
        let id = Uuid::new_v4().to_string();

        let mut sim = Self {
            cells: HashMap::new(),
            binary_pairs: Vec::new(),
            step: 0,
            planet,
            store,
            png_exporter: None,
            id,
        };

        // Initialize cells automatically
        sim.initialize_cells();
        sim
    }

    /// Initialize cells for the planet
    pub fn initialize_cells(&mut self) {
        println!("🌍 Initializing cells for simulation...");

        let mut asth_cells = Vec::new();

        let args = CellsForPlanetArgs {
            planet: self.planet.clone(),
            on_cell: |cell| {
                asth_cells.push(cell.clone());
                cell
            },
            res: crate::asthenosphere::ASTH_RES,
            joules_per_km3: JOULES_PER_KM3,
            seed: 42,
            anomaly_freq: ANOMALY_SPAWN_CHANCE * 5.0 / 4000.0,
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
    }

    /// Level cells using binary pairs (simple, fast version)
    fn level_cells(&mut self) {
        println!("⚖️  Levelling {} pairs", self.binary_pairs.len());

        // Clone pairs to avoid borrowing issues
        let pairs = self.binary_pairs.clone();
        for pair in &pairs {
            self.level_binary_pair(pair);
        }

        println!("⚖️  Completed levelling");
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
                println!("🌋 New anomaly spawned at cell {}", cell_id);
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

    /// Add initial anomalies to 0.25% of cells
    fn add_initial_anomalies(&mut self) {
        let mut rng = rand::rng();
        let cell_count = self.cells.len();
        let anomaly_count = (cell_count as f64 * 0.0025).ceil() as usize; // 0.25%

        let mut added_count = 0;
        let cell_ids: Vec<CellIndex> = self.cells.keys().cloned().collect();

        for _ in 0..anomaly_count * 10 { // Try up to 10x to account for rejections
            if added_count >= anomaly_count {
                break;
            }

            let random_index = rng.random_range(0..cell_ids.len());
            let cell_id = cell_ids[random_index];

            if let Some(cell) = self.cells.get_mut(&cell_id) {
                if cell.add_anomaly() {
                    added_count += 1;
                }
            }
        }

        println!("🌋 Added {} initial anomalies to {} cells ({:.2}%)",
                 added_count, cell_count, (added_count as f64 / cell_count as f64) * 100.0);
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

        // Run simulation steps
        for step in 1..=steps {
            self.run_step();

            if step % 50 == 0 {
                println!("⏳ Completed step {}/{}", step, steps);
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
