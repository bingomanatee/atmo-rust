use crate::asth_cell_next::AsthenosphereCellNext;
use crate::asthenosphere::{AsthenosphereCell, CellsForPlanetArgs, ASTH_RES};
use crate::binary_pair::BinaryPair;
use crate::constants::{JOULES_PER_KM3, VOLCANO_MAX_VOLUME};
use crate::convection::Convection;
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
    pub resolution: h3o::Resolution,
    pub joules_per_km3: f64,
    pub save_to_db: bool,
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
            resolution: crate::asthenosphere::ASTH_RES,
            joules_per_km3: JOULES_PER_KM3,
            save_to_db: true,
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

    pub fn with_resolution(mut self, resolution: h3o::Resolution) -> Self {
        self.resolution = resolution;
        self
    }

    pub fn with_joules_per_km3(mut self, joules_per_km3: f64) -> Self {
        self.joules_per_km3 = joules_per_km3;
        self
    }

    pub fn with_database_saving(mut self, save_to_db: bool) -> Self {
        self.save_to_db = save_to_db;
        self
    }
}

/// Fast asthenosphere simulation without RefCell complexity
pub struct SimNext {
    /// All simulation cells - direct HashMap access, no RefCell overhead
    pub cells: HashMap<CellIndex, AsthenosphereCellNext>,

    /// Binary pairs for levelling (precomputed for efficiency)
    pub binary_pairs: Vec<BinaryPair>,

    /// Global convection system for material addition/subtraction
    pub convection: Convection,

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
    pub save_to_db: bool,
}

impl SimNext {
    /// Create a new simulation with configuration props
    pub fn new(props: SimNextProps) -> Self {
        let id = Uuid::new_v4().to_string();

        let mut sim = Self {
            cells: HashMap::new(),
            binary_pairs: Vec::new(),
            convection: Convection::new(props.seed as u32, 5882), // Will be updated after initialization
            step: 0,
            planet: props.planet,
            store: props.store,
            png_exporter: None,
            id,
            visualize: props.visualize,
            vis_freq: props.vis_freq,
            debug: props.debug,
            save_to_db: props.save_to_db,
        };

        sim.initialize_cells_with_props(
            props.seed,
            props.resolution,
            props.joules_per_km3,
        );
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
    pub fn initialize_cells_with_props(
        &mut self,
        seed: u64,
        resolution: h3o::Resolution,
        joules_per_km3: f64,
    ) {
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
        };
        println!(" -------- volcano max volume = {}", VOLCANO_MAX_VOLUME);

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

        self.generate_binary_pairs();

        // Initialize convection with correct cell count
        self.convection = Convection::new(seed as u32, self.cells.len());

        println!(
            "🌍 Initialized {} cells with {} binary pairs",
            self.cells.len(),
            self.binary_pairs.len()
        );
        println!(
            "🌪️ Convection system initialized: +{:.1}% -{:.1}% (unaffected: {:.1}%, balanced: {})",
            self.convection.addition_fraction * 100.0,
            self.convection.subtraction_fraction * 100.0,
            self.convection.unaffected_percentage() * 100.0,
            self.convection.is_balanced()
        );
    }

    /// Initialize cells for the planet (backward compatibility)
    pub fn initialize_cells(&mut self) {
        self.initialize_cells_with_props(
            42,
            crate::asthenosphere::ASTH_RES,
            JOULES_PER_KM3,
        );
    }

    /// Setup visualization components
    fn setup_visualization(&mut self) {
        self.png_exporter = Some(PngExporter::new(720, 480, self.planet.clone()));

        let output_dir = "vis/sim_next";
        if let Err(e) = fs::create_dir_all(output_dir) {
            eprintln!(
                "Warning: Could not create visualization directory {}: {}",
                output_dir, e
            );
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
    
    fn back_fill(&mut self) {
        let cell_list: Vec<AsthenosphereCellNext> = self.cells.clone().into_values().collect();
        for mut cell in cell_list {
            let mut in_or_near_anomaly = cell.has_any_anomaly();
            
            for  neighbor_id in cell.clone().cell.neighbors {
                let next = self.cells.get(&neighbor_id).unwrap();
                if next.has_any_anomaly() {
                    in_or_near_anomaly = true;
                }
            }
            if !in_or_near_anomaly {
                cell.back_fill();
            } 
        }
    }

    /// Run a single simulation step
    pub fn run_step(&mut self) {
        
        self.back_fill();
        self.debug_print(&format!("🔄 Running step {}", self.step + 1));

        // 1. Level cells (equilibrate volumes/energies)
        self.level_cells();

        // 2. Cool cells
        self.cool_cells();

        // 3. Process existing anomalies and potentially add new ones
        self.process_volcanoes_and_sinkholes();
        self.try_spawn_volcanoes_and_sinkholes();

        // 4. Apply global convection (final step)
        self.apply_convection();

        // 5. Commit next state to current state
        self.commit_step();

        self.step += 1;

        // 5. Export visualization if needed
        if self.visualize && self.step % self.vis_freq == 0 {
            self.debug_print(&format!(
                "🖼️ Exporting visualization for step {}",
                self.step
            ));
            self.export_visualization();
        }
    }

    /// Export visualization of current simulation state
    fn export_visualization(&mut self) {
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
                let mut export_cell = cell.cell.clone();
                export_cell.energy_j = (export_cell.energy_j / 100.0).round() * 100.0;
                (cell_id, export_cell)
            })
            .collect();

        self.debug_print(&format!(
            "  🖌️ Rendering Voronoi image with {} cells...",
            cells_for_export.len()
        ));

        let image = self
            .png_exporter
            .as_mut()
            .unwrap()
            .render_voronoi_image_from_cells(&cells_for_export);

        self.debug_print("  💾 Saving image file...");
        let filename = format!("{}/step_{:04}.png", output_dir, self.step);
        match image.save(&filename) {
            Ok(_) => self.debug_print(&format!("  ✅ Exported visualization: {}", filename)),
            Err(e) => eprintln!(
                "  ❌ Failed to export visualization for step {}: {}",
                self.step, e
            ),
        }
    }

    /// Level cells using binary pairs (simple, fast version)
    fn level_cells(&mut self) {
        self.debug_print(&format!("⚖️ Levelling {} pairs", self.binary_pairs.len()));

        let pairs = self.binary_pairs.clone();
        for pair in &pairs {
            self.level_binary_pair(pair);
        }

        self.debug_print("⚖️ Completed levelling");
    }

    /// Level a single binary pair using conservative volume transfer
    fn level_binary_pair(&mut self, pair: &BinaryPair) {
        let (higher_id, lower_id) = {
            let vol_a = match self.cells.get(&pair.cell_a) {
                Some(cell) => cell.next_cell.volume,
                None => return,
            };
            let vol_b = match self.cells.get(&pair.cell_b) {
                Some(cell) => cell.next_cell.volume,
                None => return,
            };

            if vol_a > vol_b {
                (pair.cell_a, pair.cell_b)
            } else {
                (pair.cell_b, pair.cell_a)
            }
        };

        let (mut higher_cell, mut lower_cell) = (
            self.cells[&higher_id].clone(),
            self.cells[&lower_id].clone(),
        );

        let total_volume = higher_cell.next_cell.volume + lower_cell.next_cell.volume;
        let equilibrium_volume = total_volume / 2.0;
        let excess = higher_cell.next_cell.volume - equilibrium_volume;

        let transfer_amount = excess * 0.25;

        higher_cell
            .next_cell
            .transfer_volume(transfer_amount, &mut lower_cell.next_cell);

        self.cells.insert(higher_id, higher_cell);
        self.cells.insert(lower_id, lower_cell);
    }

    /// Cool all cells with Perlin heating from below
    fn cool_cells(&mut self) {
        for cell in self.cells.values_mut() {
            cell.cool_with_heating(&self.planet, self.step, ASTH_RES);
        }
    }

    fn process_volcanoes_and_sinkholes(&mut self) {
        let cell_ids: Vec<CellIndex> = self.cells.keys().cloned().collect();

        for cell_id in cell_ids {
            if let Some(cell) = self.cells.get_mut(&cell_id) {
                cell.process_volcanoes_and_sinkholes();
            }
        }
    }

    fn try_spawn_volcanoes_and_sinkholes(&mut self) {
        let cell_ids: Vec<CellIndex> = self.cells.keys().cloned().collect();
        let mut volcano_clusters = Vec::new();
        let mut sinkhole_clusters = Vec::new();
        let mut volcano_spawns = Vec::new();
        let mut sinkhole_spawns = Vec::new();

        for cell_id in &cell_ids {
            if let Some(cell) = self.cells.get_mut(cell_id) {
                if let Some(cluster) = cell.try_create_volcano_cluster() {
                    volcano_spawns.push((*cell_id, format!("cluster with {} volcanoes", cluster.len())));
                    volcano_clusters.push(cluster);
                } else if cell.try_add_volcano() {
                    volcano_spawns.push((*cell_id, "single volcano".to_string()));
                }

                if let Some(cluster) = cell.try_create_massive_sink_event() {
                    sinkhole_spawns.push((*cell_id, format!("massive sink affecting {} neighbors", cluster.len())));
                    sinkhole_clusters.push(cluster);
                } else if cell.try_add_sinkhole() {
                    sinkhole_spawns.push((*cell_id, "single sinkhole".to_string()));
                }
            }
        }

        for (cell_id, description) in volcano_spawns {
            self.debug_print(&format!("🌋 Volcano {} spawned at cell {}", description, cell_id));
        }

        for (cell_id, description) in sinkhole_spawns {
            self.debug_print(&format!("🕳️ Sinkhole {} spawned at cell {}", description, cell_id));
        }

        for cluster in volcano_clusters {
            for (neighbor_id, volume) in cluster {
                if let Some(neighbor_cell) = self.cells.get_mut(&neighbor_id) {
                    if !neighbor_cell.has_sinkhole() {
                        neighbor_cell.add_volcano_with_volume(volume);
                    }
                }
            }
        }

        for cluster in sinkhole_clusters {
            for (neighbor_id, volume) in cluster {
                if let Some(neighbor_cell) = self.cells.get_mut(&neighbor_id) {
                    if !neighbor_cell.has_volcano() {
                        neighbor_cell.add_sinkhole_with_volume(volume);
                    }
                }
            }
        }
    }

    /// Apply global convection to all cells
    fn apply_convection(&mut self) {
        self.debug_print("🌪️ Applying global convection");
        
        for cell in self.cells.values_mut() {
            self.convection.apply_convection(
                &mut cell.next_cell,
                &self.planet,
                ASTH_RES
            );
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
        let volumes: Vec<f64> = self
            .cells
            .values()
            .map(|cell| cell.next_cell.volume)
            .collect();
        let energies: Vec<f64> = self
            .cells
            .values()
            .map(|cell| cell.next_cell.energy_j)
            .collect();

        let total_volume: f64 = volumes.iter().sum();
        let total_energy: f64 = energies.iter().sum();
        let cell_count = volumes.len() as f64;

        let avg_volume = total_volume / cell_count;
        let avg_energy = total_energy / cell_count;

        let volume_variance: f64 = volumes
            .iter()
            .map(|v| (v - avg_volume).powi(2))
            .sum::<f64>()
            / cell_count;
        let energy_variance: f64 = energies
            .iter()
            .map(|e| (e - avg_energy).powi(2))
            .sum::<f64>()
            / cell_count;

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

        let initial_stats = self.calculate_statistics();
        println!("📊 Initial Statistics:");
        self.print_statistics(&initial_stats);

        if self.visualize {
            self.export_visualization();
        }

        for _ in 1..=steps {
            self.run_step();

            if self.step % 50 == 0 {
                println!("⏳ Completed step {}/{}", self.step, steps);
            }
        }

        let final_stats = self.calculate_statistics();
        println!("📊 Final Statistics:");
        self.print_statistics(&final_stats);

        println!("📈 Changes:");
        println!(
            "  Volume: {:.2} → {:.2} ({:+.2}%)",
            initial_stats.total_volume,
            final_stats.total_volume,
            ((final_stats.total_volume - initial_stats.total_volume) / initial_stats.total_volume)
                * 100.0
        );
        println!(
            "  Energy: {:.2e} → {:.2e} ({:+.2}%)",
            initial_stats.total_energy,
            final_stats.total_energy,
            ((final_stats.total_energy - initial_stats.total_energy) / initial_stats.total_energy)
                * 100.0
        );
        println!(
            "  Volume Std Dev: {:.2} → {:.2}",
            initial_stats.volume_std_dev, final_stats.volume_std_dev
        );
        println!(
            "  Energy Std Dev: {:.2e} → {:.2e}",
            initial_stats.energy_std_dev, final_stats.energy_std_dev
        );

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
