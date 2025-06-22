use crate::asth_cell_next::AsthenosphereCellNext;
use crate::asthenosphere::{AsthenosphereCell, CellsForPlanetArgs, ASTH_RES};
use crate::binary_pair::BinaryPair;
use crate::constants::{JOULES_PER_KM3, VOLCANO_MAX_VOLUME, LAYER_COUNT, VERTICAL_ENERGY_MIXING, USE_BATCH_TRANSFERS};
use crate::convection::{Convection, ConvectionSystem};
use crate::planet::Planet;
use crate::png_exporter::PngExporter;
use crate::rock_store::RockStore;
use crate::batch_transfer::{BatchTransferTracker, CellLayerKey, StepTimer, time_operation, TransferAccountingSystem};
use h3o::CellIndex;
use noise::NoiseFn;
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

/// Fast asthenosphere simulation without RefCell complexity, with array-based layers
pub struct SimNext {
    /// Simulation cells with array-based volume and energy per layer
    pub cells: HashMap<CellIndex, AsthenosphereCellNext>,

    /// Binary pairs for levelling (precomputed for efficiency)
    pub binary_pairs: Vec<BinaryPair>,

    /// Dynamic convection system with template interpolation (applied only to bottom layer)
    pub convection_system: ConvectionSystem,

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

    /// Performance timing
    pub step_timer: StepTimer,
}

impl SimNext {
    /// Create a new simulation with configuration props
    pub fn new(props: SimNextProps) -> Self {
        let id = Uuid::new_v4().to_string();

        let mut sim = Self {
            cells: HashMap::new(),
            binary_pairs: Vec::new(),
            convection_system: ConvectionSystem::new(props.seed as u32, 5882), // Will be updated after initialization
            step: 0,
            planet: props.planet,
            store: props.store,
            png_exporter: None,
            id,
            visualize: props.visualize,
            vis_freq: props.vis_freq,
            debug: props.debug,
            save_to_db: props.save_to_db,
            step_timer: StepTimer::new(),
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

    /// Initialize array-based cells for the planet using provided configuration
    pub fn initialize_cells_with_props(
        &mut self,
        seed: u64,
        resolution: h3o::Resolution,
        joules_per_km3: f64,
    ) {
        println!("🌍 Initializing array-based cells for simulation...");

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

        // Use the regular single-cell initialization method
        AsthenosphereCell::initial_cells_for_planet(args);
        println!("🌍 Generated {} asthenosphere cells with {} layers each", asth_cells.len(), LAYER_COUNT);

        // Create simulation cells with array-based layer data
        for asth_cell in asth_cells {
            let sim_cell = AsthenosphereCellNext::new(asth_cell.clone());
            self.cells.insert(sim_cell.cell_index(), sim_cell);
        }

        self.generate_binary_pairs();

        // Initialize convection system with correct cell count
        let total_cells = self.cells.len();
        self.convection_system = ConvectionSystem::new(seed as u32, total_cells);

        println!(
            "🌍 Initialized {} cells with {} layers each ({} total layer entries) with {} binary pairs",
            total_cells,
            LAYER_COUNT,
            total_cells * LAYER_COUNT,
            self.binary_pairs.len()
        );
        let current_template = self.convection_system.get_interpolated_template();
        println!(
            "🌪️ Convection system initialized: +{:.1}% -{:.1}% (unaffected: {:.1}%, balanced: {}, cycle: {}/{})",
            current_template.addition_fraction * 100.0,
            current_template.subtraction_fraction * 100.0,
            current_template.unaffected_percentage() * 100.0,
            current_template.is_balanced(),
            self.convection_system.step_in_cycle,
            self.convection_system.cycle_lifespan
        );
        
        // Apply jump-start convection for dramatic initial patterns (only to upper layer)
        self.jumpstart_convection();
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
        // Apply back_fill to all cells
        let cell_list: Vec<AsthenosphereCellNext> = self.cells.clone().into_values().collect();
        for mut cell in cell_list {
            let mut in_or_near_anomaly = cell.has_any_anomaly();
            
            for neighbor_id in cell.clone().cell.neighbors {
                if let Some(next) = self.cells.get(&neighbor_id) {
                    if next.has_any_anomaly() {
                        in_or_near_anomaly = true;
                        break;
                    }
                }
            }
            if !in_or_near_anomaly {
                cell.back_fill();
                self.cells.insert(cell.cell_index(), cell);
            }
        }
    }

    /// Run a single simulation step
    pub fn run_step(&mut self) {
        let (_, total_time) = time_operation(|| {
            // Set current step in accounting system
            let accounting = TransferAccountingSystem::instance();
            accounting.set_step(self.step + 1);
            
            self.back_fill();
            self.debug_print(&format!("🔄 Running step {}", self.step + 1));

            // 1. Level cells (equilibrate volumes/energies) across all layers
            let (_, leveling_time) = time_operation(|| self.level_cells());
            self.step_timer.leveling_time_ms = leveling_time;

            // 2. Vertical energy mixing between layers
            let (_, mixing_time) = time_operation(|| self.apply_vertical_energy_mixing());
            self.step_timer.mixing_time_ms = mixing_time;

            // 3. Process existing anomalies and potentially add new ones (all layers)
            self.process_volcanoes_and_sinkholes();
            self.try_spawn_volcanoes_and_sinkholes();

            // 4. Apply global convection (only to bottom layer)
            let (_, convection_time) = time_operation(|| self.apply_convection());
            self.step_timer.convection_time_ms = convection_time;

            // 5. Commit batch transfers if using batch system
            let (_, transfer_time) = time_operation(|| {
                if USE_BATCH_TRANSFERS {
                    self.commit_batch_transfers();
                }
            });
            self.step_timer.transfer_commit_time_ms = transfer_time;

            // 6. Apply quantum transfers from accounting system
            self.apply_quantum_transfers();

            // 7. Commit next state to current state
            self.commit_step();

            // 8. Cool cells after all transfers are reconciled
            self.cool_cells();
        });
        
        self.step_timer.total_step_time_ms = total_time;
        
        // Print timing stats every 10 steps when debug is enabled
        if self.debug && self.step % 10 == 0 {
            self.step_timer.print_stats(self.step + 1);
        }

        self.step += 1;

        // 8. Export visualization if needed
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

        // Convert simulation cells to AsthenosphereCell for rendering (use upper layer from arrays for visualization)
        let upper_layer_index = LAYER_COUNT - 1;
        let cells_for_export: Vec<(CellIndex, AsthenosphereCell)> = self
            .cells
            .iter()
            .map(|(&cell_id, cell)| {
                let export_cell = cell.cell.clone();
                // The PNG exporter will use the compatibility methods volume() and energy_j()
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

    /// Level a single binary pair using conservative volume transfer (across all layer arrays)
    fn level_binary_pair(&mut self, pair: &BinaryPair) {
        if USE_BATCH_TRANSFERS {
            self.level_binary_pair_batch(pair);
        } else {
            self.level_binary_pair_direct(pair);
        }
    }

    /// Direct transfer version (original implementation)
    fn level_binary_pair_direct(&mut self, pair: &BinaryPair) {
        // Get both cells
        let (mut cell_a, mut cell_b) = match (
            self.cells.get(&pair.cell_a).cloned(),
            self.cells.get(&pair.cell_b).cloned(),
        ) {
            (Some(a), Some(b)) => (a, b),
            _ => return, // Skip if either cell is missing
        };

        // Apply leveling to each layer in the arrays
        for layer_index in 0..LAYER_COUNT {
            let vol_a = cell_a.next_cell.volume_layers[layer_index];
            let vol_b = cell_b.next_cell.volume_layers[layer_index];

            let (higher_vol, lower_vol, is_a_higher) = if vol_a > vol_b {
                (vol_a, vol_b, true)
            } else {
                (vol_b, vol_a, false)
            };

            let total_volume = higher_vol + lower_vol;
            let equilibrium_volume = total_volume / 2.0;
            let excess = higher_vol - equilibrium_volume;
            let transfer_amount = excess * 0.25;

            if transfer_amount > 0.0 {
                // Calculate proportional energy transfer
                let energy_density = if higher_vol > 0.0 {
                    if is_a_higher {
                        cell_a.next_cell.energy_layers[layer_index] / vol_a
                    } else {
                        cell_b.next_cell.energy_layers[layer_index] / vol_b
                    }
                } else {
                    0.0
                };
                let energy_transfer = energy_density * transfer_amount;

                // Apply the transfer and record in accounting system
                let accounting = TransferAccountingSystem::instance();
                if is_a_higher {
                    cell_a.next_cell.volume_layers[layer_index] -= transfer_amount;
                    cell_a.next_cell.energy_layers[layer_index] -= energy_transfer;
                    cell_b.next_cell.volume_layers[layer_index] += transfer_amount;
                    cell_b.next_cell.energy_layers[layer_index] += energy_transfer;
                    
                    // Record transfer in accounting system
                    let from_key = CellLayerKey::new(pair.cell_a, layer_index);
                    let to_key = CellLayerKey::new(pair.cell_b, layer_index);
                    accounting.record_transfer(from_key, to_key, transfer_amount);
                } else {
                    cell_b.next_cell.volume_layers[layer_index] -= transfer_amount;
                    cell_b.next_cell.energy_layers[layer_index] -= energy_transfer;
                    cell_a.next_cell.volume_layers[layer_index] += transfer_amount;
                    cell_a.next_cell.energy_layers[layer_index] += energy_transfer;
                    
                    // Record transfer in accounting system
                    let from_key = CellLayerKey::new(pair.cell_b, layer_index);
                    let to_key = CellLayerKey::new(pair.cell_a, layer_index);
                    accounting.record_transfer(from_key, to_key, transfer_amount);
                }

                // Ensure non-negative values
                cell_a.next_cell.volume_layers[layer_index] = cell_a.next_cell.volume_layers[layer_index].max(0.0);
                cell_a.next_cell.energy_layers[layer_index] = cell_a.next_cell.energy_layers[layer_index].max(0.0);
                cell_b.next_cell.volume_layers[layer_index] = cell_b.next_cell.volume_layers[layer_index].max(0.0);
                cell_b.next_cell.energy_layers[layer_index] = cell_b.next_cell.energy_layers[layer_index].max(0.0);
            }
        }

        // Update the cells back in the HashMap
        self.cells.insert(pair.cell_a, cell_a);
        self.cells.insert(pair.cell_b, cell_b);
    }

    /// Batch transfer version (registers transfers instead of applying directly)
    fn level_binary_pair_batch(&mut self, pair: &BinaryPair) {
        // Get both cells for calculation only (no cloning needed)
        let (cell_a, cell_b) = match (
            self.cells.get(&pair.cell_a),
            self.cells.get(&pair.cell_b),
        ) {
            (Some(a), Some(b)) => (a, b),
            _ => return, // Skip if either cell is missing
        };

        let tracker = BatchTransferTracker::instance();

        // Calculate transfers for each layer and register them
        for layer_index in 0..LAYER_COUNT {
            let vol_a = cell_a.next_cell.volume_layers[layer_index];
            let vol_b = cell_b.next_cell.volume_layers[layer_index];

            let (higher_vol, lower_vol, is_a_higher) = if vol_a > vol_b {
                (vol_a, vol_b, true)
            } else {
                (vol_b, vol_a, false)
            };

            let total_volume = higher_vol + lower_vol;
            let equilibrium_volume = total_volume / 2.0;
            let excess = higher_vol - equilibrium_volume;
            let transfer_amount = excess * 0.25;

            if transfer_amount > 0.0 {
                // Calculate proportional energy transfer
                let energy_density = if higher_vol > 0.0 {
                    if is_a_higher {
                        cell_a.next_cell.energy_layers[layer_index] / vol_a
                    } else {
                        cell_b.next_cell.energy_layers[layer_index] / vol_b
                    }
                } else {
                    0.0
                };
                let energy_transfer = energy_density * transfer_amount;

                // Register transfer instead of applying directly
                if is_a_higher {
                    let from_key = CellLayerKey::new(pair.cell_a, layer_index);
                    let to_key = CellLayerKey::new(pair.cell_b, layer_index);
                    tracker.register_transfer(from_key, to_key, transfer_amount);
                } else {
                    let from_key = CellLayerKey::new(pair.cell_b, layer_index);
                    let to_key = CellLayerKey::new(pair.cell_a, layer_index);
                    tracker.register_transfer(from_key, to_key, transfer_amount);
                }
            }
        }
    }

    /// Cool all cells with Perlin heating from below (across all layer arrays)
    fn cool_cells(&mut self) {
        for cell in self.cells.values_mut() {
            cell.cool_with_heating(&self.planet, self.step, ASTH_RES);
        }
    }

    fn process_volcanoes_and_sinkholes(&mut self) {
        // Process volcanoes and sinkholes for all cells
        let cell_ids: Vec<CellIndex> = self.cells.keys().cloned().collect();

        for cell_id in cell_ids {
            if let Some(cell) = self.cells.get_mut(&cell_id) {
                cell.process_volcanoes_and_sinkholes();
            }
        }
    }

    fn try_spawn_volcanoes_and_sinkholes(&mut self) {
        // Apply volcano/sinkhole spawning to all cells (surface activity)
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

    /// Apply dynamic convection to bottom layer of arrays only and advance cycle
    fn apply_convection(&mut self) {
        let progress = self.convection_system.cycle_progress();
        self.debug_print(&format!("🌪️ Applying dynamic convection to bottom layer arrays (cycle {:.1}%)", progress * 100.0));
        
        // Get needed data from convection system first
        let template = self.convection_system.get_interpolated_template();
        let t = self.convection_system.step_in_cycle as f64 / self.convection_system.cycle_lifespan as f64;
        let planet = self.planet.clone();
        
        // Apply convection only to the bottom layer (layer 0) of each cell's arrays
        for cell in self.cells.values_mut() {
            Self::apply_convection_to_layer_static(&mut cell.next_cell, 0, &template, t, &self.convection_system, &planet);
        }
        
        // Advance convection system to next step
        self.convection_system.advance_step();
        
        // Check if we started a new cycle and log it
        if self.convection_system.step_in_cycle == 0 {
            let new_template = self.convection_system.get_interpolated_template();
            self.debug_print(&format!(
                "🌪️ New convection cycle started: +{:.1}% -{:.1}% (scale: {:.1}, lifespan: {})",
                new_template.addition_fraction * 100.0,
                new_template.subtraction_fraction * 100.0,
                new_template.noise_scale,
                self.convection_system.cycle_lifespan
            ));
        }
    }
    
    /// Apply jump-start convection with 10x amplification for dramatic initial patterns (bottom layer arrays only)
    fn jumpstart_convection(&mut self) {
        println!("🚀 Applying jump-start convection to bottom layer arrays (10x amplification)...");
        
        // Get needed data from convection system first
        let mut template = self.convection_system.get_interpolated_template();
        template.per_cell_addition *= 10.0;
        template.per_cell_subtraction *= 10.0;
        let t = self.convection_system.step_in_cycle as f64 / self.convection_system.cycle_lifespan as f64;
        let planet = self.planet.clone();
        
        // Apply jump-start convection only to the bottom layer (layer 0) of each cell's arrays
        for cell in self.cells.values_mut() {
            Self::apply_convection_to_layer_static(&mut cell.next_cell, 0, &template, t, &self.convection_system, &planet);
        }
        
        // Commit the jump-start changes immediately for all cells
        for cell in self.cells.values_mut() {
            cell.commit_step();
        }
        
        println!("✅ Jump-start convection applied successfully to bottom layer arrays");
    }

    /// Apply convection to a specific layer in a cell's arrays (static method to avoid borrow issues)
    fn apply_convection_to_layer_static(
        cell: &mut AsthenosphereCell,
        layer_index: usize,
        template: &crate::convection::ConvectionTemplate,
        t: f64,
        convection_system: &crate::convection::ConvectionSystem,
        planet: &crate::planet::Planet,
    ) {
        if layer_index >= LAYER_COUNT {
            return;
        }
        
        // Get cell location for perlin noise sampling
        let gc = crate::geoconverter::GeoCellConverter::new(planet.radius_km as f64, ASTH_RES);
        let location = gc.cell_to_vec3(cell.id);
        
        // Use template's noise scale
        let scaled_location = location.normalize() * template.noise_scale;
        
        // Interpolate noise values between current and next perlin
        let current_noise = convection_system.current_perlin.get(scaled_location.to_array().map(|n| n as f64));
        let next_noise = convection_system.next_perlin.get(scaled_location.to_array().map(|n| n as f64));
        let noise_value = current_noise * (1.0 - t) + next_noise * t;
        
        // Calculate thresholds for top and bottom percentiles
        let addition_threshold = 1.0 - (template.addition_fraction * 2.0);
        let subtraction_threshold = -1.0 + (template.subtraction_fraction * 2.0);
        
        if noise_value > addition_threshold {
            // Top percentile: Addition - add material at start temperature
            let volume_to_add = template.per_cell_addition;
            let energy_to_add = volume_to_add * (crate::constants::CELL_JOULES_START / crate::constants::AVG_STARTING_VOLUME_KM_3);
            
            cell.volume_layers[layer_index] += volume_to_add;
            cell.energy_layers[layer_index] += energy_to_add;
        } else if noise_value < subtraction_threshold {
            // Bottom percentile: Subtraction - remove material proportionally
            let volume_to_remove = template.per_cell_subtraction;
            let volume_to_remove = volume_to_remove.min(cell.volume_layers[layer_index] * 0.9); // Don't remove more than 90%
            
            if volume_to_remove > 0.0 && cell.volume_layers[layer_index] > 0.0 {
                // Remove energy proportionally to volume
                let energy_ratio = volume_to_remove / cell.volume_layers[layer_index];
                let energy_to_remove = cell.energy_layers[layer_index] * energy_ratio;
                
                cell.volume_layers[layer_index] -= volume_to_remove;
                cell.energy_layers[layer_index] -= energy_to_remove;
                
                // Ensure no negative values
                cell.volume_layers[layer_index] = cell.volume_layers[layer_index].max(0.0);
                cell.energy_layers[layer_index] = cell.energy_layers[layer_index].max(0.0);
            }
        }
        // Middle cells remain unaffected
    }

    /// Apply vertical mixing between layers (both volume and energy transfer) using arrays
    fn apply_vertical_energy_mixing(&mut self) {
        if LAYER_COUNT < 2 {
            return; // No mixing with single layer
        }

        self.debug_print(&format!("🔄 Applying vertical volume and energy mixing ({:.1}% exchange)", VERTICAL_ENERGY_MIXING * 100.0));

        if USE_BATCH_TRANSFERS {
            self.apply_vertical_energy_mixing_batch();
        } else {
            self.apply_vertical_energy_mixing_direct();
        }
    }

    /// Direct transfer version for vertical mixing
    fn apply_vertical_energy_mixing_direct(&mut self) {
        // Get all cell indices
        let cell_indices: Vec<CellIndex> = self.cells.keys().cloned().collect();

        for cell_index in cell_indices {
            if let Some(mut cell) = self.cells.get(&cell_index).cloned() {
                // Mix both volume and energy between adjacent layers within the same cell
                for layer_index in 0..(LAYER_COUNT - 1) {
                    let lower_layer = layer_index;
                    let upper_layer = layer_index + 1;

                    // Calculate volume and energy exchange amounts
                    let lower_volume_to_transfer = cell.next_cell.volume_layers[lower_layer] * VERTICAL_ENERGY_MIXING;
                    let upper_volume_to_transfer = cell.next_cell.volume_layers[upper_layer] * VERTICAL_ENERGY_MIXING;
                    
                    let lower_energy_to_transfer = cell.next_cell.energy_layers[lower_layer] * VERTICAL_ENERGY_MIXING;
                    let upper_energy_to_transfer = cell.next_cell.energy_layers[upper_layer] * VERTICAL_ENERGY_MIXING;

                    // Apply volume exchange (creates imbalances that leveling will fix)
                    cell.next_cell.volume_layers[lower_layer] -= lower_volume_to_transfer;
                    cell.next_cell.volume_layers[lower_layer] += upper_volume_to_transfer;

                    cell.next_cell.volume_layers[upper_layer] -= upper_volume_to_transfer;
                    cell.next_cell.volume_layers[upper_layer] += lower_volume_to_transfer;
                    
                    // Record transfers in accounting system
                    let accounting = TransferAccountingSystem::instance();
                    let lower_key = CellLayerKey::new(cell_index, lower_layer);
                    let upper_key = CellLayerKey::new(cell_index, upper_layer);
                    accounting.record_transfer(lower_key.clone(), upper_key.clone(), lower_volume_to_transfer);
                    accounting.record_transfer(upper_key, lower_key, upper_volume_to_transfer);

                    // Apply energy exchange
                    cell.next_cell.energy_layers[lower_layer] -= lower_energy_to_transfer;
                    cell.next_cell.energy_layers[lower_layer] += upper_energy_to_transfer;

                    cell.next_cell.energy_layers[upper_layer] -= upper_energy_to_transfer;
                    cell.next_cell.energy_layers[upper_layer] += lower_energy_to_transfer;

                    // Ensure non-negative values
                    cell.next_cell.volume_layers[lower_layer] = cell.next_cell.volume_layers[lower_layer].max(0.0);
                    cell.next_cell.energy_layers[lower_layer] = cell.next_cell.energy_layers[lower_layer].max(0.0);
                    cell.next_cell.volume_layers[upper_layer] = cell.next_cell.volume_layers[upper_layer].max(0.0);
                    cell.next_cell.energy_layers[upper_layer] = cell.next_cell.energy_layers[upper_layer].max(0.0);
                }

                // Update the cell back in the HashMap
                self.cells.insert(cell_index, cell);
            }
        }
    }

    /// Batch transfer version for vertical mixing
    fn apply_vertical_energy_mixing_batch(&mut self) {
        let tracker = BatchTransferTracker::instance();

        // Register transfers for all cells without cloning
        for (&cell_index, cell) in &self.cells {
            // Mix both volume and energy between adjacent layers within the same cell
            for layer_index in 0..(LAYER_COUNT - 1) {
                let lower_layer = layer_index;
                let upper_layer = layer_index + 1;

                // Calculate volume and energy exchange amounts
                let lower_volume_to_transfer = cell.next_cell.volume_layers[lower_layer] * VERTICAL_ENERGY_MIXING;
                let upper_volume_to_transfer = cell.next_cell.volume_layers[upper_layer] * VERTICAL_ENERGY_MIXING;
                
                let lower_energy_to_transfer = cell.next_cell.energy_layers[lower_layer] * VERTICAL_ENERGY_MIXING;
                let upper_energy_to_transfer = cell.next_cell.energy_layers[upper_layer] * VERTICAL_ENERGY_MIXING;

                // Register transfers instead of applying directly
                let lower_key = CellLayerKey::new(cell_index, lower_layer);
                let upper_key = CellLayerKey::new(cell_index, upper_layer);

                // Register bidirectional transfer
                tracker.register_transfer(
                    lower_key.clone(),
                    upper_key.clone(),
                    lower_volume_to_transfer
                );
                tracker.register_transfer(
                    upper_key,
                    lower_key,
                    upper_volume_to_transfer
                );
            }
        }
    }

    /// Commit next state to current state (across all cells)
    fn commit_step(&mut self) {
        for cell in self.cells.values_mut() {
            cell.commit_step();
        }
    }

    /// Commit all batch transfers to cells
    fn commit_batch_transfers(&mut self) {
        let tracker = BatchTransferTracker::instance();
        let transfers = tracker.consume_transfers();
        
        if transfers.is_empty() {
            return;
        }

        self.debug_print(&format!("🔄 Committing {} batch transfer operations", 
            transfers.values().map(|to_map| to_map.len()).sum::<usize>()));

        let accounting = TransferAccountingSystem::instance();

        // Process all transfers
        for (from_key, to_map) in transfers {
            // Sum all outgoing transfers from this source
            let total_outgoing_volume: f64 = to_map.values().sum();

            // Apply outgoing transfers (subtract from source)
            if let Some(from_cell) = self.cells.get_mut(&from_key.cell_id) {
                if from_key.layer_index < LAYER_COUNT {
                    from_cell.next_cell.volume_layers[from_key.layer_index] -= total_outgoing_volume;
                    
                    // Ensure non-negative values
                    from_cell.next_cell.volume_layers[from_key.layer_index] = 
                        from_cell.next_cell.volume_layers[from_key.layer_index].max(0.0);
                }
            }

            // Apply individual incoming transfers (add to destinations)
            for (to_key, volume) in to_map {
                // Record transfer in accounting system
                accounting.record_transfer(from_key.clone(), to_key.clone(), volume);
                
                if let Some(to_cell) = self.cells.get_mut(&to_key.cell_id) {
                    if to_key.layer_index < LAYER_COUNT {
                        to_cell.next_cell.volume_layers[to_key.layer_index] += volume;
                        // Energy will be computed during reconciliation
                        
                        // Ensure non-negative values
                        to_cell.next_cell.volume_layers[to_key.layer_index] = 
                            to_cell.next_cell.volume_layers[to_key.layer_index].max(0.0);
                    }
                }
            }
        }
    }

    /// Apply quantum transfers from accounting system - consolidate all transfers into single operations
    fn apply_quantum_transfers(&mut self) {
        let accounting = TransferAccountingSystem::instance();
        let current_step = self.step + 1; // We set the accounting step to step + 1 earlier
        let flattened_transfers = accounting.flatten_transfers_for_step(current_step);
        
        if flattened_transfers.is_empty() {
            return;
        }

        self.debug_print(&format!("🔬 Applying quantum transfers for {} source cells", flattened_transfers.len()));

        // Process each source cell+layer - only handle volume transfers
        for (from_key, to_map) in &flattened_transfers {
            if let Some(source_cell) = self.cells.get_mut(&from_key.cell_id) {
                if from_key.layer_index < LAYER_COUNT {
                    // Calculate total outgoing volume
                    let total_outgoing_volume: f64 = to_map.values().sum();
                    
                    if total_outgoing_volume <= 0.0 {
                        continue;
                    }

                    // Get current source volume
                    let source_volume = source_cell.next_cell.volume_layers[from_key.layer_index];
                    
                    // Remove volume from source (ensure we don't exceed available)
                    let actual_volume_transfer = total_outgoing_volume.min(source_volume * 0.95); // Don't transfer more than 95%
                    source_cell.next_cell.volume_layers[from_key.layer_index] -= actual_volume_transfer;
                    
                    // Ensure non-negative values
                    source_cell.next_cell.volume_layers[from_key.layer_index] = 
                        source_cell.next_cell.volume_layers[from_key.layer_index].max(0.0);
                }
            }
        }

        // Second pass: distribute volume to destinations
        for (from_key, to_map) in &flattened_transfers {
            let total_outgoing_volume: f64 = to_map.values().sum();
            
            if total_outgoing_volume <= 0.0 {
                continue;
            }

            // Distribute volume to each destination
            for (to_key, volume_to_transfer) in to_map {
                if let Some(dest_cell) = self.cells.get_mut(&to_key.cell_id) {
                    if to_key.layer_index < LAYER_COUNT {
                        dest_cell.next_cell.volume_layers[to_key.layer_index] += volume_to_transfer;
                        
                        // Ensure non-negative values
                        dest_cell.next_cell.volume_layers[to_key.layer_index] = 
                            dest_cell.next_cell.volume_layers[to_key.layer_index].max(0.0);
                    }
                }
            }
        }

        // Clear accounting records for this step
        accounting.clear_step(current_step);
    }

    /// Get total cell count
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }

    /// Get current step
    pub fn current_step(&self) -> u32 {
        self.step
    }

    /// Calculate simulation statistics (using upper layer arrays for primary stats)
    pub fn calculate_statistics(&self) -> SimulationStats {
        // Use upper layer from arrays for primary visualization stats
        let upper_layer_index = LAYER_COUNT - 1;
        let volumes: Vec<f64> = self.cells
            .values()
            .map(|cell| cell.next_cell.volume_layers[upper_layer_index])
            .collect();
        let energies: Vec<f64> = self.cells
            .values()
            .map(|cell| cell.next_cell.energy_layers[upper_layer_index])
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
