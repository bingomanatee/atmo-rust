use crate::asth_cell_next::AsthenosphereCellNext;
use crate::asthenosphere::{AsthenosphereCell, CellsForPlanetArgs, ASTH_RES};
use crate::binary_pair::BinaryPair;
use crate::constants::{JOULES_PER_KM3, VOLCANO_MAX_VOLUME, LAYER_COUNT, VERTICAL_ENERGY_MIXING};
use crate::convection::{Convection, ConvectionSystem};
use crate::planet::Planet;
use crate::png_exporter::PngExporter;
use crate::rock_store::RockStore;
use h3o::CellIndex;
use noise::NoiseFn;
use rand::Rng;
use std::collections::HashMap;
use std::fs;
use uuid::Uuid;
use rayon::prelude::*;

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
    
    /// Overall execution timing
    pub start_time: std::time::Instant,
    pub initialization_time_s: f64,
}

const LIFT_THRESHOLD: f64 = 0.9;
const LIFT_REPLACEMENT_RATE : f64= 0.33;

impl SimNext {
    /// Create a new simulation with configuration props
    pub fn new(props: SimNextProps) -> Self {
        let start_time = std::time::Instant::now();
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
            start_time,
            initialization_time_s: 0.0,
        };

        sim.initialize_cells_with_props(
            props.seed,
            props.resolution,
            props.joules_per_km3,
        );
        if sim.visualize {
            sim.setup_visualization();
        }

        // Record initialization time
        sim.initialization_time_s = start_time.elapsed().as_secs_f64();
        
        sim.debug_print(&format!("🚀 Simulation initialized in {:.3}s with {} cells", 
            sim.initialization_time_s, sim.cells.len()));

        sim
    }

    /// Helper method for debug printing
    fn debug_print(&self, message: &str) {
        if self.debug {
            println!("{}", message);
        }
    }
    
    /// Print comprehensive timing summary
    pub fn print_timing_summary(&self, steps_completed: u32) {
        let total_execution_time_s = self.start_time.elapsed().as_secs_f64();
        let simulation_time_s = total_execution_time_s - self.initialization_time_s;
        let avg_step_time_s = if steps_completed > 0 { 
            simulation_time_s / steps_completed as f64 
        } else { 
            0.0 
        };
        
        println!("\n⏱️ === EXECUTION TIMING SUMMARY ===");
        println!("📊 Initialization time: {:.3}s", self.initialization_time_s);
        println!("🔄 Steps completed: {}", steps_completed);
        println!("⚡ Average time per step: {:.3}s ({:.1}ms)", avg_step_time_s, avg_step_time_s * 1000.0);
        println!("🎯 Total simulation time: {:.3}s", simulation_time_s);
        println!("🚀 Total execution time: {:.3}s", total_execution_time_s);
        println!("===============================\n");
    }

    pub fn initialize_cells_with_props(
        &mut self,
        seed: u64,
        resolution: h3o::Resolution,
        joules_per_km3: f64,
    ) {
        self.debug_print("🌍 Initializing array-based cells for simulation...");

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

        AsthenosphereCell::initial_cells_for_planet(args);
        self.debug_print(&format!("🌍 Generated {} asthenosphere cells with {} layers each", asth_cells.len(), LAYER_COUNT));

        for asth_cell in asth_cells {
            let sim_cell = AsthenosphereCellNext::new(asth_cell.clone());
            self.cells.insert(sim_cell.cell_index(), sim_cell);
        }

        self.generate_binary_pairs();

        let total_cells = self.cells.len();
        self.convection_system = ConvectionSystem::new(seed as u32, total_cells);

        self.debug_print(&format!(
            "🌍 Initialized {} cells with {} layers each ({} total layer entries) with {} binary pairs",
            total_cells,
            LAYER_COUNT,
            total_cells * LAYER_COUNT,
            self.binary_pairs.len()
        ));
        let current_template = self.convection_system.get_interpolated_template();
        self.debug_print(&format!(
            "🌪️ Convection system initialized: +{:.1}% -{:.1}% (unaffected: {:.1}%, balanced: {}, cycle: {}/{})",
            current_template.addition_fraction * 100.0,
            current_template.subtraction_fraction * 100.0,
            current_template.unaffected_percentage() * 100.0,
            current_template.is_balanced(),
            self.convection_system.step_in_cycle,
            self.convection_system.cycle_lifespan
        ));
        
        self.jumpstart_convection();
    }

    pub fn initialize_cells(&mut self) {
        self.initialize_cells_with_props(
            42,
            crate::asthenosphere::ASTH_RES,
            JOULES_PER_KM3,
        );
    }

    fn setup_visualization(&mut self) {
        self.png_exporter = Some(PngExporter::new(720, 480, self.planet.clone()));

        let output_dir = "vis/sim_next";
        if let Err(e) = fs::create_dir_all(output_dir) {
            eprintln!(
                "Warning: Could not create visualization directory {}: {}",
                output_dir, e
            );
        } else {
            self.debug_print(&format!("📁 Created visualization directory: {}", output_dir));
        }
    }

    fn generate_binary_pairs(&mut self) {
        self.debug_print("🔗 Generating binary pairs...");

        let mut pairs_map = HashMap::new();

        for (&cell_index, cell) in &self.cells {
            for &neighbor_index in &cell.cell.neighbors {
                let pair = BinaryPair::new(cell_index, neighbor_index);
                pairs_map.insert(pair.to_string_id(), pair);
            }
        }

        self.binary_pairs = pairs_map.into_values().collect();
    }
    

    pub fn run_step(&mut self) {
        
        self.debug_print(&format!("🔄 Running step {}", self.step + 1));

        self.level_cells();
        self.apply_material_lift();
        self.cool_cells();
        self.apply_vertical_energy_mixing();
        self.process_volcanoes_and_sinkholes();
        self.try_spawn_volcanoes_and_sinkholes();
        self.apply_convection();
        self.commit_step();

        self.step += 1;

        if self.visualize && self.step % self.vis_freq == 0 {
            self.debug_print(&format!(
                "🖼️ Exporting visualization for step {}",
                self.step
            ));
            self.export_visualization();
        }
    }

    fn export_visualization(&mut self) {
        if !self.visualize || self.png_exporter.is_none() {
            return;
        }

        self.debug_print("  📁 Preparing visualization export...");
        let output_dir = "vis/sim_next";

        let cells_for_export: Vec<(CellIndex, AsthenosphereCell)> = self
            .cells
            .iter()
            .map(|(&cell_id, cell)| {
                let export_cell = cell.cell.clone();
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

    fn level_cells(&mut self) {
        self.debug_print(&format!("⚖️ Levelling {} pairs", self.binary_pairs.len()));

        let pairs = self.binary_pairs.clone();
        for pair in &pairs {
            self.level_binary_pair(pair);
        }

        self.debug_print("⚖️ Completed levelling");
    }

    fn level_binary_pair(&mut self, pair: &BinaryPair) {
        let (mut cell_a, mut cell_b) = match (
            self.cells.get(&pair.cell_a).cloned(),
            self.cells.get(&pair.cell_b).cloned(),
        ) {
            (Some(a), Some(b)) => (a, b),
            _ => return,
        };

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
            let transfer_amount = excess * 0.4;

            if transfer_amount > 0.0 {
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

                if is_a_higher {
                    cell_a.next_cell.volume_layers[layer_index] -= transfer_amount;
                    cell_a.next_cell.energy_layers[layer_index] -= energy_transfer;
                    cell_b.next_cell.volume_layers[layer_index] += transfer_amount;
                    cell_b.next_cell.energy_layers[layer_index] += energy_transfer;
                } else {
                    cell_b.next_cell.volume_layers[layer_index] -= transfer_amount;
                    cell_b.next_cell.energy_layers[layer_index] -= energy_transfer;
                    cell_a.next_cell.volume_layers[layer_index] += transfer_amount;
                    cell_a.next_cell.energy_layers[layer_index] += energy_transfer;
                }

                cell_a.next_cell.volume_layers[layer_index] = cell_a.next_cell.volume_layers[layer_index].max(0.0);
                cell_a.next_cell.energy_layers[layer_index] = cell_a.next_cell.energy_layers[layer_index].max(0.0);
                cell_b.next_cell.volume_layers[layer_index] = cell_b.next_cell.volume_layers[layer_index].max(0.0);
                cell_b.next_cell.energy_layers[layer_index] = cell_b.next_cell.energy_layers[layer_index].max(0.0);
            }
        }

        self.cells.insert(pair.cell_a, cell_a);
        self.cells.insert(pair.cell_b, cell_b);
    }

    fn cool_cells(&mut self) {
        self.cells.par_iter_mut().for_each(|(_, cell)| {
            cell.cool_with_heating(&self.planet, self.step, ASTH_RES);
        });
    }

    fn apply_material_lift(&mut self) {
        use crate::constants::AVG_STARTING_VOLUME_KM_3;
        
        self.debug_print("🔺 Applying material lift for low volume layers");
        
        self.cells.par_iter_mut().for_each(|(_, cell)| {
            let threshold_volume = AVG_STARTING_VOLUME_KM_3 * LIFT_THRESHOLD; // 90% of average starting volume
            
            // Process layers from top to bottom to create true cascading effect
            // Each layer's deficit can trigger lifting that creates deficits in layers below
            for layer_index in (0..LAYER_COUNT).rev() {
                let current_volume = cell.next_cell.volume_layers[layer_index];
                
                if current_volume < threshold_volume {
                    let deficit = threshold_volume - current_volume;
                    let lift_amount = deficit * LIFT_REPLACEMENT_RATE; // 5% of deficit
                    
                    if layer_index == 0 {
                        // Bottom layer: lift from imaginary layer below with 50% higher energy
                        let base_energy_density = if current_volume > 0.0 {
                            cell.next_cell.energy_layers[layer_index] / current_volume
                        } else {
                            crate::constants::JOULES_PER_KM3
                        };
                        let enhanced_energy_density = base_energy_density * 1.5; // 50% higher energy
                        let energy_to_add = lift_amount * enhanced_energy_density;
                        
                        cell.next_cell.volume_layers[layer_index] += lift_amount;
                        cell.next_cell.energy_layers[layer_index] += energy_to_add;
                    } else {
                        // Upper layers: lift from layer below, which may create cascading deficit
                        let lower_layer = layer_index - 1;
                        let lower_volume = cell.next_cell.volume_layers[lower_layer];
                        
                        if lower_volume > lift_amount {
                            // Calculate proportional energy transfer
                            let energy_density = if lower_volume > 0.0 {
                                cell.next_cell.energy_layers[lower_layer] / lower_volume
                            } else {
                                0.0
                            };
                            let energy_to_transfer = lift_amount * energy_density;
                            
                            // Transfer material from lower to upper layer
                            cell.next_cell.volume_layers[lower_layer] -= lift_amount;
                            cell.next_cell.energy_layers[lower_layer] -= energy_to_transfer;
                            cell.next_cell.volume_layers[layer_index] += lift_amount;
                            cell.next_cell.energy_layers[layer_index] += energy_to_transfer;
                            
                            // The removal from lower_layer may create a deficit that will be
                            // processed when we reach that layer in the loop (cascading effect)
                        }
                    }
                    
                    // Ensure non-negative values
                    cell.next_cell.volume_layers[layer_index] = cell.next_cell.volume_layers[layer_index].max(0.0);
                    cell.next_cell.energy_layers[layer_index] = cell.next_cell.energy_layers[layer_index].max(0.0);
                }
            }
            
            // Ensure all layers have non-negative values after cascading
            for layer_index in 0..LAYER_COUNT {
                cell.next_cell.volume_layers[layer_index] = cell.next_cell.volume_layers[layer_index].max(0.0);
                cell.next_cell.energy_layers[layer_index] = cell.next_cell.energy_layers[layer_index].max(0.0);
            }
        });
    }

    fn process_volcanoes_and_sinkholes(&mut self) {
        // Process volcanoes and sinkholes for all cells
        self.cells.par_iter_mut().for_each(|(_, cell)| {
            cell.process_volcanoes_and_sinkholes();
        });
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

    fn apply_convection(&mut self) {
        let progress = self.convection_system.cycle_progress();
        self.debug_print(&format!("🌪️ Applying dynamic convection to bottom layer arrays (cycle {:.1}%)", progress * 100.0));
        
        let template = self.convection_system.get_interpolated_template();
        let t = self.convection_system.step_in_cycle as f64 / self.convection_system.cycle_lifespan as f64;
        let planet = self.planet.clone();
        
        self.cells.par_iter_mut().for_each(|(_, cell)| {
            Self::apply_convection_to_layer_static(&mut cell.next_cell, 0, &template, t, &self.convection_system, &planet);
        });
        
        self.convection_system.advance_step();
        
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
    
    fn jumpstart_convection(&mut self) {
        self.debug_print("🚀 Applying jump-start convection to bottom layer arrays (10x amplification)...");
        
        let mut template = self.convection_system.get_interpolated_template();
        template.per_cell_addition *= 10.0;
        template.per_cell_subtraction *= 10.0;
        let t = self.convection_system.step_in_cycle as f64 / self.convection_system.cycle_lifespan as f64;
        let planet = self.planet.clone();
        
        self.cells.par_iter_mut().for_each(|(_, cell)| {
            Self::apply_convection_to_layer_static(&mut cell.next_cell, 0, &template, t, &self.convection_system, &planet);
        });
        
        for cell in self.cells.values_mut() {
            cell.commit_step();
        }
        
        self.debug_print("✅ Jump-start convection applied successfully to bottom layer arrays");
    }

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
        
        let gc = crate::geoconverter::GeoCellConverter::new(planet.radius_km as f64, ASTH_RES);
        let location = gc.cell_to_vec3(cell.id);
        
        let scaled_location = location.normalize() * template.noise_scale;
        
        let current_noise = convection_system.current_perlin.get(scaled_location.to_array().map(|n| n as f64));
        let next_noise = convection_system.next_perlin.get(scaled_location.to_array().map(|n| n as f64));
        let noise_value = current_noise * (1.0 - t) + next_noise * t;
        
        let addition_threshold = 1.0 - (template.addition_fraction * 2.0);
        let subtraction_threshold = -1.0 + (template.subtraction_fraction * 2.0);
        
        if noise_value > addition_threshold {
            let volume_to_add = template.per_cell_addition;
            let energy_to_add = volume_to_add * (crate::constants::CELL_JOULES_START / crate::constants::AVG_STARTING_VOLUME_KM_3);
            
            cell.volume_layers[layer_index] += volume_to_add;
            cell.energy_layers[layer_index] += energy_to_add;
        } else if noise_value < subtraction_threshold {
            let volume_to_remove = template.per_cell_subtraction;
            let volume_to_remove = volume_to_remove.min(cell.volume_layers[layer_index] * 0.9);
            
            if volume_to_remove > 0.0 && cell.volume_layers[layer_index] > 0.0 {
                let energy_ratio = volume_to_remove / cell.volume_layers[layer_index];
                let energy_to_remove = cell.energy_layers[layer_index] * energy_ratio;
                
                cell.volume_layers[layer_index] -= volume_to_remove;
                cell.energy_layers[layer_index] -= energy_to_remove;
                
                cell.volume_layers[layer_index] = cell.volume_layers[layer_index].max(0.0);
                cell.energy_layers[layer_index] = cell.energy_layers[layer_index].max(0.0);
            }
        }
    }

    fn apply_vertical_energy_mixing(&mut self) {
        if LAYER_COUNT < 2 {
            return;
        }

        self.debug_print(&format!("🔄 Applying vertical volume and energy mixing ({:.1}% exchange)", VERTICAL_ENERGY_MIXING * 100.0));

        self.cells.par_iter_mut().for_each(|(_, cell)| {
            for layer_index in 0..(LAYER_COUNT - 1) {
                let lower_layer = layer_index;
                let upper_layer = layer_index + 1;

                let lower_volume_to_transfer = cell.next_cell.volume_layers[lower_layer] * VERTICAL_ENERGY_MIXING;
                let upper_volume_to_transfer = cell.next_cell.volume_layers[upper_layer] * VERTICAL_ENERGY_MIXING;
                
                let lower_energy_to_transfer = cell.next_cell.energy_layers[lower_layer] * VERTICAL_ENERGY_MIXING;
                let upper_energy_to_transfer = cell.next_cell.energy_layers[upper_layer] * VERTICAL_ENERGY_MIXING;

                cell.next_cell.volume_layers[lower_layer] -= lower_volume_to_transfer;
                cell.next_cell.volume_layers[lower_layer] += upper_volume_to_transfer;

                cell.next_cell.volume_layers[upper_layer] -= upper_volume_to_transfer;
                cell.next_cell.volume_layers[upper_layer] += lower_volume_to_transfer;

                cell.next_cell.energy_layers[lower_layer] -= lower_energy_to_transfer;
                cell.next_cell.energy_layers[lower_layer] += upper_energy_to_transfer;

                cell.next_cell.energy_layers[upper_layer] -= upper_energy_to_transfer;
                cell.next_cell.energy_layers[upper_layer] += lower_energy_to_transfer;

                cell.next_cell.volume_layers[lower_layer] = cell.next_cell.volume_layers[lower_layer].max(0.0);
                cell.next_cell.energy_layers[lower_layer] = cell.next_cell.energy_layers[lower_layer].max(0.0);
                cell.next_cell.volume_layers[upper_layer] = cell.next_cell.volume_layers[upper_layer].max(0.0);
                cell.next_cell.energy_layers[upper_layer] = cell.next_cell.energy_layers[upper_layer].max(0.0);
            }
        });
    }

    fn commit_step(&mut self) {
        self.cells.par_iter_mut().for_each(|(_, cell)| {
            cell.commit_step();
        });
    }

    /// Get total cell count
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }

    /// Get current step
    pub fn current_step(&self) -> u32 {
        self.step
    }

    fn calculate_statistics(&self) -> SimulationStats {
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
            let step_start = std::time::Instant::now();
            self.run_step();
            let step_time_ms = step_start.elapsed().as_secs_f64() * 1000.0;

            if self.step % 50 == 0 {
                println!("⏳ Completed step {}/{} (last step: {:.1}ms)", 
                    self.step, steps, step_time_ms);
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

        // Print comprehensive timing summary
        self.print_timing_summary(steps);
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
