use crate::constants::{
    AVG_STARTING_VOLUME_KM_3, CELL_JOULES_START, GLOBAL_CONVECTION,
    CONVECTION_ADDITION_MIN, CONVECTION_ADDITION_MAX,
    CONVECTION_SUBTRACTION_MIN, CONVECTION_SUBTRACTION_MAX,
    CONVECTION_BALANCE_TOLERANCE, CONVECTION_TEMPLATE_LIFESPAN_MIN,
    CONVECTION_TEMPLATE_LIFESPAN_MAX, CONVECTION_NOISE_SCALE_MIN,
    CONVECTION_NOISE_SCALE_MAX
};
use crate::asthenosphere::AsthenosphereCell;
use crate::planet::Planet;
use crate::geoconverter::GeoCellConverter;
use h3o::{CellIndex, Resolution};
use glam::Vec3;
use noise::{NoiseFn, Perlin};
use rand::Rng;

/// Template defining convection pattern parameters
#[derive(Clone, Debug)]
pub struct ConvectionTemplate {
    pub seed: u32,
    pub addition_fraction: f64,
    pub subtraction_fraction: f64,
    pub noise_scale: f32,
    pub per_cell_addition: f64,
    pub per_cell_subtraction: f64,
}

impl ConvectionTemplate {
    /// Generate a new random convection template
    pub fn generate_random(seed: u32, total_cells: usize) -> Self {
        let mut rng = rand::rng();
        
        // Generate addition and subtraction fractions within specified ranges
        let addition_fraction = rng.random_range(CONVECTION_ADDITION_MIN..CONVECTION_ADDITION_MAX);
        let subtraction_fraction = Self::balanced_subtraction(addition_fraction, &mut rng);
        
        // Random noise scale for variety in convection patterns
        let noise_scale = rng.random_range(CONVECTION_NOISE_SCALE_MIN..CONVECTION_NOISE_SCALE_MAX);
        
        // Calculate global asthenosphere mass (30% of total)
        let global_mass = total_cells as f64 * AVG_STARTING_VOLUME_KM_3 * GLOBAL_CONVECTION;
        
        // Calculate per-cell amounts
        let per_cell_addition = (global_mass * addition_fraction) / total_cells as f64;
        let per_cell_subtraction = (global_mass * subtraction_fraction) / total_cells as f64;
        
        Self {
            seed,
            addition_fraction,
            subtraction_fraction,
            noise_scale,
            per_cell_addition,
            per_cell_subtraction,
        }
    }
    
    /// Generate a subtraction fraction that's within 3% of the addition fraction
    fn balanced_subtraction(addition_fraction: f64, rng: &mut rand::rngs::ThreadRng) -> f64 {
        let tolerance = CONVECTION_BALANCE_TOLERANCE;
        let min_sub = (addition_fraction * (1.0 - tolerance)).max(CONVECTION_SUBTRACTION_MIN);
        let max_sub = (addition_fraction * (1.0 + tolerance)).min(CONVECTION_SUBTRACTION_MAX);
        
        rng.random_range(min_sub..max_sub)
    }
    
    /// Check if the template is balanced (addition ≈ subtraction)
    pub fn is_balanced(&self) -> bool {
        let diff = (self.addition_fraction - self.subtraction_fraction).abs();
        diff <= CONVECTION_BALANCE_TOLERANCE * self.addition_fraction
    }
    
    /// Get the percentage of cells that remain unaffected by convection
    pub fn unaffected_percentage(&self) -> f64 {
        1.0 - (self.addition_fraction + self.subtraction_fraction)
    }
    
    /// Interpolate between templates with expanded range at extremes and power scaling
    /// At t=0.0 and t=1.0: full range (0% to 100%)
    /// At t=0.5: compressed range (10% to 90%)
    /// Power scaling of 1.5 makes transitions favor one template or the other
    pub fn lerp(&self, other: &ConvectionTemplate, t: f64) -> ConvectionTemplate {
        let t = t.clamp(0.0, 1.0);
        
        // Apply power scaling to make transitions more extreme
        // This makes the interpolation spend more time near the extremes (0 or 1)
        let power_scaled_t = if t <= 0.5 {
            (t * 2.0).powf(1.5) * 0.5
        } else {
            1.0 - ((1.0 - t) * 2.0).powf(1.5) * 0.5
        };
        
        // Calculate expansion factor based on distance from midpoint
        // At t=0.0 or t=1.0: expansion_factor = 1.0 (full range)
        // At t=0.5: expansion_factor = 0.8 (10% to 90% range)
        let distance_from_midpoint = (power_scaled_t - 0.5).abs() * 2.0; // 0.0 at midpoint, 1.0 at extremes
        let expansion_factor = 0.8 + (distance_from_midpoint * 0.2); // 0.8 to 1.0 range
        let compression_offset = (1.0 - expansion_factor) * 0.5; // Center the compressed range
        
        // Apply the expansion curve to the interpolation parameter
        let adjusted_t = compression_offset + (power_scaled_t * expansion_factor);
        let adjusted_inv_t = 1.0 - adjusted_t;
        
        ConvectionTemplate {
            seed: if t < 0.5 { self.seed } else { other.seed }, // Switch seed at midpoint
            addition_fraction: self.addition_fraction * adjusted_inv_t + other.addition_fraction * adjusted_t,
            subtraction_fraction: self.subtraction_fraction * adjusted_inv_t + other.subtraction_fraction * adjusted_t,
            noise_scale: self.noise_scale * adjusted_inv_t as f32 + other.noise_scale * adjusted_t as f32,
            per_cell_addition: self.per_cell_addition * adjusted_inv_t + other.per_cell_addition * adjusted_t,
            per_cell_subtraction: self.per_cell_subtraction * adjusted_inv_t + other.per_cell_subtraction * adjusted_t,
        }
    }
}

/// Convection system that interpolates between templates over time
pub struct ConvectionSystem {
    pub current_template: ConvectionTemplate,
    pub next_template: ConvectionTemplate,
    pub current_perlin: Perlin,
    pub next_perlin: Perlin,
    pub step_in_cycle: u32,
    pub cycle_lifespan: u32,
    pub total_cells: usize,
}

impl ConvectionSystem {
    /// Create a new convection system with initial templates
    pub fn new(seed: u32, total_cells: usize) -> Self {
        let mut rng = rand::rng();
        
        // Generate initial current and next templates
        let current_template = ConvectionTemplate::generate_random(seed, total_cells);
        let next_seed = rng.random_range(1..u32::MAX);
        let next_template = ConvectionTemplate::generate_random(next_seed, total_cells);
        
        // Random cycle lifespan
        let cycle_lifespan = rng.random_range(CONVECTION_TEMPLATE_LIFESPAN_MIN..=CONVECTION_TEMPLATE_LIFESPAN_MAX);
        
        Self {
            current_perlin: Perlin::new(current_template.seed),
            next_perlin: Perlin::new(next_template.seed),
            current_template,
            next_template,
            step_in_cycle: 0,
            cycle_lifespan,
            total_cells,
        }
    }
    
    /// Get the current interpolated template based on cycle progress
    pub fn get_interpolated_template(&self) -> ConvectionTemplate {
        let t = self.step_in_cycle as f64 / self.cycle_lifespan as f64;
        self.current_template.lerp(&self.next_template, t)
    }
    
    /// Advance to the next step, potentially creating a new cycle
    pub fn advance_step(&mut self) {
        self.step_in_cycle += 1;
        
        // Check if we need to start a new cycle
        if self.step_in_cycle >= self.cycle_lifespan {
            self.start_new_cycle();
        }
    }
    
    /// Start a new convection cycle with new target template
    fn start_new_cycle(&mut self) {
        let mut rng = rand::rng();
        
        // Current "next" becomes new "current"
        self.current_template = self.next_template.clone();
        self.current_perlin = self.next_perlin;
        
        // Generate new "next" template
        let next_seed = rng.random_range(1..u32::MAX);
        self.next_template = ConvectionTemplate::generate_random(next_seed, self.total_cells);
        self.next_perlin = Perlin::new(self.next_template.seed);
        
        // Reset cycle
        self.step_in_cycle = 0;
        self.cycle_lifespan = rng.random_range(CONVECTION_TEMPLATE_LIFESPAN_MIN..=CONVECTION_TEMPLATE_LIFESPAN_MAX);
    }
    
    /// Apply convection to a single asthenosphere cell using interpolated templates
    pub fn apply_convection(&self, cell: &mut AsthenosphereCell, planet: &Planet, resolution: Resolution) {
        // Get interpolated template for current step
        let template = self.get_interpolated_template();
        
        // Get cell location for perlin noise sampling
        let gc = GeoCellConverter::new(planet.radius_km as f64, resolution);
        let location = gc.cell_to_vec3(cell.id);
        
        // Use template's noise scale
        let scaled_location = location.normalize() * template.noise_scale;
        
        // Interpolate noise values between current and next perlin
        let t = self.step_in_cycle as f64 / self.cycle_lifespan as f64;
        let current_noise = self.current_perlin.get(scaled_location.to_array().map(|n| n as f64));
        let next_noise = self.next_perlin.get(scaled_location.to_array().map(|n| n as f64));
        let noise_value = current_noise * (1.0 - t) + next_noise * t;
        
        // Calculate thresholds for top and bottom percentiles
        let addition_threshold = 1.0 - (template.addition_fraction * 2.0);
        let subtraction_threshold = -1.0 + (template.subtraction_fraction * 2.0);
        
        // NOTE: This method is deprecated - the new SimNext uses array-based convection
        // This is kept for backward compatibility but should not be used with array-based cells
        unimplemented!("Old convection system is not compatible with array-based cells. Use SimNext::apply_convection_to_layer_static instead.");
        // Middle cells remain unaffected
    }
    
    /// Get total material added by current interpolated template
    pub fn total_material_added(&self) -> f64 {
        let template = self.get_interpolated_template();
        self.total_cells as f64 * template.addition_fraction * template.per_cell_addition
    }
    
    /// Get total material removed by current interpolated template
    pub fn total_material_removed(&self) -> f64 {
        let template = self.get_interpolated_template();
        self.total_cells as f64 * template.subtraction_fraction * template.per_cell_subtraction
    }
    
    /// Check if the current interpolated template is balanced
    pub fn is_balanced(&self) -> bool {
        self.get_interpolated_template().is_balanced()
    }
    
    /// Get the percentage of cells that remain unaffected by current template
    pub fn unaffected_percentage(&self) -> f64 {
        self.get_interpolated_template().unaffected_percentage()
    }
    
    /// Get current cycle progress (0.0 to 1.0)
    pub fn cycle_progress(&self) -> f64 {
        self.step_in_cycle as f64 / self.cycle_lifespan as f64
    }
    
    /// Apply jump-start convection with amplified effect for initial dramatic patterns
    pub fn apply_jumpstart_convection(&self, cell: &mut AsthenosphereCell, planet: &Planet, resolution: Resolution, amplification: f64) {
        // Get interpolated template for current step
        let mut template = self.get_interpolated_template();
        
        // Amplify the per-cell amounts for jump-start
        template.per_cell_addition *= amplification;
        template.per_cell_subtraction *= amplification;
        
        // Get cell location for perlin noise sampling
        let gc = GeoCellConverter::new(planet.radius_km as f64, resolution);
        let location = gc.cell_to_vec3(cell.id);
        
        // Use template's noise scale
        let scaled_location = location.normalize() * template.noise_scale;
        
        // Interpolate noise values between current and next perlin
        let t = self.step_in_cycle as f64 / self.cycle_lifespan as f64;
        let current_noise = self.current_perlin.get(scaled_location.to_array().map(|n| n as f64));
        let next_noise = self.next_perlin.get(scaled_location.to_array().map(|n| n as f64));
        let noise_value = current_noise * (1.0 - t) + next_noise * t;
        
        // Calculate thresholds for top and bottom percentiles
        let addition_threshold = 1.0 - (template.addition_fraction * 2.0);
        let subtraction_threshold = -1.0 + (template.subtraction_fraction * 2.0);
        
        // NOTE: This method is deprecated - the new SimNext uses array-based convection
        // This is kept for backward compatibility but should not be used with array-based cells
        unimplemented!("Old jumpstart convection system is not compatible with array-based cells. Use SimNext::apply_convection_to_layer_static instead.");
        // Middle cells remain unaffected
    }
}

pub struct Convection {
    pub seed: u32,
    pub addition_fraction: f64,
    pub subtraction_fraction: f64,
    pub perlin: Perlin,
    pub per_cell_addition: f64,
    pub per_cell_subtraction: f64,
}

impl Convection {
    pub fn new(seed: u32, total_cells: usize) -> Self {
        let mut rng = rand::rng();
        
        // Generate addition and subtraction fractions within specified ranges
        let addition_fraction = rng.random_range(CONVECTION_ADDITION_MIN..CONVECTION_ADDITION_MAX);
        let subtraction_fraction = Self::balanced_subtraction(addition_fraction, &mut rng);
        
        // Calculate global asthenosphere mass (30% of total)
        let global_mass = total_cells as f64 * AVG_STARTING_VOLUME_KM_3 * GLOBAL_CONVECTION;
        
        // Calculate per-cell amounts
        let per_cell_addition = (global_mass * addition_fraction) / total_cells as f64;
        let per_cell_subtraction = (global_mass * subtraction_fraction) / total_cells as f64;
        
        let perlin = Perlin::new(seed);
        
        Self {
            seed,
            addition_fraction,
            subtraction_fraction,
            perlin,
            per_cell_addition,
            per_cell_subtraction,
        }
    }
    
    /// Generate a subtraction fraction that's within 3% of the addition fraction
    fn balanced_subtraction(addition_fraction: f64, rng: &mut rand::rngs::ThreadRng) -> f64 {
        let tolerance = CONVECTION_BALANCE_TOLERANCE;
        let min_sub = (addition_fraction * (1.0 - tolerance)).max(CONVECTION_SUBTRACTION_MIN);
        let max_sub = (addition_fraction * (1.0 + tolerance)).min(CONVECTION_SUBTRACTION_MAX);
        
        rng.random_range(min_sub..max_sub)
    }
    
    /// Apply convection to a single asthenosphere cell
    /// Only affects cells in the extreme percentiles of the perlin noise curve
    pub fn apply_convection(&self, cell: &mut AsthenosphereCell, planet: &Planet, resolution: Resolution) {
        // Get cell location for perlin noise sampling
        let gc = GeoCellConverter::new(planet.radius_km as f64, resolution);
        let location = gc.cell_to_vec3(cell.id);
        
        // Use same scale as asthenosphere initialization for consistency
        let noise_scale = 5.0;
        let scaled_location = location.normalize() * noise_scale;
        let noise_value = self.perlin.get(scaled_location.to_array().map(|n| n as f64));
        
        // Calculate thresholds for top and bottom percentiles
        // addition_fraction determines what percentage of cells get material added (top percentile)
        // subtraction_fraction determines what percentage get material removed (bottom percentile)
        let addition_threshold = 1.0 - (self.addition_fraction * 2.0); // Convert to threshold in [-1, 1] range
        let subtraction_threshold = -1.0 + (self.subtraction_fraction * 2.0); // Convert to threshold in [-1, 1] range
        
        // NOTE: This method is deprecated - the new SimNext uses array-based convection
        // This is kept for backward compatibility but should not be used with array-based cells
        unimplemented!("Old convection system is not compatible with array-based cells. Use SimNext::apply_convection_to_layer_static instead.");
        // Middle 60-80% of cells remain unaffected
    }
    
    /// Get total material added by this convection system
    pub fn total_material_added(&self, total_cells: usize) -> f64 {
        // Only addition_fraction percentage of cells get material added
        total_cells as f64 * self.addition_fraction * self.per_cell_addition
    }
    
    /// Get total material removed by this convection system
    pub fn total_material_removed(&self, total_cells: usize) -> f64 {
        // Only subtraction_fraction percentage of cells get material removed
        total_cells as f64 * self.subtraction_fraction * self.per_cell_subtraction
    }
    
    /// Check if the convection system is balanced (addition ≈ subtraction)
    pub fn is_balanced(&self) -> bool {
        let diff = (self.addition_fraction - self.subtraction_fraction).abs();
        diff <= CONVECTION_BALANCE_TOLERANCE * self.addition_fraction
    }
    
    /// Get the percentage of cells that remain unaffected by convection
    pub fn unaffected_percentage(&self) -> f64 {
        1.0 - (self.addition_fraction + self.subtraction_fraction)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::EARTH;
    use h3o::Resolution;

    #[test]
    fn test_convection_creation() {
        let convection = Convection::new(42, 1000);
        
        assert!(convection.addition_fraction >= CONVECTION_ADDITION_MIN);
        assert!(convection.addition_fraction <= CONVECTION_ADDITION_MAX);
        assert!(convection.subtraction_fraction >= CONVECTION_SUBTRACTION_MIN);
        assert!(convection.subtraction_fraction <= CONVECTION_SUBTRACTION_MAX);
        assert!(convection.is_balanced());
    }
    
    #[test]
    fn test_convection_balance() {
        let convection = Convection::new(123, 1000);
        let total_added = convection.total_material_added(1000);
        let total_removed = convection.total_material_removed(1000);
        
        // Should be approximately balanced
        let ratio = total_added / total_removed;
        assert!((ratio - 1.0).abs() < 0.1); // Within 10%
    }
    
    #[test]
    fn test_apply_convection() {
        let convection = Convection::new(42, 1000);
        let planet = EARTH.clone();
        
        let mut cell = AsthenosphereCell {
            id: CellIndex::base_cells().next().unwrap(),
            neighbors: vec![],
            step: 0,
            volume_layers: vec![AVG_STARTING_VOLUME_KM_3; crate::constants::LAYER_COUNT],
            energy_layers: vec![CELL_JOULES_START; crate::constants::LAYER_COUNT],
            volcano_volume: 0.0,
            sinkhole_volume: 0.0,
        };
        
        let original_volume = cell.volume_layers[0];
        let _original_energy = cell.energy_layers[0];
        
        convection.apply_convection(&mut cell, &planet, Resolution::Two);
        
        // Energy should be non-negative for all layers
        for layer_idx in 0..crate::constants::LAYER_COUNT {
            assert!(cell.energy_layers[layer_idx] >= 0.0);
            assert!(cell.volume_layers[layer_idx] >= 0.0);
        }
        
        // Note: Volume may or may not change depending on perlin noise value for this specific cell
        // Most cells (60-80%) should remain unchanged
    }
    
    #[test]
    fn test_unaffected_percentage() {
        let convection = Convection::new(123, 1000);
        let unaffected = convection.unaffected_percentage();
        
        // Should be between 60-80% unaffected
        assert!(unaffected >= 0.6);
        assert!(unaffected <= 0.8);
        
        // Should equal 1 minus the sum of addition and subtraction fractions
        let expected = 1.0 - (convection.addition_fraction + convection.subtraction_fraction);
        assert!((unaffected - expected).abs() < 1e-10);
    }

    #[test]
    fn test_convection_template_creation() {
        let template = ConvectionTemplate::generate_random(42, 1000);
        
        assert!(template.addition_fraction >= CONVECTION_ADDITION_MIN);
        assert!(template.addition_fraction <= CONVECTION_ADDITION_MAX);
        assert!(template.subtraction_fraction >= CONVECTION_SUBTRACTION_MIN);
        assert!(template.subtraction_fraction <= CONVECTION_SUBTRACTION_MAX);
        assert!(template.noise_scale >= CONVECTION_NOISE_SCALE_MIN);
        assert!(template.noise_scale <= CONVECTION_NOISE_SCALE_MAX);
        assert!(template.is_balanced());
    }
    
    #[test]
    fn test_template_interpolation() {
        let template1 = ConvectionTemplate::generate_random(42, 1000);
        let template2 = ConvectionTemplate::generate_random(123, 1000);
        
        // Test interpolation at midpoint (should be closer to center but not exactly 50%)
        let mid_template = template1.lerp(&template2, 0.5);
        // With expansion curve, midpoint is compressed to 10%-90% range
        // So the interpolated value should be closer to the center than simple linear interpolation
        let simple_mid = (template1.addition_fraction + template2.addition_fraction) / 2.0;
        let diff_from_simple = (mid_template.addition_fraction - simple_mid).abs();
        // The compressed interpolation should be closer to the average than linear interpolation
        assert!(diff_from_simple < 0.1); // Should be reasonably close but not exact
        
        // Test boundary conditions - these should still be exact
        let start_template = template1.lerp(&template2, 0.0);
        assert!((start_template.addition_fraction - template1.addition_fraction).abs() < 1e-10);
        
        let end_template = template1.lerp(&template2, 1.0);
        assert!((end_template.addition_fraction - template2.addition_fraction).abs() < 1e-10);
        
        // Test that values at extremes (0.25, 0.75) are more different from center than midpoint
        let quarter_template = template1.lerp(&template2, 0.25);
        let three_quarter_template = template1.lerp(&template2, 0.75);
        
        let quarter_diff = (quarter_template.addition_fraction - simple_mid).abs();
        let three_quarter_diff = (three_quarter_template.addition_fraction - simple_mid).abs();
        let mid_diff = (mid_template.addition_fraction - simple_mid).abs();
        
        // Values at 0.25 and 0.75 should be further from center than midpoint due to expansion
        assert!(quarter_diff > mid_diff || three_quarter_diff > mid_diff);
    }
    
    #[test]
    fn test_convection_system_creation() {
        let system = ConvectionSystem::new(42, 1000);
        
        assert!(system.cycle_lifespan >= CONVECTION_TEMPLATE_LIFESPAN_MIN);
        assert!(system.cycle_lifespan <= CONVECTION_TEMPLATE_LIFESPAN_MAX);
        assert_eq!(system.step_in_cycle, 0);
        assert_eq!(system.total_cells, 1000);
        assert!(system.is_balanced());
    }
    
    #[test]
    fn test_convection_system_step_advance() {
        let mut system = ConvectionSystem::new(42, 1000);
        let initial_template = system.get_interpolated_template();
        let initial_cycle_lifespan = system.cycle_lifespan;
        
        // Advance one step
        system.advance_step();
        assert_eq!(system.step_in_cycle, 1);
        
        // Template should be slightly different due to interpolation
        let new_template = system.get_interpolated_template();
        // Since we're only 1 step in, the difference should be small
        let progress = 1.0 / initial_cycle_lifespan as f64;
        assert!(progress < 0.1); // Should be less than 10% progress
        
        // Progress should be non-zero
        assert!(system.cycle_progress() > 0.0);
        assert!(system.cycle_progress() <= 1.0);
    }
    
    #[test]
    fn test_convection_system_cycle_completion() {
        let mut system = ConvectionSystem::new(42, 100); // Small total_cells for faster test
        let initial_cycle_lifespan = system.cycle_lifespan;
        let initial_current_seed = system.current_template.seed;
        let initial_next_seed = system.next_template.seed;
        
        // Advance to just before cycle completion
        for _ in 0..(initial_cycle_lifespan - 1) {
            system.advance_step();
        }
        assert_eq!(system.step_in_cycle, initial_cycle_lifespan - 1);
        
        // One more step should trigger new cycle
        system.advance_step();
        assert_eq!(system.step_in_cycle, 0); // Reset to 0
        
        // The old "next" should become new "current"
        assert_eq!(system.current_template.seed, initial_next_seed);
        
        // Should have a new "next" template
        assert_ne!(system.next_template.seed, initial_next_seed);
        assert_ne!(system.next_template.seed, initial_current_seed);
    }
}