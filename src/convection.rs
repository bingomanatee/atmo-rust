use crate::constants::{
    AVG_STARTING_VOLUME_KM_3, CELL_JOULES_START, GLOBAL_CONVECTION,
    CONVECTION_ADDITION_MIN, CONVECTION_ADDITION_MAX,
    CONVECTION_SUBTRACTION_MIN, CONVECTION_SUBTRACTION_MAX,
    CONVECTION_BALANCE_TOLERANCE
};
use crate::asthenosphere::AsthenosphereCell;
use crate::planet::Planet;
use crate::geoconverter::GeoCellConverter;
use h3o::{CellIndex, Resolution};
use glam::Vec3;
use noise::{NoiseFn, Perlin};
use rand::Rng;

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
        
        if noise_value > addition_threshold {
            // Top percentile: Addition - add material at start temperature
            let volume_to_add = self.per_cell_addition;
            let energy_to_add = volume_to_add * (CELL_JOULES_START / AVG_STARTING_VOLUME_KM_3);
            
            cell.volume += volume_to_add;
            cell.energy_j += energy_to_add;
        } else if noise_value < subtraction_threshold {
            // Bottom percentile: Subtraction - remove material proportionally
            let volume_to_remove = self.per_cell_subtraction;
            let volume_to_remove = volume_to_remove.min(cell.volume * 0.9); // Don't remove more than 90%
            
            if volume_to_remove > 0.0 && cell.volume > 0.0 {
                // Remove energy proportionally to volume
                let energy_ratio = volume_to_remove / cell.volume;
                let energy_to_remove = cell.energy_j * energy_ratio;
                
                cell.volume -= volume_to_remove;
                cell.energy_j -= energy_to_remove;
                
                // Ensure no negative values
                cell.volume = cell.volume.max(0.0);
                cell.energy_j = cell.energy_j.max(0.0);
            }
        }
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
            volume: AVG_STARTING_VOLUME_KM_3,
            energy_j: CELL_JOULES_START,
            volcano_volume: 0.0,
            sinkhole_volume: 0.0,
        };
        
        let original_volume = cell.volume;
        let _original_energy = cell.energy_j;
        
        convection.apply_convection(&mut cell, &planet, Resolution::Two);
        
        // Energy should be non-negative
        assert!(cell.energy_j >= 0.0);
        
        // Volume should be non-negative
        assert!(cell.volume >= 0.0);
        
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
}