use noise::{NoiseFn, Perlin};
use glam::Vec3;

/// A component for handling Perlin noise-based scaling for asthenosphere cells
pub struct PerlinScaler {
    perlin: Perlin,
    noise_scale: f64,
    pub volume_config: ScalingConfig,
    pub energy_config: ScalingConfig,
}

#[derive(Clone, Debug)]
pub struct ScalingConfig {
    /// Multiplier applied to location coordinates for noise sampling
    pub coordinate_multiplier: f64,
    /// Power applied to the raw noise value for exponential scaling
    pub noise_power: f64,
    /// Curve shape: 0.0 = linear, 1.0 = quadratic, 2.0 = cubic, etc.
    pub curve_power: f64,
    /// Minimum scaling factor (e.g., 0.8 means minimum 80% of base value)
    pub min_scale: f64,
    /// Maximum scaling factor (e.g., 1.2 means maximum 120% of base value)
    pub max_scale: f64,
}

impl ScalingConfig {
    /// Create a new scaling config with builder pattern
    pub fn new() -> ScalingConfigBuilder {
        ScalingConfigBuilder::default()
    }
    
    /// Create a conservative scaling config with minimal variation
    pub fn conservative() -> Self {
        Self {
            coordinate_multiplier: 1.0,
            noise_power: 1.5,
            curve_power: 1.0, // quadratic curve
            min_scale: 0.9,   // 90% minimum
            max_scale: 1.1,   // 110% maximum
        }
    }
    
    /// Create a moderate scaling config 
    pub fn moderate() -> Self {
        Self {
            coordinate_multiplier: 1.0,
            noise_power: 2.0,
            curve_power: 1.5, // cubic curve
            min_scale: 0.8,   // 80% minimum
            max_scale: 1.2,   // 120% maximum
        }
    }
    
    /// Create an extreme scaling config with high variation
    pub fn extreme() -> Self {
        Self {
            coordinate_multiplier: 1.7,
            noise_power: 1.5,
            curve_power: 2.0, // quartic curve
            min_scale: 0.6,   // 60% minimum
            max_scale: 1.4,   // 140% maximum
        }
    }
    
    /// Create a linear scaling config for testing/comparison
    pub fn linear(min_scale: f64, max_scale: f64) -> Self {
        Self {
            coordinate_multiplier: 1.0,
            noise_power: 1.0,
            curve_power: 0.0, // linear curve
            min_scale,
            max_scale,
        }
    }
}

/// Builder for ScalingConfig with clear parameter names
#[derive(Debug)]
pub struct ScalingConfigBuilder {
    pub coordinate_multiplier: f64,
    pub noise_power: f64,
    pub curve_power: f64,
    pub min_scale: f64,
    pub max_scale: f64,
}

impl Default for ScalingConfigBuilder {
    fn default() -> Self {
        Self {
            coordinate_multiplier: 1.0,
            noise_power: 2.0,
            curve_power: 1.0,
            min_scale: 0.8,
            max_scale: 1.2,
        }
    }
}

impl ScalingConfigBuilder {
    /// Set the coordinate multiplier (affects noise frequency)
    /// Higher values = more variation across space
    pub fn coordinate_multiplier(mut self, value: f64) -> Self {
        self.coordinate_multiplier = value;
        self
    }
    
    /// Set the noise power (shapes the raw Perlin noise)
    /// Higher values = more extreme noise variations
    pub fn noise_power(mut self, value: f64) -> Self {
        self.noise_power = value;
        self
    }
    
    /// Set the curve power (shapes the final 0..1 distribution)
    /// 0.0 = linear, 1.0 = quadratic, 2.0 = cubic, etc.
    /// < 1.0 = more values near max, > 1.0 = more values near min
    pub fn curve_power(mut self, value: f64) -> Self {
        self.curve_power = value;
        self
    }
    
    /// Set the minimum scaling factor
    pub fn min_scale(mut self, value: f64) -> Self {
        self.min_scale = value;
        self
    }
    
    /// Set the maximum scaling factor
    pub fn max_scale(mut self, value: f64) -> Self {
        self.max_scale = value;
        self
    }
    
    /// Build the final ScalingConfig
    pub fn build(self) -> ScalingConfig {
        ScalingConfig {
            coordinate_multiplier: self.coordinate_multiplier,
            noise_power: self.noise_power,
            curve_power: self.curve_power,
            min_scale: self.min_scale,
            max_scale: self.max_scale,
        }
    }
}

impl PerlinScaler {
    pub fn new(seed: u32, noise_scale: f64, volume_config: ScalingConfig, energy_config: ScalingConfig) -> Self {
        Self {
            perlin: Perlin::new(seed),
            noise_scale,
            volume_config,
            energy_config,
        }
    }
    
    /// Create a default scaler with current asthenosphere settings
    pub fn default_asthenosphere(seed: u32) -> Self {
        let volume_config = ScalingConfig::new()
            .coordinate_multiplier(1.0)
            .noise_power(2.0)
            .curve_power(1.0)
            .min_scale(0.8)
            .max_scale(1.0)
            .build();
            
        let energy_config = ScalingConfig::new()
            .coordinate_multiplier(1.0)
            .noise_power(1.5)
            .curve_power(1.0)
            .min_scale(0.6)
            .max_scale(1.0)
            .build();
        
        Self::new(seed, 5.0, volume_config, energy_config)
    }
    
    /// Get volume scaling factor for a given 3D location
    pub fn get_volume_scale(&self, location: Vec3) -> f64 {
        self.get_scale_for_config(location, &self.volume_config)
    }
    
    /// Get energy scaling factor for a given 3D location  
    pub fn get_energy_scale(&self, location: Vec3) -> f64 {
        self.get_scale_for_config(location, &self.energy_config)
    }
    
    /// Internal method to compute scaling based on config
    /// This works in normalized 0..1 space until the final mapping
    fn get_scale_for_config(&self, location: Vec3, config: &ScalingConfig) -> f64 {
        // Step 1: Get raw Perlin noise [-1, 1]
        let scaled_location = location.normalize() * self.noise_scale as f32;
        let sample_coords = [
            (scaled_location.x * config.coordinate_multiplier as f32) as f64,
            (scaled_location.y * config.coordinate_multiplier as f32) as f64, 
            (scaled_location.z * config.coordinate_multiplier as f32) as f64,
        ];
        let raw_noise = self.perlin.get(sample_coords); // [-1, 1]
        
        // Step 2: Apply noise power while preserving sign and normalize to [0, 1]
        let shaped_noise = if raw_noise >= 0.0 {
            raw_noise.powf(config.noise_power)
        } else {
            -((-raw_noise).powf(config.noise_power))
        };
        let normalized_noise = (shaped_noise + 1.0) / 2.0; // [0, 1]
        
        // Step 3: Apply curve shaping in [0, 1] space
        let curved_value = if config.curve_power == 0.0 {
            normalized_noise // Linear (no curve)
        } else {
            // Apply curve power to shape the distribution
            // For values [0, 1], this creates different curve shapes:
            // curve_power < 1.0: more values near 1.0 (convex)
            // curve_power = 1.0: linear
            // curve_power > 1.0: more values near 0.0 (concave)
            normalized_noise.powf(config.curve_power)
        };
        
        // Step 4: Map from [0, 1] to [min_scale, max_scale]
        let scale_range = config.max_scale - config.min_scale;
        config.min_scale + curved_value * scale_range
    }
    
    /// Get the raw normalized value [0, 1] before final scaling (for analysis)
    pub fn get_normalized_value(&self, location: Vec3, config: &ScalingConfig) -> f64 {
        let scaled_location = location.normalize() * self.noise_scale as f32;
        let sample_coords = [
            (scaled_location.x * config.coordinate_multiplier as f32) as f64,
            (scaled_location.y * config.coordinate_multiplier as f32) as f64, 
            (scaled_location.z * config.coordinate_multiplier as f32) as f64,
        ];
        let raw_noise = self.perlin.get(sample_coords);
        let shaped_noise = if raw_noise >= 0.0 {
            raw_noise.powf(config.noise_power)
        } else {
            -((-raw_noise).powf(config.noise_power))
        };
        let normalized_noise = (shaped_noise + 1.0) / 2.0;
        
        if config.curve_power == 0.0 {
            normalized_noise
        } else {
            normalized_noise.powf(config.curve_power)
        }
    }
    
    /// Get both volume and energy scales at once for efficiency
    pub fn get_scales(&self, location: Vec3) -> (f64, f64) {
        (self.get_volume_scale(location), self.get_energy_scale(location))
    }
    
    /// Debug method to print scaling statistics for a set of sample locations
    pub fn analyze_scaling(&self, sample_locations: &[Vec3]) -> ScalingStats {
        let mut volume_scales = Vec::new();
        let mut energy_scales = Vec::new();
        
        for &location in sample_locations {
            volume_scales.push(self.get_volume_scale(location));
            energy_scales.push(self.get_energy_scale(location));
        }
        
        ScalingStats {
            volume_stats: compute_stats(&volume_scales),
            energy_stats: compute_stats(&energy_scales),
            sample_count: sample_locations.len(),
        }
    }
}

#[derive(Debug)]
pub struct ScalingStats {
    pub volume_stats: Stats,
    pub energy_stats: Stats,
    pub sample_count: usize,
}

#[derive(Debug)]
pub struct Stats {
    pub min: f64,
    pub max: f64,
    pub mean: f64,
    pub std_dev: f64,
}

fn compute_stats(values: &[f64]) -> Stats {
    let min = values.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let max = values.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    
    let variance = values.iter()
        .map(|&x| (x - mean).powi(2))
        .sum::<f64>() / values.len() as f64;
    let std_dev = variance.sqrt();
    
    Stats { min, max, mean, std_dev }
}

impl std::fmt::Display for ScalingStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Perlin Scaling Analysis ({} samples):", self.sample_count)?;
        writeln!(f, "Volume Scaling:")?;
        writeln!(f, "  Range: {:.3} to {:.3}", self.volume_stats.min, self.volume_stats.max)?;
        writeln!(f, "  Mean: {:.3} ± {:.3}", self.volume_stats.mean, self.volume_stats.std_dev)?;
        writeln!(f, "Energy Scaling:")?;
        writeln!(f, "  Range: {:.3} to {:.3}", self.energy_stats.min, self.energy_stats.max)?;
        writeln!(f, "  Mean: {:.3} ± {:.3}", self.energy_stats.mean, self.energy_stats.std_dev)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_perlin_scaler_basic() {
        let scaler = PerlinScaler::default_asthenosphere(42);
        
        let location = Vec3::new(1.0, 0.0, 0.0);
        let volume_scale = scaler.get_volume_scale(location);
        let energy_scale = scaler.get_energy_scale(location);
        
        // Should be within expected ranges
        assert!(volume_scale >= 0.6 && volume_scale <= 1.2);
        assert!(energy_scale >= 0.4 && energy_scale <= 1.2);
    }
    
    #[test]
    fn test_scaling_configs() {
        let conservative = ScalingConfig::conservative();
        let moderate = ScalingConfig::moderate(); 
        let extreme = ScalingConfig::extreme();
        
        assert!(conservative.min_scale > moderate.min_scale);
        assert!(moderate.min_scale > extreme.min_scale);
    }
    
    #[test]
    fn test_builder_pattern() {
        let config = ScalingConfig::new()
            .coordinate_multiplier(2.0)
            .noise_power(1.5)
            .curve_power(2.0)
            .min_scale(0.5)
            .max_scale(1.5)
            .build();
            
        assert_eq!(config.coordinate_multiplier, 2.0);
        assert_eq!(config.noise_power, 1.5);
        assert_eq!(config.curve_power, 2.0);
        assert_eq!(config.min_scale, 0.5);
        assert_eq!(config.max_scale, 1.5);
    }
}