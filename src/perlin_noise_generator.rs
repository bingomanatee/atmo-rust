use noise::{NoiseFn, Perlin};
use glam::Vec3;

/// Simple Perlin noise generator with just two intuitive controls
pub struct PerlinNoiseGenerator {
    perlin: Perlin,
    /// Detail level: how much spatial variation (0.1 = very smooth, 10.0 = very detailed)
    pub detail: f64,
    /// Sharpness: how extreme the values get (0.5 = gentle curves, 3.0 = sharp peaks/valleys)
    pub sharpness: f64,
}

impl PerlinNoiseGenerator {
    /// Create a new noise generator
    /// - seed: random seed for reproducible noise
    /// - detail: spatial frequency (0.1 = smooth, 10.0 = detailed)
    /// - sharpness: value distribution (0.5 = gentle, 3.0 = sharp)
    pub fn new(seed: u32, detail: f64, sharpness: f64) -> Self {
        Self {
            perlin: Perlin::new(seed),
            detail,
            sharpness,
        }
    }
    
    /// Generate noise value from -1 to 1 for a normalized 3D coordinate
    /// Input should be a normalized coordinate (like from Vec3::normalize())
    pub fn sample(&self, normalized_coord: Vec3) -> f64 {
        // Apply detail scaling to the coordinate
        let sample_coords = [
            (normalized_coord.x * self.detail as f32) as f64,
            (normalized_coord.y * self.detail as f32) as f64,
            (normalized_coord.z * self.detail as f32) as f64,
        ];
        
        // Get raw Perlin noise [-1, 1]
        let raw_noise = self.perlin.get(sample_coords);
        
        // Apply sharpness by using power function while preserving sign
        if raw_noise >= 0.0 {
            raw_noise.powf(1.0 / self.sharpness).min(1.0)
        } else {
            -((-raw_noise).powf(1.0 / self.sharpness).min(1.0))
        }
    }
    
    /// Sample multiple locations efficiently
    pub fn sample_batch(&self, coords: &[Vec3]) -> Vec<f64> {
        coords.iter().map(|&coord| self.sample(coord)).collect()
    }
    
    /// Analyze the distribution of values for a set of sample points
    pub fn analyze(&self, sample_coords: &[Vec3]) -> NoiseStats {
        let values: Vec<f64> = sample_coords.iter().map(|&coord| self.sample(coord)).collect();
        
        let min = values.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max = values.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        
        let variance = values.iter()
            .map(|&x| (x - mean).powi(2))
            .sum::<f64>() / values.len() as f64;
        let std_dev = variance.sqrt();
        
        NoiseStats {
            min,
            max,
            mean,
            std_dev,
            sample_count: values.len(),
            detail: self.detail,
            sharpness: self.sharpness,
        }
    }
}

#[derive(Debug)]
pub struct NoiseStats {
    pub min: f64,
    pub max: f64,
    pub mean: f64,
    pub std_dev: f64,
    pub sample_count: usize,
    pub detail: f64,
    pub sharpness: f64,
}

impl std::fmt::Display for NoiseStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Noise Analysis ({} samples):", self.sample_count)?;
        writeln!(f, "  Settings: detail={:.1}, sharpness={:.1}", self.detail, self.sharpness)?;
        writeln!(f, "  Range: {:.3} to {:.3}", self.min, self.max)?;
        writeln!(f, "  Mean: {:.3} ± {:.3}", self.mean, self.std_dev)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_noise_generation() {
        let generator = PerlinNoiseGenerator::new(42, 1.0, 1.0);
        let coord = Vec3::new(1.0, 0.0, 0.0).normalize();
        let value = generator.sample(coord);
        
        // Should be in range [-1, 1]
        assert!(value >= -1.0 && value <= 1.0);
    }
    
    #[test]
    fn test_different_settings() {
        let smooth = PerlinNoiseGenerator::new(42, 0.5, 0.8);  // Low detail, gentle
        let sharp = PerlinNoiseGenerator::new(42, 5.0, 2.0);   // High detail, sharp
        
        let coord = Vec3::new(1.0, 0.5, 0.0).normalize();
        
        let smooth_val = smooth.sample(coord);
        let sharp_val = sharp.sample(coord);
        
        assert!(smooth_val >= -1.0 && smooth_val <= 1.0);
        assert!(sharp_val >= -1.0 && sharp_val <= 1.0);
    }
}