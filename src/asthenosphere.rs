use crate::constants::{
    AVG_STARTING_VOLUME_KM_3, JOULES_PER_KM3, MAX_SUNK_PERCENT, MAX_SUNK_TEMP, MIN_SUNK_TEMP,
    LAYER_COUNT, ENERGY_INCREASE_PER_LAYER,
};
use crate::geoconverter::GeoCellConverter;
use crate::h3_utils::H3Utils;
use crate::planet::Planet;
use glam::Vec3;
use h3o::{CellIndex, Resolution};
use noise::{NoiseFn, Perlin};
use serde::{Deserialize, Serialize};

/**
some assumptions: the meaningful volume of the asthenosphere we track is a 10km
at l2, the area is 426 km^2, the volume of the average cell is 4260 km^3
*/
#[derive(Serialize, Deserialize, Clone, Debug)]
pub(crate) struct AsthenosphereCell {
    pub id: CellIndex,
    pub neighbors: Vec<CellIndex>,
    pub step: u32,
    pub volume_layers: Vec<f64>,   // Volume at each layer depth (c. 2 MIO KM3 per layer)
    pub energy_layers: Vec<f64>,   // Energy at each layer depth (Joules)
    pub volcano_volume: f64, // Additional volume from volcanic activity (surface only)
    pub sinkhole_volume: f64, // Volume removed by sinkhole activity (surface only)
}

impl Default for AsthenosphereCell {
    fn default() -> Self {
        let cell = CellIndex::base_cells().last().expect("no last cell");
        AsthenosphereCell {
            id: cell,
            volume_layers: vec![0.0; LAYER_COUNT],
            energy_layers: vec![0.0; LAYER_COUNT],
            neighbors: Vec::new(),
            step: 0,
            volcano_volume: 0.0,
            sinkhole_volume: 0.0,
        }
    }
}

pub const ASTH_RES: Resolution = Resolution::Two;

#[derive(Clone)]
pub struct CellsForPlanetArgs<F>
where
    F: FnMut(AsthenosphereCell) -> AsthenosphereCell,
{
    pub planet: Planet,
    pub on_cell: F,
    pub res: Resolution,
    pub joules_per_km3: f64,
    pub seed: u32,
}

#[derive(Clone)]
struct AsthenosphereCellParams {
    cell: CellIndex,
    volume_kg_3: f64,
    energy_per_volume: f64,
}

impl AsthenosphereCell {
    /// Get surface layer volume (top layer only - for visualization)
    pub fn volume(&self) -> f64 {
        let surface_layer = LAYER_COUNT - 1;
        self.volume_layers[surface_layer]
    }
    
    /// Get surface layer energy (top layer only - for visualization)
    pub fn energy_j(&self) -> f64 {
        let surface_layer = LAYER_COUNT - 1;
        self.energy_layers[surface_layer]
    }
    
    /// Apply smooth noise scaling to raw Perlin noise value
    /// Uses gentle exponential curve and maps to a specified range
    fn smooth_noise_scale(noise_value: f64, exponent: f64, min_scale: f64, max_scale: f64) -> f64 {
        let curved_noise = if noise_value > 0.0 {
            noise_value.powf(exponent)
        } else {
            -((-noise_value).powf(exponent))
        };
        
        // Map [-1,1] to [min_scale, max_scale]
        let range = max_scale - min_scale;
        min_scale + (curved_noise + 1.0) * range / 2.0
    }

    /// Calculate initial volume for a cell based on noise
    fn initial_volume(noise: &Perlin, location: Vec3, noise_scale: f32) -> f64 {
        let scaled_location = location.normalize() * noise_scale;
        let noise_val = noise.get(scaled_location.to_array().map(|n| n as f64));

        // Gentler variation: 0.7x to 1.3x the base volume
        let random_scale = Self::smooth_noise_scale(noise_val, 0.8, 0.9, 1.1);
        AVG_STARTING_VOLUME_KM_3 * random_scale
    }

    /// Calculate initial energy per volume for a cell based on independent noise
    fn initial_energy_j(noise: &Perlin, location: Vec3, noise_scale: f32, base_energy_per_volume: f64) -> f64 {
        let scaled_location = location.normalize() * noise_scale;
        let perlin_value = noise.get([
            scaled_location.x as f64, 
            scaled_location.y as f64, 
            scaled_location.z as f64
        ]);
        
        // Gentler energy variation: 0.8x to 1.2x base energy per volume
        let energy_scale = Self::smooth_noise_scale(perlin_value, 0.3, 0.9, 1.1);
        base_energy_per_volume * energy_scale
    }

    /**
        this utility method creates a series of Asthenosphere cells for multiple layers
        given a planet and resolution and feed it back to the callback;
        it is a utility for multi-layer asth_sim.
    */
    pub fn initial_cells_for_planet_multilayer<F>(args: CellsForPlanetArgs<F>)
    where
        F: FnMut(AsthenosphereCell) -> AsthenosphereCell,
    {
        let CellsForPlanetArgs {
            planet,
            mut on_cell,
            res,
            joules_per_km3: base_energy_per_volume,
            seed,
        } = args;
        let resolution = res;
        let gc = GeoCellConverter::new(planet.radius_km as f64, resolution);
        
        // Create separate Perlin noise generators for each layer
        let volume_noise = Perlin::new(seed);
        let energy_noise = Perlin::new(seed + 1000); // Offset seed for energy variation
        
        let mut index = 0u32;
        H3Utils::iter_at(res, |l2_cell| {
            let location = gc.cell_to_vec3(l2_cell);
            let noise_scale = 5.0;
            
            // Create cells for each layer
            for level in 0..LAYER_COUNT {
                // Calculate volume (same for all layers)
                let volume = Self::initial_volume(&volume_noise, location, noise_scale);
                
                // Calculate energy per volume with layer-based increase
                let layer_energy_multiplier = 1.0 + (level as f64 * ENERGY_INCREASE_PER_LAYER);
                let layer_base_energy = base_energy_per_volume * layer_energy_multiplier;
                let adjusted_energy_per_volume = Self::initial_energy_j(&energy_noise, location, noise_scale, layer_base_energy);
                
                let asth_cell = AsthenosphereCell::at_cell(AsthenosphereCellParams {
                    cell: l2_cell,
                    volume_kg_3: volume,
                    energy_per_volume: adjusted_energy_per_volume,
                });

                on_cell(asth_cell);
            }
            index += 1;
        });
    }

    /**
        this utility method creates a series of Asthenosphere cells for each cell in an h3 grid
        given a planet and resolution and feed it back to the callback;
        it is a utility for asth_sim.
    */
    pub fn initial_cells_for_planet<F>(args: CellsForPlanetArgs<F>)
    where
        F: FnMut(AsthenosphereCell) -> AsthenosphereCell,
    {
        let CellsForPlanetArgs {
            planet,
            mut on_cell,
            res,
            joules_per_km3: energy_per_volume,
            seed,
        } = args;
        let resolution = res;
        let gc = GeoCellConverter::new(planet.radius_km as f64, resolution);
        let per = Perlin::new(seed);
        let mut index = 0u32;
        H3Utils::iter_at(res, |l2_cell| {
            let location = gc.cell_to_vec3(l2_cell);
            let noise_scale = 5.0;
            
            // Calculate volume and energy using utility methods
            let volume = Self::initial_volume(&per, location, noise_scale);
            let adjusted_energy_per_volume = Self::initial_energy_j(&per, location, noise_scale, energy_per_volume);
            
            let asth_cell = AsthenosphereCell::at_cell(AsthenosphereCellParams {
                cell: l2_cell,
                volume_kg_3: volume,
                energy_per_volume: adjusted_energy_per_volume,
            });

            on_cell(asth_cell);
            index += 1;
        });
    }

    pub fn sunk_volume(&self) -> f64 {
        // Use surface layer (top layer) for sinkhole calculations
        let surface_layer = LAYER_COUNT - 1;
        if self.energy_layers[surface_layer] >= MAX_SUNK_TEMP {
            return 0.0;
        }
        let effective_energy: f64 = self.energy_layers[surface_layer].max(MIN_SUNK_TEMP);
        let percent = MAX_SUNK_PERCENT * effective_energy / (MAX_SUNK_TEMP / MIN_SUNK_TEMP);
        self.volume_layers[surface_layer] * percent
    }

    pub fn at_cell(params: AsthenosphereCellParams) -> AsthenosphereCell {
        let mut volume_layers = vec![0.0; LAYER_COUNT];
        let mut energy_layers = vec![0.0; LAYER_COUNT];
        
        // Initialize arrays with layer-specific values
        for layer_idx in 0..LAYER_COUNT {
            volume_layers[layer_idx] = params.volume_kg_3;
            // Apply energy increase per layer (layer 0 = base, layer 1 = 1.5x base, etc.)
            let layer_energy_multiplier = 1.0 + (layer_idx as f64 * ENERGY_INCREASE_PER_LAYER);
            energy_layers[layer_idx] = params.energy_per_volume * params.volume_kg_3 * layer_energy_multiplier;
        }
        
        AsthenosphereCell {
            id: params.cell,
            volume_layers,
            energy_layers,
            step: 0,
            neighbors: params
                .cell
                .grid_disk::<Vec<_>>(1)
                .into_iter()
                .filter(|&c| c != params.cell)
                .collect::<Vec<CellIndex>>(),
            volcano_volume: 0.0,
            sinkhole_volume: 0.0,
        }
    }

    pub fn location(&self, planet: &Planet, res: Option<Resolution>) -> Vec3 {
        let resolution = res.unwrap_or(ASTH_RES);
        let gc = GeoCellConverter::new(planet.radius_km as f64, resolution);
        gc.cell_to_vec3(self.id)
    }

    /// Move volume and proportional energy to another cell (conservative binary operation)
    /// Returns true if transfer occurred, false if no transfer (invalid volume or insufficient source)
    /// NOTE: This method works with surface layer only for backward compatibility
    pub fn transfer_volume(&mut self, volume: f64, other: &mut AsthenosphereCell) -> bool {
        let surface_layer = LAYER_COUNT - 1;
        
        // Validation: volume must be positive
        if volume <= 0.0 {
            return false; // Silent return for zero/negative volume
        }

        // Validation: source must have sufficient volume
        if self.volume_layers[surface_layer] < volume {
            return false; // Cannot transfer more than available
        }

        // Validation: source would not go below zero (additional safety check)
        let remaining_volume = self.volume_layers[surface_layer] - volume;
        if remaining_volume < 0.0 {
            return false; // Safety check against floating point errors
        }

        // Calculate proportional energy transfer
        let energy_to_transfer = if self.volume_layers[surface_layer] > 0.0 {
            // Transfer energy proportional to volume being moved
            let energy_density = self.energy_layers[surface_layer] / self.volume_layers[surface_layer];
            energy_density * volume
        } else {
            0.0 // No energy if no volume
        };

        // Validation: source must have sufficient energy
        if self.energy_layers[surface_layer] < energy_to_transfer {
            return false; // Cannot transfer more energy than available
        }

        // Perform conservative binary transfer on surface layer
        // Source loses volume and energy
        self.volume_layers[surface_layer] -= volume;
        self.energy_layers[surface_layer] -= energy_to_transfer;

        // Target gains volume and energy
        other.volume_layers[surface_layer] += volume;
        other.energy_layers[surface_layer] += energy_to_transfer;

        // Ensure no negative values due to floating point precision
        self.volume_layers[surface_layer] = self.volume_layers[surface_layer].max(0.0);
        self.energy_layers[surface_layer] = self.energy_layers[surface_layer].max(0.0);
        other.volume_layers[surface_layer] = other.volume_layers[surface_layer].max(0.0);
        other.energy_layers[surface_layer] = other.energy_layers[surface_layer].max(0.0);

        true // Transfer successful
    }

    /// Move a percentage of volume to another cell
    /// Returns true if transfer occurred, false if no transfer
    pub fn transfer_volume_fraction(&mut self, fraction: f64, other: &mut AsthenosphereCell) -> bool {
        if fraction <= 0.0 || fraction > 1.0 {
            return false; // Invalid fraction
        }

        let surface_layer = LAYER_COUNT - 1;
        let volume_to_move = self.volume_layers[surface_layer] * fraction;
        self.transfer_volume(volume_to_move, other)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{CELL_JOULES_START, EARTH};
    use h3o::Resolution;

    #[test]
    fn test_asthenosphere_cell_init() {
        let planet = EARTH.clone();

        let mut cell_count = 0;

        AsthenosphereCell::initial_cells_for_planet(CellsForPlanetArgs {
            planet: planet.clone(),
            joules_per_km3: 200.0,
            on_cell: |asth_cell| {
                assert_eq!(asth_cell.id.resolution(), Resolution::Two);
                assert!(!asth_cell.neighbors.is_empty(), "Cell has no neighbors");
                assert!(
                    !asth_cell.neighbors.contains(&asth_cell.id),
                    "Cell neighbors include itself"
                );
                assert!(asth_cell.volume_layers[0] > 0.0, "Cell volume layer 0 should be positive");
                assert!(asth_cell.energy_layers[0] > CELL_JOULES_START / 2.0, "Cell energy layer 0 should be above minimum");

                cell_count += 1;
                return asth_cell;
            },
            res: Resolution::Two,
            seed: 42,
        });

        assert_eq!(cell_count, 5882, "wrong cell count");
    }
}
