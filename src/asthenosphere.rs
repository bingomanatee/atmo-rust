use crate::constants::{
    AVG_STARTING_VOLUME_KM_3, JOULES_PER_KM3, MAX_SUNK_PERCENT, MAX_SUNK_TEMP, MIN_SUNK_TEMP,
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
    pub volume: f64,   // c. 2 MIO KM3 per
    pub energy_j: f64, // Temp in kelvin; each unit of volume adds 2073, cools by 20 per MIO years
    pub anomaly_volume: f64, // Additional volume from anomalies (flows)
}

impl Default for AsthenosphereCell {
    fn default() -> Self {
        let cell = CellIndex::base_cells().last().expect("no last cell");
        AsthenosphereCell {
            id: cell,
            volume: 0.0,
            energy_j: 0.0,
            neighbors: Vec::new(),
            step: 0,
            anomaly_volume: 0.0,
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
    pub anomaly_freq: f64, // Frequency of initial anomalies (0.0 to 1.0)
}

#[derive(Clone)]
struct AsthenosphereCellParams {
    cell: CellIndex,
    volume_kg_3: f64,
    energy_per_volume: f64,
}

impl AsthenosphereCell {
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
            let noise_scale = 25.0;
            let scaled_location = location.normalize() * noise_scale as f32;
            let noise_val = per.get(scaled_location.to_array().map(|n| n as f64));

            let exponential_noise = if noise_val > 0.0 {
                noise_val.powf(2.0) // Less extreme exponent for smoother variation
            } else {
                -((-noise_val).powf(2.0))
            };

            // Much more extreme variation: 0.3x to 2.5x the base volume
            let random_scale = 0.3 + (exponential_noise + 1.0) * 0.7; // Maps [-1,1] to [0.3, 2.5]
            let volume = AVG_STARTING_VOLUME_KM_3 * random_scale as f64;
            
            // Add independent energy variation using different noise parameters
            let energy_noise_val = per.get([
                (scaled_location.x * 1.7) as f64, 
                (scaled_location.y * 1.7) as f64, 
                (scaled_location.z * 1.7) as f64
            ]);
            
            let energy_exponential_noise = if energy_noise_val > 0.0 {
                energy_noise_val.powf(1.5)
            } else {
                -((-energy_noise_val).powf(1.5))
            };
            
            // Independent energy variation: 0.2x to 3.0x base energy per volume
            let energy_scale = 0.2 + (energy_exponential_noise + 1.0) * 0.8; // Maps [-1,1] to [0.2, 3.0]
            let adjusted_energy_per_volume = energy_per_volume * energy_scale;
            
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
        if self.energy_j as f64 >= MAX_SUNK_TEMP {
            return 0.0;
        }
        let effective_energy: f64 = (self.energy_j as f64).max(MIN_SUNK_TEMP);
        let percent = MAX_SUNK_PERCENT * effective_energy / (MAX_SUNK_TEMP / MIN_SUNK_TEMP);
        self.volume * percent
    }

    pub fn at_cell(params: AsthenosphereCellParams) -> AsthenosphereCell {
        AsthenosphereCell {
            id: params.cell,
            volume: params.volume_kg_3,
            energy_j: params.energy_per_volume * params.volume_kg_3,
            step: 0,
            neighbors: params
                .cell
                .grid_disk::<Vec<_>>(1)
                .into_iter()
                .filter(|&c| c != params.cell)
                .collect::<Vec<CellIndex>>(),
            anomaly_volume: 0.0,
        }
    }

    pub fn location(&self, planet: &Planet, res: Option<Resolution>) -> Vec3 {
        let resolution = res.unwrap_or(ASTH_RES);
        let gc = GeoCellConverter::new(planet.radius_km as f64, resolution);
        gc.cell_to_vec3(self.id)
    }

    /// Move volume and proportional energy to another cell (conservative binary operation)
    /// Returns true if transfer occurred, false if no transfer (invalid volume or insufficient source)
    pub fn transfer_volume(&mut self, volume: f64, other: &mut AsthenosphereCell) -> bool {
        // Validation: volume must be positive
        if volume <= 0.0 {
            return false; // Silent return for zero/negative volume
        }

        // Validation: source must have sufficient volume
        if self.volume < volume {
            return false; // Cannot transfer more than available
        }

        // Validation: source would not go below zero (additional safety check)
        let remaining_volume = self.volume - volume;
        if remaining_volume < 0.0 {
            return false; // Safety check against floating point errors
        }

        // Calculate proportional energy transfer
        let energy_to_transfer = if self.volume > 0.0 {
            // Transfer energy proportional to volume being moved
            let energy_density = self.energy_j / self.volume;
            energy_density * volume
        } else {
            0.0 // No energy if no volume
        };

        // Validation: source must have sufficient energy
        if self.energy_j < energy_to_transfer {
            return false; // Cannot transfer more energy than available
        }

        // Perform conservative binary transfer
        // Source loses volume and energy
        self.volume -= volume;
        self.energy_j -= energy_to_transfer;

        // Target gains volume and energy
        other.volume += volume;
        other.energy_j += energy_to_transfer;

        // Ensure no negative values due to floating point precision
        self.volume = self.volume.max(0.0);
        self.energy_j = self.energy_j.max(0.0);
        other.volume = other.volume.max(0.0);
        other.energy_j = other.energy_j.max(0.0);

        true // Transfer successful
    }

    /// Move a percentage of volume to another cell
    /// Returns true if transfer occurred, false if no transfer
    pub fn transfer_volume_fraction(&mut self, fraction: f64, other: &mut AsthenosphereCell) -> bool {
        if fraction <= 0.0 || fraction > 1.0 {
            return false; // Invalid fraction
        }

        let volume_to_move = self.volume * fraction;
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
                assert!(asth_cell.volume > 0.0, "Cell volume should be positive");
                assert!(asth_cell.energy_j > CELL_JOULES_START / 2.0);

                cell_count += 1;
                return asth_cell;
            },
            res: Resolution::Two,
            seed: 42
        });

        assert_eq!(cell_count, 5882, "wrong cell count");
    }
}
