use crate::geoconverter::GeoCellConverter;
use crate::h30_utils::H3Utils;
use crate::planet::Planet;
use glam::Vec3;
use h3o::{CellIndex, Resolution};
use noise::{NoiseFn, Perlin};
use serde::{Deserialize, Serialize};
use crate::asth_constants::{AVG_STARTING_VOLUME, K_PER_VOLUME, MAX_SUNK_PERCENT, MAX_SUNK_TEMP, MIN_SUNK_TEMP};

/**
some assumptions: the meaningful volume of the asthenosphere we track is a 10km
at l2, the area is 426 km^2, the volume of the average cell is 4260 km^3
*/
#[derive(Serialize, Deserialize, Clone, Debug)]
pub(crate) struct AsthenosphereCell {
    pub cell: CellIndex,
    pub neighbors: Vec<CellIndex>,
    pub step: u32,
    pub volume: f64,   // c. 2 MIO KM3 per
    pub energy_k: f64, // Temp in kelvin; each unit of volume adds 2073, cools by 20 per MIO years
}

impl Default for AsthenosphereCell {
    fn default() -> Self {
        let cell= CellIndex::base_cells().last().expect("no last cell");
        AsthenosphereCell {
            cell,
            volume: 0.0,
            energy_k: 0.0,
            neighbors: Vec::new(),
            // initialize other fields with sensible defaults
            step: 0,
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
    pub res: Option<Resolution>,
    pub energy_per_volume: f64,
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
    pub fn cells_for_planet<F>(args: CellsForPlanetArgs<F>)
    where
        F: FnMut(AsthenosphereCell) -> AsthenosphereCell,
    {
        let CellsForPlanetArgs {
            planet,
            mut on_cell,
            res,
            energy_per_volume,
        } = args;
        let resolution = res.unwrap_or(ASTH_RES);
        let gc = GeoCellConverter::new(planet.radius_km as f64, resolution);
        // Use a different seed for initial generation to avoid conflicts with step-based seeds
        let per = Perlin::new(42);
        let mut index = 0u32;
        H3Utils::iter_at(ASTH_RES, |l2_cell| {
            let location = gc.cell_to_vec3(l2_cell);
            let [x, y, z] = location.normalize().to_array();

            // Enhanced Perlin noise with much higher frequency for very detailed initial patterns
            let noise_scale = 25.0; // Much higher scale for very detailed initial patterns
            let scaled_x = x as f64 * noise_scale;
            let scaled_y = y as f64 * noise_scale;
            let scaled_z = z as f64 * noise_scale;
            let noise_val = per.get([scaled_x, scaled_y, scaled_z]);

            // Use exponential function to create more extreme spikes
            // This creates sharp peaks and valleys instead of smooth gradients
            let exponential_noise = if noise_val > 0.0 {
                noise_val.powf(3.0) // Cube positive values for sharp peaks
            } else {
                -((-noise_val).powf(3.0)) // Cube negative values for sharp valleys
            };

            let random_scale = 1.0 + (exponential_noise / 6.0); // Slightly reduced range: ~0.83...1.17
            let volume = AVG_STARTING_VOLUME * random_scale as f64;
            let asth_cell = AsthenosphereCell::at_cell(AsthenosphereCellParams {
                cell: l2_cell,
                volume_kg_3: volume,
                energy_per_volume,
            });

            on_cell(asth_cell);
            index += 1;
        });
    }

    pub fn sunk_volume(&self) -> f64 {
        if self.energy_k as f64 >= MAX_SUNK_TEMP {
            return 0.0;
        }
        let effective_energy: f64 = (self.energy_k as f64).max(MIN_SUNK_TEMP);
        let percent = MAX_SUNK_PERCENT * effective_energy / (MAX_SUNK_TEMP / MIN_SUNK_TEMP);
        self.volume * percent
    }

    pub fn at_cell(params: AsthenosphereCellParams) -> AsthenosphereCell {
        AsthenosphereCell {
            cell: params.cell,
            volume: params.volume_kg_3,
            energy_k: K_PER_VOLUME * params.volume_kg_3, // should sum to around 2000k at start
            step: 0,
            neighbors: params
                .cell
                .grid_disk::<Vec<_>>(1)
                .into_iter()
                .filter(|&c| c != params.cell)
                .collect::<Vec<CellIndex>>(),
        }
    }

    /// Dynamically compute the location vector for this cell using the given planet and resolution.
    pub fn location(&self, planet: &Planet, res: Option<Resolution>) -> Vec3 {
        let resolution = res.unwrap_or(ASTH_RES);
        let gc = GeoCellConverter::new(planet.radius_km as f64, resolution);
        gc.cell_to_vec3(self.cell)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::planet::EARTH;
    use h3o::Resolution;
    use crate::asth_constants::CELL_ENERGY_START;

    #[test]
    fn test_asthenosphere_cell_init() {
        // Use the predefined Earth planet constant
        let planet = EARTH.clone();

        let mut cell_count = 0;

        AsthenosphereCell::cells_for_planet(CellsForPlanetArgs {
            planet: planet.clone(),
            energy_per_volume: 200.0, // average cell starts with
            on_cell: |asth_cell| {
                // Assertions inside the callback
                assert_eq!(asth_cell.cell.resolution(), Resolution::Two);
                assert!(!asth_cell.neighbors.is_empty(), "Cell has no neighbors");
                assert!(
                    !asth_cell.neighbors.contains(&asth_cell.cell),
                    "Cell neighbors include itself"
                );
                assert!(asth_cell.volume > 0.0, "Cell volume should be positive");
                assert!(asth_cell.energy_k > CELL_ENERGY_START/2.0);

                cell_count += 1;
                return asth_cell
            },
            res: None,
        });

        assert_eq!(cell_count, 5882, "wrong cell count");
    }
}
