use crate::constants::{
    AVG_STARTING_VOLUME_KM_3, JOULES_PER_KM3, MAX_SUNK_PERCENT, MAX_SUNK_TEMP, MIN_SUNK_TEMP,
};
use crate::geoconverter::GeoCellConverter;
use crate::h3_utils::H3Utils;
use crate::planet::Planet;
use glam::Vec3;
use h3o::{CellIndex, Resolution};
use serde::{Deserialize, Serialize};
use crate::perlin_noise_generator::PerlinNoiseGenerator;

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
    pub anomaly_energy: f64, // Additional energy from anomalies (flows)
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
            anomaly_energy: 0.0,
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
        
        // Create simple noise generators for volume and energy variation
        let volume_noise = PerlinNoiseGenerator::new(seed, 3.0, 1.0);      // Moderate detail and sharpness
        let energy_noise = PerlinNoiseGenerator::new(seed + 1, 4.0, 1.5);  // Slightly different settings
        
        let mut index = 0u32;
        H3Utils::iter_at(res, |l2_cell| {
            let location = gc.cell_to_vec3(l2_cell).normalize();
            
            // Get noise values [-1, 1] and map to ranges centered around 1.0
            let volume_noise_val = volume_noise.sample(location);
            let energy_noise_val = energy_noise.sample(location);
            
            // Map volume noise [-1, 1] to [0.8, 1.2] centered around 1.0 (80% to 120% of base)
            let volume_scale = 1.0 + volume_noise_val * 0.2; // Maps [-1,1] to [0.8, 1.2]
            
            // Map energy noise [-1, 1] to [0.7, 1.3] centered around 1.0 (70% to 130% of base)  
            let energy_scale = 1.0 + energy_noise_val * 0.3; // Maps [-1,1] to [0.7, 1.3]
            
            let volume = AVG_STARTING_VOLUME_KM_3 * volume_scale;
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
            anomaly_energy: 0.0,
            anomaly_volume: 0.0,
        }
    }

    pub fn location(&self, planet: &Planet, res: Option<Resolution>) -> Vec3 {
        let resolution = res.unwrap_or(ASTH_RES);
        let gc = GeoCellConverter::new(planet.radius_km as f64, resolution);
        gc.cell_to_vec3(self.id)
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
