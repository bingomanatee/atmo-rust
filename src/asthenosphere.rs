use crate::geoconverter::GeoCellConverter;
use crate::h30_utils::H3Utils;
use crate::planet::Planet;
use glam::Vec3;
use h3o::{CellIndex, Resolution};
use noise::{NoiseFn, Perlin};
use serde::{Deserialize, Serialize};
use crate::asth_constants::{AVG_STARTING_VOLUME, K_PER_VOLUME, MAX_SUNK_PERCENT, MAX_SUNK_TEMP, MIN_SUNK_TEMP};

/**
some assumptions: the meaningfful volume of the aesthenosphere we track is a 10km
at l2, the area is 426 km^2, the volume of the average cell is 4260 km^3
*/
#[derive(Serialize, Deserialize, Clone, Debug)]
pub(crate) struct AsthenosphereCell {
    pub cell: CellIndex,
    pub location: Vec3,
    pub neighbors: Vec<CellIndex>,
    pub step: u32,
    pub volume: f64,   // c. 2 MIO KM3 per
    pub energy_k: f64, // Temp in kelvin; each unit of volume adds 2073, cools by 20 per MIO years
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
    location: Vec3,
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
        let per = Perlin::new(100);
        let mut index = 0u32;
        H3Utils::iter_at(ASTH_RES, |l2_cell| {
            let location = gc.cell_to_vec3(l2_cell);
            let [x, y, z] = location.normalize().to_array();
            let noise_val = per.get([x as f64, y as f64, z as f64]);
            let random_scale = 1.0 +  (noise_val  / 4.0); // 0.5...1.5
            let count = AVG_STARTING_VOLUME * random_scale as f64;
            let asth_cell = AsthenosphereCell::at_cell(AsthenosphereCellParams {
                cell: l2_cell,
                location,
                volume_kg_3: count,
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
            energy_k: K_PER_VOLUME * params.volume_kg_3, // average starting cound is 10
            // avg starting energy is 2000K
            location: params.location,
            step: 0,
            neighbors: params
                .cell
                .grid_disk::<Vec<_>>(1)
                .into_iter()
                .filter(|&c| c != params.cell)
                .collect::<Vec<CellIndex>>(),
        }
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

// see notes in asthenosphere.rs
