use crate::geoconverter::GeoCellConverter;
use crate::h30_utils::H3Utils;
use crate::planet::Planet;
use glam::Vec3;
use h3o::{CellIndex, Resolution};
use noise::{NoiseFn, Perlin};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Clone)]
pub(crate) struct AsthenosphereCell {
    pub cell: CellIndex,
    pub location: Vec3,
    pub neighbors: Vec<CellIndex>,
    pub step: u32,
    pub volume: f32,   // c. 2 MIO KM3 per
    pub energy_k: u32, // Temp in kelvin; each unit of volume adds 2073, cools by 20 per MIO years
}

pub const ASTH_RES: Resolution = Resolution::Two;

#[derive(Clone)]
pub struct CellsForPlanetArgs<F>
where
    F: FnMut(AsthenosphereCell),
{
    pub planet: Planet,
    pub on_cell: F,
    pub res: Option<Resolution>,
}

impl AsthenosphereCell {
    /**
        this utility method creates a series of Asthenosphere cells for each cell in an h3 grid
        given a planet and resolution and feed it back to the callback;
        it is a utility for asth_sim.
    */
    pub fn cells_for_planet<F>(args: CellsForPlanetArgs<F>)
    where
        F: FnMut(AsthenosphereCell),
    {
        let CellsForPlanetArgs {
            planet,
            mut on_cell,
            res,
        } = args;
        let resolution = res.unwrap_or(ASTH_RES);
        let gc = GeoCellConverter::new(planet.radius_km as f64, resolution);
        let per = Perlin::new(100);
        H3Utils::iter_at(ASTH_RES, |l2_cell| {
            let location = gc.cell_to_vec3(l2_cell);
            let [x, y, z] = location.normalize().to_array();
            let noise_val = per.get([x as f64, y as f64, z as f64]);
            on_cell(AsthenosphereCell::at_cell(
                l2_cell,
                location,
                (10.0 + (noise_val * 10.0)) as f32,
            ));
        });
    }

    pub fn sunk_volume(&self) -> f32 {
        if self.energy_k as f32 >= MAX_SUNK_TEMP {
            return 0.0;
        }
        let effective_energy: f32 = (self.energy_k as f32).max(MIN_SUNK_TEMP);
        let percent = MAX_SUNK_PERCENT * effective_energy / MAX_SUNK_TEMP;
        self.volume * percent
    }

    pub fn at_cell(cell: CellIndex, location: Vec3, count: f32) -> AsthenosphereCell {
        AsthenosphereCell {
            cell,
            volume: count,
            energy_k: (16000.0 * count as f32) as u32,
            location,
            step: 0,
            neighbors: cell
                .grid_disk::<Vec<_>>(1)
                .into_iter()
                .filter(|&c| c != cell)
                .collect::<Vec<CellIndex>>(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::planet::EARTH;
    use h3o::Resolution;

    #[test]
    fn test_asthenosphere_cell_init() {
        // Use the predefined Earth planet constant
        let planet = EARTH.clone();

        let mut cell_count = 0;

        AsthenosphereCell::cells_for_planet(CellsForPlanetArgs {
            planet: planet.clone(),
            on_cell: |cell| {
                // Assertions inside the callback
            assert_eq!(cell.cell.resolution(), Resolution::Two);
            assert!(!cell.neighbors.is_empty(), "Cell has no neighbors");
            assert!(
                !cell.neighbors.contains(&cell.cell),
                "Cell neighbors include itself"
            );
            assert!(cell.volume > 0.0, "Cell volume should be positive");
            assert!(cell.energy_k > 10000);

                cell_count += 1;
            },
            res: None,
        });

        assert_eq!(cell_count, 5882, "wrong cell count");
    }
}

const MAX_SUNK_PERCENT: f32 = 0.15;
const MAX_SUNK_TEMP: f32 = 2000.0;
const MIN_SUNK_TEMP: f32 = 1000.0;
