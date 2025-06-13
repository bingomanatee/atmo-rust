use crate::planet::{Planet, PlanetParams};
use serde::{Deserialize, Serialize};
use uuid::Uuid;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Sim {
    pub id: Uuid,
    pub planet_id: Option<Uuid>,
    pub sim_id: Uuid,
}

pub struct SimPlanetParams {
    pub radius: i32,
    pub mantle_density_gcm3: Option<f64>,
}

impl SimPlanetParams {
    pub fn into_planet_params(self, sim_id: Uuid) -> PlanetParams {
        PlanetParams {
            sim_id,
            radius: self.radius,
            mantle_density_gcm3: self.mantle_density_gcm3,
            plate_ids: None
        }
    }
}

impl Sim {
    pub fn new() -> Self {
        Sim {
            id: Uuid::new_v4(),
            planet_id: None,
            sim_id: Uuid::new_v4(),
        }
    }

    pub fn generate_planet(&mut self, params: SimPlanetParams) -> Planet {
        if self.planet_id.is_some() {
            panic!("Sim already has a planet id");
        }

        let new_params = params.into_planet_params(self.id);

        let planet = Planet::new(new_params);
        self.planet_id = Option::from(planet.id);
        planet
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::planet::EARTH_RADIUS_KM;
    use uuid::Uuid;

    #[test]
    fn test_sim_new() {
        let sim = Sim::new();
        assert_ne!(sim.id, Uuid::nil());
        assert!(sim.planet_id.is_none());
    }

    #[test]
    fn test_generate_planet_success() {
        let mut sim = Sim::new();
        let planet = sim.generate_planet(SimPlanetParams {
            radius: EARTH_RADIUS_KM,
            mantle_density_gcm3: None
        });

        // The planet's owner id should be the sim's id
        assert_eq!(planet.sim_id, sim.id);
        // The planet's radius should be as specified
        assert_eq!(planet.radius_km, EARTH_RADIUS_KM);
        // The sim's planet_id should now be set to the planet's id
        assert_eq!(sim.planet_id, Some(planet.id));
    }

    #[test]
    #[should_panic(expected = "Sim already has a planet id")]
    fn test_generate_planet_already_has_planet() {
        let mut sim = Sim::new();
        sim.generate_planet(SimPlanetParams {
            radius: EARTH_RADIUS_KM / 2,
            mantle_density_gcm3: None
        });
        sim.generate_planet(SimPlanetParams {
            radius: EARTH_RADIUS_KM * 2,
            mantle_density_gcm3: None
        });
    }
}
