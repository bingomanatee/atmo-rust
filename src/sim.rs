use serde::{Serialize, Deserialize};
use uuid::Uuid;
use crate::planet::{Planet, PlanetParams};

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Sim {
    pub id: Uuid,
    pub planet_id: Option<Uuid>,
}

impl Sim {
    pub fn new() -> Self {
        Sim {
            id: Uuid::new_v4(),
            planet_id: None,
        }
    }

    pub fn generate_planet(&mut self, radius: i32) -> Planet {
        if self.planet_id.is_some() {
            panic!("Sim already has a planet id");
        }
        let planet = Planet::new(PlanetParams {
            sim_id: self.id,
            radius,
            mantle_density_gcm3: None
        });
        self.planet_id = Option::from(planet.id);
        planet
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use uuid::Uuid;
    use crate::planet::EARTH_RADIUS_KM;

    #[test]
    fn test_sim_new() {
        let sim = Sim::new();
        assert_ne!(sim.id, Uuid::nil());
        assert!(sim.planet_id.is_none());
    }

    #[test]
    fn test_generate_planet_success() {
        let mut sim = Sim::new();
        let planet = sim.generate_planet(EARTH_RADIUS_KM);

        // The planet's owner id should be the sim's id
        assert_eq!(planet.sim_id, sim.id);
        // The planet's radius should be as specified
        assert_eq!(planet.radius, EARTH_RADIUS_KM);
        // The sim's planet_id should now be set to the planet's id
        assert_eq!(sim.planet_id, Some(planet.id));
    }

    #[test]
    #[should_panic(expected = "Sim already has a planet id")]
    fn test_generate_planet_already_has_planet() {
        let mut sim = Sim::new();
        sim.generate_planet(EARTH_RADIUS_KM /2);
        sim.generate_planet(EARTH_RADIUS_KM * 2);
    }
}
