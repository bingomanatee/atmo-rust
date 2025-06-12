use crate::planet::{Planet, PlanetParams};
use crate::rock_store::RockStore;
use crate::sim::{Sim, SimPlanetParams};
use uuid::Uuid;

pub struct SimManager {
    sim_id: Uuid,
    store: RockStore,
}

pub enum SimManagerParams {
    Create {
        db_path: String,
        planet_config: SimPlanetParams,
    },
    Load {
        db_path: String,
        sim_id: Uuid,
    },
}

impl SimManager {
    pub fn new(params: SimManagerParams) -> Self {
        match params {
            SimManagerParams::Load { db_path, sim_id } => {
                let mut store =
                    RockStore::open(&db_path).expect("failed to open RocksDB for SimManager");
                let sim = store
                    .get_sim(&sim_id)
                    .expect("DB error loading Sim")
                    .expect("No Sim found with given id");

                SimManager {
                    sim_id: sim.id,
                    store,
                }
            }
            SimManagerParams::Create {
                db_path,
                planet_config,
            } => {
                println!("-------- SimManager new ---------");
                let mut store =
                    RockStore::open(&db_path).expect("failed to open RocksDB for SimManager");
                let mut sim = Sim::new();
                let planet = sim.generate_planet(planet_config);
                store.put_planet(&planet).expect("failed to save Planet");
                store.put_sim(&sim).expect("failed to save Sim");

                println!("saved planet {} for sim {}", planet.id, sim.id);
                SimManager {
                    sim_id: sim.id,
                    store,
                }
            }
        }
    }
    pub fn sim(&self) -> Result<Sim, String> {
        match self.store.get_sim(&self.sim_id) {
            Err(e) => Err(format!("DB error loading Sim: {}", e)),
            Ok(None) => Err(format!("No Sim found with id {}", self.sim_id)),
            Ok(Some(sim)) => Ok(sim),
        }
    }

    pub fn planet(&self) -> Result<Planet, String> {
        let sim = self.sim()?;
        let pid = sim
            .planet_id
            .ok_or_else(|| "No Planet associated with this Sim".to_string())?;

        match self.store.get_planet(&pid) {
            Err(e) => Err(format!("DB error loading Planet: {}", e)),
            Ok(None) => Err(format!("No Planet found with id {}", pid)),
            Ok(Some(planet)) => Ok(planet),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::sim::{Sim, SimPlanetParams};
    use crate::sim_manager::{SimManager, SimManagerParams};
    use tempfile::tempdir;

    fn make_test_planet_config() -> SimPlanetParams {
        SimPlanetParams {
            radius: 42,
            mantle_density_gcm3: Some(3.3),
        }
    }

    #[test]
    fn sim_and_planet_are_persisted_in_store() {
        let dir = tempdir().expect("create temp dir");
        let db_path = dir.path().join("testdb2");
        let db_path_str = db_path.to_str().unwrap().to_string();

        let params = SimManagerParams::Create {
            db_path: db_path_str.clone(),
            planet_config: make_test_planet_config(),
        };
        let mut manager = SimManager::new(params);

        // part two: see that the stores have matching data

        let sim_id = manager.sim_id;
        println!(" laoding planet for sim {}", sim_id);
        let sim = manager.sim().expect("Failed to load sim");
        println!(" (sim has planet_id: {} )", sim.planet_id.unwrap());

        let planet = manager.planet().expect("Failed to laod planet");
        let planet_id = planet.id;

        let loaded_sim = manager
            .store
            .get_sim(&sim_id)
            .expect("DB error loading sim")
            .expect("Sim not found in store");
        assert_eq!(loaded_sim.id, sim_id);

        let loaded_planet = manager
            .store
            .get_planet(&planet_id)
            .expect("DB error loading planet")
            .expect("Planet not found in store");
        assert_eq!(loaded_planet.id, planet_id);

        assert_eq!(loaded_planet.sim_id, sim_id);
        assert_eq!(manager.sim().unwrap().planet_id, Some(planet_id));
    }

    #[test]
    fn load_existing_sim_and_planet() {
        let dir = tempdir().unwrap();
        let db = dir.path().join("db").to_string_lossy().to_string();

        let create_params = SimManagerParams::Create {
            db_path: db.clone(),
            planet_config: make_test_planet_config(),
        };
        let existing_id;
        let existing_pid;
        {
            let mgr0 = SimManager::new(create_params);
            existing_id = mgr0.sim_id;
            existing_pid = mgr0.planet().unwrap().id;
            let load_params = SimManagerParams::Load {
                db_path: db,
                sim_id: existing_id,
            };
        }

        let dir2 = tempdir().unwrap();
        let db2 = dir.path().join("db").to_string_lossy().to_string();

        let mgr1 = SimManager::new(SimManagerParams::Load {
            db_path: db2,
            sim_id: existing_id,
        });
        assert_eq!(mgr1.sim_id, existing_id);
        assert_eq!(mgr1.planet().unwrap().id, existing_pid);
    }
}
