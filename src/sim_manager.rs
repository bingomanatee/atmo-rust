use crate::planet::{Planet, PlanetParams};
use crate::plate_generator::{GenerateRadiiParams, PartialPlateGenConfig, PlateGenerator};
use crate::rock_store::RockStore;
use crate::sim::{Sim, SimPlanetParams};
use rocksdb::LogLevel::Error;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use uuid::Uuid;
use crate::plate::Plate;

#[derive(Debug)]
pub enum SimManagerError {
    RocksDb(rocksdb::Error),
    Bincode(bincode::Error),
}

impl std::fmt::Display for SimManagerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SimManagerError::RocksDb(e) => write!(f, "RocksDB error: {}", e),
            SimManagerError::Bincode(e) => write!(f, "Bincode error: {}", e),
        }
    }
}

impl std::error::Error for SimManagerError {}

impl From<rocksdb::Error> for SimManagerError {
    fn from(e: rocksdb::Error) -> Self {
        SimManagerError::RocksDb(e)
    }
}

impl From<bincode::Error> for SimManagerError {
    fn from(e: bincode::Error) -> Self {
        SimManagerError::Bincode(e)
    }
}

#[derive(Serialize, Deserialize, Clone)]
pub struct SimExportData {
    pub sim: Sim,
    pub planet: Planet,
    pub plates: Vec<Plate>,
}

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
                let store =
                    RockStore::open(&db_path).expect("failed to open RocksDB for SimManager");
                let sim = store
                    .get_sim(&sim_id)
                    .expect("failed to load Sim")
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
                let mut store =
                    RockStore::open(&db_path).expect("failed to open RocksDB for SimManager");
                let mut sim = Sim::new();
                let planet = sim.generate_planet(planet_config);
                store.put_planet(&planet).expect("failed to save Planet");
                store.put_sim(&sim).expect("failed to save Sim");

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

    pub fn each_plate<F, E>(&self, mut callback: F) -> Result<(), E>
    where
        F: FnMut(Plate, Box<[u8]>) -> Result<(), E>,
        E: From<rocksdb::Error> + From<bincode::Error>,
    {
        self.store.each_plate(callback)
    }

    /**
    note: this method was formed when there was fewer plates.
    Going forward avoid using this method
    - prefer each_plate as it pulls data one at a time from memory.
    */
    pub fn plates(&self) -> Result<Vec<Plate>, String> {
        let mut plates = Vec::new();
        self.each_plate::<_, SimManagerError>(|plate, _| {
            plates.push(plate);
            Ok(())
        })
        .map_err(|e| format!("Failed to iterate plates: {}", e))?;
        Ok(plates)
    }

    pub fn make_plates(&self, target_coverage: f64) {
        // @deprecated - going forward we wil use the asthenosphere to generate plates
        let planet = self.planet().expect("Failed to load planet");

        let mut generator = PlateGenerator::new(
            PartialPlateGenConfig {
                target_coverage: Some(target_coverage),
                power_law_exponent: Some(1.3),
                min_density: None,
                max_density: None,
                min_thickness_km: None,
                max_thickness_km: None,
                variation_factor: None,
                max_plate_radius_radians: None,
            },
            &planet,
        );
        let radii = generator.generate_radii(GenerateRadiiParams {
            target_coverage: target_coverage as f64,
            min_radius: (planet.radius_km / 20),
            max_radius: (planet.radius_km as f64 * 0.8) as i32,
            exponent: 0.3,
        });

        let mut plate_ids: HashSet<Uuid> = HashSet::new();
        for radius in &radii {
            // @TODO: put in a different index
            let plate = generator.generate_one(*radius, planet.id);
            let _ = self.store.put_plate(&plate);
            plate_ids.insert(plate.id);
        }

        let mut planet = self.planet().expect("cannot retrieve planet"); // reloading for integrity
        self.store.put_planet(&planet).expect("Cannot put planet");
    }

    /// Export simulation data to JSON format
    pub fn export_data(&self) -> Result<SimExportData, String> {
        let sim = self.sim()?;
        let planet = self.planet()?;
        let mut plates = Vec::new();

        /**
                note - this system now writes all plates for all _STEPS_ in a single collection
                there may be a LOT More data than is expected in this loop!
        */
        self.each_plate::<_, SimManagerError>(|plate, _| {
            plates.push(plate);
            Ok(())
        })
        .map_err(|e| format!("Failed to iterate plates: {}", e))?;

        Ok(SimExportData {
            sim,
            planet,
            plates,
        })
    }

    /// Save simulation data to a JSON file
    pub fn save_to_json<P: AsRef<Path>>(&self, path: P) -> Result<(), String> {
        let export_data = self.export_data()?;
        let json_string = serde_json::to_string_pretty(&export_data)
            .map_err(|e| format!("Failed to serialize to JSON: {}", e))?;

        let mut file = File::create(path).map_err(|e| format!("Failed to create file: {}", e))?;

        file.write_all(json_string.as_bytes())
            .map_err(|e| format!("Failed to write to file: {}", e))?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::sim::SimPlanetParams;
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
        let sim = manager.sim().expect("Failed to load sim");

        let planet = manager.planet().expect("Failed to load planet");
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
            let _ = SimManagerParams::Load {
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

    #[test]
    fn make_plates_creates_and_stores_plates() {
        let dir = tempdir().expect("create temp dir");
        let db_path = dir.path().join("testdb_plates");
        let db_path_str = db_path.to_str().unwrap().to_string();

        let params = SimManagerParams::Create {
            db_path: db_path_str.clone(),
            planet_config: make_test_planet_config(),
        };
        let manager = SimManager::new(params);

        // Call make_plates with a target coverage
        let target_coverage = 0.5;
        manager.make_plates(target_coverage);

        // Load the planet to get its id
        let planet = manager.planet().expect("Failed to load planet");
        let planet_id = planet.id;

        // Retrieve plates from the store for this planet
        let plates = manager.plates().expect("Failed to get plates for planet");

        // Assert that some plates were created and stored
        assert!(!plates.is_empty(), "No plates were created");

        // Optionally, check that plate radii are within expected bounds
        for plate in plates {
            assert_eq!(plate.planet_id, planet_id);
            assert!(plate.radius_km > 0);
        }
    }

    #[test]
    fn save_to_json_saves_sim_data() {
        let dir = tempdir().expect("create temp dir");
        let db_path = dir.path().join("testdb_json");
        let db_path_str = db_path.to_str().unwrap().to_string();

        let params = SimManagerParams::Create {
            db_path: db_path_str.clone(),
            planet_config: make_test_planet_config(),
        };
        let manager = SimManager::new(params);

        let output_path = dir
            .path()
            .join("sim_data.json")
            .to_str()
            .unwrap()
            .to_string();
        manager
            .save_to_json(&output_path)
            .expect("Failed to save sim data");
    }
}
