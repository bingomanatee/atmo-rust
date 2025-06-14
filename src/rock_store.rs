use rocksdb::{DB, Options};
use uuid::Uuid;
use bincode;
use crate::planet::Planet;
use crate::plate::{Plate};
use crate::sim::Sim;

pub struct RockStore {
    db: DB,
}

impl RockStore {
    pub fn open(path: &str) -> Result<Self, rocksdb::Error> {
        let mut opts = Options::default();
        opts.create_if_missing(true);
        let db = DB::open(&opts, path)?;
        Ok(RockStore { db })
    }

    // Helper to generate keys with prefixes
    fn key_sim(id: &Uuid) -> Vec<u8> {
        let mut key = b"sim:".to_vec();
        key.extend(id.as_bytes());
        key
    }

    fn key_planet(id: &Uuid) -> Vec<u8> {
        let mut key = b"planet:".to_vec();
        key.extend(id.as_bytes());
        key
    }

    fn key_plate(id: &Uuid) -> Vec<u8> {
        let mut key = b"plate:".to_vec();
        key.extend(id.as_bytes());
        key
    }

    // Save Sim
    pub fn put_sim(&self, sim: &Sim) -> Result<(), rocksdb::Error> {
        let key = Self::key_sim(&sim.id);
        let value = bincode::serialize(sim).unwrap();
        self.db.put(key, value)
    }

    // Get Sim
    pub fn get_sim(&self, id: &Uuid) -> Result<Option<Sim>, rocksdb::Error> {
        let key = Self::key_sim(id);
        match self.db.get(key)? {
            Some(value) => {
                let sim: Sim = bincode::deserialize(&value).unwrap();
                Ok(Some(sim))
            }
            None => Ok(None),
        }
    }

    // Save Planet
    pub fn put_planet(&self, planet: &Planet) -> Result<(), rocksdb::Error> {
        let key = Self::key_planet(&planet.id);
        let value = bincode::serialize(planet).unwrap();
        self.db.put(key, value)
    }

    // Get Planet
    pub fn get_planet(&self, id: &Uuid) -> Result<Option<Planet>, rocksdb::Error> {
        let key = Self::key_planet(id);
        match self.db.get(key)? {
            Some(value) => {
                let planet: Planet = bincode::deserialize(&value).unwrap();
                Ok(Some(planet))
            }
            None => Ok(None),
        }
    }

    // Save Plate
    pub fn put_plate(&self, plate: &Plate) -> Result<(), rocksdb::Error> {
        let key = Self::key_plate(&plate.id);
        let value = bincode::serialize(plate).unwrap();
        self.db.put(key, value)
    }

    // Get Plate
    pub fn get_plate(&self, id: &Uuid) -> Result<Option<Plate>, rocksdb::Error> {
        let key = Self::key_plate(id);
        match self.db.get(key)? {
            Some(value) => {
                let plate: Plate = bincode::deserialize(&value).unwrap();
                Ok(Some(plate))
            }
            None => Ok(None),
        }
    }
}
