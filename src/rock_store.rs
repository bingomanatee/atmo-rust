use rocksdb::{DB, Options, WriteBatch};
use uuid::Uuid;
use bincode;
use h3o::CellIndex;
use crate::asthenosphere::AsthenosphereCell;
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
    fn key_asth(id: &CellIndex, step: u32) -> Vec<u8> {
        let mut key = b"asth:".to_vec();
        key.extend_from_slice(id.to_string().as_bytes());
        key.extend_from_slice(b":");
        key.extend_from_slice(step.to_string().as_bytes());
        key
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

    // Get AsthenosphereCell
    pub fn get_asth(&self, cell: CellIndex, step: u32) -> Result<Option<AsthenosphereCell>, rocksdb::Error> {
        let key = Self::key_asth(&cell, step);
        match self.db.get(key)? {
            Some(value) => {
                let a: AsthenosphereCell = bincode::deserialize(&value).unwrap();
                Ok(Some(a))
            }
            None => Ok(None),
        }
    }
    // Save AsthenosphereCell
    pub fn put_asth(&self, a: &AsthenosphereCell) -> Result<(), rocksdb::Error> {
        let key = Self::key_asth(&a.cell, a.step);
        let value = bincode::serialize(a).unwrap();
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

    pub fn each_plate<F, E>(&self, mut callback: F) -> Result<(), E>
    where
        F: FnMut(Plate, Box<[u8]>) -> Result<(), E>,
        E: From<rocksdb::Error> + From<bincode::Error>,
    {
        let prefix = b"plate:";
        let iter = self.db.prefix_iterator(prefix);

        for item in iter {
            let (key, value) = item.map_err(E::from)?;
            if !key.starts_with(prefix) {
                break;
            }
            let plate: Plate = bincode::deserialize(&value).map_err(E::from)?;
            callback(plate, key)?;
        }
        Ok(())
    }

    pub fn put_asth_batch(&self, cells: &[AsthenosphereCell]) {
        
      let mut batch = WriteBatch::default();
        for cell in cells {
            let key = RockStore::key_asth(&cell.cell, cell.step);
           let value = bincode::serialize(cell).expect("serialize failed");
           // self.put_asth( &cell).unwrap();
            batch.put(key, value);
        }
        self.db.write(batch);
    }
}
