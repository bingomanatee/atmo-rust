use glam::Vec3;
use serde::{Serialize, Deserialize};
use uuid::Uuid;

pub const PLATELET_RESOLUTION: u8 = 3;

#[derive(Serialize, Deserialize, Clone)]
pub struct Plate {
    pub id: Uuid,
    pub center: Vec3,
    pub radius_km: i32,
    pub thickness_km: i32,
    pub density: f64,
    pub planet_id: Uuid,
}

#[derive(Clone)]
pub struct PlateParams {
    pub center: Vec3,
    pub radius_km: i32,
    pub thickness_km: i32,
    pub density: f64,
    pub planet_id: Uuid,
}

impl Plate {
    /// Creates a new Plate instance with the given parameters.
    pub fn new(params: PlateParams) -> Self {

        let center = params.center;
        let radius_km = params.radius_km;
        let density = params.density;
        let thickness_km = params.thickness_km;
        let planet_id = params.planet_id;

        if center.x == 0.0 && center.y == 0.0 && center.z == 0.0 {
            panic!("plate cannot be at origin")
        }

        if radius_km <= 0 {
            panic!("radius must be a positive number")
        }

        if density <= 0.0 {
            panic!("density must be a positive number")
        }

        if thickness_km <= 0 {
            panic!("thickness must be a positive number")
        }

        Self {
            id: Uuid::new_v4(),
            center,
            radius_km,
            thickness_km,
            density,
            planet_id,
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Platelet {
    pub id: Uuid,
    pub plate_id: Option<u64>,
    pub position_start: Option<Vec3>,
    pub velocity_start: Option<Vec3>,
    pub radius_km: i32,
    pub h3_cell: u64, // CellIndex
    pub neighbors: [Option<Uuid>; 6], // exactly 6 slots, some may be None
}

