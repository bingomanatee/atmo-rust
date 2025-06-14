use glam::Vec3;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use h3o::{CellIndex, Resolution};
use uuid::Uuid;
use crate::geoconverter::GeoCellConverter;

pub const PLATE_RESOLUTION: Resolution = Resolution::Two;

#[derive(Serialize, Deserialize, Clone)]
pub struct Plate {
    pub id: Uuid,
    pub center: Vec3,
    pub radius_km: i32,
    pub thickness_km: i32,
    pub density: f64,
    pub planet_id: Uuid,
    pub h3o_l2_cell: CellIndex, // starting cell
}

#[derive(Clone)]
pub struct PlateParams {
    pub center: Vec3,
    pub radius_km: i32,
    pub thickness_km: i32,
    pub density: f64,
    pub planet_id: Uuid,
    pub planet_radius_km: i32,
}

impl Plate {
    /// Creates a new Plate instance with the given parameters.
    pub fn new(params: PlateParams) -> Self {
        let center = params.center;
        let radius_km = params.radius_km;
        let density = params.density;
        let thickness_km = params.thickness_km;
        let planet_id = params.planet_id;
        let planet_radius_km = params.planet_radius_km;

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
        
        let gc = GeoCellConverter::new(planet_radius_km as f64, PLATE_RESOLUTION);
        let location = gc.randomize_point_on_sphere(center);
        Self {
            id: Uuid::new_v4(),
            center: location,
            radius_km,
            thickness_km,
            density,
            planet_id,
            h3o_l2_cell: gc.vec3_to_cell(location).unwrap()
        }
    }
}