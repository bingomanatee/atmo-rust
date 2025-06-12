use serde::{Deserialize, Serialize};
use uuid::Uuid;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Planet {
    pub id: Uuid,
    pub sim_id: Uuid,
    pub radius: i32,
    pub mantle_density_gcm3: f64,
}

pub const EARTH_ID: Uuid = Uuid::from_u128(0x1234567890abcdef1234567890abcdef);
pub const EARTH_SIM_ID: Uuid = Uuid::from_u128(0xfedcba9876543210fedcba9876543210);

pub const EARTH: Planet = Planet {
    radius: EARTH_RADIUS_KM,
    id: EARTH_ID,
    sim_id: EARTH_SIM_ID,
    mantle_density_gcm3: RHO_EARTH,
};

pub struct PlanetParams {
    pub sim_id: Uuid,
    pub radius: i32,
    pub mantle_density_gcm3: Option<f64>,
}

impl Planet {
    pub fn new(params: PlanetParams) -> Self {
        if params.radius <= 0{
            panic!("planet must have a positive radius");
        }

        Planet {
            id: Uuid::new_v4(),
            sim_id: params.sim_id,
            radius: params.radius,
            mantle_density_gcm3: params
                .mantle_density_gcm3
                .unwrap_or_else(|| estimate_mantle_density(params.radius, None)),
        }
    }

    /// Returns the surface area of the planet in square kilometers.
    pub fn surface_area_km2(&self) -> f32 {
        (4.0 * std::f64::consts::PI * (self.radius as f32).powi(2) as f64 )as f32
    }
}

use crate::vary::vary_within_range;

pub const EARTH_RADIUS_KM: i32 = 6372;
pub const RHO_EARTH: f64 = 4.5; // g/cm³

pub fn estimate_mantle_density(planet_radius_km: i32, exp: Option<f64>) -> f64 {
    if planet_radius_km < EARTH_RADIUS_KM / 10 {
        panic!(
            "estimate_mantle_density:  radii must be at least {}; is {}",
            EARTH_RADIUS_KM / 10, planet_radius_km
        )
    }
    let scale_factor = ((planet_radius_km / EARTH_RADIUS_KM) as f64).powf(exp.unwrap_or(0.25));
    let base_density = RHO_EARTH * scale_factor;

    // 2) Define a ±10% band around that base density
    let min_density = base_density * 0.9;
    let max_density = base_density * 1.1;

    vary_within_range(min_density, max_density, 1.0)
}
