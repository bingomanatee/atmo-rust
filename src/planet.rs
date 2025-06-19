use std::collections::HashSet;
use once_cell::sync::Lazy;
use serde::{Deserialize, Serialize};
use uuid::Uuid;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Planet {
    pub id: Uuid,
    pub sim_id: Uuid,
    pub radius_km: i32,
    pub mantle_density_gcm3: f64
}

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
            radius_km: params.radius,
            mantle_density_gcm3: params
                .mantle_density_gcm3
                .unwrap_or_else(|| estimate_mantle_density(params.radius, None)),
        }
    }

    /// Returns the surface area of the planet in square kilometers.
    pub fn surface_area_km2(&self) -> f64 {
        (4.0 * std::f64::consts::PI * (self.radius_km as f64).powi(2) as f64 )as f64
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
