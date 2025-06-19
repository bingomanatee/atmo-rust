use uuid::Uuid;
use num_traits::real::Real;
use once_cell::sync::Lazy;
use crate::planet::Planet;

pub const EARTH_ID: Uuid = Uuid::from_u128(0x1234567890abcdef1234567890abcdef);
pub const EARTH_SIM_ID: Uuid = Uuid::from_u128(0xfedcba9876543210fedcba9876543210);

pub const CELL_JOULES_START: f64 = 5.78e24;
pub const CELL_JOULES_EQUILIBRIUM: f64 = 5.13e24;
pub const JOULES_TO_KELVIN:f64= 3.7879e-7;
pub const AE_SPAN: f64 = CELL_JOULES_START - CELL_JOULES_EQUILIBRIUM;

pub const MAX_SUNK_PERCENT: f64 = 0.0000025;
pub const MAX_SUNK_TEMP: f64 = 1875.0;
pub const MIN_SUNK_TEMP: f64 = 1700.0;

pub const AVG_STARTING_VOLUME_KM_3: f64 = 1217830.0;
pub const JOULES_PER_KM3: f64 = CELL_JOULES_START / AVG_STARTING_VOLUME_KM_3;

pub const AVG_VOLUME_TO_ADD: f64 = 0.0325;
pub const LEVEL_AMT: f64 = 1.0; // Percentage for volume levelling between cells
pub const STANDARD_STEPS: u32 = 500;

// Anomaly (flow) constants - reduced magnitude to not overwhelm cooling
pub const ANOMALY_SPAWN_CHANCE: f64 = 0.15; // 15% chance per cycle
pub const ANOMALY_DECAY_RATE: f64 = 0.03; // 3% decay per step
pub const ANOMALY_ENERGY_AMOUNT: f64 = CELL_JOULES_START * 1.0; // 2% of starting energy (was 20%)
pub const ANOMALY_VOLUME_AMOUNT: f64 = AVG_STARTING_VOLUME_KM_3 * 0.1; // 1% of starting volume (was 10%)

pub const EARTH_RADIUS_KM: i32 = 6372;
pub const RHO_EARTH: f64 = 4.5; // g/cm³

pub static EARTH: Planet = Planet {
    id: EARTH_ID,
    sim_id: EARTH_SIM_ID,
    radius_km: EARTH_RADIUS_KM,
    mantle_density_gcm3: RHO_EARTH,
};

pub static COOLING_RATE: Lazy<f64> = Lazy::new(|| {
    let mid_target = (CELL_JOULES_EQUILIBRIUM + MAX_SUNK_TEMP)/ 2.0;
    (mid_target/ CELL_JOULES_START).powf(1.0 / (STANDARD_STEPS as f64) )
});
