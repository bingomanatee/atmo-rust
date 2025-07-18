use uuid::Uuid;
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
pub const LEVEL_AMT: f64 = 0.2; // 100% equilibration rate, adjusted for ~6 neighbors per cell in binary pair leveller
pub const STANDARD_STEPS: u32 = 500;
pub const HEIGHT_PER_ASTH_LAYER : f64 = 10.0; // km
pub const EARTH_RADIUS_KM: i32 = 6372;
pub const RHO_EARTH: f64 = 4.5; // g/cm³
pub const VOLCANO_BIAS: f64 = 0.05;
pub static EARTH: Planet = Planet {
    id: EARTH_ID,
    sim_id: EARTH_SIM_ID,
    radius_km: EARTH_RADIUS_KM,
    mantle_density_gcm3: RHO_EARTH,
};

pub static COOLING_RATE: Lazy<f64> = Lazy::new(|| {
    let mid_target = CELL_JOULES_EQUILIBRIUM ;
    (mid_target/ CELL_JOULES_START).powf(  1.1/ (STANDARD_STEPS as f64) )
});

// === Volcano Properties ===
pub const VOLCANO_MIN_VOLUME: f64 = AVG_STARTING_VOLUME_KM_3 * 0.005; // ~6,089 km³
pub const VOLCANO_MAX_VOLUME: f64 = AVG_STARTING_VOLUME_KM_3 * 0.03;  // ~24,356 km³
pub const VOLCANO_JOULES_PER_VOLUME: f64 = JOULES_PER_KM3 * 1.5;

// === Sinkhole Properties ===
pub const SINKHOLE_MIN_VOLUME_TO_REMOVE: f64 = AVG_STARTING_VOLUME_KM_3 * 0.0001; // ~1,218 km³
pub const SINKHOLE_MAX_VOLUME_TO_REMOVE: f64 = AVG_STARTING_VOLUME_KM_3 * 0.001; // ~6,089 km³

// These are per-cell, per-step spawn probabilities.
pub const VOLCANO_CHANCE: f64 = 0.0005;       // 1 in 2000 cells per step
pub const SINKHOLE_CHANCE: f64 = 0.005;        // 1 in 100 cells per step

// === Decay Rates ===
// These control how fast anomalies fade each step (multiplicative decay).
pub const VOLCANO_DECAY_RATE: f64 = 0.666;
pub const SINKHOLE_DECAY_RATE: f64 = 0.85;


// ===== Convection =====
pub const GLOBAL_CONVECTION: f64 = 0.3;
pub const CONVECTION_ADDITION_MIN: f64 = 0.1;
pub const CONVECTION_ADDITION_MAX: f64 = 0.2;
pub const CONVECTION_SUBTRACTION_MIN: f64 = 0.1;
pub const CONVECTION_SUBTRACTION_MAX: f64 = 0.2;
pub const CONVECTION_BALANCE_TOLERANCE: f64 = 0.03; // 3% tolerance between addition and subtraction

// ===== Convection Templates =====
pub const CONVECTION_TEMPLATE_LIFESPAN_MIN: u32 = 100;
pub const CONVECTION_TEMPLATE_LIFESPAN_MAX: u32 = 300;
pub const CONVECTION_NOISE_SCALE_MIN: f32 = 3.0;
pub const CONVECTION_NOISE_SCALE_MAX: f32 = 8.0;

// ======= Layers 

pub const LAYER_COUNT: usize = 2;
pub const ENERGY_INCREASE_PER_LAYER: f64 = 0.5; // each layer is 50% more joules based on the starting joules abov
// eg the second layer has  1 1/2 as many joules per cell stariting the next layer has 2x , then 2.5x, etc. 
// note - we are not computing additional density per layer as at this scale it only increaeses by < 10% so the added math burden is 
// not going to create added realism.
pub const VERTICAL_ENERGY_MIXING : f64 = 0.2; // that is, each turn, 
// 20% of the energy from the layer below goes up 
// and 20% of the energy from the layer above goes down.

// ======= Lithosphere 
pub const DENSITY_KG_M3: f64 = 3300.0;
pub const SPECIFIC_HEAT_KJ_KG_C: f64 = 1.0; // or 1000 J/kg°C
pub const LITHOSPHERE_CREATION_TEMP_C: f64 = 850.0;
pub const LITHOSPHERE_CREATION_JOULES_PER_KM3 : f64 = LITHOSPHERE_CREATION_TEMP_C * SPECIFIC_HEAT_KJ_KG_C * DENSITY_KG_M3;
