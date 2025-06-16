use num_traits::real::Real;
use once_cell::sync::Lazy;

pub const MAX_SUNK_PERCENT: f64 = 0.0000025;
pub const MAX_SUNK_TEMP: f64 = 1875.0;
pub const MIN_SUNK_TEMP: f64 = 1700.0;
pub const AVG_STARTING_VOLUME: f64 = 4260.0;
pub const CELL_ENERGY_START: f64 = 2000.0;
pub const CELL_ENERGY_EQUILIBRIUM: f64 = 1700.0;
pub const AE_SPAN: f64 = CELL_ENERGY_START - CELL_ENERGY_EQUILIBRIUM;
pub const K_PER_VOLUME: f64 = CELL_ENERGY_START/AVG_STARTING_VOLUME;
pub const AVG_VOLUME_TO_ADD: f64 = 0.0325;
pub const STANDARD_STEPS: u32 = 500;

pub static COOLING_RATE: Lazy<f64> = Lazy::new(|| {
    let mid_target = (CELL_ENERGY_EQUILIBRIUM + MAX_SUNK_TEMP)/ 2.0;
    (mid_target/CELL_ENERGY_START).powf(1.0 / (STANDARD_STEPS as f64) )
});
