use crate::planet::Planet;
use rand::Rng;

pub fn sample_power_law(min: i32, max: i32, exp: f64) -> f64 {
    let mut rng = rand::rng();
    let u: f64 = rng.random();
    let inv = 1.0 / (1.0 - exp);
    let a = (min as f64).powf(1.0 - exp);
    let b = (max as f64).powf(1.0 - exp);
    (u * (b - a) + a).powf(inv)
}

pub fn rad_to_area_int(radius_km: i32) -> i32 {
    (std::f64::consts::PI * radius_km.pow(2) as f64).round() as i32
}