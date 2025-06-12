use rand::Rng;
use num_traits::{Float, FromPrimitive};

pub fn vary_within_range<T>(min: T, max: T, variation: T) -> T
where
    T: Float + FromPrimitive,
{
    let mut rng = rand::rng();
    let base: T = (min + max) /T::from_f64(2.0).unwrap();
    let range = (max - min) * variation;
    let offset = T::from_f64(rng.random_range(-1.0..=1.0)).unwrap() * range;
    (base + offset).clamp(min, max)
}