use std::collections::HashSet;
use crate::helpers::{rad_to_area_int, sample_power_law};
use crate::planet::Planet;
use crate::deprecated::plate::{Plate, PlateParams};
use crate::vary::vary_within_range;
use glam::Vec3;
use rand::Rng;
use std::f64::consts::FRAC_PI_6;
use uuid::Uuid;
use crate::geoconverter::GeoCellConverter;
use crate::h3_utils::PointSampler;

pub struct PlateGeneratorConfig {
    pub target_coverage: f64,
    pub power_law_exponent: f64,
    pub min_density: f64,
    pub max_density: f64,
    pub min_thickness_km: i32,
    pub max_thickness_km: i32,
    pub variation_factor: f64,
    pub max_plate_radius: Option<f64>,
}

pub struct PartialPlateGenConfig {
    pub target_coverage: Option<f64>,
    pub power_law_exponent: Option<f64>,
    pub min_density: Option<f64>,
    pub max_density: Option<f64>,
    pub min_thickness_km: Option<i32>,
    pub max_thickness_km: Option<i32>,
    pub variation_factor: Option<f64>,
    pub max_plate_radius_radians: Option<f64>,
}

impl PlateGeneratorConfig {
    pub fn from_partial(partial: PartialPlateGenConfig) -> PlateGeneratorConfig {
        let variation_factor = partial.variation_factor.unwrap_or(0.1);

        if variation_factor < 0.0 || variation_factor > 0.2 {
            panic!("Variation factor must be between 0 and 0.2");
        }
        Self {
            target_coverage: partial.target_coverage.unwrap_or(0.7),
            power_law_exponent: partial.power_law_exponent.unwrap_or(3.0),
            min_density: partial.min_density.unwrap_or(2.6),
            max_density: partial.max_density.unwrap_or(3.3),
            min_thickness_km: partial.min_thickness_km.unwrap_or(7),
            max_thickness_km: partial.max_thickness_km.unwrap_or(35),
            variation_factor: partial.variation_factor.unwrap_or(0.1),
            max_plate_radius: Option::from(partial.max_plate_radius_radians.unwrap_or(FRAC_PI_6)),
        }
    }
}

pub struct GenerateRadiiParams {
    pub target_coverage: f64,
    pub min_radius: i32,
    pub max_radius: i32,
    pub exponent: f64,
}

pub struct PlateGenerator<'a> {
    config: PlateGeneratorConfig,
    planet: &'a Planet,
    point_sampler: PointSampler,
}


impl<'a> PlateGenerator<'a>  {

    pub fn new(partial_plate_gen_config: PartialPlateGenConfig, planet: &'a Planet) -> Self {
        let config = PlateGeneratorConfig::from_partial(partial_plate_gen_config);
        let point_sampler = PointSampler::new(planet.radius_km);
        Self { config, planet, point_sampler }
    }
    pub fn surface_area_km2(&self) -> f64 {
        self.planet.surface_area_km2()
    }

    pub fn random_planet_point(&mut self) -> Vec3 {
        self.point_sampler.random_point_on_planet()
    }

    pub fn generate_one(&mut self, radius_km: i32, planet_id: Uuid) -> Plate {
        let center = self.random_planet_point();

        // thickness: vary based on radius
        let t_variation = self.config.variation_factor;
        let base_thickness = self.scaled_thickness(radius_km as f64) ;
        let max_thickness = base_thickness as f64 * (1.0 + t_variation);
        let min_thickness = base_thickness as f64 * (1.0 - t_variation);

        let thickness_km =
            vary_within_range(min_thickness, max_thickness, self.config.variation_factor)
                .clamp(
                    self.config.min_thickness_km as f64,
                    self.config.max_thickness_km as f64,
                )
                .floor();

        // density: vary based on radius
        let base_density = self.scaled_density(radius_km);
        let max_density = base_density * (1.0 + t_variation);
        let min_density = base_density * (1.0 - t_variation);
        let density = vary_within_range(min_density, max_density, self.config.variation_factor)
            .clamp(self.config.min_density, self.config.max_density);



        Plate::new(
            PlateParams {
                center,
                radius_km: radius_km as i32,
                thickness_km: thickness_km as i32,
                density,
                planet_id,
                planet_radius_km: self.planet.radius_km
            }
        )
    }

    // Convert max_plate_radius (radians) to kilometers on the planet surface
    fn max_plate_radius_km(&self) -> f64 {
        let max_rad = self.config.max_plate_radius.unwrap_or(FRAC_PI_6);
        max_rad * (self.planet.radius_km as f64)
    }

    fn scaled_thickness(&self, radius_km: f64) -> i32 {
        let config = &self.config;
        let t_min = config.min_thickness_km as f64;
        let t_max = config.max_thickness_km as f64;

        let max_radius_km = self.max_plate_radius_km();

        let ratio = 1.0 - (radius_km / max_radius_km).clamp(0.0, 1.0); // inverse scaling
        let thickness = (t_min + (t_max - t_min) * ratio).clamp(
            config.min_thickness_km as f64,
            config.max_thickness_km as f64,
        );
        thickness.round() as i32
    }

    fn scaled_density(&self, radius_km: i32) -> f64 {
        let config = &self.config;
        let d_min = config.min_density;
        let d_max = config.max_density;

        let max_radius_km = self.max_plate_radius_km();

        let ratio = (radius_km as f64 / max_radius_km).clamp(0.0, 1.0); // direct scaling
        (d_min + (d_max - d_min) * ratio).clamp(config.min_density, config.max_density)
    }
    pub fn generate_radii(&self, params: GenerateRadiiParams) -> Vec<i32> {
        let target_coverage = params.target_coverage;
        let min_radius = params.min_radius;
        let max_radius = params.max_radius;
        let exponent = params.exponent;
        let planet = &self.planet;

        println!("generate_radii between {} and {}", min_radius, max_radius);
        assert!(min_radius < max_radius, "min rad must be < max_rad");
        assert!(
            max_radius < planet.radius_km,
            "radii cannot be > planet radius {}",
            planet.radius_km
        );

        let mut radii: Vec<i32> = Vec::new();
        let mut status = area_status(&radii, &planet, target_coverage);
        while status.actual_area < status.target_area {
            let r = sample_power_law(min_radius, max_radius, exponent);
            radii.push(r.round() as i32);
            status = area_status(&radii, &planet, target_coverage);
        }

        // sort descending
        radii.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let status = area_status(&radii, planet, target_coverage);

        if status.actual_coverage > status.target_coverage {
            loop {
                if radii.is_empty() {
                    break;
                }
                let status = area_status(&radii, &planet, target_coverage);

                if status.actual_area <= status.target_area {
                    break;
                }

                let smallest_area = rad_to_area_int(radii[0]);
                let overage = status.actual_area - status.target_area;
                if smallest_area > overage {
                    break;
                }

                for (idx, &val) in radii.iter().enumerate() {
                    if rad_to_area_int(val) > overage {
                        let mut randomizer = rand::rng();
                        let index_to_delete = randomizer.random_range(0..idx);
                        radii.remove(index_to_delete);
                        break;
                    }
                }
            }
        }

        radii
    }
}

fn find_index_exceeding_area(radii: &Vec<i32>, total_area: i32) -> Option<usize> {
    radii
        .iter()
        .scan(0, |acc, rad| {
            *acc += rad_to_area_int(*rad);
            Some(*acc)
        })
        .position(|cum| cum > total_area)
}

// Helper to compute total spherical cap area of current radii
fn total_area(radii: &Vec<i32>) -> i32 {
    radii.iter().map(|&r| rad_to_area_int(r)).sum()
}

fn area_status(radii: &Vec<i32>, planet: &Planet, target_coverage: f64) -> AreaStatus {
    let planet_surface = planet.surface_area_km2();
    let target_area = (planet_surface * target_coverage) as i32;

    let actual_area = total_area(radii);
    let actual_coverage = actual_area as f64 / planet_surface;
    /*    println!(
        "TARGET: coverage {}, area: {} \n ACTUAL: coverage: {}, area: {}\nOVERAGE: area {}, ratio {}",
        target_coverage,
        target_area,
        actual_coverage,
        actual_area,
        actual_area - target_area,
        actual_area as f64 / target_area as f64
    );*/

    AreaStatus {
        planet_surface,

        target_coverage,
        target_area,

        actual_coverage,
        actual_area,
    }
}
struct AreaStatus {
    planet_surface: f64,

    target_coverage: f64,
    target_area: i32,

    actual_coverage: f64,
    actual_area: i32,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::planet::{EARTH_RADIUS_KM};
    use uuid::Uuid;
    use crate::constants::EARTH;

    #[test]
    fn test_generate_one_with_full_config() {
        let full_config = PartialPlateGenConfig {
            target_coverage: Some(0.8),
            power_law_exponent: Some(2.5),
            min_density: Some(2.7),
            max_density: Some(3.1),
            min_thickness_km: Some(10),
            max_thickness_km: Some(30),
            variation_factor: Some(0.1),
            max_plate_radius_radians: Some(FRAC_PI_6),
        };

        let mut generator = PlateGenerator::new(full_config, &EARTH);

        let radius_km = 100;
        let planet_id = Uuid::new_v4();
        let plate = generator.generate_one(radius_km, planet_id);

        assert_eq!(plate.radius_km, radius_km as i32);
        assert!(plate.thickness_km >= generator.config.min_thickness_km);
        assert!(plate.thickness_km <= generator.config.max_thickness_km);

        assert!((plate.density) >= generator.config.min_density);
        assert!((plate.density) <= generator.config.max_density);
        assert_eq!(plate.planet_id, planet_id);
    }

    #[test]
    fn test_generate_one_with_partial_config_defaults() {
        // Empty partial config to use all defaults
        let partial_config = PartialPlateGenConfig {
            target_coverage: None,
            power_law_exponent: None,
            min_density: None,
            max_density: None,
            min_thickness_km: None,
            max_thickness_km: None,
            variation_factor: None,
            max_plate_radius_radians: None,
        };

        let mut generator = PlateGenerator::new(partial_config, &EARTH);

        let radius_km = 50;
        let planet_id = Uuid::new_v4();
        let plate = generator.generate_one(radius_km, planet_id);

        // Check radius
        assert_eq!(plate.radius_km, radius_km as i32);

        // Thickness within default min/max
        assert!(plate.thickness_km >= generator.config.min_thickness_km);
        assert!(plate.thickness_km <= generator.config.max_thickness_km);

        // Density within default min/max
        assert!((plate.density) >= generator.config.min_density);
        assert!((plate.density) <= generator.config.max_density);

        // Check that planet_id matches
        assert_eq!(plate.planet_id, planet_id);
    }

    #[test]
    fn test_generate_radii() {
        let target_coverage = 0.5;
        let min_radius = 200;
        let max_radius = EARTH_RADIUS_KM / 3;
        let exponent = 3.0;
        let partial_config = PartialPlateGenConfig {
            target_coverage: None,
            power_law_exponent: None,
            min_density: None,
            max_density: None,
            min_thickness_km: None,
            max_thickness_km: None,
            variation_factor: None,
            max_plate_radius_radians: None,
        };
        let generator = PlateGenerator::new(partial_config, &EARTH);

        let radii = generator.generate_radii(GenerateRadiiParams {
            target_coverage,
            min_radius,
            max_radius,
            exponent,
        });

        // Check that all radii are within bounds
        for &r in &radii {
            assert!(r >= min_radius);
            assert!(r <= max_radius);
        }

        let total = total_area(&radii);
        let surface_area = EARTH.surface_area_km2();
        let target_area = (surface_area * target_coverage).round() as i32;

        println!("total_area: {}, target_area: {}", total, target_area);
        assert!(total >= target_area);
        assert!(total as f64 <= target_area as f64 * 1.5);
    }
}
