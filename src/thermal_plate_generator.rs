use std::f64::consts::PI;

#[derive(Debug)]
pub struct Platelet {
    pub radius_km: f64,     // in kilometers
    pub thickness_km: f64,  // in kilometers
    pub mass_kg: f64,
    pub density_kg_m3: f64,
}

/// Create an oceanic platelet using pressure in GPa
/// Returns None if pressure is too low to generate crust
pub fn create_platelet_from_pressure(pressure_gpa: f64) -> Option<Platelet> {
    let threshold_gpa = 3.1;
    if pressure_gpa <= threshold_gpa {
        return None;
    }

    // Pressure surplus
    let surplus_gpa = pressure_gpa - threshold_gpa;

    // Assume 100_000 kg/m² of crust mass added per GPa
    let mass_per_m2 = surplus_gpa * 100_000.0; // kg/m²

    // Radius in km (approx for H3 L2)
    let radius_km = 900.0;
    let radius_m = radius_km * 1_000.0;

    // Area in m²
    let area_m2 = PI * radius_m * radius_m;

    // Mass in kg
    let mass = mass_per_m2 * area_m2;

    // Oceanic crust density in kg/m³
    let density = 3_000.0;

    // Volume = mass / density, in m³
    let volume_m3 = mass / density;

    // Thickness in meters, then km
    let thickness_km = (volume_m3 / area_m2) / 1_000.0;

    Some(Platelet {
        radius_km,
        thickness_km,
        mass_kg: mass,
        density_kg_m3: density,
    })
}