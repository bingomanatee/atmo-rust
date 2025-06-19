use atmo_rust::perlin_noise_generator::PerlinNoiseGenerator;
use atmo_rust::constants::EARTH;
use atmo_rust::geoconverter::GeoCellConverter;
use atmo_rust::h3_utils::H3Utils;
use atmo_rust::asthenosphere::ASTH_RES;
use glam::Vec3;

fn main() {
    println!("🌊 Simple Perlin Noise Generator Demo");
    println!("=====================================");
    println!();
    
    // Create different noise generators to compare
    let generators = [
        ("Smooth & Gentle", PerlinNoiseGenerator::new(42, 0.5, 0.8)),
        ("Moderate", PerlinNoiseGenerator::new(42, 1.0, 1.0)),
        ("Detailed & Sharp", PerlinNoiseGenerator::new(42, 3.0, 2.0)),
        ("Very Detailed", PerlinNoiseGenerator::new(42, 8.0, 1.5)),
        ("Sharp Peaks", PerlinNoiseGenerator::new(42, 1.5, 3.0)),
    ];
    
    println!("📊 Generator Settings:");
    for (name, generator) in &generators {
        println!("  {}: detail={:.1}, sharpness={:.1}", name, generator.detail, generator.sharpness);
    }
    println!();
    
    // Sample from real H3 cell locations
    println!("🌍 Sampling from actual H3 cell locations...");
    let sample_coords = collect_sample_coordinates(100);
    
    // Analyze each generator
    for (name, generator) in &generators {
        println!("📈 {} Analysis:", name);
        let stats = generator.analyze(&sample_coords);
        print!("{}", stats);
        println!();
    }
    
    // Show specific examples at test locations
    println!("🎯 Sample Values at Specific Locations:");
    let test_locations = [
        Vec3::new(1.0, 0.0, 0.0).normalize(),
        Vec3::new(0.0, 1.0, 0.0).normalize(),
        Vec3::new(0.0, 0.0, 1.0).normalize(),
        Vec3::new(1.0, 1.0, 1.0).normalize(),
    ];
    
    for (i, coord) in test_locations.iter().enumerate() {
        println!("Location {}: ({:.3}, {:.3}, {:.3})", i+1, coord.x, coord.y, coord.z);
        for (name, generator) in &generators {
            let value = generator.sample(*coord);
            println!("  {}: {:.3}", name, value);
        }
        println!();
    }
    
    // Show how to map noise to custom ranges centered around 1.0
    println!("🔧 Example: Mapping to Ranges Centered Around 1.0");
    let noise_gen = PerlinNoiseGenerator::new(42, 1.0, 1.0);
    let test_coord = Vec3::new(1.0, 0.5, 0.0).normalize();
    let noise_val = noise_gen.sample(test_coord); // [-1, 1]
    
    println!("Raw noise value: {:.3}", noise_val);
    
    // Map to volume scale [0.8, 1.2] centered around 1.0 (±20%)
    let volume_scale = 1.0 + noise_val * 0.2;
    println!("Volume scale: {:.3}x (80% to 120%, centered at 1.0)", volume_scale);
    
    // Map to energy scale [0.7, 1.3] centered around 1.0 (±30%)  
    let energy_scale = 1.0 + noise_val * 0.3;
    println!("Energy scale: {:.3}x (70% to 130%, centered at 1.0)", energy_scale);
    
    // Map to symmetric range [0.5, 1.5] centered around 1.0 (±50%)
    let wide_scale = 1.0 + noise_val * 0.5;
    println!("Wide scale: {:.3}x (50% to 150%, centered at 1.0)", wide_scale);
    
    println!();
    println!("💡 Tuning Tips:");
    println!("  • detail: 0.1=very smooth, 1.0=moderate, 10.0=very detailed");
    println!("  • sharpness: 0.5=gentle curves, 1.0=moderate, 3.0=sharp peaks");
    println!("  • Use different seeds for independent noise patterns");
    println!("  • Map [-1,1] output to any range you need in your code");
}

fn collect_sample_coordinates(max_samples: usize) -> Vec<Vec3> {
    let planet = EARTH.clone();
    let gc = GeoCellConverter::new(planet.radius_km as f64, ASTH_RES);
    let mut coords = Vec::new();
    let mut count = 0;
    
    H3Utils::iter_at(ASTH_RES, |cell_index| {
        if count < max_samples {
            let location = gc.cell_to_vec3(cell_index).normalize();
            coords.push(location);
            count += 1;
        }
    });
    
    coords
}