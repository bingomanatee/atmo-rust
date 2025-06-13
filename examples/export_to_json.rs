use atmo_rust::sim_manager::{SimManager, SimManagerParams};
use atmo_rust::sim::{SimPlanetParams};
use std::env;
use std::path::Path;
use tempfile::tempdir;

fn main() {
    let args: Vec<String> = env::args().collect();
    
    // Default output path
    let output_path = if args.len() > 1 {
        args[1].clone()
    } else {
        "public/simulation_data.json".to_string()
    };

    println!("Generating simulation data...");

    // Create a temporary directory for the database
    let temp_dir = tempdir().expect("Failed to create temp directory");
    let db_path = temp_dir.path().join("sim_db");
    let db_path_str = db_path.to_str().unwrap().to_string();

    // Create simulation with Earth-like parameters
    let planet_config = SimPlanetParams {
        radius: 6372, // Earth radius in km
        mantle_density_gcm3: Some(4.5), // Earth-like density
    };

    let params = SimManagerParams::Create {
        db_path: db_path_str,
        planet_config,
    };

    let manager = SimManager::new(params);
    
    // Generate plates with 60% coverage
    println!("Generating tectonic plates...");
    manager.make_plates(0.6);

    // Ensure output directory exists
    if let Some(parent) = Path::new(&output_path).parent() {
        std::fs::create_dir_all(parent).expect("Failed to create output directory");
    }

    // Export to JSON
    println!("Exporting to JSON: {}", output_path);
    match manager.save_to_json(&output_path) {
        Ok(()) => {
            println!("Successfully exported simulation data to {}", output_path);
            
            // Print some stats
            let export_data = manager.export_data().unwrap();
            println!("Exported data contains:");
            println!("  - Planet radius: {} km", export_data.planet.radius_km);
            println!("  - Planet density: {:.2} g/cm³", export_data.planet.mantle_density_gcm3);
            println!("  - Number of plates: {}", export_data.plates.len());
            
            let total_area: f32 = export_data.plates.iter()
                .map(|plate| std::f32::consts::PI * (plate.radius_km as f32).powi(2))
                .sum();
            let planet_surface_area = export_data.planet.surface_area_km2();
            let coverage = (total_area / planet_surface_area) * 100.0;
            println!("  - Plate coverage: {:.1}%", coverage);
        }
        Err(e) => {
            eprintln!("Failed to export simulation data: {}", e);
            std::process::exit(1);
        }
    }
}
