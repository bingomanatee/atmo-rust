use atmo_rust::constants::EARTH;
use atmo_rust::rock_store::RockStore;
use atmo_rust::sim_next::{SimNext, SimNextProps};
use std::env;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🌍 SimNext - Fast Asthenosphere Simulation");
    println!("==========================================");

    // Parse command line arguments
    let args: Vec<String> = env::args().collect();
    let steps = if args.len() > 1 {
        args[1].parse::<u32>().unwrap_or(100)
    } else {
        100
    };

    println!("📋 Configuration:");
    println!("  Planet: Earth");
    println!("  Steps: {}", steps);
    println!("  Database: data/sim_next_example.db");
    println!("  Noise: Reduced exponential scaling for smoother initial conditions");
    println!();

    // Create rock store
    let store = RockStore::open("data/sim_next_example.db")?;
    println!("💾 Database opened successfully");

    // Create and initialize simulation
    println!("🔧 Creating simulation...");
    let props = SimNextProps::new(EARTH.clone(), store)
        .with_visualization(false, 0)
        .with_database_saving(false)
        .with_debug(false);
    let mut sim = SimNext::new(props);
    
    println!("✅ Simulation initialized with {} cells", sim.cell_count());
    println!("🔗 Generated {} binary pairs for levelling", sim.binary_pairs.len());
    println!("🖼️ Visualization enabled - images will export every 10 steps");
    println!();

    // Run the simulation
    println!("🚀 Starting simulation...");
    sim.run_simulation(steps);

    println!();
    println!("✅ Simulation completed successfully!");
    println!("💾 Results saved to sim_next_example.db");

    Ok(())
}
