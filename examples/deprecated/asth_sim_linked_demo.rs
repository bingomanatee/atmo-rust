use atmo_rust::asth_sim_linked::{AsthSimLinked, ASLParams};
use atmo_rust::constants::EARTH;

fn main() {
    println!("🌍 AsthSimLinked Demo - Asthenosphere Simulation with Anomalies");
    println!("================================================================");
    
    // Configure the simulation
    let config = ASLParams {
        planet: EARTH.clone(),
        steps: 500, // Total simulation steps
        store_path: String::from("./data/asth_sim_linked_demo"),
        visualize: true,
        vis_freq: 20, // Export visualization every 5 steps
        debug: false, // Set to true for detailed progress output
    };
    
    // Create and initialize the simulation
    println!("🔧 Initializing simulation...");
    let mut sim = AsthSimLinked::new(config);
    
    println!("📊 Simulation initialized with {} cells", sim.cell_count());
    println!("🎯 Features enabled:");
    println!("   • Volume leveling between neighboring cells");
    println!("   • Energy cooling over time");  
    println!("   • Anomaly spawning (15% chance per cycle)");
    println!("   • Anomaly decay (3% per step)");
    println!("   • PNG visualization export every 5 steps");
    println!();
    
    // Run the simulation
    println!("🚀 Starting simulation...");
    sim.run_simulation(500);
    
    println!();
    println!("🎉 Demo completed!");
    println!("📁 Check the 'vis/asth_sim_linked/' folder for PNG visualizations");
    println!("📊 Check the './data/asth_sim_linked_demo/' folder for saved cell data");
    println!();
    println!("🔍 What to look for in the visualizations:");
    println!("   • White/Yellow areas = Very high energy (anomaly hotspots)");
    println!("   • Red areas = High energy (hot regions)");
    println!("   • Purple areas = Medium energy (cooling areas)");
    println!("   • Black areas = Low energy (cold/cooled regions)");
    println!("   • Brightness = Volume levels (bright=high volume, dark=low volume)");
    println!("   • Watch anomaly hotspots appear as white/yellow and fade through red→purple→black");
}