use atmo_rust::constants::{EARTH, LAYER_COUNT};
use atmo_rust::rock_store::RockStore;
use atmo_rust::sim_next::{SimNext, SimNextProps};
use tempfile;

fn main() {
    println!("🧪 Testing multi-layer asthenosphere simulation...");
    
    // Create a temporary database
    let temp_dir = tempfile::tempdir().expect("Failed to create temp dir");
    let db_path = temp_dir.path().join("test_db");
    
    let store = RockStore::open(db_path.to_str().unwrap()).expect("Failed to create RockStore");
    
    let props = SimNextProps::new(EARTH.clone(), store)
        .with_debug(true)
        .with_visualization(false, 1);
    
    println!("🌍 Creating simulation with {} layers...", LAYER_COUNT);
    let mut sim = SimNext::new(props);
    
    println!("📊 Initial statistics:");
    println!("  Total cells: {}", sim.cell_count());
    println!("  Cells per layer: {}", sim.cell_count() / LAYER_COUNT);
    
    // Check that both layers have cells
    for layer_idx in 0..LAYER_COUNT {
        let layer_cell_count = sim.cell_layers[layer_idx].len();
        println!("  Layer {}: {} cells", layer_idx, layer_cell_count);
    }
    
    println!("🔄 Running one simulation step...");
    sim.run_step();
    
    println!("✅ Test completed successfully!");
    println!("📈 Current step: {}", sim.current_step());
}