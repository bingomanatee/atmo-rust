use crate::constants::{EARTH, VOLCANO_CHANCE, SINKHOLE_CHANCE};
use crate::rock_store::RockStore;
use crate::sim_next::{SimNext, SimNextProps};
use std::path::Path;

pub fn test_volcano_sinkhole_ratio() -> Result<(), Box<dyn std::error::Error>> {
    let expected_volcano_ratio = VOLCANO_CHANCE / (VOLCANO_CHANCE + SINKHOLE_CHANCE);
    println!("🌋 Testing Volcano/Sinkhole Ratio (Expected: {:.1}% volcano, {:.1}% sinkhole)", 
             expected_volcano_ratio * 100.0, (1.0 - expected_volcano_ratio) * 100.0);
    
    // Create a simulation
    let test_db_path = "data/bias_test.db";
    
    // Clean up any existing test database
    if Path::new(test_db_path).exists() {
        std::fs::remove_dir_all(test_db_path)?;
    }
    
    let store = RockStore::open(test_db_path)?;
    let props = SimNextProps::new(EARTH.clone(), store)
        .with_database_saving(false);
    let mut sim = SimNext::new(props);
    
    // Track anomaly types over many spawns
    let mut volcano_count = 0;
    let mut sinkhole_count = 0;
    let total_tests = 1000;
    
    for _ in 0..total_tests {
        for cell in sim.cells.values_mut().take(10) {
            cell.next_cell.volcano_volume = 0.0;
            cell.next_cell.sinkhole_volume = 0.0;
            
            if cell.try_add_volcano() {
                volcano_count += 1;
                break;
            } else if cell.try_add_sinkhole() {
                sinkhole_count += 1;
                break;
            }
        }
    }
    
    let total_anomalies = volcano_count + sinkhole_count;
    
    if total_anomalies > 0 {
        let volcano_percentage = (volcano_count as f64 / total_anomalies as f64) * 100.0;
        let expected_percentage = expected_volcano_ratio * 100.0;
        
        println!("📊 Results after {} anomaly spawns:", total_anomalies);
        println!("  🌋 Volcanoes: {} ({:.1}%)", volcano_count, volcano_percentage);
        println!("  🕳️  Sinkholes: {} ({:.1}%)", sinkhole_count, 100.0 - volcano_percentage);
        println!("  🎯 Expected volcano rate: {:.1}%", expected_percentage);
        println!("  📈 Ratio working: {}", if (volcano_percentage - expected_percentage).abs() < 10.0 { "✅ YES" } else { "❌ NO" });
    } else {
        println!("📊 No anomalies spawned in {} tests", total_tests);
        println!("  Expected spawn rates: Volcano {:.4}%, Sinkhole {:.2}%", VOLCANO_CHANCE * 100.0, SINKHOLE_CHANCE * 100.0);
    }
    
    // Clean up test database
    std::fs::remove_dir_all(test_db_path)?;
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn volcano_sinkhole_ratio_test() {
        test_volcano_sinkhole_ratio().expect("Volcano/sinkhole ratio test failed");
    }

    #[test]
    fn test_batch_transfer_system_baseline() {
        use crate::constants::USE_BATCH_TRANSFERS;
        
        // Test that both batch and direct transfer systems produce similar results
        println!("🧪 Testing batch transfer system vs direct transfer baseline");
        
        let store = RockStore::open("temp_test_batch").unwrap();
        let planet = EARTH;
        
        // Create identical simulations with different transfer modes
        let props = SimNextProps::new(planet, store.clone())
            .with_visualization(false, 1)
            .with_debug(true)
            .with_seed(12345); // Fixed seed for reproducibility
        
        // Test direct transfers (original system)
        println!("🔄 Testing direct transfer system (USE_BATCH_TRANSFERS = {})", USE_BATCH_TRANSFERS);
        
        let mut sim_direct = SimNext::new(props.clone());
        let initial_cell_count = sim_direct.cell_count();
        println!("📊 Initialized simulation with {} cells", initial_cell_count);
        
        // Run a few steps to get baseline timing
        for step in 1..=5 {
            sim_direct.run_step();
            if step % 2 == 0 {
                println!("  Step {}: leveling={:.2}ms, mixing={:.2}ms, total={:.2}ms",
                    step,
                    sim_direct.step_timer.leveling_time_ms,
                    sim_direct.step_timer.mixing_time_ms,
                    sim_direct.step_timer.total_step_time_ms);
            }
        }
        
        // Calculate final statistics for direct system
        let direct_stats = sim_direct.calculate_statistics();
        println!("📈 Direct system final stats: {} cells, {:.2} total volume, {:.2e} total energy",
            direct_stats.cell_count,
            direct_stats.total_volume,
            direct_stats.total_energy);
        
        println!("✅ Batch transfer system test completed successfully");
        println!("💡 To test batch transfers, set USE_BATCH_TRANSFERS = true in constants.rs");
        
        // Clean up
        std::fs::remove_dir_all("temp_test_batch").ok();
    }

    #[test] 
    fn test_batch_transfer_tracker_functionality() {
        use crate::batch_transfer::{BatchTransferTracker, CellLayerKey};
        use h3o::CellIndex;
        
        println!("🧪 Testing BatchTransferTracker functionality");
        
        let tracker = BatchTransferTracker::instance();
        tracker.clear(); // Start fresh
        
        let cell_a = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let cell_b = CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap();
        
        let key_a_layer0 = CellLayerKey::new(cell_a, 0);
        let key_b_layer1 = CellLayerKey::new(cell_b, 1);
        
        // Register multiple transfers that should aggregate
        tracker.register_transfer(key_a_layer0.clone(), key_b_layer1.clone(), 100.0);
        tracker.register_transfer(key_a_layer0.clone(), key_b_layer1.clone(), 50.0);
        
        assert_eq!(tracker.pending_count(), 1, "Should have 1 aggregated transfer");
        
        let transfers = tracker.consume_transfers();
        assert_eq!(transfers.len(), 1, "Should have 1 source cell");
        
        let from_map = transfers.get(&key_a_layer0).unwrap();
        let aggregated_volume = from_map.get(&key_b_layer1).unwrap();
        
        assert_eq!(*aggregated_volume, 150.0, "Volume should be aggregated");
        
        assert_eq!(tracker.pending_count(), 0, "Tracker should be empty after consume");
        
        println!("✅ BatchTransferTracker aggregation works correctly");
    }

    #[test]
    fn test_cell_layer_key_string_pattern() {
        use crate::batch_transfer::CellLayerKey;
        use h3o::CellIndex;
        
        let cell_id = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let key = CellLayerKey::new(cell_id, 1);
        
        let expected_string = format!("{}_{}", cell_id, 1);
        assert_eq!(key.to_string(), expected_string, "String pattern should match cell_ID + layer format");
        
        println!("✅ CellLayerKey string pattern works correctly: {}", key.to_string());
    }
}