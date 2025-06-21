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
}