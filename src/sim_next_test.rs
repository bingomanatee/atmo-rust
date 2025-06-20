use crate::constants::EARTH;
use crate::rock_store::RockStore;
use crate::sim_next::{SimNext, SimNextProps};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sim_next_creation() {
        let planet = EARTH.clone();
        let store = RockStore::open("data/test_sim_next.db").expect("Failed to create store");

        let props = SimNextProps::new(planet, store);
        let sim = SimNext::new(props);
        
        assert!(sim.cell_count() > 0, "Should have cells after initialization");
        assert!(sim.binary_pairs.len() > 0, "Should have binary pairs after initialization");
        assert_eq!(sim.current_step(), 0, "Should start at step 0");
        
        println!("✅ SimNext created with {} cells and {} pairs",
                 sim.cell_count(), sim.binary_pairs.len());
    }

    #[test]
    fn test_sim_next_step() {
        let planet = EARTH.clone();
        let store = RockStore::open("data/test_sim_next_step.db").expect("Failed to create store");

        let props = SimNextProps::new(planet, store);
        let mut sim = SimNext::new(props);

        let initial_step = sim.current_step();
        let initial_cell_count = sim.cell_count();

        // Run one step
        sim.run_step();

        assert_eq!(sim.current_step(), initial_step + 1, "Step should increment");
        assert_eq!(sim.cell_count(), initial_cell_count, "Cell count should remain same");

        println!("✅ SimNext step completed successfully");
    }

    #[test]
    fn test_sim_next_multiple_steps() {
        let planet = EARTH.clone();
        let store = RockStore::open("data/test_sim_next_multi.db").expect("Failed to create store");

        let props = SimNextProps::new(planet, store);
        let mut sim = SimNext::new(props);

        let steps_to_run = 5;

        for i in 0..steps_to_run {
            println!("Running step {}", i);
            sim.run_step();
        }

        assert_eq!(sim.current_step(), steps_to_run, "Should complete all steps");

        println!("✅ SimNext completed {} steps successfully", steps_to_run);
    }

    #[test]
    fn test_binary_pair_uniqueness() {
        let planet = EARTH.clone();
        let store = RockStore::open("data/test_binary_pairs.db").expect("Failed to create store");

        let props = SimNextProps::new(planet, store);
        let sim = SimNext::new(props);

        // Check that all binary pairs are unique
        let mut pair_ids = std::collections::HashSet::new();
        let mut duplicate_count = 0;

        for pair in &sim.binary_pairs {
            let pair_id = pair.to_string_id();
            if !pair_ids.insert(pair_id.clone()) {
                duplicate_count += 1;
                println!("Duplicate pair found: {}", pair_id);
            }
        }

        assert_eq!(duplicate_count, 0, "Should have no duplicate binary pairs");
        assert_eq!(pair_ids.len(), sim.binary_pairs.len(), "All pairs should be unique");

        println!("✅ All {} binary pairs are unique", sim.binary_pairs.len());
    }

    #[test]
    fn test_cell_volume_conservation() {
        let planet = EARTH.clone();
        let store = RockStore::open("data/test_conservation.db").expect("Failed to create store");

        let props = SimNextProps::new(planet, store);
        let mut sim = SimNext::new(props);
        
        // Calculate initial total volume and energy
        let initial_volume: f64 = sim.cells.values().map(|cell| cell.next_cell.volume).sum();
        let initial_energy: f64 = sim.cells.values().map(|cell| cell.next_cell.energy_j).sum();
        
        println!("Initial totals - Volume: {:.2}, Energy: {:.2e}", initial_volume, initial_energy);
        
        // Run a few steps
        for i in 0..3 {
            sim.run_step();
            
            let current_volume: f64 = sim.cells.values().map(|cell| cell.next_cell.volume).sum();
            let current_energy: f64 = sim.cells.values().map(|cell| cell.next_cell.energy_j).sum();
            
            println!("Step {} totals - Volume: {:.2}, Energy: {:.2e}", 
                     i + 1, current_volume, current_energy);
            
            // Allow for small floating point differences and random material changes
            let volume_diff_percent = ((current_volume - initial_volume).abs() / initial_volume) * 100.0;
            let energy_diff_percent = ((current_energy - initial_energy).abs() / initial_energy) * 100.0;
            
            // Should be close to conserved (within 5% due to random material changes)
            assert!(volume_diff_percent < 5.0, 
                    "Volume should be approximately conserved (diff: {:.2}%)", volume_diff_percent);
            assert!(energy_diff_percent < 10.0, 
                    "Energy should be approximately conserved (diff: {:.2}%)", energy_diff_percent);
        }
        
        println!("✅ Volume and energy approximately conserved through simulation steps");
    }

    #[test]
    fn test_cell_data_access() {
        let planet = EARTH.clone();
        let store = RockStore::open("data/test_cell_access.db").expect("Failed to create store");

        let props = SimNextProps::new(planet, store);
        let sim = SimNext::new(props);

        // Test accessing cell data
        for (cell_id, cell) in sim.cells.iter().take(5) {
            assert!(cell.next_cell.volume > 0.0, "Cell volume should be positive");
            assert!(cell.next_cell.energy_j > 0.0, "Cell energy should be positive");
            assert!(!cell.cell.neighbors.is_empty(), "Cell should have neighbors");
            
            println!("Cell {}: Volume={:.2}, Energy={:.2e}, Neighbors={}", 
                     cell_id, cell.next_cell.volume, cell.next_cell.energy_j, cell.cell.neighbors.len());
        }

        println!("✅ Cell data access working correctly");
    }
}