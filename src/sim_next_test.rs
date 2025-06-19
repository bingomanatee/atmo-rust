use crate::constants::EARTH;
use crate::rock_store::RockStore;
use crate::sim_next::SimNext;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sim_next_creation() {
        let planet = EARTH.clone();
        let store = RockStore::open("test_sim_next.db").expect("Failed to create store");

        let mut sim = SimNext::new(planet, store);
        
        // Test initialization
        sim.initialize_cells();
        
        assert!(sim.cell_count() > 0, "Should have cells after initialization");
        assert!(sim.binary_pairs.len() > 0, "Should have binary pairs after initialization");
        assert_eq!(sim.current_step(), 0, "Should start at step 0");
        
        println!("✅ SimNext created with {} cells and {} pairs",
                 sim.cell_count(), sim.binary_pairs.len());
    }

    #[test]
    fn test_sim_next_step() {
        let planet = EARTH.clone();
        let store = RockStore::open("test_sim_next_step.db").expect("Failed to create store");

        let mut sim = SimNext::new(planet, store);
        sim.initialize_cells();

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
        let store = RockStore::open("test_sim_next_multi.db").expect("Failed to create store");

        let mut sim = SimNext::new(planet, store);
        sim.initialize_cells();

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
        let store = RockStore::open("test_binary_pairs.db").expect("Failed to create store");

        let mut sim = SimNext::new(planet, store);
        sim.initialize_cells();

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
        let store = RockStore::open("test_conservation.db").expect("Failed to create store");

        let mut sim = SimNext::new(planet, store);
        sim.initialize_cells();
        
        // Calculate initial total volume and energy
        let initial_volume: f64 = sim.cells.values().map(|cell| cell.next_volume()).sum();
        let initial_energy: f64 = sim.cells.values().map(|cell| cell.next_energy()).sum();
        
        println!("Initial totals - Volume: {:.2}, Energy: {:.2e}", initial_volume, initial_energy);
        
        // Run a few steps
        for i in 0..3 {
            sim.run_step();
            
            let current_volume: f64 = sim.cells.values().map(|cell| cell.next_volume()).sum();
            let current_energy: f64 = sim.cells.values().map(|cell| cell.next_energy()).sum();
            
            println!("Step {} totals - Volume: {:.2}, Energy: {:.2e}", 
                     i + 1, current_volume, current_energy);
            
            // Allow for small floating point differences and random material changes
            let volume_diff_percent = ((current_volume - initial_volume).abs() / initial_volume) * 100.0;
            let energy_diff_percent = ((current_energy - initial_energy).abs() / initial_energy) * 100.0;
            
            // Should be close to conserved (within 5% due to random material changes)
            assert!(volume_diff_percent < 5.0, 
                    "Volume should be approximately conserved (diff: {:.2}%)", volume_diff_percent);
            assert!(energy_diff_percent < 5.0, 
                    "Energy should be approximately conserved (diff: {:.2}%)", energy_diff_percent);
        }
        
        println!("✅ Volume and energy approximately conserved through simulation steps");
    }

    #[test]
    fn test_conservative_volume_transfers() {
        let planet = EARTH.clone();
        let store = RockStore::open("test_transfers.db").expect("Failed to create store");

        let mut sim = SimNext::new(planet, store);
        sim.initialize_cells();

        // Get two cells for testing
        let cell_ids: Vec<_> = sim.cells.keys().take(2).cloned().collect();
        if cell_ids.len() < 2 {
            panic!("Need at least 2 cells for transfer test");
        }

        let cell_a_id = cell_ids[0];
        let cell_b_id = cell_ids[1];

        // Record initial state
        let initial_vol_a = sim.cells[&cell_a_id].next_volume();
        let initial_vol_b = sim.cells[&cell_b_id].next_volume();
        let initial_energy_a = sim.cells[&cell_a_id].next_energy();
        let initial_energy_b = sim.cells[&cell_b_id].next_energy();

        let total_initial_volume = initial_vol_a + initial_vol_b;
        let total_initial_energy = initial_energy_a + initial_energy_b;

        // Perform a direct transfer using our new method
        let transfer_volume = initial_vol_a * 0.1; // Transfer 10%
        sim.transfer_volume_between_cells(&cell_a_id, &cell_b_id, transfer_volume);

        // Check results
        let final_vol_a = sim.cells[&cell_a_id].next_volume();
        let final_vol_b = sim.cells[&cell_b_id].next_volume();
        let final_energy_a = sim.cells[&cell_a_id].next_energy();
        let final_energy_b = sim.cells[&cell_b_id].next_energy();

        let total_final_volume = final_vol_a + final_vol_b;
        let total_final_energy = final_energy_a + final_energy_b;

        // Verify conservation
        assert!((total_initial_volume - total_final_volume).abs() < 1e-10,
                "Volume should be perfectly conserved");
        assert!((total_initial_energy - total_final_energy).abs() < 1e6,
                "Energy should be approximately conserved");

        // Verify transfer occurred
        assert!(final_vol_a < initial_vol_a, "Source cell should lose volume");
        assert!(final_vol_b > initial_vol_b, "Target cell should gain volume");
        assert!(final_vol_a >= 0.0, "Source volume should never go negative");
        assert!(final_energy_a >= 0.0, "Source energy should never go negative");

        println!("✅ Conservative volume transfers working correctly");
        println!("   Volume conserved: {:.2} -> {:.2}", total_initial_volume, total_final_volume);
        println!("   Energy conserved: {:.2e} -> {:.2e}", total_initial_energy, total_final_energy);
    }

    #[test]
    fn test_zero_volume_safety() {
        let planet = EARTH.clone();
        let store = RockStore::open("test_zero_volume.db").expect("Failed to create store");

        let mut sim = SimpleSimulation::new(planet, store);
        sim.initialize_cells();

        // Get a cell and try to transfer zero volume
        let cell_id = *sim.cells.keys().next().expect("Need at least one cell");
        let other_id = *sim.cells.keys().nth(1).expect("Need at least two cells");

        let initial_volume = sim.cells[&cell_id].next_volume();

        // Try to transfer zero volume (should be ignored)
        sim.transfer_volume_between_cells(&cell_id, &other_id, 0.0);

        let final_volume = sim.cells[&cell_id].next_volume();

        assert_eq!(initial_volume, final_volume, "Zero volume transfer should be ignored");

        // Try to transfer negative volume (should be ignored)
        sim.transfer_volume_between_cells(&cell_id, &other_id, -100.0);

        let final_volume_2 = sim.cells[&cell_id].next_volume();

        assert_eq!(initial_volume, final_volume_2, "Negative volume transfer should be ignored");

        println!("✅ Zero and negative volume transfers safely ignored");
    }
}
