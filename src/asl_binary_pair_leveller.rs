use crate::asthenosphere_linked::AsthenosphereCellLinked;
use crate::constants::LEVEL_AMT;
use h3o::CellIndex;
use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::rc::Rc;

/// Binary pair for levelling - represents two cells that will equalize towards each other
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct BinaryPair {
    pub cell_a: CellIndex,
    pub cell_b: CellIndex,
}

impl BinaryPair {
    pub fn new(cell_a: CellIndex, cell_b: CellIndex) -> Self {
        // Ensure consistent ordering to avoid duplicate pairs (A,B) and (B,A)
        if cell_a < cell_b {
            Self { cell_a, cell_b }
        } else {
            Self { cell_a: cell_b, cell_b: cell_a }
        }
    }
}

/// Binary pair leveller that equalizes pairs of cells towards equilibrium
pub struct AslBinaryPairLeveller {
    pairs: Vec<BinaryPair>,
}

impl AslBinaryPairLeveller {
    /// Create a new binary pair leveller with pre-computed pairs
    pub fn new(pairs: Vec<BinaryPair>) -> Self {
        Self { pairs }
    }

    /// Generate binary pairs from a collection of cells
    /// Each cell is paired with its immediate neighbors
    pub fn generate_pairs_from_cells(
        cells: &HashMap<CellIndex, Rc<RefCell<AsthenosphereCellLinked>>>,
    ) -> Vec<BinaryPair> {
        let mut pairs = HashSet::new();
        
        for &cell_index in cells.keys() {
            // Get neighbors for this cell using grid_disk
            let neighbors: Vec<CellIndex> = cell_index
                .grid_disk::<Vec<_>>(1)
                .into_iter()
                .filter(|&c| c != cell_index) // Exclude self
                .filter(|neighbor| cells.contains_key(neighbor)) // Only include existing cells
                .collect();

            // Create binary pairs with each neighbor
            for neighbor in neighbors {
                let pair = BinaryPair::new(cell_index, neighbor);
                pairs.insert(pair);
            }
        }

        pairs.into_iter().collect()
    }

    /// Level all binary pairs towards equilibrium
    pub fn level_pairs(
        &self,
        cells: &HashMap<CellIndex, Rc<RefCell<AsthenosphereCellLinked>>>,
    ) {
        for pair in &self.pairs {
            self.level_binary_pair(pair, cells);
        }
    }

    /// Level a single binary pair towards equilibrium
    fn level_binary_pair(
        &self,
        pair: &BinaryPair,
        cells: &HashMap<CellIndex, Rc<RefCell<AsthenosphereCellLinked>>>,
    ) {
        // Get both cells
        let cell_a = match cells.get(&pair.cell_a) {
            Some(cell) => cell,
            None => return, // Cell doesn't exist, skip
        };
        let cell_b = match cells.get(&pair.cell_b) {
            Some(cell) => cell,
            None => return, // Cell doesn't exist, skip
        };

        // Get current volumes and energies
        let (volume_a, energy_a) = {
            let borrowed_a = cell_a.borrow();
            (borrowed_a.cell.volume, borrowed_a.cell.energy_j)
        };

        let (volume_b, energy_b) = {
            let borrowed_b = cell_b.borrow();
            (borrowed_b.cell.volume, borrowed_b.cell.energy_j)
        };

        // Calculate equilibrium volumes (mean of the two)
        let total_volume = volume_a + volume_b;
        let equilibrium_volume = total_volume / 2.0;

        // Calculate how much to transfer (based on LEVEL_AMT)
        // Account for the fact that each cell participates in ~6 binary pairs
        // So we need to scale down the transfer to avoid over-correction
        let volume_diff_a = equilibrium_volume - volume_a;
        let neighbor_count_adjustment = 3.0; // because the transfer should be from the part of the 
        // cell that is nearest the data we are only transferring a portion of the volume 
        let volume_to_transfer = (volume_diff_a * LEVEL_AMT) / neighbor_count_adjustment;

        // Only transfer if there's a meaningful difference
        if volume_to_transfer.abs() < 1e-10 {
            return;
        }

        // Calculate energy transfer as a fraction of the source cell's total energy
        // energy_to_transfer = (volume_moved / source_volume) * source_energy
        let energy_to_transfer = if volume_to_transfer > 0.0 {
            // A is gaining volume from B, so energy comes from B
            // Transfer fraction: volume_to_transfer / volume_b
            if volume_b > 0.0 {
                let fraction = volume_to_transfer / volume_b;
                fraction * energy_b
            } else {
                0.0
            }
        } else {
            // A is losing volume to B, so energy comes from A
            // Transfer fraction: |volume_to_transfer| / volume_a
            if volume_a > 0.0 {
                let fraction = volume_to_transfer.abs() / volume_a;
                -(fraction * energy_a) // Negative because A is losing energy
            } else {
                0.0
            }
        };

        // Perform the direct 1:1 transfer
        self.execute_binary_transfer(
            cell_a,
            cell_b,
            volume_to_transfer,
            energy_to_transfer,
        );
    }

    /// Execute the direct 1:1 transfer between two cells
    fn execute_binary_transfer(
        &self,
        cell_a: &Rc<RefCell<AsthenosphereCellLinked>>,
        cell_b: &Rc<RefCell<AsthenosphereCellLinked>>,
        volume_to_transfer: f64,
        energy_to_transfer: f64,
    ) {
        // Borrow both cells mutably for the transfer
        let cell_a_borrow = cell_a.borrow_mut();
        let cell_b_borrow = cell_b.borrow_mut();

        if let (Some(next_a), Some(next_b)) = (&cell_a_borrow.next, &cell_b_borrow.next) {
            let mut next_a_ref = next_a.borrow_mut();
            let mut next_b_ref = next_b.borrow_mut();

            // Direct 1:1 transfer: what leaves one cell exactly enters the other
            next_a_ref.cell.volume += volume_to_transfer;
            next_a_ref.cell.energy_j += energy_to_transfer;

            next_b_ref.cell.volume -= volume_to_transfer;
            next_b_ref.cell.energy_j -= energy_to_transfer;

            // Safety checks to prevent negative values
            if next_a_ref.cell.volume < 0.0 {
                next_a_ref.cell.volume = 0.0;
            }
            if next_a_ref.cell.energy_j < 0.0 {
                next_a_ref.cell.energy_j = 0.0;
            }
            if next_b_ref.cell.volume < 0.0 {
                next_b_ref.cell.volume = 0.0;
            }
            if next_b_ref.cell.energy_j < 0.0 {
                next_b_ref.cell.energy_j = 0.0;
            }
        }
    }

    /// Get the number of binary pairs
    pub fn pair_count(&self) -> usize {
        self.pairs.len()
    }

    /// Get a reference to the pairs (for debugging/inspection)
    pub fn pairs(&self) -> &[BinaryPair] {
        &self.pairs
    }
}

/// Convenience function to create and use binary pair levelling
pub fn level_cells_binary_pairs(
    cells: &HashMap<CellIndex, Rc<RefCell<AsthenosphereCellLinked>>>,
    pairs: &[BinaryPair],
) {
    let leveller = AslBinaryPairLeveller::new(pairs.to_vec());
    leveller.level_pairs(cells);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::asthenosphere::AsthenosphereCell;
    use crate::constants::AVG_STARTING_VOLUME_KM_3;

    #[test]
    fn test_binary_pair_ordering() {
        let cell_a = h3o::CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let cell_b = h3o::CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap();

        let pair1 = BinaryPair::new(cell_a, cell_b);
        let pair2 = BinaryPair::new(cell_b, cell_a);

        // Should be the same regardless of order
        assert_eq!(pair1, pair2);
    }

    #[test]
    fn test_binary_pair_conservation() {
        // Create test cells with different volumes
        let cell_a_data = AsthenosphereCell {
            id: h3o::CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap(),
            step: 1,
            volume: AVG_STARTING_VOLUME_KM_3 * 2.0, // High volume
            energy_j: 1000.0,
            neighbors: vec![],
            anomaly_energy: 0.0,
            anomaly_volume: 0.0,
        };

        let cell_b_data = AsthenosphereCell {
            id: h3o::CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap(),
            step: 1,
            volume: AVG_STARTING_VOLUME_KM_3 * 0.5, // Low volume
            energy_j: 200.0,
            neighbors: vec![],
            anomaly_energy: 0.0,
            anomaly_volume: 0.0,
        };

        // Create linked cells
        let cell_a_linked = Rc::new(RefCell::new(
            AsthenosphereCellLinked::new(cell_a_data.clone())
        ));
        let cell_b_linked = Rc::new(RefCell::new(
            AsthenosphereCellLinked::new(cell_b_data.clone())
        ));

        // Add next cells
        let cell_a_next = AsthenosphereCellLinked::add(&cell_a_linked);
        let cell_b_next = AsthenosphereCellLinked::add(&cell_b_linked);

        // Create cells map
        let mut cells = HashMap::new();
        cells.insert(cell_a_data.id, cell_a_linked);
        cells.insert(cell_b_data.id, cell_b_linked);

        // Calculate initial totals
        let initial_total_volume = cell_a_data.volume + cell_b_data.volume;
        let initial_total_energy = cell_a_data.energy_j + cell_b_data.energy_j;

        // Create binary pair and level
        let pair = BinaryPair::new(cell_a_data.id, cell_b_data.id);
        let leveller = AslBinaryPairLeveller::new(vec![pair]);
        leveller.level_pairs(&cells);

        // Calculate final totals
        let final_volume_a = cell_a_next.borrow().cell.volume;
        let final_energy_a = cell_a_next.borrow().cell.energy_j;
        let final_volume_b = cell_b_next.borrow().cell.volume;
        let final_energy_b = cell_b_next.borrow().cell.energy_j;

        let final_total_volume = final_volume_a + final_volume_b;
        let final_total_energy = final_energy_a + final_energy_b;

        // Verify conservation
        assert!((initial_total_volume - final_total_volume).abs() < 1e-10, 
               "Volume not conserved: {} -> {}", initial_total_volume, final_total_volume);
        assert!((initial_total_energy - final_total_energy).abs() < 1e-10, 
               "Energy not conserved: {} -> {}", initial_total_energy, final_total_energy);

        // Verify that volumes moved towards equilibrium
        let equilibrium = initial_total_volume / 2.0;
        assert!(final_volume_a < cell_a_data.volume, "High volume cell should have lost volume");
        assert!(final_volume_b > cell_b_data.volume, "Low volume cell should have gained volume");
        
        // Verify they moved towards equilibrium but didn't overshoot (with LEVEL_AMT < 1.0)
        assert!(final_volume_a > equilibrium, "Should not overshoot equilibrium");
        assert!(final_volume_b < equilibrium, "Should not overshoot equilibrium");
    }
}
