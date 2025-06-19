use crate::asthenosphere_linked::AsthenosphereCellLinked;
use crate::constants::LEVEL_AMT;
use h3o::CellIndex;
use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::rc::Rc;

/// Binary pair for levelling - represents two cells that will equalize towards each other
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct BinaryPair {
    pub cell_a: CellIndex,
    pub cell_b: CellIndex,
}

impl BinaryPair {
    pub fn new(cell_a: CellIndex, cell_b: CellIndex) -> Self {
        // Ensure consistent ordering: first id should be higher than second to avoid redundant pairs
        if cell_a > cell_b {
            Self { cell_a, cell_b }
        } else {
            Self { cell_a: cell_b, cell_b: cell_a }
        }
    }

    /// Create a unique string ID for this binary pair
    pub fn to_string_id(&self) -> String {
        format!("{}:{}", self.cell_a, self.cell_b)
    }

    /// Create a BinaryPair from a string ID (for debugging/testing)
    pub fn from_string_id(id: &str) -> Result<Self, String> {
        let parts: Vec<&str> = id.split(':').collect();
        if parts.len() != 2 {
            return Err(format!("Invalid pair ID format: {}", id));
        }

        let cell_a = parts[0].parse::<u64>()
            .map_err(|_| format!("Invalid cell_a in pair ID: {}", parts[0]))?;
        let cell_b = parts[1].parse::<u64>()
            .map_err(|_| format!("Invalid cell_b in pair ID: {}", parts[1]))?;

        let cell_a_index = CellIndex::try_from(cell_a)
            .map_err(|_| format!("Invalid CellIndex for cell_a: {}", cell_a))?;
        let cell_b_index = CellIndex::try_from(cell_b)
            .map_err(|_| format!("Invalid CellIndex for cell_b: {}", cell_b))?;

        Ok(BinaryPair::new(cell_a_index, cell_b_index))
    }
}

impl fmt::Display for BinaryPair {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}", self.cell_a, self.cell_b)
    }
}

/// Binary pair leveller that equalizes pairs of cells towards equilibrium
pub struct AslBinaryPairLeveller {
    pairs: Vec<BinaryPair>,
    // Efficient neighbor lookup: cell_id -> set of neighbor_ids
    neighbor_map: HashMap<CellIndex, HashSet<CellIndex>>,
}

impl AslBinaryPairLeveller {
    /// Create a new binary pair leveller with pre-computed pairs
    pub fn new(pairs: Vec<BinaryPair>) -> Self {
        // Build efficient neighbor lookup map
        let mut neighbor_map: HashMap<CellIndex, HashSet<CellIndex>> = HashMap::new();

        for pair in &pairs {
            neighbor_map.entry(pair.cell_a).or_insert_with(HashSet::new).insert(pair.cell_b);
            neighbor_map.entry(pair.cell_b).or_insert_with(HashSet::new).insert(pair.cell_a);
        }

        println!("🔗 Built neighbor map for {} cells", neighbor_map.len());

        Self { pairs, neighbor_map }
    }

    /// Generate binary pairs from a collection of cells
    /// Each cell is paired with its immediate neighbors, ensuring no duplicate pairs
    /// Uses string IDs to guarantee absolute uniqueness
    pub fn generate_pairs_from_cells(
        cells: &HashMap<CellIndex, Rc<RefCell<AsthenosphereCellLinked>>>,
    ) -> Vec<BinaryPair> {
        println!("🔍 Starting binary pair generation for {} cells", cells.len());
        let mut pair_map: HashMap<String, BinaryPair> = HashMap::new();
        let mut processed_cells = 0;

        for &cell_index in cells.keys() {
            processed_cells += 1;
            if processed_cells % 100 == 0 {
                println!("🔍 Processed {} cells, found {} pairs so far", processed_cells, pair_map.len());
            }

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
                let pair_id = pair.to_string_id();

                // Only insert if we haven't seen this pair before
                if !pair_map.contains_key(&pair_id) {
                    pair_map.insert(pair_id, pair);
                }
            }
        }

        println!("🔍 Finished binary pair generation: {} pairs from {} cells", pair_map.len(), cells.len());
        pair_map.into_values().collect()
    }

    /// Level all binary pairs towards equilibrium
    pub fn level_pairs(
        &self,
        cells: &HashMap<CellIndex, Rc<RefCell<AsthenosphereCellLinked>>>,
    ) {
        println!("⚖️  Starting binary pair levelling for {} pairs", self.pairs.len());
        let mut processed = 0;
        for pair in &self.pairs {
            processed += 1;
            if processed % 5000 == 0 {
                println!("⚖️  Processed {} pairs so far", processed);
            }
            self.level_binary_pair(pair, cells);
        }
        println!("⚖️  Completed binary pair levelling for {} pairs", processed);
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

    /// Hole balancing: find neighbors within 2 steps that are significantly different
    /// and move volume between them to achieve 75% equilibration
    pub fn hole_balance(
        &self,
        cells: &HashMap<CellIndex, Rc<RefCell<AsthenosphereCellLinked>>>,
    ) {
        println!("🕳️  Starting hole balancing for {} cells", cells.len());

        // Build volume/energy cache for fast lookups
        println!("🕳️  Building volume cache...");
        let mut volume_cache: HashMap<CellIndex, (f64, f64)> = HashMap::new();
        for (&cell_id, cell) in cells {
            let borrowed = cell.borrow();
            if let Some(next) = &borrowed.next {
                let next_ref = next.borrow();
                volume_cache.insert(cell_id, (next_ref.cell.volume, next_ref.cell.energy_j));
            }
        }
        println!("🕳️  Built volume cache for {} cells", volume_cache.len());

        let mut processed = 0;
        for &cell_id in cells.keys() {
            processed += 1;
            if processed % 1000 == 0 {
                println!("🕳️  Hole balanced {} cells so far", processed);
            }
            self.hole_balance_cell_cached(cell_id, cells, &mut volume_cache);
        }
        println!("🕳️  Completed hole balancing for {} cells", processed);
    }

    /// Hole balance a single cell with its 2-step neighbors (EFFICIENT VERSION)
    fn hole_balance_cell(
        &self,
        cell_id: CellIndex,
        cells: &HashMap<CellIndex, Rc<RefCell<AsthenosphereCellLinked>>>,
    ) {
        // EFFICIENT: Get 1-step neighbors using O(1) HashMap lookup
        let one_step_neighbors = match self.neighbor_map.get(&cell_id) {
            Some(neighbors) => neighbors,
            None => return, // No neighbors found
        };

        // EFFICIENT: Find 2-step neighbors using O(neighbors) lookup instead of O(pairs)
        let mut two_step_neighbors = HashSet::new();
        for &neighbor in one_step_neighbors {
            if let Some(neighbor_neighbors) = self.neighbor_map.get(&neighbor) {
                for &two_step_neighbor in neighbor_neighbors {
                    // Exclude self and 1-step neighbors
                    if two_step_neighbor != cell_id && !one_step_neighbors.contains(&two_step_neighbor) {
                        two_step_neighbors.insert(two_step_neighbor);
                    }
                }
            }
        }

        if two_step_neighbors.is_empty() {
            return;
        }

        // Get current cell's next volume (after binary levelling)
        let current_volume = {
            let cell = match cells.get(&cell_id) {
                Some(cell) => cell,
                None => return,
            };
            let borrowed = cell.borrow();
            if let Some(next) = &borrowed.next {
                next.borrow().cell.volume
            } else {
                return;
            }
        };

        // Find the most different neighbor within 2 steps
        let mut best_neighbor: Option<CellIndex> = None;
        let mut best_difference = 0.0;

        for &neighbor_id in &two_step_neighbors {
            let neighbor_volume = {
                let neighbor = match cells.get(&neighbor_id) {
                    Some(neighbor) => neighbor,
                    None => continue,
                };
                let borrowed = neighbor.borrow();
                if let Some(next) = &borrowed.next {
                    next.borrow().cell.volume
                } else {
                    continue;
                }
            };

            // Calculate percentage difference
            let avg_volume = (current_volume + neighbor_volume) / 2.0;
            if avg_volume > 0.0 {
                let percent_diff = ((current_volume - neighbor_volume).abs() / avg_volume) * 100.0;

                // Must be more than 2% different
                if percent_diff > 2.0 && percent_diff > best_difference {
                    best_difference = percent_diff;
                    best_neighbor = Some(neighbor_id);
                }
            }
        }

        // If we found a suitable neighbor, perform hole balancing
        if let Some(neighbor_id) = best_neighbor {
            self.execute_hole_balance(cell_id, neighbor_id, cells);
        }
    }

    /// Hole balance a single cell with its 2-step neighbors (CACHED VERSION)
    fn hole_balance_cell_cached(
        &self,
        cell_id: CellIndex,
        cells: &HashMap<CellIndex, Rc<RefCell<AsthenosphereCellLinked>>>,
        volume_cache: &mut HashMap<CellIndex, (f64, f64)>,
    ) {
        // EFFICIENT: Get 1-step neighbors using O(1) HashMap lookup
        let one_step_neighbors = match self.neighbor_map.get(&cell_id) {
            Some(neighbors) => neighbors,
            None => return, // No neighbors found
        };

        // EFFICIENT: Find 2-step neighbors using O(neighbors) lookup instead of O(pairs)
        let mut two_step_neighbors = HashSet::new();
        for &neighbor in one_step_neighbors {
            if let Some(neighbor_neighbors) = self.neighbor_map.get(&neighbor) {
                for &two_step_neighbor in neighbor_neighbors {
                    // Exclude self and 1-step neighbors
                    if two_step_neighbor != cell_id && !one_step_neighbors.contains(&two_step_neighbor) {
                        two_step_neighbors.insert(two_step_neighbor);
                    }
                }
            }
        }

        if two_step_neighbors.is_empty() {
            return;
        }

        // FAST: Get current cell's volume from cache
        let current_volume = match volume_cache.get(&cell_id) {
            Some(&(volume, _)) => volume,
            None => return,
        };

        // Find the most different neighbor within 2 steps using cache
        let mut best_neighbor: Option<CellIndex> = None;
        let mut best_difference = 0.0;

        for &neighbor_id in &two_step_neighbors {
            let neighbor_volume = match volume_cache.get(&neighbor_id) {
                Some(&(volume, _)) => volume,
                None => continue,
            };

            // Calculate percentage difference
            let avg_volume = (current_volume + neighbor_volume) / 2.0;
            if avg_volume > 0.0 {
                let percent_diff = ((current_volume - neighbor_volume).abs() / avg_volume) * 100.0;

                // Must be more than 2% different
                if percent_diff > 2.0 && percent_diff > best_difference {
                    best_difference = percent_diff;
                    best_neighbor = Some(neighbor_id);
                }
            }
        }

        // If we found a suitable neighbor, perform hole balancing
        if let Some(neighbor_id) = best_neighbor {
            self.execute_hole_balance_cached(cell_id, neighbor_id, cells, volume_cache);
        }
    }

    /// Execute hole balancing between two cells (75% equilibration)
    fn execute_hole_balance(
        &self,
        cell_a_id: CellIndex,
        cell_b_id: CellIndex,
        cells: &HashMap<CellIndex, Rc<RefCell<AsthenosphereCellLinked>>>,
    ) {
        let cell_a = match cells.get(&cell_a_id) {
            Some(cell) => cell,
            None => return,
        };
        let cell_b = match cells.get(&cell_b_id) {
            Some(cell) => cell,
            None => return,
        };

        // Get current volumes and energies from next cells (after binary levelling)
        let (volume_a, energy_a) = {
            let borrowed_a = cell_a.borrow();
            if let Some(next) = &borrowed_a.next {
                let next_ref = next.borrow();
                (next_ref.cell.volume, next_ref.cell.energy_j)
            } else {
                return;
            }
        };

        let (volume_b, energy_b) = {
            let borrowed_b = cell_b.borrow();
            if let Some(next) = &borrowed_b.next {
                let next_ref = next.borrow();
                (next_ref.cell.volume, next_ref.cell.energy_j)
            } else {
                return;
            }
        };

        // Calculate 75% equilibration
        let total_volume = volume_a + volume_b;
        let equilibrium_volume = total_volume / 2.0;
        let volume_diff_a = equilibrium_volume - volume_a;
        let volume_to_transfer = volume_diff_a * 0.75; // 75% equilibration

        // Only transfer if there's a meaningful difference
        if volume_to_transfer.abs() < 1e-10 {
            return;
        }

        // Calculate energy transfer using fractional approach
        let energy_to_transfer = if volume_to_transfer > 0.0 {
            // A is gaining volume from B
            if volume_b > 0.0 {
                let fraction = volume_to_transfer / volume_b;
                fraction * energy_b
            } else {
                0.0
            }
        } else {
            // A is losing volume to B
            if volume_a > 0.0 {
                let fraction = volume_to_transfer.abs() / volume_a;
                -(fraction * energy_a)
            } else {
                0.0
            }
        };

        // Execute the transfer on next cells
        self.execute_hole_balance_transfer(
            cell_a,
            cell_b,
            volume_to_transfer,
            energy_to_transfer,
        );
    }

    /// Execute the hole balance transfer between two cells
    fn execute_hole_balance_transfer(
        &self,
        cell_a: &Rc<RefCell<AsthenosphereCellLinked>>,
        cell_b: &Rc<RefCell<AsthenosphereCellLinked>>,
        volume_to_transfer: f64,
        energy_to_transfer: f64,
    ) {
        let cell_a_borrow = cell_a.borrow();
        let cell_b_borrow = cell_b.borrow();

        if let (Some(next_a), Some(next_b)) = (&cell_a_borrow.next, &cell_b_borrow.next) {
            let mut next_a_ref = next_a.borrow_mut();
            let mut next_b_ref = next_b.borrow_mut();

            // Direct 1:1 transfer for hole balancing
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

    /// Execute hole balancing between two cells (75% equilibration) - CACHED VERSION
    fn execute_hole_balance_cached(
        &self,
        cell_a_id: CellIndex,
        cell_b_id: CellIndex,
        cells: &HashMap<CellIndex, Rc<RefCell<AsthenosphereCellLinked>>>,
        volume_cache: &mut HashMap<CellIndex, (f64, f64)>,
    ) {
        // FAST: Get volumes and energies from cache
        let (volume_a, energy_a) = match volume_cache.get(&cell_a_id) {
            Some(&(vol, energy)) => (vol, energy),
            None => return,
        };

        let (volume_b, energy_b) = match volume_cache.get(&cell_b_id) {
            Some(&(vol, energy)) => (vol, energy),
            None => return,
        };

        // Calculate 75% equilibration
        let total_volume = volume_a + volume_b;
        let equilibrium_volume = total_volume / 2.0;
        let volume_diff_a = equilibrium_volume - volume_a;
        let volume_to_transfer = volume_diff_a * 0.75; // 75% equilibration

        // Only transfer if there's a meaningful difference
        if volume_to_transfer.abs() < 1e-10 {
            return;
        }

        // Calculate energy transfer using fractional approach
        let energy_to_transfer = if volume_to_transfer > 0.0 {
            // A is gaining volume from B
            if volume_b > 0.0 {
                let fraction = volume_to_transfer / volume_b;
                fraction * energy_b
            } else {
                0.0
            }
        } else {
            // A is losing volume to B
            if volume_a > 0.0 {
                let fraction = volume_to_transfer.abs() / volume_a;
                -(fraction * energy_a)
            } else {
                0.0
            }
        };

        // Calculate new values
        let new_volume_a = volume_a + volume_to_transfer;
        let new_energy_a = energy_a + energy_to_transfer;
        let new_volume_b = volume_b - volume_to_transfer;
        let new_energy_b = energy_b - energy_to_transfer;

        // Safety checks to prevent negative values
        let final_volume_a = if new_volume_a < 0.0 { 0.0 } else { new_volume_a };
        let final_energy_a = if new_energy_a < 0.0 { 0.0 } else { new_energy_a };
        let final_volume_b = if new_volume_b < 0.0 { 0.0 } else { new_volume_b };
        let final_energy_b = if new_energy_b < 0.0 { 0.0 } else { new_energy_b };

        // Update cache
        volume_cache.insert(cell_a_id, (final_volume_a, final_energy_a));
        volume_cache.insert(cell_b_id, (final_volume_b, final_energy_b));

        // Update actual cells
        if let Some(cell_a) = cells.get(&cell_a_id) {
            let cell_a_borrow = cell_a.borrow();
            if let Some(next_a) = &cell_a_borrow.next {
                let mut next_a_ref = next_a.borrow_mut();
                next_a_ref.cell.volume = final_volume_a;
                next_a_ref.cell.energy_j = final_energy_a;
            }
        }

        if let Some(cell_b) = cells.get(&cell_b_id) {
            let cell_b_borrow = cell_b.borrow();
            if let Some(next_b) = &cell_b_borrow.next {
                let mut next_b_ref = next_b.borrow_mut();
                next_b_ref.cell.volume = final_volume_b;
                next_b_ref.cell.energy_j = final_energy_b;
            }
        }
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
