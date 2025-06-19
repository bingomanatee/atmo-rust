use crate::asthenosphere_linked::AsthenosphereCellLinked;
use crate::constants::LEVEL_AMT;
use h3o::CellIndex;
use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

/// Direct transfer leveller that moves volume and energy directly from center to neighbors
/// ensuring perfect conservation with no net loss
pub struct AslDirectTransferMapper<'a> {
    root_cell: &'a Rc<RefCell<AsthenosphereCellLinked>>,
    downstream_cells: Vec<&'a Rc<RefCell<AsthenosphereCellLinked>>>,
}

impl<'a> AslDirectTransferMapper<'a> {
    pub fn new(
        root_cell: &'a Rc<RefCell<AsthenosphereCellLinked>>,
        downstream_cells: Vec<&'a Rc<RefCell<AsthenosphereCellLinked>>>,
    ) -> Self {
        Self {
            root_cell,
            downstream_cells,
        }
    }

    /// Construct a direct transfer mapper from AsthSimLinked and a specific CellIndex
    pub fn from_sim_and_cell(
        sim: &'a crate::asth_sim_linked::AsthSimLinked,
        cell_index: CellIndex,
    ) -> Option<Self> {
        // Get the root cell from the simulation
        let root_cell = sim.cells.get(&cell_index)?;

        // Get neighbors for this cell using grid_disk (same as in asthenosphere.rs)
        let neighbors: Vec<CellIndex> = cell_index
            .grid_disk::<Vec<_>>(1)
            .into_iter()
            .filter(|&c| c != cell_index) // Exclude self
            .collect();

        // Find neighbor cells that exist in the simulation
        let downstream_cells: Vec<&Rc<RefCell<AsthenosphereCellLinked>>> = neighbors
            .iter()
            .filter_map(|neighbor_index| sim.cells.get(neighbor_index))
            .collect();

        // Only create mapper if we have neighbors
        if downstream_cells.is_empty() {
            None
        } else {
            Some(Self {
                root_cell,
                downstream_cells,
            })
        }
    }

    /// Level cells by directly transferring volume and energy from center to neighbors
    /// This ensures perfect conservation - what leaves the center exactly equals what neighbors gain
    pub fn level_cells(&self) {
        if self.downstream_cells.is_empty() {
            return; // No neighbors to level with
        }

        // Get current state of all cells
        let root_borrow = self.root_cell.borrow();
        let root_volume = root_borrow.cell.volume;
        let root_energy = root_borrow.cell.energy_j;
        drop(root_borrow);

        // Collect neighbor volumes to determine who needs volume
        let neighbor_volumes: Vec<f64> = self.downstream_cells.iter()
            .map(|cell| cell.borrow().cell.volume)
            .collect();

        // Calculate equilibrium volume (mean of all cells)
        let total_volume = root_volume + neighbor_volumes.iter().sum::<f64>();
        let cell_count = 1.0 + neighbor_volumes.len() as f64;
        let equilibrium_volume = total_volume / cell_count;

        // Only transfer if root has excess volume
        if root_volume <= equilibrium_volume {
            return; // Root doesn't have excess to give
        }

        // Calculate how much volume root should give up
        let volume_excess = root_volume - equilibrium_volume;
        let volume_to_transfer = volume_excess * LEVEL_AMT;

        // Only transfer if we have meaningful volume to move
        if volume_to_transfer <= 0.0 {
            return;
        }

        // Calculate energy density of root cell
        let root_energy_density = if root_volume > 0.0 {
            root_energy / root_volume
        } else {
            0.0
        };

        // Find neighbors that need volume (below equilibrium)
        let mut needy_neighbors: Vec<(usize, f64)> = Vec::new();
        let mut total_need = 0.0;

        for (i, &neighbor_volume) in neighbor_volumes.iter().enumerate() {
            if neighbor_volume < equilibrium_volume {
                let need = equilibrium_volume - neighbor_volume;
                needy_neighbors.push((i, need));
                total_need += need;
            }
        }

        if needy_neighbors.is_empty() || total_need <= 0.0 {
            return; // No neighbors need volume
        }

        // Limit transfer to what we can actually give
        let actual_volume_to_transfer = volume_to_transfer.min(total_need);
        let energy_to_transfer = root_energy_density * actual_volume_to_transfer;

        // Perform the direct transfers
        self.execute_direct_transfers(
            actual_volume_to_transfer,
            energy_to_transfer,
            &needy_neighbors,
            total_need,
        );
    }

    /// Execute direct 1:1 transfers from root to neighbors
    /// Each unit of volume/energy taken from root goes directly to a specific neighbor
    fn execute_direct_transfers(
        &self,
        total_volume_to_transfer: f64,
        _total_energy_to_transfer: f64, // Not used in 1:1 transfer
        needy_neighbors: &[(usize, f64)],
        total_need: f64,
    ) {
        // Calculate energy density for 1:1 transfer
        let root_borrow = self.root_cell.borrow();
        let root_volume = root_borrow.cell.volume;
        let root_energy_density = if root_volume > 0.0 {
            root_borrow.cell.energy_j / root_volume
        } else {
            0.0
        };
        drop(root_borrow);

        // Perform direct 1:1 transfers to each neighbor
        for &(neighbor_index, need) in needy_neighbors {
            let neighbor_cell = &self.downstream_cells[neighbor_index];

            // Calculate this neighbor's share based on their need
            let share_fraction = need / total_need;
            let volume_for_neighbor = total_volume_to_transfer * share_fraction;
            let energy_for_neighbor = volume_for_neighbor * root_energy_density; // Direct 1:1 calculation

            // DIRECT 1:1 TRANSFER: Remove from root and add to neighbor in one atomic operation
            {
                let root_borrow = self.root_cell.borrow_mut();
                let neighbor_borrow = neighbor_cell.borrow_mut();

                if let (Some(root_next), Some(neighbor_next)) = (&root_borrow.next, &neighbor_borrow.next) {
                    let mut root_next_ref = root_next.borrow_mut();
                    let mut neighbor_next_ref = neighbor_next.borrow_mut();

                    // Atomic 1:1 transfer - what leaves root exactly equals what enters neighbor
                    root_next_ref.cell.volume -= volume_for_neighbor;
                    root_next_ref.cell.energy_j -= energy_for_neighbor;

                    neighbor_next_ref.cell.volume += volume_for_neighbor;
                    neighbor_next_ref.cell.energy_j += energy_for_neighbor;

                    // Safety checks
                    if root_next_ref.cell.volume < 0.0 {
                        root_next_ref.cell.volume = 0.0;
                    }
                    if root_next_ref.cell.energy_j < 0.0 {
                        root_next_ref.cell.energy_j = 0.0;
                    }
                }
            }
        }
    }
}

/// Convenience function to level a collection of cells using direct transfer
pub fn level_cells_direct_transfer(
    cells: &HashMap<CellIndex, Rc<RefCell<AsthenosphereCellLinked>>>,
) {
    for (cell_index, root_cell) in cells {
        // Get neighbors for this cell using grid_disk (like in asthenosphere.rs)
        let neighbors: Vec<CellIndex> = cell_index
            .grid_disk::<Vec<_>>(1)
            .into_iter()
            .filter(|&c| c != *cell_index) // Exclude self
            .collect();

        // Find neighbor cells that exist in our collection
        let downstream_cells: Vec<&Rc<RefCell<AsthenosphereCellLinked>>> = neighbors
            .iter()
            .filter_map(|neighbor_index| cells.get(neighbor_index))
            .collect();

        if !downstream_cells.is_empty() {
            let mapper = AslDirectTransferMapper::new(root_cell, downstream_cells);
            mapper.level_cells();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::asthenosphere::AsthenosphereCell;
    use crate::constants::AVG_STARTING_VOLUME_KM_3;

    #[test]
    fn test_direct_transfer_conservation() {
        // Create test cells with different volumes
        let root_cell = AsthenosphereCell {
            id: h3o::CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap(),
            step: 1,
            volume: AVG_STARTING_VOLUME_KM_3 * 2.0, // High volume
            energy_j: 1000.0,
            neighbors: vec![],
            anomaly_energy: 0.0,
            anomaly_volume: 0.0,
        };

        let neighbor1_cell = AsthenosphereCell {
            id: h3o::CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap(),
            step: 1,
            volume: AVG_STARTING_VOLUME_KM_3 * 0.5, // Low volume
            energy_j: 200.0,
            neighbors: vec![],
            anomaly_energy: 0.0,
            anomaly_volume: 0.0,
        };

        let neighbor2_cell = AsthenosphereCell {
            id: h3o::CellIndex::try_from(0x8a1fb46622effff_u64).unwrap(),
            step: 1,
            volume: AVG_STARTING_VOLUME_KM_3 * 0.8, // Medium volume
            energy_j: 300.0,
            neighbors: vec![],
            anomaly_energy: 0.0,
            anomaly_volume: 0.0,
        };

        // Create linked cells
        let root_linked = Rc::new(RefCell::new(
            AsthenosphereCellLinked::new(root_cell.clone())
        ));
        let neighbor1_linked = Rc::new(RefCell::new(
            AsthenosphereCellLinked::new(neighbor1_cell.clone())
        ));
        let neighbor2_linked = Rc::new(RefCell::new(
            AsthenosphereCellLinked::new(neighbor2_cell.clone())
        ));

        // Add next cells for all
        let root_next = AsthenosphereCellLinked::add(&root_linked);
        let neighbor1_next = AsthenosphereCellLinked::add(&neighbor1_linked);
        let neighbor2_next = AsthenosphereCellLinked::add(&neighbor2_linked);

        // Calculate initial totals
        let initial_total_volume = root_cell.volume + neighbor1_cell.volume + neighbor2_cell.volume;
        let initial_total_energy = root_cell.energy_j + neighbor1_cell.energy_j + neighbor2_cell.energy_j;

        // Perform levelling
        let downstream_cells = vec![&neighbor1_linked, &neighbor2_linked];
        let mapper = AslDirectTransferMapper::new(&root_linked, downstream_cells);
        mapper.level_cells();

        // Calculate final totals
        let final_root_volume = root_next.borrow().cell.volume;
        let final_root_energy = root_next.borrow().cell.energy_j;
        let final_neighbor1_volume = neighbor1_next.borrow().cell.volume;
        let final_neighbor1_energy = neighbor1_next.borrow().cell.energy_j;
        let final_neighbor2_volume = neighbor2_next.borrow().cell.volume;
        let final_neighbor2_energy = neighbor2_next.borrow().cell.energy_j;

        let final_total_volume = final_root_volume + final_neighbor1_volume + final_neighbor2_volume;
        let final_total_energy = final_root_energy + final_neighbor1_energy + final_neighbor2_energy;

        // Verify conservation
        assert!((initial_total_volume - final_total_volume).abs() < 10.0, 
               "Volume not conserved: {} -> {}", initial_total_volume, final_total_volume);
        assert!((initial_total_energy - final_total_energy).abs() < 10.0, 
               "Energy not conserved: {} -> {}", initial_total_energy, final_total_energy);

        // Verify that high-volume cell gave up volume
        assert!(final_root_volume < root_cell.volume, "Root should have lost volume");
        
        // Verify that low-volume cells gained volume
        assert!(final_neighbor1_volume > neighbor1_cell.volume, "Neighbor1 should have gained volume");
        assert!(final_neighbor2_volume > neighbor2_cell.volume, "Neighbor2 should have gained volume");
    }
}
