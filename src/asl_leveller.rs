use crate::asth_sim_linked::AsthSimLinked;
use crate::asthenosphere::AsthenosphereCell;
use crate::asthenosphere_linked::AsthenosphereCellLinked;
use h3o::CellIndex;
use std::cell::{Ref, RefCell};
use std::rc::Rc;

pub struct AslMapper<'a> {
    map: &'a AsthSimLinked,
    id: &'a CellIndex,
}

impl<'a> AslMapper<'a> {
    pub fn new(map_input: &'a AsthSimLinked, id: &'a CellIndex) -> Self {
        Self { map: map_input, id }
    }

    /// Returns the Rc<RefCell> of the current cell.
    /// The caller is responsible for borrowing from it.
    pub fn cell_rc(&self) -> Rc<RefCell<AsthenosphereCellLinked>> {
        self.map
            .cells
            .get(self.id)
            .expect("Cell not found for given id")
            .clone()
    }

    /// Returns the Rc<RefCell> of the next cell if it exists.
    pub fn next_cell_rc(&self) -> Option<Rc<RefCell<AsthenosphereCellLinked>>> {
        let current_rc = self.cell_rc();
        let current_ref = current_rc.borrow();
        current_ref.next.as_ref().map(|next| next.clone())
    }

    /// Returns the neighbors of the current cell.
    pub fn neighbors(&self) -> Vec<CellIndex> {
        let cell_rc = self.cell_rc();
        let cell_ref = cell_rc.borrow();
        cell_ref.cell.neighbors.clone()
    }

    /// Returns downstream neighbor cells (those with volume less than root).
    pub fn downstream_neighbor_cells(&self) -> Vec<AsthenosphereCell> {
        let root_rc = self.cell_rc();
        let root_ref = root_rc.borrow();
        let cell = &root_ref.cell;
        let root_volume = cell.volume;

        cell.neighbors
            .iter()
            .filter_map(|neighbor_id| {
                let neighbor_rc = self
                    .map
                    .cells
                    .get(neighbor_id)
                    .expect("Neighbor cell not found")
                    .clone();
                let neighbor_ref = neighbor_rc.borrow();
                if neighbor_ref.cell.volume < root_volume {
                    Some(neighbor_ref.cell.clone())
                } else {
                    None
                }
            })
            .collect()
    }

    /// Calculate equilibrium volume as mean of root and downstream neighbors.
    pub fn equilibrium(&self) -> f64 {
        let root_rc = self.cell_rc();
        let root_ref = root_rc.borrow();
        let root_volume = root_ref.cell.volume;

        let downstream_cells = self.downstream_neighbor_cells();

        let total_volume: f64 =
            downstream_cells.iter().map(|c| c.volume).sum::<f64>() + root_volume;
        let count = downstream_cells.len() as f64 + 1.0;

        if count > 0.0 {
            total_volume / count
        } else {
            root_volume
        }
    }

    /// Level the root cell volume towards equilibrium.
    pub fn level(&self) {
        use crate::constants::LEVEL_AMT;

        let root_rc = self.cell_rc();
        let root_ref = root_rc.borrow();
        let root_volume = root_ref.cell.volume;
        let root_energy = root_ref.cell.energy_j;
        let root_energy_per_volume = if root_volume > 0.0 {
            root_energy / root_volume
        } else {
            0.0
        };

        // Check if next cell exists
        if root_ref.next.is_none() {
            return; // Skip if no next cell exists
        }

        // Get reference to the next cell
        let next_rc = root_ref.next.as_ref().unwrap().clone();
        drop(root_ref); // Release borrow before proceeding

        let equilibrium_volume = self.equilibrium();

        let volume_diff = equilibrium_volume - root_volume;
        let volume_to_move = volume_diff * LEVEL_AMT;

        if volume_to_move.abs() < 1e-12 {
            return;
        }

        let downstream_neighbors = self.downstream_neighbor_cells();

        let total_gap: f64 = downstream_neighbors
            .iter()
            .map(|neighbor| (root_volume - neighbor.volume).max(0.0))
            .sum();

        // Update next cell
        let mut next_ref = next_rc.borrow_mut();

        next_ref.cell.volume += volume_to_move;
        if root_volume > 0.0 {
            let energy_change = root_energy * volume_to_move / root_volume;
            next_ref.cell.energy_j += energy_change;
        }
        drop(next_ref); // Release mutable borrow before borrowing neighbors

        for neighbor_cell in downstream_neighbors {
            let neighbor_id = &neighbor_cell.id;
            let neighbor_rc = self
                .map
                .cells
                .get(neighbor_id)
                .expect("Neighbor cell not found")
                .clone();

            let neighbor_ref = neighbor_rc.borrow();

            // Skip if neighbor has no next cell
            if neighbor_ref.next.is_none() {
                continue;
            }

            // Get reference to the neighbor's next cell
            let neighbor_next_rc = neighbor_ref.next.as_ref().unwrap().clone();
            drop(neighbor_ref); // Release borrow before mutable borrow

            let gap = (root_volume - neighbor_cell.volume).max(0.0);
            if gap <= 0.0 || total_gap <= 0.0 {
                continue;
            }
            let mut neighbor_next_ref = neighbor_next_rc.borrow_mut();

            let fraction_of_moved_volume = gap / total_gap;
            let neighbor_volume_add = volume_to_move * fraction_of_moved_volume;
            neighbor_next_ref.cell.volume += neighbor_volume_add;
            neighbor_next_ref.cell.energy_j += root_energy_per_volume * neighbor_volume_add;
        }
    }
}
