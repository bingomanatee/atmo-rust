use crate::asthenosphere::AsthenosphereCell;
use crate::constants::{ANOMALY_DECAY_RATE, ANOMALY_VOLUME_AMOUNT, COOLING_RATE, JOULES_PER_KM3};
use h3o::CellIndex;
use rand::Rng;

/// Simulation cell with current state and next state
/// Uses full AsthenosphereCell for both states for maximum extensibility
#[derive(Clone, Debug)]
pub struct AsthenosphereCellNext {
    /// Current committed asthenosphere cell state
    pub cell: AsthenosphereCell,

    /// Next step asthenosphere cell state (working state)
    pub next_cell: AsthenosphereCell,
}

impl AsthenosphereCellNext {
    /// Create a new simulation cell from an asthenosphere cell
    pub fn new(cell: AsthenosphereCell) -> Self {
        Self {
            next_cell: cell.clone(),
            cell,
        }
    }



    /// Commit next state to current state (call at end of simulation step)
    pub fn commit_step(&mut self) {
        self.cell = self.next_cell.clone();
        self.next_cell.step = self.cell.step + 1;
    }

    /// Get cell index
    pub fn cell_index(&self) -> CellIndex {
        self.cell.id
    }

    /// Move volume and proportional energy to another cell (conservative binary operation)
    /// Returns true if transfer occurred, false if no transfer (invalid volume or insufficient source)
    pub fn move_volume_to(&mut self, other: &mut AsthenosphereCellNext, volume: f64) -> bool {
        self.next_cell.transfer_volume(volume, &mut other.next_cell)
    }

    /// Move a percentage of volume to another cell
    /// Returns true if transfer occurred, false if no transfer
    pub fn move_volume_fraction_to(&mut self, other: &mut AsthenosphereCellNext, fraction: f64) -> bool {
        self.next_cell.transfer_volume_fraction(fraction, &mut other.next_cell)
    }

    /// Cool the cell using the global cooling rate
    pub fn cool(&mut self) {
        self.next_cell.energy_j *= *COOLING_RATE;
    }

    /// Try to add a new anomaly to this cell
    /// Returns false if cell already has an anomaly (no new anomaly added)
    /// Returns true if a new anomaly was successfully added
    pub fn add_anomaly(&mut self) -> bool {
        // Return false if cell already has an anomaly
        if self.has_anomaly() {
            return false;
        }

        // Generate random anomaly parameters using constants
        let mut rng = rand::rng();
        let is_volcano = rng.random_bool(); // 50/50 chance volcano vs sinkhole
        let intensity_factor = rng.random_range(0.1..=1.0); // 10% to 100% of base amount

        let volume_change = if is_volcano {
            ANOMALY_VOLUME_AMOUNT * intensity_factor
        } else {
            -ANOMALY_VOLUME_AMOUNT * intensity_factor
        };

        // Add the anomaly (energy will be derived when volume is applied)
        self.next_cell.anomaly_volume = volume_change;

        true
    }

    /// Apply and decay existing anomaly effects
    pub fn process_anomaly(&mut self) {
        if self.next_cell.anomaly_volume.abs() < 1e-10 {
            return; // No anomaly to process
        }

        // Apply the anomaly effect to volume and derive energy automatically
        let volume_change = self.next_cell.anomaly_volume;
        self.next_cell.volume = (self.next_cell.volume + volume_change).max(0.0);
        self.next_cell.energy_j = (self.next_cell.energy_j + volume_change * JOULES_PER_KM3).max(0.0);

        // Decay the anomaly
        self.next_cell.anomaly_volume *= (1.0 - ANOMALY_DECAY_RATE);

        // Clear anomaly if it's become negligible
        if self.next_cell.anomaly_volume.abs() < 1e-10 {
            self.next_cell.anomaly_volume = 0.0;
        }
    }

    /// Check if this cell has an active anomaly
    pub fn has_anomaly(&self) -> bool {
        self.next_cell.anomaly_volume.abs() > 1e-10
    }

    /// Apply and decay existing anomaly effects
    pub fn process_anomaly(&mut self) {
        if self.next_cell.anomaly_volume.abs() < 1e-10 {
            return; // No anomaly to process
        }

        // Apply the anomaly effect to volume and derive energy automatically
        let volume_change = self.next_cell.anomaly_volume;
        self.next_cell.volume = (self.next_cell.volume + volume_change).max(0.0);
        self.next_cell.energy_j = (self.next_cell.energy_j + volume_change * JOULES_PER_KM3).max(0.0);

        // Decay the anomaly
        self.next_cell.anomaly_volume *= (1.0 - ANOMALY_DECAY_RATE);

        // Clear anomaly if it's become negligible
        if self.next_cell.anomaly_volume.abs() < 1e-10 {
            self.next_cell.anomaly_volume = 0.0;
        }
    }

    /// Add anomaly to this cell and degraded strength to neighbors
    /// Returns list of neighbor cell IDs that received degraded anomalies
    pub fn add_anomaly_with_neighbors(&mut self) -> (bool, Vec<CellIndex>) {
        // Try to add anomaly to this cell first
        let added_to_self = self.add_anomaly();
        if !added_to_self {
            return (false, vec![]); // Cell already has anomaly
        }

        // Get the anomaly strength that was just added
        let base_strength = self.next_cell.anomaly_volume;
        let degraded_strength = base_strength * 0.3; // 30% strength for neighbors

        // Return the neighbor IDs (caller will handle adding to neighbors)
        (true, self.cell.neighbors.clone())
    }
}

impl From<AsthenosphereCell> for AsthenosphereCellNext {
    fn from(cell: AsthenosphereCell) -> Self {
        Self::new(cell)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use h3o::CellIndex;

    fn create_test_cell(cell_index: CellIndex, volume: f64, energy_j: f64) -> AsthenosphereCellNext {
        let cell = AsthenosphereCell {
            id: cell_index,
            volume,
            energy_j,
            step: 0,
            neighbors: vec![],
            anomaly_volume: 0.0,
        };
        AsthenosphereCellNext::new(cell)
    }

    #[test]
    fn test_simulation_cell_creation() {
        let cell_index = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let sim_cell = create_test_cell(cell_index, 1000.0, 2e12);

        assert_eq!(sim_cell.cell.volume, 1000.0);
        assert_eq!(sim_cell.cell.energy_j, 2e12);
        assert_eq!(sim_cell.next_cell.volume, 1000.0);
        assert_eq!(sim_cell.next_cell.energy_j, 2e12);
        assert_eq!(sim_cell.cell_index(), cell_index);
    }

    #[test]
    fn test_next_state_updates() {
        let cell_index = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let mut sim_cell = create_test_cell(cell_index, 1000.0, 2e12);

        sim_cell.next_cell.volume = 1100.0;
        sim_cell.next_cell.energy_j = 2.2e12;

        assert_eq!(sim_cell.cell.volume, 1000.0); // Current unchanged
        assert_eq!(sim_cell.next_cell.volume, 1100.0);    // Next updated
        assert_eq!(sim_cell.cell.energy_j, 2e12);   // Current unchanged
        assert_eq!(sim_cell.next_cell.energy_j, 2.2e12);    // Next updated
    }

    #[test]
    fn test_commit_step() {
        let cell_index = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let mut sim_cell = create_test_cell(cell_index, 1000.0, 2e12);

        sim_cell.next_cell.volume = 1100.0;
        sim_cell.next_cell.energy_j = 2.2e12;
        sim_cell.next_cell.step = 5;
        sim_cell.commit_step();

        assert_eq!(sim_cell.cell.volume, 1100.0); // Current now updated
        assert_eq!(sim_cell.next_cell.volume, 1100.0);    // Next same as current
        assert_eq!(sim_cell.cell.energy_j, 2.2e12); // Current now updated
        assert_eq!(sim_cell.next_cell.energy_j, 2.2e12);    // Next same as current
        assert_eq!(sim_cell.cell.step, 5);
        assert_eq!(sim_cell.next_cell.step, 6); // Next step incremented
    }

    #[test]
    fn test_move_volume_to_basic() {
        let cell_index_a = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let cell_index_b = CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap();

        let mut cell_a = create_test_cell(cell_index_a, 1000.0, 2e12);
        let mut cell_b = create_test_cell(cell_index_b, 500.0, 1e12);

        // Record initial totals for conservation check
        let initial_total_volume = cell_a.next_cell.volume + cell_b.next_cell.volume;
        let initial_total_energy = cell_a.next_cell.energy_j + cell_b.next_cell.energy_j;

        // Transfer 100 km³ from A to B
        let success = cell_a.move_volume_to(&mut cell_b, 100.0);

        assert!(success, "Transfer should succeed");
        assert_eq!(cell_a.next_cell.volume, 900.0); // 1000 - 100
        assert_eq!(cell_b.next_cell.volume, 600.0); // 500 + 100

        // Energy should transfer proportionally
        // A had 2e12 J / 1000 km³ = 2e9 J/km³ density
        // 100 km³ should transfer 100 * 2e9 = 2e11 J
        assert_eq!(cell_a.next_cell.energy_j, 1.8e12); // 2e12 - 2e11
        assert_eq!(cell_b.next_cell.energy_j, 1.2e12); // 1e12 + 2e11

        // Conservation check: totals should be within 0.1% of initial
        let final_total_volume = cell_a.next_cell.volume + cell_b.next_cell.volume;
        let final_total_energy = cell_a.next_cell.energy_j + cell_b.next_cell.energy_j;

        let volume_diff_percent = ((final_total_volume - initial_total_volume).abs() / initial_total_volume) * 100.0;
        let energy_diff_percent = ((final_total_energy - initial_total_energy).abs() / initial_total_energy) * 100.0;

        assert!(volume_diff_percent < 0.1, "Volume conservation failed: {:.3}% difference", volume_diff_percent);
        assert!(energy_diff_percent < 0.1, "Energy conservation failed: {:.3}% difference", energy_diff_percent);
    }

    #[test]
    fn test_move_volume_to_insufficient_volume() {
        let cell_index_a = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let cell_index_b = CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap();

        let mut cell_a = create_test_cell(cell_index_a, 100.0, 2e11);
        let mut cell_b = create_test_cell(cell_index_b, 500.0, 1e12);

        // Try to transfer more than available
        let success = cell_a.move_volume_to(&mut cell_b, 200.0);

        assert!(!success, "Transfer should fail");
        assert_eq!(cell_a.next_cell.volume, 100.0); // Unchanged
        assert_eq!(cell_b.next_cell.volume, 500.0); // Unchanged
        assert_eq!(cell_a.next_cell.energy_j, 2e11); // Unchanged
        assert_eq!(cell_b.next_cell.energy_j, 1e12); // Unchanged
    }

    #[test]
    fn test_move_volume_to_zero_volume() {
        let cell_index_a = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let cell_index_b = CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap();

        let mut cell_a = create_test_cell(cell_index_a, 1000.0, 2e12);
        let mut cell_b = create_test_cell(cell_index_b, 500.0, 1e12);

        // Try to transfer zero volume
        let success = cell_a.move_volume_to(&mut cell_b, 0.0);

        assert!(!success, "Transfer should fail for zero volume");
        assert_eq!(cell_a.next_cell.volume, 1000.0); // Unchanged
        assert_eq!(cell_b.next_cell.volume, 500.0); // Unchanged
    }

    #[test]
    fn test_move_volume_to_negative_volume() {
        let cell_index_a = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let cell_index_b = CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap();

        let mut cell_a = create_test_cell(cell_index_a, 1000.0, 2e12);
        let mut cell_b = create_test_cell(cell_index_b, 500.0, 1e12);

        // Try to transfer negative volume
        let success = cell_a.move_volume_to(&mut cell_b, -100.0);

        assert!(!success, "Transfer should fail for negative volume");
        assert_eq!(cell_a.next_cell.volume, 1000.0); // Unchanged
        assert_eq!(cell_b.next_cell.volume, 500.0); // Unchanged
    }

    #[test]
    fn test_move_volume_fraction_to() {
        let cell_index_a = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let cell_index_b = CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap();

        let mut cell_a = create_test_cell(cell_index_a, 1000.0, 2e12);
        let mut cell_b = create_test_cell(cell_index_b, 500.0, 1e12);

        // Record initial totals for conservation check
        let initial_total_volume = cell_a.next_cell.volume + cell_b.next_cell.volume;
        let initial_total_energy = cell_a.next_cell.energy_j + cell_b.next_cell.energy_j;

        // Transfer 10% of A's volume to B
        let success = cell_a.move_volume_fraction_to(&mut cell_b, 0.1);

        assert!(success, "Fraction transfer should succeed");
        assert_eq!(cell_a.next_cell.volume, 900.0); // 1000 - 100 (10%)
        assert_eq!(cell_b.next_cell.volume, 600.0); // 500 + 100

        // Energy should transfer proportionally
        assert_eq!(cell_a.next_cell.energy_j, 1.8e12); // 90% of original
        assert_eq!(cell_b.next_cell.energy_j, 1.2e12); // Original + transferred

        // Conservation check: totals should be within 0.1% of initial
        let final_total_volume = cell_a.next_cell.volume + cell_b.next_cell.volume;
        let final_total_energy = cell_a.next_cell.energy_j + cell_b.next_cell.energy_j;

        let volume_diff_percent = ((final_total_volume - initial_total_volume).abs() / initial_total_volume) * 100.0;
        let energy_diff_percent = ((final_total_energy - initial_total_energy).abs() / initial_total_energy) * 100.0;

        assert!(volume_diff_percent < 0.1, "Volume conservation failed: {:.3}% difference", volume_diff_percent);
        assert!(energy_diff_percent < 0.1, "Energy conservation failed: {:.3}% difference", energy_diff_percent);
    }

    #[test]
    fn test_conservation_in_transfer() {
        let cell_index_a = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let cell_index_b = CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap();

        let mut cell_a = SimulationCell::from_cell_index(cell_index_a, 1000.0, 2e12);
        let mut cell_b = SimulationCell::from_cell_index(cell_index_b, 500.0, 1e12);

        // Record initial totals
        let initial_volume = cell_a.next_volume() + cell_b.next_volume();
        let initial_energy = cell_a.next_energy() + cell_b.next_energy();

        // Transfer volume
        cell_a.move_volume_to(&mut cell_b, 200.0);

        // Check conservation
        let final_volume = cell_a.next_volume() + cell_b.next_volume();
        let final_energy = cell_a.next_energy() + cell_b.next_energy();

        assert!((initial_volume - final_volume).abs() < 1e-10, "Volume should be conserved");
        assert!((initial_energy - final_energy).abs() < 1e6, "Energy should be conserved");
    }
}
