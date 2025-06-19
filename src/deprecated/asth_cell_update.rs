use crate::constants::CELL_JOULES_START;
use crate::asthenosphere::AsthenosphereCell;

#[derive(Clone, Debug)]
pub struct AsthCellUpdate {
    pub cell: AsthenosphereCell,
    pub energy_k: f64,
    pub volume: f64,
}
impl AsthCellUpdate {
    pub fn new(cell: AsthenosphereCell) -> AsthCellUpdate {
        AsthCellUpdate {
            cell: cell.clone(),
            energy_k: cell.energy_j,
            volume: cell.volume,
        }
    }

    pub fn is_dirty(&self) -> bool {
         self.energy_k != self.cell.energy_j || self.volume != self.cell.volume
    }
    
    pub fn add_volume(&mut self, volume: f64) ->&mut AsthCellUpdate {
        if volume == 0.0 { return self; }
        if volume < 0.0 { return self.remove(volume); }
        
        let new_volume = self.volume + volume;
        let new_volume_energy = CELL_JOULES_START * volume;
        let old_volume_energy = self.energy_k * self.volume;
        let new_energy = (new_volume_energy + old_volume_energy) / new_volume;

        self.volume = new_volume;
        self.energy_k = new_energy;
        self
    }
    
    pub(crate) fn update(&self) -> AsthenosphereCell {
        AsthenosphereCell {
            volume: self.volume,
            energy_j: self.energy_k,
            ..self.cell.clone()
        }
    }

    pub fn add_by(&mut self, fraction: f64) {
        if fraction < 0.0 {
            return;
        }
        self.add_volume(self.cell.volume * fraction);
    }

    pub fn remove(&mut self, volume: f64) -> &mut AsthCellUpdate {
        let new_volume = self.volume - volume;
        let energy_k = self.energy_k * new_volume / self.volume;

        self.volume = new_volume;
        self.energy_k = energy_k;
        self
    }

    fn remove_by(&mut self, fraction: f64) {
        if fraction < 0.0 {
            return;
        }
        self.remove(self.cell.volume * fraction);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::asthenosphere::AsthenosphereCell;
    use crate::constants::CELL_JOULES_START;

    // Helper to create a default AsthenosphereCell for testing
    fn default_cell() -> AsthenosphereCell {
        AsthenosphereCell {
            volume: 10.0,
            energy_j: 5.0,
            ..Default::default() // assuming AsthenosphereCell implements Default
        }
    }

    #[test]
    fn test_new_and_is_dirty() {
        let cell = default_cell();
        let update = AsthCellUpdate::new(cell.clone());

        // Initially, energy_k and volume should match cell's
        assert_eq!(update.energy_k, cell.energy_j);
        assert_eq!(update.volume, cell.volume);

        // is_dirty should be false initially
        assert!(!update.is_dirty());

        // Modify energy_k to make it dirty
        let mut update2 = update.clone();
        update2.energy_k += 1.0;
        assert!(update2.is_dirty());

        // Modify volume to make it dirty
        let mut update3 = update.clone();
        update3.volume += 1.0;
        assert!(update3.is_dirty());
    }

    #[test]
    fn test_add() {
        let cell = default_cell();
        let mut update = AsthCellUpdate::new(cell);

        let add_volume = 5.0;
        update.add_volume(add_volume);

        // Volume should increase by add_volume
        assert_eq!(update.volume, 10.0 + add_volume);

        // Energy should be weighted average of old and new volume energies
        let old_energy_total = 5.0 * 10.0;
        let new_energy_total = CELL_JOULES_START * add_volume;
        let expected_energy = (old_energy_total + new_energy_total) / (10.0 + add_volume);
        assert!((update.energy_k - expected_energy).abs() < 1e-12);
    }

    #[test]
    fn test_remove() {
        let cell = default_cell();
        let mut update = AsthCellUpdate::new(cell);

        let remove_volume = 4.0;
        update.remove(remove_volume);

        // Volume should decrease by remove_volume
        assert_eq!(update.volume, 10.0 - remove_volume);

        // Energy should scale proportionally
        let expected_energy = 5.0 * (10.0 - remove_volume) / 10.0;
        assert!((update.energy_k - expected_energy).abs() < 1e-12);
    }

    #[test]
    fn test_add_by() {
        let cell = default_cell();
        let mut update = AsthCellUpdate::new(cell.clone());

        let fraction = 0.5;
        update.add_by(fraction);

        let added_volume = cell.volume * fraction;
        assert_eq!(update.volume, cell.volume + added_volume);

        // Energy calculation same as add
        let old_energy_total = cell.energy_j * cell.volume;
        let new_energy_total = CELL_JOULES_START * added_volume;
        let expected_energy = (old_energy_total + new_energy_total) / (cell.volume + added_volume);
        assert!((update.energy_k - expected_energy).abs() < 1e-12);
    }

    #[test]
    fn test_add_by_negative_fraction() {
        let cell = default_cell();
        let mut update = AsthCellUpdate::new(cell.clone());

        let fraction = -0.1;
        update.add_by(fraction);

        // Should not change anything
        assert_eq!(update.volume, cell.volume);
        assert_eq!(update.energy_k, cell.energy_j);
    }

    #[test]
    fn test_remove_by() {
        let cell = default_cell();
        let mut update = AsthCellUpdate::new(cell.clone());

        let fraction = 0.3;
        update.remove_by(fraction);

        let removed_volume = cell.volume * fraction;
        assert_eq!(update.volume, cell.volume - removed_volume);

        let expected_energy = cell.energy_j * (cell.volume - removed_volume) / cell.volume;
        assert!((update.energy_k - expected_energy).abs() < 1e-12);
    }

    #[test]
    fn test_remove_by_negative_fraction() {
        let cell = default_cell();
        let mut update = AsthCellUpdate::new(cell.clone());

        let fraction = -0.2;
        update.remove_by(fraction);

        // Should not change anything
        assert_eq!(update.volume, cell.volume);
        assert_eq!(update.energy_k, cell.energy_j);
    }

    #[test]
    fn test_update_method() {
        let cell = default_cell();
        let mut update = AsthCellUpdate::new(cell.clone());

        update.add_volume(5.0);
        let updated_cell = update.update();

        assert_eq!(updated_cell.volume, update.volume);
        assert_eq!(updated_cell.energy_j, update.energy_k);

        // Other fields should be same as original cell
        // Assuming AsthenosphereCell has other fields, they should be equal
        // Here we just check that the id or some field is same if applicable
    }
}
