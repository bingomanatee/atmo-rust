use crate::asthenosphere::AsthenosphereCell;
use crate::constants::{
    AVG_STARTING_VOLUME_KM_3, BACK_FILL_LEVEL, BACK_FILL_RATE, CELL_JOULES_EQUILIBRIUM,
    COOLING_RATE, JOULES_PER_KM3, LAYER_COUNT, SINKHOLE_CHANCE, SINKHOLE_DECAY_RATE,
    SINKHOLE_MAX_VOLUME_TO_REMOVE, SINKHOLE_MIN_VOLUME_TO_REMOVE, VOLCANO_BIAS, VOLCANO_CHANCE,
    VOLCANO_DECAY_RATE, VOLCANO_JOULES_PER_VOLUME, VOLCANO_MAX_VOLUME,
    VOLCANO_MIN_VOLUME, VOLCANO_ROTATION_SKEW_PERCENTAGE,
};
use crate::geoconverter::GeoCellConverter;
use crate::planet::Planet;
use h3o::CellIndex;
use noise::{NoiseFn, Perlin};
use rand::Rng;

#[derive(Clone, Debug)]
pub struct AsthenosphereCellNext {
    pub cell: AsthenosphereCell,
    pub next_cell: AsthenosphereCell,
}

impl AsthenosphereCellNext {
    pub fn new(cell: AsthenosphereCell) -> Self {
        Self {
            next_cell: cell.clone(),
            cell,
        }
    }

    pub fn commit_step(&mut self) {
        self.cell = self.next_cell.clone();
        self.next_cell.step = self.cell.step + 1;
    }

    pub fn cell_index(&self) -> CellIndex {
        self.cell.id
    }

    pub fn move_volume_to(&mut self, other: &mut AsthenosphereCellNext, volume: f64) -> bool {
        self.next_cell.transfer_volume(volume, &mut other.next_cell)
    }

    pub fn move_volume_fraction_to(
        &mut self,
        other: &mut AsthenosphereCellNext,
        fraction: f64,
    ) -> bool {
        self.next_cell
            .transfer_volume_fraction(fraction, &mut other.next_cell)
    }

    pub fn cool_with_heating(&mut self, planet: &Planet, step: u32, resolution: h3o::Resolution) {
        // Apply cooling to all layers
        for layer_idx in 0..LAYER_COUNT {
            self.next_cell.energy_layers[layer_idx] *= *COOLING_RATE;
        }
        //  let heating_energy = Self::calculate_perlin_heating(self.cell_index(), planet, step, resolution);
        //  self.next_cell.energy_layers[0] += heating_energy; // Apply heating to bottom layer
    }

    pub fn back_fill(&mut self) {
        if self.has_any_anomaly() {
            return;
        }

        let trigger_rate = AVG_STARTING_VOLUME_KM_3 * BACK_FILL_LEVEL;
        // Apply back-fill to all layers
        for layer_idx in 0..LAYER_COUNT {
            if self.cell.volume_layers[layer_idx] < trigger_rate {
                let diff = trigger_rate - self.cell.volume_layers[layer_idx];
                let back_fill = diff * BACK_FILL_RATE;
                self.next_cell.volume_layers[layer_idx] += back_fill;
                self.next_cell.energy_layers[layer_idx] += back_fill * JOULES_PER_KM3;
            }
        }
    }

    fn calculate_perlin_heating(
        cell_index: CellIndex,
        planet: &Planet,
        step: u32,
        resolution: h3o::Resolution,
    ) -> f64 {
        let noise_seed = (step % 1000) as u32;
        let noise = Perlin::new(noise_seed);
        let gc = GeoCellConverter::new(planet.radius_km as f64, resolution);
        let location = gc.cell_to_vec3(cell_index);
        let noise_scale = 10.0;
        let scaled_location = location.normalize() * noise_scale;
        let noise_value = noise.get([
            scaled_location.x as f64,
            scaled_location.y as f64,
            scaled_location.z as f64,
        ]);

        // lower the temperature a bit to allow for excess energy aded by volcano bias
        let base_heating = (1.0 - *COOLING_RATE) * CELL_JOULES_EQUILIBRIUM * (1.0 - VOLCANO_BIAS);
        let variation = 1.0 + (noise_value * 0.02);
        base_heating * variation
    }

    pub fn cool(&mut self) {
        // Apply cooling to all layers
        for layer_idx in 0..LAYER_COUNT {
            self.next_cell.energy_layers[layer_idx] *= *COOLING_RATE;
        }
    }

    pub fn try_add_volcano(&mut self) -> bool {
        if self.has_volcano() {
            return false;
        }

        let mut rng = rand::rng();
        if rng.random::<f64>() > VOLCANO_CHANCE {
            return false;
        }

        let volume = rng.random_range(VOLCANO_MIN_VOLUME..=VOLCANO_MAX_VOLUME);
        self.next_cell.volcano_volume = volume;
        true
    }

    pub fn try_add_sinkhole(&mut self) -> bool {
        if self.has_sinkhole() {
            return false;
        }

        let mut rng = rand::rng();
        if rng.random::<f64>() > SINKHOLE_CHANCE {
            return false;
        }

        let volume_to_remove =
            rng.random_range(SINKHOLE_MIN_VOLUME_TO_REMOVE..=SINKHOLE_MAX_VOLUME_TO_REMOVE);
        self.next_cell.sinkhole_volume = volume_to_remove;
        true
    }

    pub fn add_volcano_with_volume(&mut self, volume: f64) -> bool {
        if self.has_volcano() {
            return false;
        }
        self.next_cell.volcano_volume = volume;
        true
    }

    pub fn add_sinkhole_with_volume(&mut self, volume: f64) -> bool {
        if self.has_sinkhole() {
            return false;
        }
        self.next_cell.sinkhole_volume = volume;
        true
    }

    pub fn has_volcano(&self) -> bool {
        self.next_cell.volcano_volume > 1e-10
    }

    pub fn has_sinkhole(&self) -> bool {
        self.next_cell.sinkhole_volume > 1e-10
    }

    pub fn has_any_anomaly(&self) -> bool {
        self.has_volcano() || self.has_sinkhole()
    }

    pub fn process_volcanoes_and_sinkholes(&mut self) -> Option<(f64, f64)> {
        let surface_layer = LAYER_COUNT - 1; // Apply volcano/sinkhole effects to surface layer
        let mut rotation_skew_data = None;
        
        if self.has_volcano() {
            let volcano_volume = self.next_cell.volcano_volume;
            
            // Calculate the amount that gets skewed from rotation (30%)
            let skewed_volume = volcano_volume * VOLCANO_ROTATION_SKEW_PERCENTAGE;
            let skewed_energy = skewed_volume * VOLCANO_JOULES_PER_VOLUME;
            
            // Add the remaining 70% normally to the surface layer
            let remaining_volume = volcano_volume - skewed_volume;
            let remaining_energy = remaining_volume * VOLCANO_JOULES_PER_VOLUME;
            
            self.next_cell.volume_layers[surface_layer] = (self.next_cell.volume_layers[surface_layer] + remaining_volume).max(0.0);
            self.next_cell.energy_layers[surface_layer] = (self.next_cell.energy_layers[surface_layer] + remaining_energy).max(0.0);

            // Store the skewed portion for rotation processing
            if skewed_volume > 0.0 {
                rotation_skew_data = Some((skewed_volume, skewed_energy));
            }

            self.next_cell.volcano_volume *= VOLCANO_DECAY_RATE;
            if self.next_cell.volcano_volume < 10.0 {
                self.next_cell.volcano_volume = 0.0;
            }
        }

        if self.has_sinkhole() {
            let sinkhole_volume = self.next_cell.sinkhole_volume;
            self.next_cell.volume_layers[surface_layer] = (self.next_cell.volume_layers[surface_layer] - sinkhole_volume).max(0.0);
            self.next_cell.energy_layers[surface_layer] =
                (self.next_cell.energy_layers[surface_layer] - sinkhole_volume * JOULES_PER_KM3).max(0.0);

            self.next_cell.sinkhole_volume *= SINKHOLE_DECAY_RATE;
            if self.next_cell.sinkhole_volume < 1e-10 {
                self.next_cell.sinkhole_volume = 0.0;
            }
        }
        
        rotation_skew_data
    }

    pub fn add_volcano_cluster(&mut self, volume: f64) {
        let surface_layer = LAYER_COUNT - 1;
        self.next_cell.volume_layers[surface_layer] = (self.next_cell.volume_layers[surface_layer] + volume).max(0.0);
        self.next_cell.energy_layers[surface_layer] = (self.next_cell.energy_layers[surface_layer] + volume * JOULES_PER_KM3).max(0.0);
    }

    pub fn add_sinkhole_cluster(&mut self, volume: f64) {
        let surface_layer = LAYER_COUNT - 1;
        self.next_cell.volume_layers[surface_layer] = (self.next_cell.volume_layers[surface_layer] - volume).max(0.0);
        self.next_cell.energy_layers[surface_layer] = (self.next_cell.energy_layers[surface_layer] - volume * JOULES_PER_KM3).max(0.0);
    }

    pub fn try_create_volcano_cluster(&mut self) -> Option<Vec<(CellIndex, f64)>> {
        if !self.try_add_volcano() {
            return None;
        }

        let mut rng = rand::rng();
        if rng.random::<f64>() > 0.1 {
            return None;
        }

        let base_volume = self.next_cell.volcano_volume;
        let cluster_size = rng.random_range(2..=8);
        let mut neighbors_within_4: Vec<CellIndex> = Vec::new();

        for &neighbor in &self.cell.neighbors {
            neighbors_within_4.push(neighbor);
            if neighbors_within_4.len() >= 4 {
                break;
            }
        }

        let mut cluster_volcanoes = Vec::new();
        for i in 0..cluster_size.min(neighbors_within_4.len()) {
            let neighbor_id = neighbors_within_4[i];
            let strength = rng.random_range(0.0..=0.8);
            let volcano_volume = base_volume * strength;
            cluster_volcanoes.push((neighbor_id, volcano_volume));
        }

        Some(cluster_volcanoes)
    }

    pub fn try_create_massive_sink_event(&mut self) -> Option<Vec<(CellIndex, f64)>> {
        if !self.try_add_sinkhole() {
            return None;
        }

        let mut rng = rand::rng();
        if rng.random::<f64>() > 0.02 {
            return None;
        }

        let base_volume = self.next_cell.sinkhole_volume;
        let mut neighbor_sinkholes = Vec::new();

        for &neighbor_id in &self.cell.neighbors {
            let sinkhole_volume = base_volume * 0.5;
            neighbor_sinkholes.push((neighbor_id, sinkhole_volume));
        }

        Some(neighbor_sinkholes)
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
    use crate::h3_utils::H3Utils;
    use h3o::CellIndex;

    fn create_test_cell(
        cell_index: CellIndex,
        volume: f64,
        energy_j: f64,
    ) -> AsthenosphereCellNext {
        let cell = AsthenosphereCell {
            id: cell_index,
            volume_layers: vec![volume; LAYER_COUNT],
            energy_layers: vec![energy_j; LAYER_COUNT],
            step: 0,
            neighbors: vec![],
            volcano_volume: 0.0,
            sinkhole_volume: 0.0,
        };
        AsthenosphereCellNext::new(cell)
    }

    #[test]
    fn test_simulation_cell_creation() {
        let cell_index = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let sim_cell = create_test_cell(cell_index, 1000.0, 2e12);

        assert_eq!(sim_cell.cell.volume_layers[0], 1000.0);
        assert_eq!(sim_cell.cell.energy_layers[0], 2e12);
        assert_eq!(sim_cell.next_cell.volume_layers[0], 1000.0);
        assert_eq!(sim_cell.next_cell.energy_layers[0], 2e12);
        assert_eq!(sim_cell.cell_index(), cell_index);
    }

    #[test]
    fn test_next_state_updates() {
        let cell_index = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let mut sim_cell = create_test_cell(cell_index, 1000.0, 2e12);

        sim_cell.next_cell.volume_layers[0] = 1100.0;
        sim_cell.next_cell.energy_layers[0] = 2.2e12;

        assert_eq!(sim_cell.cell.volume_layers[0], 1000.0); // Current unchanged
        assert_eq!(sim_cell.next_cell.volume_layers[0], 1100.0); // Next updated
        assert_eq!(sim_cell.cell.energy_layers[0], 2e12); // Current unchanged
        assert_eq!(sim_cell.next_cell.energy_layers[0], 2.2e12); // Next updated
    }

    #[test]
    fn test_commit_step() {
        let cell_index = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let mut sim_cell = create_test_cell(cell_index, 1000.0, 2e12);

        sim_cell.next_cell.volume_layers[0] = 1100.0;
        sim_cell.next_cell.energy_layers[0] = 2.2e12;
        sim_cell.next_cell.step = 5;
        sim_cell.commit_step();

        assert_eq!(sim_cell.cell.volume_layers[0], 1100.0); // Current now updated
        assert_eq!(sim_cell.next_cell.volume_layers[0], 1100.0); // Next same as current
        assert_eq!(sim_cell.cell.energy_layers[0], 2.2e12); // Current now updated
        assert_eq!(sim_cell.next_cell.energy_layers[0], 2.2e12); // Next same as current
        assert_eq!(sim_cell.cell.step, 5);
        assert_eq!(sim_cell.next_cell.step, 6); // Next step incremented
    }

    #[test]
    fn test_move_volume_to_basic() {
        let cell_index_a = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let cell_index_b = CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap();

        let mut cell_a = create_test_cell(cell_index_a, 1000.0, 2e12);
        let mut cell_b = create_test_cell(cell_index_b, 500.0, 1e12);

        // Record initial totals for conservation check (using surface layer)
        let surface_layer = LAYER_COUNT - 1;
        let initial_total_volume = cell_a.next_cell.volume_layers[surface_layer] + cell_b.next_cell.volume_layers[surface_layer];
        let initial_total_energy = cell_a.next_cell.energy_layers[surface_layer] + cell_b.next_cell.energy_layers[surface_layer];

        // Transfer 100 km³ from A to B
        let success = cell_a.move_volume_to(&mut cell_b, 100.0);

        assert!(success, "Transfer should succeed");
        assert_eq!(cell_a.next_cell.volume_layers[surface_layer], 900.0); // 1000 - 100
        assert_eq!(cell_b.next_cell.volume_layers[surface_layer], 600.0); // 500 + 100

        // Energy should transfer proportionally
        // A had 2e12 J / 1000 km³ = 2e9 J/km³ density
        // 100 km³ should transfer 100 * 2e9 = 2e11 J
        assert_eq!(cell_a.next_cell.energy_layers[surface_layer], 1.8e12); // 2e12 - 2e11
        assert_eq!(cell_b.next_cell.energy_layers[surface_layer], 1.2e12); // 1e12 + 2e11

        // Conservation check: totals should be within 0.1% of initial
        let final_total_volume = cell_a.next_cell.volume_layers[surface_layer] + cell_b.next_cell.volume_layers[surface_layer];
        let final_total_energy = cell_a.next_cell.energy_layers[surface_layer] + cell_b.next_cell.energy_layers[surface_layer];

        let volume_diff_percent =
            ((final_total_volume - initial_total_volume).abs() / initial_total_volume) * 100.0;
        let energy_diff_percent =
            ((final_total_energy - initial_total_energy).abs() / initial_total_energy) * 100.0;

        assert!(
            volume_diff_percent < 0.1,
            "Volume conservation failed: {:.3}% difference",
            volume_diff_percent
        );
        assert!(
            energy_diff_percent < 0.1,
            "Energy conservation failed: {:.3}% difference",
            energy_diff_percent
        );
    }

    #[test]
    fn test_move_volume_to_insufficient_volume() {
        let cell_index_a = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let cell_index_b = CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap();

        let mut cell_a = create_test_cell(cell_index_a, 100.0, 2e11);
        let mut cell_b = create_test_cell(cell_index_b, 500.0, 1e12);

        // Try to transfer more than available
        let success = cell_a.move_volume_to(&mut cell_b, 200.0);

        let surface_layer = LAYER_COUNT - 1;
        assert!(!success, "Transfer should fail");
        assert_eq!(cell_a.next_cell.volume_layers[surface_layer], 100.0); // Unchanged
        assert_eq!(cell_b.next_cell.volume_layers[surface_layer], 500.0); // Unchanged
        assert_eq!(cell_a.next_cell.energy_layers[surface_layer], 2e11); // Unchanged
        assert_eq!(cell_b.next_cell.energy_layers[surface_layer], 1e12); // Unchanged
    }

    #[test]
    fn test_move_volume_to_zero_volume() {
        let cell_index_a = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let cell_index_b = CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap();

        let mut cell_a = create_test_cell(cell_index_a, 1000.0, 2e12);
        let mut cell_b = create_test_cell(cell_index_b, 500.0, 1e12);

        // Try to transfer zero volume
        let success = cell_a.move_volume_to(&mut cell_b, 0.0);

        let surface_layer = LAYER_COUNT - 1;
        assert!(!success, "Transfer should fail for zero volume");
        assert_eq!(cell_a.next_cell.volume_layers[surface_layer], 1000.0); // Unchanged
        assert_eq!(cell_b.next_cell.volume_layers[surface_layer], 500.0); // Unchanged
    }

    #[test]
    fn test_move_volume_to_negative_volume() {
        let cell_index_a = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let cell_index_b = CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap();

        let mut cell_a = create_test_cell(cell_index_a, 1000.0, 2e12);
        let mut cell_b = create_test_cell(cell_index_b, 500.0, 1e12);

        // Try to transfer negative volume
        let success = cell_a.move_volume_to(&mut cell_b, -100.0);

        let surface_layer = LAYER_COUNT - 1;
        assert!(!success, "Transfer should fail for negative volume");
        assert_eq!(cell_a.next_cell.volume_layers[surface_layer], 1000.0); // Unchanged
        assert_eq!(cell_b.next_cell.volume_layers[surface_layer], 500.0); // Unchanged
    }

    #[test]
    fn test_move_volume_fraction_to() {
        let cell_index_a = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let cell_index_b = CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap();

        let mut cell_a = create_test_cell(cell_index_a, 1000.0, 2e12);
        let mut cell_b = create_test_cell(cell_index_b, 500.0, 1e12);

        // Record initial totals for conservation check
        let surface_layer = LAYER_COUNT - 1;
        let initial_total_volume = cell_a.next_cell.volume_layers[surface_layer] + cell_b.next_cell.volume_layers[surface_layer];
        let initial_total_energy = cell_a.next_cell.energy_layers[surface_layer] + cell_b.next_cell.energy_layers[surface_layer];

        // Transfer 10% of A's volume to B
        let success = cell_a.move_volume_fraction_to(&mut cell_b, 0.1);

        assert!(success, "Fraction transfer should succeed");
        assert_eq!(cell_a.next_cell.volume_layers[surface_layer], 900.0); // 1000 - 100 (10%)
        assert_eq!(cell_b.next_cell.volume_layers[surface_layer], 600.0); // 500 + 100

        // Energy should transfer proportionally
        assert_eq!(cell_a.next_cell.energy_layers[surface_layer], 1.8e12); // 90% of original
        assert_eq!(cell_b.next_cell.energy_layers[surface_layer], 1.2e12); // Original + transferred

        // Conservation check: totals should be within 0.1% of initial
        let final_total_volume = cell_a.next_cell.volume_layers[surface_layer] + cell_b.next_cell.volume_layers[surface_layer];
        let final_total_energy = cell_a.next_cell.energy_layers[surface_layer] + cell_b.next_cell.energy_layers[surface_layer];

        let volume_diff_percent =
            ((final_total_volume - initial_total_volume).abs() / initial_total_volume) * 100.0;
        let energy_diff_percent =
            ((final_total_energy - initial_total_energy).abs() / initial_total_energy) * 100.0;

        assert!(
            volume_diff_percent < 0.1,
            "Volume conservation failed: {:.3}% difference",
            volume_diff_percent
        );
        assert!(
            energy_diff_percent < 0.1,
            "Energy conservation failed: {:.3}% difference",
            energy_diff_percent
        );
    }

    #[test]
    fn test_conservation_in_transfer() {
        let cell_index_a = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let cell_index_b = CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap();

        let mut cell_a = AsthenosphereCellNext::new(AsthenosphereCell {
            id: cell_index_a,
            neighbors: H3Utils::neighbors_for(cell_index_a),
            step: 0,
            volume_layers: vec![1000.0; LAYER_COUNT],
            energy_layers: vec![2e12; LAYER_COUNT],
            volcano_volume: 0.0,
            sinkhole_volume: 0.0,
        });
        let mut cell_b = AsthenosphereCellNext::new(AsthenosphereCell {
            id: cell_index_b,
            neighbors: H3Utils::neighbors_for(cell_index_b),
            step: 0,
            volume_layers: vec![500.0; LAYER_COUNT],
            energy_layers: vec![1e12; LAYER_COUNT],
            volcano_volume: 0.0,
            sinkhole_volume: 0.0,
        });

        // Record initial totals (surface layer)
        let surface_layer = LAYER_COUNT - 1;
        let initial_volume = cell_a.next_cell.volume_layers[surface_layer] + cell_b.next_cell.volume_layers[surface_layer];
        let initial_energy = cell_a.next_cell.energy_layers[surface_layer] + cell_b.next_cell.energy_layers[surface_layer];

        // Transfer volume
        cell_a.move_volume_to(&mut cell_b, 200.0);

        // Check conservation
        let final_volume = cell_a.next_cell.volume_layers[surface_layer] + cell_b.next_cell.volume_layers[surface_layer];
        let final_energy = cell_a.next_cell.energy_layers[surface_layer] + cell_b.next_cell.energy_layers[surface_layer];

        assert!(
            (initial_volume - final_volume).abs() < 1e-10,
            "Volume should be conserved"
        );
        assert!(
            (initial_energy - final_energy).abs() < 1e6,
            "Energy should be conserved"
        );
    }
}
