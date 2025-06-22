use crate::constants::{CELL_JOULES_START, LAYER_COUNT};
use crate::asthenosphere::AsthenosphereCell;

impl AsthenosphereCell {
    pub fn add_volume(&self, volume: f64) -> AsthenosphereCell {
        if volume == 0.0 {
            return self.clone();
        }
        if volume < 0.0 {
            return self.remove_volume(-volume);
        }

        let mut new_cell = self.clone();
        // Add volume to surface layer (top layer)
        let surface_layer = LAYER_COUNT - 1;
        let new_volume = self.volume_layers[surface_layer] + volume;
        let new_volume_energy = CELL_JOULES_START * volume;
        let old_volume_energy = self.energy_layers[surface_layer] * self.volume_layers[surface_layer];
        let new_energy = (new_volume_energy + old_volume_energy) / new_volume;
        
        new_cell.volume_layers[surface_layer] = new_volume;
        new_cell.energy_layers[surface_layer] = new_energy;
        new_cell
    }

    pub fn remove_volume(&self, volume: f64) -> AsthenosphereCell {
        let mut new_cell = self.clone();
        // Remove volume from surface layer (top layer)
        let surface_layer = LAYER_COUNT - 1;
        let new_volume = self.volume_layers[surface_layer] - volume;
        let energy_k = self.energy_layers[surface_layer] * new_volume / self.volume_layers[surface_layer];
        
        new_cell.volume_layers[surface_layer] = new_volume;
        new_cell.energy_layers[surface_layer] = energy_k;
        new_cell
    }

    pub fn add_by_fraction(&self, fraction: f64) -> AsthenosphereCell {
        if fraction < 0.0 {
            return self.clone();
        }
        let surface_layer = LAYER_COUNT - 1;
        self.add_volume(self.volume_layers[surface_layer] * fraction)
    }

    pub fn remove_by_fraction(&self, fraction: f64) -> AsthenosphereCell {
        if fraction < 0.0 {
            return self.clone();
        }
        let surface_layer = LAYER_COUNT - 1;
        self.remove_volume(self.volume_layers[surface_layer] * fraction)
    }
}
