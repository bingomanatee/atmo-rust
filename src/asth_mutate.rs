use crate::constants::CELL_JOULES_START;
use crate::asthenosphere::AsthenosphereCell;

impl AsthenosphereCell {
    pub fn add_volume(&self, volume: f64) -> AsthenosphereCell {
        if volume == 0.0 {
            return self.clone();
        }
        if volume < 0.0 {
            return self.remove_volume(-volume);
        }

        let new_volume = self.volume + volume;
        let new_volume_energy = CELL_JOULES_START * volume;
        let old_volume_energy = self.energy_j * self.volume;
        let new_energy = (new_volume_energy + old_volume_energy) / new_volume;

        AsthenosphereCell {
            volume: new_volume,
            energy_j: new_energy,
            ..self.clone()
        }
    }

    pub fn remove_volume(&self, volume: f64) -> AsthenosphereCell {
        let new_volume = self.volume - volume;
        let energy_k = self.energy_j * new_volume / self.volume;

        AsthenosphereCell {
            volume: new_volume,
            energy_j: energy_k,
            ..self.clone()
        }
    }

    pub fn add_by_fraction(&self, fraction: f64) -> AsthenosphereCell {
        if fraction < 0.0 {
            return self.clone();
        }
        self.add_volume(self.volume * fraction)
    }

    pub fn remove_by_fraction(&self, fraction: f64) -> AsthenosphereCell {
        if fraction < 0.0 {
            return self.clone();
        }
        self.remove_volume(self.volume * fraction)
    }
}
