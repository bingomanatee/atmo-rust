use image::{ImageBuffer, Rgb, RgbImage};
use std::collections::HashMap;
use std::path::Path;
use crate::asthenosphere::AsthenosphereCell;
use crate::rock_store::RockStore;
use crate::asthenosphere::ASTH_RES;
use crate::planet::Planet;
use crate::asth_constants::{AVG_STARTING_VOLUME, CELL_ENERGY_START, CELL_ENERGY_EQUILIBRIUM};
use h3o::{CellIndex, LatLng};

pub struct PngExporter {
    width: u32,
    height: u32,
    planet: Planet,
    pixel_to_cell_map: Option<HashMap<(u32, u32), CellIndex>>,
}

impl PngExporter {
    pub fn new(width: u32, height: u32, planet: Planet) -> Self {
        Self { 
            width, 
            height, 
            planet,
            pixel_to_cell_map: None,
        }
    }

    /// Export asthenosphere simulation as PNG images every 10 steps using Voronoi pattern
    pub fn export_asthenosphere_pngs<P: AsRef<Path>>(
        &mut self,
        store: &RockStore,
        start_step: u32,
        end_step: u32,
        step_interval: u32,
        output_dir: P,
    ) -> Result<(), Box<dyn std::error::Error>> {
        // Create output directory if it doesn't exist
        std::fs::create_dir_all(&output_dir)?;
        
        // Generate pixel-to-cell mapping on first run
        if self.pixel_to_cell_map.is_none() {
            println!("Generating Voronoi pixel-to-cell mapping...");
            self.generate_pixel_to_cell_map(store, start_step)?;
        }

        let mut step = start_step;
        while step <= end_step {
            println!("Rendering step {} as PNG...", step);
            let cells = self.collect_cells_for_step(store, step)?;
            let image = self.render_voronoi_image(&cells);
            
            let filename = format!("asthenosphere_step_{:04}.png", step);
            let filepath = output_dir.as_ref().join(filename);
            image.save(filepath)?;
            
            step += step_interval;
        }

        Ok(())
    }

    fn generate_pixel_to_cell_map(&mut self, store: &RockStore, step: u32) -> Result<(), Box<dyn std::error::Error>> {
        let cells = self.collect_cells_for_step(store, step)?;
        let mut pixel_map = HashMap::new();
        
        // Convert all cells to pixel coordinates
        let cell_positions: Vec<(CellIndex, (i32, i32))> = cells.iter()
            .map(|(cell_index, _)| (*cell_index, self.cell_to_pixel(*cell_index)))
            .collect();

        // For each pixel, find the nearest cell using Voronoi logic
        for y in 0..self.height {
            for x in 0..self.width {
                let mut min_distance = f64::MAX;
                let mut nearest_cell = cell_positions[0].0; // Default to first cell
                
                for (cell_index, (cell_x, cell_y)) in &cell_positions {
                    let dx = x as i32 - cell_x;
                    let dy = y as i32 - cell_y;
                    let distance = ((dx * dx + dy * dy) as f64).sqrt();
                    
                    if distance < min_distance {
                        min_distance = distance;
                        nearest_cell = *cell_index;
                    }
                }
                
                pixel_map.insert((x, y), nearest_cell);
            }
            
            // Progress indicator
            if y % 50 == 0 {
                println!("Mapping progress: {}/{}", y, self.height);
            }
        }
        
        self.pixel_to_cell_map = Some(pixel_map);
        println!("Voronoi mapping complete!");
        Ok(())
    }

    fn collect_cells_for_step(&self, store: &RockStore, step: u32) -> Result<Vec<(CellIndex, AsthenosphereCell)>, Box<dyn std::error::Error>> {
        let mut cells = Vec::new();
        
        crate::h30_utils::H3Utils::iter_at(ASTH_RES, |cell_index| {
            if let Ok(Some(cell)) = store.get_asth(cell_index, step) {
                cells.push((cell_index, cell));
            }
        });
        
        Ok(cells)
    }

    fn render_voronoi_image(&self, cells: &[(CellIndex, AsthenosphereCell)]) -> RgbImage {
        let mut img = ImageBuffer::new(self.width, self.height);

        // Create a lookup map for quick cell data access
        let cell_data: HashMap<CellIndex, &AsthenosphereCell> = cells.iter()
            .map(|(index, cell)| (*index, cell))
            .collect();

        // Use the pre-computed pixel-to-cell mapping
        if let Some(pixel_map) = &self.pixel_to_cell_map {
            for y in 0..self.height {
                for x in 0..self.width {
                    if let Some(cell_index) = pixel_map.get(&(x, y)) {
                        if let Some(cell) = cell_data.get(cell_index) {
                            let color = self.cell_to_color(cell);
                            img.put_pixel(x, y, color);
                        } else {
                            // Cell not found in current step data, use dark red to indicate missing
                            img.put_pixel(x, y, Rgb([50, 0, 0]));
                        }
                    } else {
                        // No mapping found, use black
                        img.put_pixel(x, y, Rgb([0, 0, 0]));
                    }
                }
            }
        } else {
            // Fallback: fill with black if no mapping exists
            for pixel in img.pixels_mut() {
                *pixel = Rgb([0, 0, 0]);
            }
        }

        img
    }

    fn cell_to_pixel(&self, cell_index: CellIndex) -> (i32, i32) {
        // Convert H3 cell to lat/lon
        let lat_lon = LatLng::from(cell_index);
        let lat = lat_lon.lat_radians();
        let lon = lat_lon.lng_radians();
        
        // Simple equirectangular projection
        let x = ((lon + std::f64::consts::PI) / (2.0 * std::f64::consts::PI) * self.width as f64) as i32;
        let y = ((std::f64::consts::PI / 2.0 - lat) / std::f64::consts::PI * self.height as f64) as i32;
        
        (x, y)
    }

    fn cell_to_color(&self, cell: &AsthenosphereCell) -> Rgb<u8> {

        // Gray scale based on volume: 75-125% of AVG_STARTING_VOLUME
        let volume_min = AVG_STARTING_VOLUME * 0.75;
        let volume_max = AVG_STARTING_VOLUME * 1.25;
        let volume_normalized = ((cell.volume - volume_min) / (volume_max - volume_min)).clamp(0.0, 1.0);

        // Hue based on energy: 3000K (hot/red) to 1000K (cold/blue)
        let temp_hot = 3000.0;
        let temp_cold = 1000.0;
        let energy_normalized = ((cell.energy_k - temp_cold) / (temp_hot - temp_cold)).clamp(0.0, 1.0);

        // Interpolate between blue (cold) and red (hot)
        let red = (energy_normalized * 255.0) as u8;
        let blue = ((1.0 - energy_normalized) * 255.0) as u8;
        let green = 0u8; // Keep green at 0 for clear red-blue gradient

        // Use volume to modulate brightness, but ensure minimum visibility
        // Map volume_normalized from 0-1 to 0.3-1.0 so cells are always visible
        let intensity = 0.3 + (volume_normalized * 0.7);

        Rgb([
            (red as f64 * intensity) as u8,
            (green as f64 * intensity) as u8,
            (blue as f64 * intensity) as u8,
        ])
    }
}
