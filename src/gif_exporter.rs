use gif::{Encoder, Frame, Repeat};
use image::{ImageBuffer, Rgb, RgbImage};
use std::fs::File;
use std::path::Path;
use crate::asthenosphere::AsthenosphereCell;
use crate::rock_store::RockStore;
use crate::asthenosphere::ASTH_RES;
use crate::planet::Planet;
use crate::constants::{AVG_STARTING_VOLUME_KM_3, CELL_JOULES_START, CELL_JOULES_EQUILIBRIUM};
use h3o::{CellIndex, LatLng};

pub struct GifExporter {
    width: u16,
    height: u16,
    planet: Planet,
}

impl GifExporter {
    pub fn new(width: u16, height: u16, planet: Planet) -> Self {
        Self { width, height, planet }
    }

    /// Export asthenosphere simulation as animated GIF with spots for each H3 cell
    /// Gray scale: 75-125% of AVG_STARTING_VOLUME
    /// Hue: Red at CELL_ENERGY_START+ to Blue at CELL_ENERGY_EQUILIBRIUM
    pub fn export_asthenosphere_animation<P: AsRef<Path>>(
        &self,
        store: &RockStore,
        start_step: u32,
        end_step: u32,
        step_interval: u32,
        output_path: P,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut file = File::create(output_path)?;
        let mut encoder = Encoder::new(&mut file, self.width, self.height, &[])?;
        encoder.set_repeat(Repeat::Infinite)?;

        // Sample every step_interval steps
        let mut step = start_step;
        while step <= end_step {
            println!("Rendering step {} for GIF...", step);
            let cells = self.collect_cells_for_step(store, step)?;
            let frame_image = self.visualize_cells_as_spots(&cells);

            // Convert to indexed color for GIF
            let (_palette, indexed_data) = self.rgb_to_indexed(&frame_image);

            let mut frame = Frame::from_indexed_pixels(
                self.width,
                self.height,
                &indexed_data,
                None
            );
            frame.delay = 20; // 0.2 seconds per frame

            encoder.write_frame(&frame)?;

            step += step_interval;
        }

        Ok(())
    }

    fn collect_cells_for_step(&self, store: &RockStore, step: u32) -> Result<Vec<(CellIndex, AsthenosphereCell)>, Box<dyn std::error::Error>> {
        let mut cells = Vec::new();
        
        crate::h3_utils::H3Utils::iter_at(ASTH_RES, |cell_index| {
            if let Ok(Some(cell)) = store.get_asth(cell_index, step) {
                cells.push((cell_index, cell));
            }
        });
        
        Ok(cells)
    }

    fn visualize_cells_as_spots(&self, cells: &[(CellIndex, AsthenosphereCell)]) -> RgbImage {
        let mut img = ImageBuffer::new(self.width as u32, self.height as u32);
        
        // Fill with black background
        for pixel in img.pixels_mut() {
            *pixel = Rgb([0, 0, 0]);
        }

        // Draw each cell as a colored spot
        for (cell_index, cell) in cells {
            let (x, y) = self.cell_to_pixel(*cell_index);
            let color = self.cell_to_color(cell);

            // Draw a larger circle (8x8 pixels) for each cell to make them visible
            self.draw_spot(&mut img, x, y, 4, color);
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
        let volume_min = AVG_STARTING_VOLUME_KM_3 * 0.75;
        let volume_max = AVG_STARTING_VOLUME_KM_3 * 1.25;
        let volume_normalized = ((cell.volume - volume_min) / (volume_max - volume_min)).clamp(0.0, 1.0);

        // Hue based on energy: Red at CELL_ENERGY_START+ to Blue at CELL_ENERGY_EQUILIBRIUM
        let energy_normalized = ((cell.energy_j - CELL_JOULES_EQUILIBRIUM) / (CELL_JOULES_START - CELL_JOULES_EQUILIBRIUM)).clamp(0.0, 1.0);

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

    fn draw_spot(&self, img: &mut RgbImage, center_x: i32, center_y: i32, radius: i32, color: Rgb<u8>) {
        for dy in -radius..=radius {
            for dx in -radius..=radius {
                let x = center_x + dx;
                let y = center_y + dy;
                
                // Check if within image bounds
                if x >= 0 && x < self.width as i32 && y >= 0 && y < self.height as i32 {
                    // Check if within circle
                    if dx * dx + dy * dy <= radius * radius {
                        img.put_pixel(x as u32, y as u32, color);
                    }
                }
            }
        }
    }

    fn rgb_to_indexed(&self, img: &RgbImage) -> (Vec<u8>, Vec<u8>) {
        // Create a simple 256-color palette
        let mut palette = Vec::new();
        for i in 0..256 {
            let gray = i as u8;
            palette.extend_from_slice(&[gray, gray, gray]);
        }
        
        // Convert image to indexed colors (simple grayscale quantization)
        let mut indexed_data = Vec::new();
        for pixel in img.pixels() {
            let [r, g, b] = pixel.0;
            // Convert to grayscale
            let gray = (0.299 * r as f64 + 0.587 * g as f64 + 0.114 * b as f64) as u8;
            indexed_data.push(gray);
        }
        
        (palette, indexed_data)
    }
}
