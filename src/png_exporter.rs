use image::{ImageBuffer, Rgb, RgbImage};
use std::collections::HashMap;
use std::path::Path;
use crate::asthenosphere::AsthenosphereCell;
use crate::rock_store::RockStore;
use crate::asthenosphere::ASTH_RES;
use crate::planet::Planet;
use crate::constants::{AVG_STARTING_VOLUME_KM_3, CELL_JOULES_START, CELL_JOULES_EQUILIBRIUM};
use h3o::{CellIndex, LatLng};
#[derive(Clone, Copy, Debug)]
struct ColorPoint {
    value: f64,    // normalized energy value between 0.0 and 1.0
    color: Rgb<u8>,
}
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
        
        crate::h3_utils::H3Utils::iter_at(ASTH_RES, |cell_index| {
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

        // Expand volume range to capture more variation: 50-150% of AVG_STARTING_VOLUME
        let volume_min = AVG_STARTING_VOLUME_KM_3 * 0.5;
        let volume_max = AVG_STARTING_VOLUME_KM_3 * 1.5;
        let volume_normalized = ((cell.volume - volume_min) / (volume_max - volume_min)).clamp(0.0, 1.0);

        // Energy range from 0 to 8e24 - much wider spectrum for better color distribution
        let energy_min = 0.0e24;
        let energy_max = 8.0e24;
        let energy_normalized = ((cell.energy_j - energy_min) / (energy_max - energy_min)).clamp(0.0, 1.0);

        let rgb = self.energy_to_color(energy_normalized);
        let (red, green, blue) = (rgb.0[0], rgb.0[1], rgb.0[2]);

        // Make low-volume cells much darker - map volume from 0-1 to 0.1-1.0
        // This creates stronger contrast between high and low volume areas
        let intensity = 0.1 + (volume_normalized * 0.9);

        Rgb([
            (red as f64 * intensity) as u8,
            (green as f64 * intensity) as u8,
            (blue as f64 * intensity) as u8,
        ])
    }
    
    /// Render a Voronoi image directly from cell data without needing to read from store
    pub fn render_voronoi_image_from_cells(&mut self, cells: &[(CellIndex, AsthenosphereCell)]) -> RgbImage {
        // Generate pixel-to-cell mapping if not already done
        if self.pixel_to_cell_map.is_none() {
            println!("Generating Voronoi pixel-to-cell mapping...");
            self.generate_pixel_to_cell_map_from_cells(cells);
        }
        
        let mut image = self.render_voronoi_image(cells);
        self.add_temperature_histogram(&mut image, cells);
        self.add_color_legend(&mut image);
        image
    }
    
    /// Generate pixel-to-cell mapping from cell data instead of reading from store
    fn generate_pixel_to_cell_map_from_cells(&mut self, cells: &[(CellIndex, AsthenosphereCell)]) {
        println!("      🗺️  Generating mapping for {} pixels and {} cells...", self.width * self.height, cells.len());
        let mut pixel_map = HashMap::new();
        
        // Pre-compute cell positions and organize into spatial grid (10 degree regions)
        let cell_positions: Vec<(CellIndex, f64, f64)> = cells.iter()
            .map(|(cell_index, _)| {
                let cell_lat_lon = LatLng::from(*cell_index);
                (*cell_index, cell_lat_lon.lat_radians(), cell_lat_lon.lng_radians())
            })
            .collect();
        
        // Create spatial grid - divide world into 10-degree regions
        let grid_size = (10.0_f64).to_radians(); // 10 degrees in radians
        let mut spatial_grid: HashMap<(i32, i32), Vec<(CellIndex, f64, f64)>> = HashMap::new();
        
        println!("      📦 Organizing cells into spatial grid...");
        for (cell_index, lat, lon) in &cell_positions {
            let grid_lat = (lat / grid_size).floor() as i32;
            let grid_lon = (lon / grid_size).floor() as i32;
            spatial_grid.entry((grid_lat, grid_lon))
                .or_insert_with(Vec::new)
                .push((*cell_index, *lat, *lon));
        }
        
        println!("      📊 Created spatial grid with {} regions", spatial_grid.len());
        
        let total_pixels = self.width * self.height;
        let mut processed = 0;
        
        // Process pixels in 2x2 blocks to reduce computations by 4x
        for y in (0..self.height).step_by(2) {
            for x in (0..self.width).step_by(2) {
                let mut closest_cell = None;
                let mut min_distance = f64::INFINITY;
                
                // Convert pixel to approximate lat/lon (use center of 2x2 block)
                let lon = ((x + 1) as f64 / self.width as f64) * 2.0 * std::f64::consts::PI - std::f64::consts::PI;
                let lat = std::f64::consts::PI / 2.0 - ((y + 1) as f64 / self.height as f64) * std::f64::consts::PI;
                
                // Find which grid region this pixel is in
                let pixel_grid_lat = (lat / grid_size).floor() as i32;
                let pixel_grid_lon = (lon / grid_size).floor() as i32;
                
                // Check current region and neighboring regions (+/- 1 in each direction)
                for grid_lat_offset in -1..=1 {
                    for grid_lon_offset in -1..=1 {
                        let check_grid_lat = pixel_grid_lat + grid_lat_offset;
                        let check_grid_lon = pixel_grid_lon + grid_lon_offset;
                        
                        if let Some(region_cells) = spatial_grid.get(&(check_grid_lat, check_grid_lon)) {
                            // Only check cells in this nearby region
                            for (cell_index, cell_lat, cell_lon) in region_cells {
                                let distance = ((lat - cell_lat).powi(2) + (lon - cell_lon).powi(2)).sqrt();
                                
                                if distance < min_distance {
                                    min_distance = distance;
                                    closest_cell = Some(*cell_index);
                                }
                            }
                        }
                    }
                }
                
                // Apply the same cell to all 4 pixels in the 2x2 block
                if let Some(cell_index) = closest_cell {
                    for block_y in y..=(y + 1).min(self.height - 1) {
                        for block_x in x..=(x + 1).min(self.width - 1) {
                            pixel_map.insert((block_x, block_y), cell_index);
                        }
                    }
                }
                
                processed += 4;
                if processed % 120000 == 0 {
                    println!("      📊 Mapping progress: {}/{} pixels ({:.1}%)", processed, total_pixels, (processed as f64 / total_pixels as f64) * 100.0);
                }
            }
        }
        
        self.pixel_to_cell_map = Some(pixel_map);
        println!("      ✅ Voronoi mapping complete!");
    }
    
    /// Add a temperature distribution histogram to the bottom half of the image
    fn add_temperature_histogram(&self, image: &mut RgbImage, cells: &[(CellIndex, AsthenosphereCell)]) {
        let hist_height = self.height / 2;
        let hist_start_y = self.height - hist_height;
        let padding = 20;
        
        // Make graph only 1/3 as wide as total image, centered
        let chart_width = self.width / 3;
        let chart_x = (self.width - chart_width) / 2; // Center horizontally
        let chart_y = hist_start_y + padding;
        let chart_height = hist_height - (2 * padding);
        
        // No background - keep transparent/existing image background
        
        // Draw black border around histogram area
        for x in chart_x..(chart_x + chart_width) {
            image.put_pixel(x, chart_y, Rgb([0, 0, 0])); // Top border
            image.put_pixel(x, chart_y + chart_height - 1, Rgb([0, 0, 0])); // Bottom border
        }
        for y in chart_y..(chart_y + chart_height) {
            image.put_pixel(chart_x, y, Rgb([0, 0, 0])); // Left border
            image.put_pixel(chart_x + chart_width - 1, y, Rgb([0, 0, 0])); // Right border
        }
        
        // Use wider range for better distribution visualization (0-8e+24)
        let temp_min = 0.0e24;
        let temp_max = 8.0e24;
        let temp_range = temp_max - temp_min;
        
        // Calculate mean energy for display
        let total_energy: f64 = cells.iter().map(|(_, cell)| cell.energy_j).sum();
        let mean_energy = total_energy / cells.len() as f64;
        
        // Count cells in 100 temperature bins for higher resolution
        let mut bins = [0u32; 100];
        for (_, cell) in cells {
            let bin_index = ((cell.energy_j - temp_min) / temp_range * 100.0).floor() as usize;
            let bin_index = bin_index.min(99); // Clamp to valid range
            bins[bin_index] += 1;
        }
        
        // Find max count for scaling
        let max_count = *bins.iter().max().unwrap_or(&1);
        if max_count == 0 { return; }
        
        // Draw thin line bars (2% width)
        let line_width = (chart_width as f32 * 0.02) as u32; // 2% of chart width
        for (i, &count) in bins.iter().enumerate() {
            if count == 0 { continue; }
            
            let bar_height = (count as f32 / max_count as f32 * (chart_height - 20) as f32) as u32;
            let bar_x = chart_x + (i as f32 / 100.0 * chart_width as f32) as u32;
            let bar_start_y = chart_y + chart_height - bar_height - 10;
            
            // Calculate color for this temperature bin using the wide energy range
            let temp_value = temp_min + (i as f64 + 0.5) / 100.0 * temp_range;
            let energy_normalized = ((temp_value - 0.0e24) / (8.0e24 - 0.0e24)).clamp(0.0, 1.0);
            let bar_color = self.energy_to_color(energy_normalized);
            
            // Draw thin vertical line
            for x in bar_x..(bar_x + line_width).min(self.width) {
                for y in bar_start_y..(bar_start_y + bar_height).min(self.height) {
                    image.put_pixel(x, y, bar_color);
                }
            }
        }
        
        // Add text showing mean energy (simple bitmap text)
        let text_y = chart_y + chart_height + 10;
        let mean_text = format!("Mean Energy: {:.2e} J", mean_energy);
        self.draw_simple_text(image, &mean_text, 10, text_y);
    }
    
    /// Draw simple bitmap text (basic implementation)
    fn draw_simple_text(&self, image: &mut RgbImage, text: &str, start_x: u32, start_y: u32) {
        // Simple 3x5 pixel font for basic characters
        let font_map = self.get_simple_font_map();
        
        let mut current_x = start_x;
        for ch in text.chars() {
            if let Some(bitmap) = font_map.get(&ch) {
                for (row, &byte) in bitmap.iter().enumerate() {
                    for col in 0..8 {
                        if byte & (1 << (7 - col)) != 0 {
                            let x = current_x + col;
                            let y = start_y + row as u32;
                            if x < self.width && y < self.height {
                                image.put_pixel(x, y, Rgb([0, 0, 0])); // Black text
                            }
                        }
                    }
                }
                current_x += 6; // Character width + spacing
            } else {
                current_x += 3; // Space for unknown characters
            }
        }
    }
    
    /// Get a simple bitmap font (very basic)
    fn get_simple_font_map(&self) -> std::collections::HashMap<char, [u8; 8]> {
        let mut font = std::collections::HashMap::new();
        
        // Simple 8x8 bitmaps for essential characters
        font.insert('M', [0x00, 0x7C, 0x82, 0x82, 0x92, 0x82, 0x82, 0x00]);
        font.insert('e', [0x00, 0x00, 0x3C, 0x42, 0x7E, 0x40, 0x3C, 0x00]);
        font.insert('a', [0x00, 0x00, 0x3C, 0x02, 0x3E, 0x42, 0x3E, 0x00]);
        font.insert('n', [0x00, 0x00, 0x5C, 0x62, 0x42, 0x42, 0x42, 0x00]);
        font.insert(' ', [0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]);
        font.insert('E', [0x00, 0x7E, 0x40, 0x40, 0x7C, 0x40, 0x7E, 0x00]);
        font.insert('r', [0x00, 0x00, 0x5C, 0x62, 0x40, 0x40, 0x40, 0x00]);
        font.insert('g', [0x00, 0x00, 0x3E, 0x42, 0x42, 0x3E, 0x02, 0x3C]);
        font.insert('y', [0x00, 0x00, 0x42, 0x42, 0x42, 0x3E, 0x02, 0x3C]);
        font.insert(':', [0x00, 0x00, 0x18, 0x18, 0x00, 0x18, 0x18, 0x00]);
        font.insert('J', [0x00, 0x0E, 0x04, 0x04, 0x04, 0x44, 0x38, 0x00]);
        font.insert('o', [0x00, 0x00, 0x3C, 0x42, 0x42, 0x42, 0x3C, 0x00]);
        font.insert('u', [0x00, 0x00, 0x42, 0x42, 0x42, 0x46, 0x3A, 0x00]);
        font.insert('l', [0x00, 0x60, 0x20, 0x20, 0x20, 0x20, 0x70, 0x00]);
        font.insert('s', [0x00, 0x00, 0x3C, 0x40, 0x38, 0x04, 0x78, 0x00]);
        font.insert('.', [0x00, 0x00, 0x00, 0x00, 0x00, 0x18, 0x18, 0x00]);
        font.insert('+', [0x00, 0x00, 0x10, 0x10, 0x7C, 0x10, 0x10, 0x00]);
        font.insert('-', [0x00, 0x00, 0x00, 0x00, 0x7C, 0x00, 0x00, 0x00]);
        
        // Numbers
        font.insert('0', [0x00, 0x3C, 0x46, 0x4A, 0x52, 0x62, 0x3C, 0x00]);
        font.insert('1', [0x00, 0x18, 0x28, 0x08, 0x08, 0x08, 0x3E, 0x00]);
        font.insert('2', [0x00, 0x3C, 0x42, 0x02, 0x3C, 0x40, 0x7E, 0x00]);
        font.insert('3', [0x00, 0x3C, 0x42, 0x0C, 0x02, 0x42, 0x3C, 0x00]);
        font.insert('4', [0x00, 0x08, 0x18, 0x28, 0x48, 0x7E, 0x08, 0x00]);
        font.insert('5', [0x00, 0x7E, 0x40, 0x7C, 0x02, 0x42, 0x3C, 0x00]);
        font.insert('6', [0x00, 0x3C, 0x40, 0x7C, 0x42, 0x42, 0x3C, 0x00]);
        font.insert('7', [0x00, 0x7E, 0x02, 0x04, 0x08, 0x10, 0x10, 0x00]);
        font.insert('8', [0x00, 0x3C, 0x42, 0x3C, 0x42, 0x42, 0x3C, 0x00]);
        font.insert('9', [0x00, 0x3C, 0x42, 0x42, 0x3E, 0x02, 0x3C, 0x00]);
        
        font
    }
    
    /// Add a small color legend on the right side
    fn add_color_legend(&self, image: &mut RgbImage) {
        let legend_width = 30;
        let legend_height = 200;
        let legend_x = self.width - legend_width - 10; // 10px from right edge
        let legend_y = 20; // 20px from top
        
        // Draw legend background (light gray)
        for x in legend_x..(legend_x + legend_width) {
            for y in legend_y..(legend_y + legend_height) {
                image.put_pixel(x, y, Rgb([220, 220, 220]));
            }
        }
        
        // Draw color gradient from top (hot) to bottom (cold)
        for i in 0..legend_height {
            let energy_normalized = 1.0 - (i as f64 / legend_height as f64); // Top = hot (1.0), bottom = cold (0.0)
            let color = self.energy_to_color(energy_normalized);
            
            // Draw color bar (leave 2px margin on each side)
            for x in (legend_x + 2)..(legend_x + legend_width - 2) {
                let y = legend_y + i;
                image.put_pixel(x, y, color);
            }
        }
        
        // Draw black border around legend
        for x in legend_x..(legend_x + legend_width) {
            image.put_pixel(x, legend_y, Rgb([0, 0, 0])); // Top
            image.put_pixel(x, legend_y + legend_height - 1, Rgb([0, 0, 0])); // Bottom
        }
        for y in legend_y..(legend_y + legend_height) {
            image.put_pixel(legend_x, y, Rgb([0, 0, 0])); // Left
            image.put_pixel(legend_x + legend_width - 1, y, Rgb([0, 0, 0])); // Right
        }
        
        // Add temperature labels showing energy ranges for each color region
        let hot_y = legend_y - 8;
        let cold_y = legend_y + legend_height + 2;
        
        // Top label (hottest)
        self.draw_simple_text(image, "8e24", legend_x - 10, hot_y);
        
        // Intermediate labels for color transitions
        let white_yellow_y = legend_y + (legend_height as f64 * 0.2) as u32; // 0.8 normalized = 6.4e24
        let yellow_orange_y = legend_y + (legend_height as f64 * 0.4) as u32; // 0.6 normalized = 4.8e24  
        let orange_red_y = legend_y + (legend_height as f64 * 0.6) as u32;   // 0.4 normalized = 3.2e24
        let red_purple_y = legend_y + (legend_height as f64 * 0.8) as u32;   // 0.2 normalized = 1.6e24
        
        self.draw_simple_text(image, "6.4", legend_x + legend_width + 2, white_yellow_y);
        self.draw_simple_text(image, "4.8", legend_x + legend_width + 2, yellow_orange_y);  
        self.draw_simple_text(image, "3.2", legend_x + legend_width + 2, orange_red_y);
        self.draw_simple_text(image, "1.6", legend_x + legend_width + 2, red_purple_y);
        
        // Bottom label (coldest)
        self.draw_simple_text(image, "0e24", legend_x - 10, cold_y);
    }
    fn lerp_u8(a: u8, b: u8, t: f64) -> u8 {
        (a as f64 + (b as f64 - a as f64) * t).round() as u8
    }
    /// Extract color calculation for reuse in histogram
    /// Linear interpolation between two colors
    fn lerp_color(c1: Rgb<u8>, c2: Rgb<u8>, t: f64) -> Rgb<u8> {
        let r = Self::lerp_u8(c1[0], c2[0], t);
        let g = Self::lerp_u8(c1[1], c2[1], t);
        let b = Self::lerp_u8(c1[2], c2[2], t);
        Rgb([r, g, b])
    }

    /// Map normalized energy value [0.0, 1.0] to color by interpolating between defined color points
    fn energy_to_color(&self, energy_normalized: f64) -> Rgb<u8> {
        // Define gradient points (value must be in ascending order)
        // You can adjust these colors as you like
        let gradient = [
            ColorPoint { value: 0.0, color: Rgb([128, 0, 255]) },   // Purple
            ColorPoint { value: 0.2, color: Rgb([255, 0, 0]) },     // Red
            ColorPoint { value: 0.4, color: Rgb([255, 150, 0]) },   // Orange
            ColorPoint { value: 0.6, color: Rgb([255, 255, 0]) },   // Yellow
            ColorPoint { value: 0.8, color: Rgb([255, 255, 255]) }, // White
            // Optionally add more points for smoother gradient
        ];

        // Clamp input
        let val = energy_normalized.clamp(0.0, 1.0);

        // If val matches exactly a point, return its color
        for point in &gradient {
            if (val - point.value).abs() < 1e-8 {
                return point.color;
            }
        }

        // Find two points between which val lies
        let mut lower = &gradient[0];
        let mut upper = &gradient[gradient.len() - 1];

        for window in gradient.windows(2) {
            if val >= window[0].value && val <= window[1].value {
                lower = &window[0];
                upper = &window[1];
                break;
            }
        }

        // Compute interpolation factor t in [0,1]
        let t = (val - lower.value) / (upper.value - lower.value);

        // Interpolate color
        Self::lerp_color(lower.color, upper.color, t)
    }
}
