use image::{ImageBuffer, Rgb, RgbImage};
use std::collections::HashMap;
use crate::asthenosphere::AsthenosphereCell;
use crate::planet::Planet;
use crate::constants::AVG_STARTING_VOLUME_KM_3;
use h3o::{CellIndex, LatLng};
use once_cell::sync::Lazy;
// Energy range constants
const ENERGY_MIN: f64 = 0.0e24;
const ENERGY_MAX: f64 = 8.0e24;
const COLOR_TABLE_SIZE: usize = 1000;

#[derive(Clone, Copy, Debug)]
struct ColorPoint {
    value: f64,    // normalized energy value between 0.0 and 1.0
    color: Rgb<u8>,
}

// Precomputed color lookup table - computed once at startup
static COLOR_LOOKUP_TABLE: Lazy<Vec<Rgb<u8>>> = Lazy::new(|| {
    let mut color_table = Vec::with_capacity(COLOR_TABLE_SIZE);

    for i in 0..COLOR_TABLE_SIZE {
        let energy_normalized = i as f64 / (COLOR_TABLE_SIZE - 1) as f64;
        let color = compute_energy_to_color_static(energy_normalized);
        color_table.push(color);
    }

    color_table
});

// Static version of energy_to_color for use in const context
fn compute_energy_to_color_static(energy_normalized: f64) -> Rgb<u8> {
    // Define gradient points for magma field appearance (value must be in ascending order)
    let gradient = [
        ColorPoint { value: 0.0, color: Rgb([0, 0, 0]) },       // Black (coldest)
        ColorPoint { value: 0.15, color: Rgb([51, 0, 128]) },     // Black at 3.0e24 J
        ColorPoint { value: 0.375, color: Rgb([128, 0, 25]) },     // Black at 3.0e24 J
        ColorPoint { value: 0.6, color: Rgb([255, 0, 0]) },     // Red where yellow was (4.8e24 J)
        ColorPoint { value: 0.8, color: Rgb([255, 255, 0]) },   // Yellow (hotter)
        ColorPoint { value: 1.0, color: Rgb([255, 255, 255]) }, // White (hottest)
    ];

    // Clamp input
    let val = if energy_normalized < 0.0 { 0.0 } else if energy_normalized > 1.0 { 1.0 } else { energy_normalized };

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
    lerp_color_static(lower.color, upper.color, t)
}

// Static version of lerp_color
fn lerp_color_static(c1: Rgb<u8>, c2: Rgb<u8>, t: f64) -> Rgb<u8> {
    let r = lerp_u8_static(c1[0], c2[0], t);
    let g = lerp_u8_static(c1[1], c2[1], t);
    let b = lerp_u8_static(c1[2], c2[2], t);
    Rgb([r, g, b])
}

// Static version of lerp_u8
fn lerp_u8_static(a: u8, b: u8, t: f64) -> u8 {
    (a as f64 + (b as f64 - a as f64) * t).round() as u8
}
pub struct PngExporter {
    width: u32,
    height: u32,
    planet: Planet,
    pixel_to_cell_map: Option<HashMap<(u32, u32), CellIndex>>,
    spatial_grid: Option<HashMap<(i32, i32), Vec<(CellIndex, f64, f64)>>>,
}

impl PngExporter {
    pub fn new(width: u32, height: u32, planet: Planet) -> Self {
        Self {
            width,
            height,
            planet,
            pixel_to_cell_map: None,
            spatial_grid: None,
        }
    }





    fn cell_to_color(&self, cell: &AsthenosphereCell) -> Rgb<u8> {
        // Expand volume range to capture more variation: 50-150% of AVG_STARTING_VOLUME
        let volume_min = AVG_STARTING_VOLUME_KM_3 * 0.5;
        let volume_max = AVG_STARTING_VOLUME_KM_3 * 1.5;
        let volume_normalized = ((cell.volume - volume_min) / (volume_max - volume_min)).clamp(0.0, 1.0);

        // Use precomputed global color lookup table for fast color retrieval
        let energy_normalized = ((cell.energy_j - ENERGY_MIN) / (ENERGY_MAX - ENERGY_MIN)).clamp(0.0, 1.0);

        // Fast lookup: convert normalized energy to table index
        let table_index = (energy_normalized * (COLOR_TABLE_SIZE - 1) as f64).round() as usize;
        let table_index = table_index.min(COLOR_TABLE_SIZE - 1);
        let rgb = COLOR_LOOKUP_TABLE[table_index];

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

        let mut image = self.render_voronoi_image_optimized(cells);
        self.add_temperature_histogram(&mut image, cells);
        self.add_color_legend(&mut image);
        image
    }

    /// Optimized Voronoi image rendering using precomputed pixel-to-cell mapping
    fn render_voronoi_image_optimized(&self, cells: &[(CellIndex, AsthenosphereCell)]) -> RgbImage {
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

    /// Generate pixel-to-cell mapping from cell data instead of reading from store
    fn generate_pixel_to_cell_map_from_cells(&mut self, cells: &[(CellIndex, AsthenosphereCell)]) {
        println!("      🗺️  Generating mapping for {} pixels and {} cells...", self.width * self.height, cells.len());
        let mut pixel_map = HashMap::with_capacity((self.width * self.height / 4) as usize); // Pre-allocate for 2x2 blocks

        // Use cached spatial grid or create it if not available
        let spatial_grid = if let Some(ref grid) = self.spatial_grid {
            println!("      ♻️  Using cached spatial grid with {} regions", grid.len());
            grid
        } else {
            println!("      📦 Creating and caching spatial grid...");

            // Pre-compute cell positions and organize into spatial grid (5 degree regions for finer granularity)
            let cell_positions: Vec<(CellIndex, f64, f64)> = cells.iter()
                .map(|(cell_index, _)| {
                    let cell_lat_lon = LatLng::from(*cell_index);
                    (*cell_index, cell_lat_lon.lat_radians(), cell_lat_lon.lng_radians())
                })
                .collect();

            // Create spatial grid - divide world into 5-degree regions for better locality
            let grid_size = (5.0_f64).to_radians(); // 5 degrees in radians (smaller = more efficient)
            let mut new_spatial_grid: HashMap<(i32, i32), Vec<(CellIndex, f64, f64)>> = HashMap::new();

            for (cell_index, lat, lon) in &cell_positions {
                let grid_lat = (lat / grid_size).floor() as i32;
                let grid_lon = (lon / grid_size).floor() as i32;
                new_spatial_grid.entry((grid_lat, grid_lon))
                    .or_insert_with(Vec::new)
                    .push((*cell_index, *lat, *lon));
            }

            println!("      📊 Created spatial grid with {} regions (avg {:.1} cells/region)",
                    new_spatial_grid.len(),
                    cells.len() as f64 / new_spatial_grid.len() as f64);

            // Cache the spatial grid for future use
            self.spatial_grid = Some(new_spatial_grid);
            self.spatial_grid.as_ref().unwrap()
        };

        let grid_size = (5.0_f64).to_radians(); // Must match the grid size used above
        
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

                // Check current region first (most likely to contain closest cell)
                if let Some(region_cells) = spatial_grid.get(&(pixel_grid_lat, pixel_grid_lon)) {
                    for (cell_index, cell_lat, cell_lon) in region_cells {
                        // Use squared distance to avoid expensive sqrt() calls
                        let distance_squared = (lat - cell_lat).powi(2) + (lon - cell_lon).powi(2);

                        if distance_squared < min_distance * min_distance {
                            min_distance = distance_squared.sqrt();
                            closest_cell = Some(*cell_index);
                        }
                    }
                }

                // Only check neighboring regions if we haven't found a very close cell
                let search_threshold = (2.0_f64).to_radians(); // ~2 degrees
                if min_distance > search_threshold {
                    // Check neighboring regions (+/- 1 in each direction)
                    for grid_lat_offset in -1..=1 {
                        for grid_lon_offset in -1..=1 {
                            if grid_lat_offset == 0 && grid_lon_offset == 0 { continue; } // Skip center (already checked)

                            let check_grid_lat = pixel_grid_lat + grid_lat_offset;
                            let check_grid_lon = pixel_grid_lon + grid_lon_offset;

                            if let Some(region_cells) = spatial_grid.get(&(check_grid_lat, check_grid_lon)) {
                                for (cell_index, cell_lat, cell_lon) in region_cells {
                                    let distance_squared = (lat - cell_lat).powi(2) + (lon - cell_lon).powi(2);

                                    if distance_squared < min_distance * min_distance {
                                        min_distance = distance_squared.sqrt();
                                        closest_cell = Some(*cell_index);
                                    }
                                }
                            }
                        }
                    }
                }
                
                // Fallback: if no cell found in spatial grid, search all regions
                if closest_cell.is_none() {
                    for region_cells in spatial_grid.values() {
                        for (cell_index, cell_lat, cell_lon) in region_cells {
                            let distance_squared = (lat - cell_lat).powi(2) + (lon - cell_lon).powi(2);
                            if distance_squared < min_distance * min_distance {
                                min_distance = distance_squared.sqrt();
                                closest_cell = Some(*cell_index);
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
        // No border - clean histogram appearance
        
        // Use wider range for better distribution visualization
        let temp_min = ENERGY_MIN;
        let temp_max = ENERGY_MAX;
        let temp_range = temp_max - temp_min;
        
        // Calculate mean energy for display
        let total_energy: f64 = cells.iter().map(|(_, cell)| cell.energy_j).sum();
        let mean_energy = total_energy / cells.len() as f64;

        // Calculate standard deviations for volume and energy
        let total_volume: f64 = cells.iter().map(|(_, cell)| cell.volume).sum();
        let mean_volume = total_volume / cells.len() as f64;

        // Calculate standard deviation for volume
        let volume_variance: f64 = cells.iter()
            .map(|(_, cell)| {
                let diff = cell.volume - mean_volume;
                diff * diff
            })
            .sum::<f64>() / cells.len() as f64;
        let std_volume = volume_variance.sqrt();

        // Calculate standard deviation for energy
        let energy_variance: f64 = cells.iter()
            .map(|(_, cell)| {
                let diff = cell.energy_j - mean_energy;
                diff * diff
            })
            .sum::<f64>() / cells.len() as f64;
        let std_energy = energy_variance.sqrt();
        
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
            let energy_normalized = ((temp_value - ENERGY_MIN) / (ENERGY_MAX - ENERGY_MIN)).clamp(0.0, 1.0);

            // Fast lookup using precomputed global table
            let table_index = (energy_normalized * (COLOR_TABLE_SIZE - 1) as f64).round() as usize;
            let table_index = table_index.min(COLOR_TABLE_SIZE - 1);
            let bar_color = COLOR_LOOKUP_TABLE[table_index];
            
            // Draw thin vertical line
            for x in bar_x..(bar_x + line_width).min(self.width) {
                for y in bar_start_y..(bar_start_y + bar_height).min(self.height) {
                    image.put_pixel(x, y, bar_color);
                }
            }
        }
        
        // Add text showing mean energy (simple bitmap text) - 2x bigger, moved up to avoid clipping
        let text_y = chart_y + chart_height - 25; // Move up into the chart area
        let mean_text = format!("Mean Energy: {:.2e} J", mean_energy);
        self.draw_simple_text_2x(image, &mean_text, 10, text_y);

        // Add standard deviation statistics in the right corner - 2x bigger
        let std_text_x = self.width - 300; // Position from right edge
        let std_volume_y = 20; // Top right corner
        let std_energy_y = std_volume_y + 20; // Below std volume

        let std_volume_text = format!("Std Volume: {:.2e} km3", std_volume);
        let std_energy_text = format!("Std Energy: {:.2e} J", std_energy);

        self.draw_simple_text_2x(image, &std_volume_text, std_text_x, std_volume_y);
        self.draw_simple_text_2x(image, &std_energy_text, std_text_x, std_energy_y);
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

    /// Draw simple bitmap text at 2x size (double resolution)
    fn draw_simple_text_2x(&self, image: &mut RgbImage, text: &str, start_x: u32, start_y: u32) {
        // Simple 3x5 pixel font for basic characters, rendered at 2x scale
        let font_map = self.get_simple_font_map();

        let mut current_x = start_x;
        for ch in text.chars() {
            if let Some(bitmap) = font_map.get(&ch) {
                for (row, &byte) in bitmap.iter().enumerate() {
                    for col in 0..8 {
                        if byte & (1 << (7 - col)) != 0 {
                            // Draw 2x2 pixel block for each original pixel
                            for dy in 0..2 {
                                for dx in 0..2 {
                                    let x = current_x + (col * 2) + dx;
                                    let y = start_y + (row as u32 * 2) + dy;
                                    if x < self.width && y < self.height {
                                        image.put_pixel(x, y, Rgb([0, 0, 0])); // Black text
                                    }
                                }
                            }
                        }
                    }
                }
                current_x += 12; // Character width + spacing (2x original)
            } else {
                current_x += 6; // Space for unknown characters (2x original)
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

            // Fast lookup using precomputed global table
            let table_index = (energy_normalized * (COLOR_TABLE_SIZE - 1) as f64).round() as usize;
            let table_index = table_index.min(COLOR_TABLE_SIZE - 1);
            let color = COLOR_LOOKUP_TABLE[table_index];

            // Draw color bar (leave 2px margin on each side)
            for x in (legend_x + 2)..(legend_x + legend_width - 2) {
                let y = legend_y + i;
                image.put_pixel(x, y, color);
            }
        }

        // No border - clean legend appearance
        
        // Add temperature labels showing energy ranges for each color region - 2x bigger and on left side
        let hot_y = legend_y - 16; // More space for 2x text
        let cold_y = legend_y + legend_height + 4; // More space for 2x text

        // Top label (hottest) - moved to left side and 2x bigger
        self.draw_simple_text_2x(image, "8e24", legend_x - 50, hot_y);

        // Intermediate labels for color transitions - moved to left side and 2x bigger
        let white_y = legend_y + (legend_height as f64 * 0.0) as u32;     // 1.0 normalized = 8.0e24 (white)
        let yellow_y = legend_y + (legend_height as f64 * 0.2) as u32;    // 0.8 normalized = 6.4e24 (yellow)
        let red_y = legend_y + (legend_height as f64 * 0.4) as u32;       // 0.6 normalized = 4.8e24 (red)
        let black_transition_y = legend_y + (legend_height as f64 * 0.625) as u32; // 0.375 normalized = 3.0e24 (black)

        self.draw_simple_text_2x(image, "6.4", legend_x - 40, yellow_y);
        self.draw_simple_text_2x(image, "4.8", legend_x - 40, red_y);
        self.draw_simple_text_2x(image, "3.0", legend_x - 40, black_transition_y);

        // Bottom label (coldest) - moved to left side and 2x bigger
        self.draw_simple_text_2x(image, "0e24", legend_x - 50, cold_y);
    }
}
