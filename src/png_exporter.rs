use image::{ImageBuffer, Rgb, RgbImage};
use std::collections::HashMap;
use std::path::Path;
use crate::asthenosphere::AsthenosphereCell;
use crate::rock_store::RockStore;
use crate::asthenosphere::ASTH_RES;
use crate::planet::Planet;
use crate::constants::{AVG_STARTING_VOLUME_KM_3, CELL_JOULES_START, CELL_JOULES_EQUILIBRIUM};
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

        // Energy range from 5.0e24 (cold/blue) to 5.8e24 (hot/red)
        let energy_min = 5.0e24;
        let energy_max = 5.8e24;
        let energy_normalized = ((cell.energy_j - energy_min) / (energy_max - energy_min)).clamp(0.0, 1.0);

        // Color gradient: white -> yellow -> red -> purple -> black
        let (red, green, blue) = if energy_normalized >= 0.75 {
            // White to yellow (hottest: 0.75-1.0)
            let t = (energy_normalized - 0.75) / 0.25;
            let intensity = 255.0 * (0.5 + t * 0.5); // 50% to 100% intensity
            (intensity as u8, intensity as u8, (255.0 * (1.0 - t)) as u8)
        } else if energy_normalized >= 0.5 {
            // Yellow to red (hot: 0.5-0.75) 
            let t = (energy_normalized - 0.5) / 0.25;
            (255, (255.0 * (1.0 - t)) as u8, 0)
        } else if energy_normalized >= 0.25 {
            // Red to purple (cool: 0.25-0.5)
            let t = (energy_normalized - 0.25) / 0.25;
            (255, 0, (255.0 * t) as u8)
        } else {
            // Purple to black (coldest: 0.0-0.25)
            let t = energy_normalized / 0.25;
            ((128.0 * t) as u8, 0, (255.0 * t) as u8)
        };

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
        
        self.render_voronoi_image(cells)
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
}
