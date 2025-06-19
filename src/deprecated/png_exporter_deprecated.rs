use image::{ImageBuffer, Rgb, RgbImage};
use std::collections::HashMap;
use std::path::Path;
use crate::asthenosphere::AsthenosphereCell;
use crate::rock_store::RockStore;
use crate::asthenosphere::ASTH_RES;
use crate::planet::Planet;
use h3o::{CellIndex, LatLng};

/// Deprecated PNG exporter functions that use the old, less efficient methods
pub struct DeprecatedPngExporter {
    width: u32,
    height: u32,
    planet: Planet,
    pixel_to_cell_map: Option<HashMap<(u32, u32), CellIndex>>,
}

impl DeprecatedPngExporter {
    pub fn new(width: u32, height: u32, planet: Planet) -> Self {
        Self {
            width,
            height,
            planet,
            pixel_to_cell_map: None,
        }
    }

    /// @deprecated - Use render_voronoi_image_from_cells instead for better performance
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

    /// @deprecated - Use generate_pixel_to_cell_map_from_cells for better performance with spatial grid
    /// Generate pixel-to-cell mapping by reading from store (inefficient)
    fn generate_pixel_to_cell_map(&mut self, store: &RockStore, step: u32) -> Result<(), Box<dyn std::error::Error>> {
        let cells = self.collect_cells_for_step(store, step)?;
        let mut pixel_map = HashMap::new();
        
        // Convert all cells to pixel coordinates
        let cell_positions: Vec<(CellIndex, (i32, i32))> = cells.iter()
            .map(|(cell_index, _)| (*cell_index, self.cell_to_pixel(*cell_index)))
            .collect();

        // For each pixel, find the nearest cell using Voronoi logic (O(n*m) complexity - very slow!)
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

    /// @deprecated - Use the optimized version with precomputed color table
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
                            let color = self.cell_to_color_deprecated(cell);
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

    /// @deprecated - Use the optimized version with precomputed color lookup table
    fn cell_to_color_deprecated(&self, cell: &AsthenosphereCell) -> Rgb<u8> {
        use crate::constants::AVG_STARTING_VOLUME_KM_3;
        
        // Expand volume range to capture more variation: 50-150% of AVG_STARTING_VOLUME
        let volume_min = AVG_STARTING_VOLUME_KM_3 * 0.5;
        let volume_max = AVG_STARTING_VOLUME_KM_3 * 1.5;
        let volume_normalized = ((cell.volume - volume_min) / (volume_max - volume_min)).clamp(0.0, 1.0);

        // Energy range from 0 to 8e24 - much wider spectrum for better color distribution
        let energy_min = 0.0e24;
        let energy_max = 8.0e24;
        let energy_normalized = ((cell.energy_j - energy_min) / (energy_max - energy_min)).clamp(0.0, 1.0);

        let rgb = self.energy_to_color_deprecated(energy_normalized);
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

    /// @deprecated - Use precomputed color lookup table for better performance
    /// Map normalized energy value [0.0, 1.0] to color by interpolating between defined color points
    fn energy_to_color_deprecated(&self, energy_normalized: f64) -> Rgb<u8> {
        #[derive(Clone, Copy, Debug)]
        struct ColorPoint {
            value: f64,    // normalized energy value between 0.0 and 1.0
            color: Rgb<u8>,
        }

        // Define gradient points (value must be in ascending order)
        let gradient = [
            ColorPoint { value: 0.0, color: Rgb([128, 0, 255]) },   // Purple
            ColorPoint { value: 0.2, color: Rgb([255, 0, 0]) },     // Red
            ColorPoint { value: 0.4, color: Rgb([255, 150, 0]) },   // Orange
            ColorPoint { value: 0.6, color: Rgb([255, 255, 0]) },   // Yellow
            ColorPoint { value: 0.8, color: Rgb([255, 255, 255]) }, // White
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
        Self::lerp_color_deprecated(lower.color, upper.color, t)
    }

    /// Linear interpolation between two colors
    fn lerp_color_deprecated(c1: Rgb<u8>, c2: Rgb<u8>, t: f64) -> Rgb<u8> {
        let r = Self::lerp_u8_deprecated(c1[0], c2[0], t);
        let g = Self::lerp_u8_deprecated(c1[1], c2[1], t);
        let b = Self::lerp_u8_deprecated(c1[2], c2[2], t);
        Rgb([r, g, b])
    }

    fn lerp_u8_deprecated(a: u8, b: u8, t: f64) -> u8 {
        (a as f64 + (b as f64 - a as f64) * t).round() as u8
    }
}
