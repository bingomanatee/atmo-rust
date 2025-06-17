use gif::{Encoder, Frame, Repeat};
use image::{ImageBuffer, Rgb, RgbImage};
use std::fs::File;
use std::collections::HashMap;

// Mock data structures for demonstration
#[derive(Clone)]
struct AsthenosphereCell {
    volume: f64,
    energy_k: f64,
    x: f64,
    y: f64,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create sample asthenosphere data for multiple time steps
    let simulation_data = create_sample_data();
    
    // Create animated GIF
    let mut file = File::create("asthenosphere_animation.gif")?;
    let mut encoder = Encoder::new(&mut file, 400, 400, &[])?;
    encoder.set_repeat(Repeat::Infinite)?;

    for step in 0..simulation_data.len() {
        let frame_image = visualize_asthenosphere_step(&simulation_data[step], 400, 400);
        
        // Convert to indexed color
        let (palette, indexed_data) = create_heat_map_palette(&frame_image);
        
        let mut frame = Frame::from_indexed_pixels(400, 400, &indexed_data, Some(&palette));
        frame.delay = 20; // 20/100 seconds = 0.2 seconds per frame
        
        encoder.write_frame(&frame)?;
    }

    println!("Asthenosphere animation created: asthenosphere_animation.gif");
    Ok(())
}

fn create_sample_data() -> Vec<Vec<AsthenosphereCell>> {
    let mut steps = Vec::new();
    
    // Create 20 time steps
    for step in 0..20 {
        let mut cells = Vec::new();
        
        // Create a grid of cells
        for i in 0..20 {
            for j in 0..20 {
                let x = i as f64 / 20.0;
                let y = j as f64 / 20.0;
                
                // Simulate some wave-like energy pattern
                let time_factor = step as f64 * 0.3;
                let energy = 1500.0 + 500.0 * ((x * 10.0 + time_factor).sin() * (y * 10.0 + time_factor).cos());
                let volume = 4000.0 + 500.0 * ((x * 8.0 - time_factor).cos() * (y * 8.0 - time_factor).sin());
                
                cells.push(AsthenosphereCell {
                    volume,
                    energy_k: energy,
                    x,
                    y,
                });
            }
        }
        steps.push(cells);
    }
    
    steps
}

fn visualize_asthenosphere_step(cells: &[AsthenosphereCell], width: u32, height: u32) -> RgbImage {
    let mut img = ImageBuffer::new(width, height);
    
    // Find min/max energy for normalization
    let min_energy = cells.iter().map(|c| c.energy_k).fold(f64::INFINITY, f64::min);
    let max_energy = cells.iter().map(|c| c.energy_k).fold(f64::NEG_INFINITY, f64::max);
    
    for (x, y, pixel) in img.enumerate_pixels_mut() {
        let norm_x = x as f64 / width as f64;
        let norm_y = y as f64 / height as f64;
        
        // Find closest cell (simple nearest neighbor)
        let closest_cell = cells.iter()
            .min_by(|a, b| {
                let dist_a = ((a.x - norm_x).powi(2) + (a.y - norm_y).powi(2)).sqrt();
                let dist_b = ((b.x - norm_x).powi(2) + (b.y - norm_y).powi(2)).sqrt();
                dist_a.partial_cmp(&dist_b).unwrap()
            })
            .unwrap();
        
        // Normalize energy to 0-1 range
        let normalized_energy = (closest_cell.energy_k - min_energy) / (max_energy - min_energy);
        
        // Create heat map color (blue = cold, red = hot)
        let color = energy_to_color(normalized_energy);
        *pixel = color;
    }
    
    img
}

fn energy_to_color(normalized_energy: f64) -> Rgb<u8> {
    // Clamp to 0-1 range
    let energy = normalized_energy.max(0.0).min(1.0);
    
    if energy < 0.5 {
        // Blue to green
        let t = energy * 2.0;
        Rgb([0, (255.0 * t) as u8, (255.0 * (1.0 - t)) as u8])
    } else {
        // Green to red
        let t = (energy - 0.5) * 2.0;
        Rgb([(255.0 * t) as u8, (255.0 * (1.0 - t)) as u8, 0])
    }
}

fn create_heat_map_palette(img: &RgbImage) -> (Vec<u8>, Vec<u8>) {
    // Create a palette with 256 colors for heat map
    let mut palette = Vec::new();
    for i in 0..256 {
        let normalized = i as f64 / 255.0;
        let color = energy_to_color(normalized);
        palette.extend_from_slice(&[color.0[0], color.0[1], color.0[2]]);
    }
    
    // Convert image to indexed colors
    let mut indexed_data = Vec::new();
    for pixel in img.pixels() {
        let [r, g, b] = pixel.0;
        
        // Find closest palette color (simplified)
        let mut best_index = 0;
        let mut best_distance = f64::INFINITY;
        
        for (i, chunk) in palette.chunks(3).enumerate() {
            let pr = chunk[0] as f64;
            let pg = chunk[1] as f64;
            let pb = chunk[2] as f64;
            
            let distance = ((r as f64 - pr).powi(2) + (g as f64 - pg).powi(2) + (b as f64 - pb).powi(2)).sqrt();
            if distance < best_distance {
                best_distance = distance;
                best_index = i;
            }
        }
        
        indexed_data.push(best_index as u8);
    }
    
    (palette, indexed_data)
}
