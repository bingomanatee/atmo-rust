use gif::{Encoder, Frame, Repeat};
use image::{ImageBuffer, Rgb, RgbImage};
use std::fs::File;
use std::borrow::Cow;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create a new GIF file
    let mut file = File::create("animation.gif")?;
    let mut encoder = Encoder::new(&mut file, 200, 200, &[])?;
    encoder.set_repeat(Repeat::Infinite)?;

    // Create 30 frames for a simple animation
    for frame_num in 0..30 {
        let frame_image = create_frame(frame_num, 200, 200);
        
        // Convert RGB image to indexed color (required for GIF)
        let (palette, indexed_data) = rgb_to_indexed(&frame_image);
        
        let mut frame = Frame::from_indexed_pixels(200, 200, &indexed_data, Some(&palette));
        frame.delay = 10; // 10/100 seconds = 0.1 seconds per frame
        
        encoder.write_frame(&frame)?;
    }

    println!("Animated GIF created: animation.gif");
    Ok(())
}

fn create_frame(frame_num: u32, width: u32, height: u32) -> RgbImage {
    let mut img = ImageBuffer::new(width, height);
    
    // Create a simple animation: a moving circle
    let center_x = (frame_num as f32 * 6.0) % (width as f32);
    let center_y = height as f32 / 2.0;
    let radius = 20.0;
    
    for (x, y, pixel) in img.enumerate_pixels_mut() {
        let dx = x as f32 - center_x;
        let dy = y as f32 - center_y;
        let distance = (dx * dx + dy * dy).sqrt();
        
        if distance <= radius {
            // Red circle
            *pixel = Rgb([255, 0, 0]);
        } else {
            // White background
            *pixel = Rgb([255, 255, 255]);
        }
    }
    
    img
}

fn rgb_to_indexed(img: &RgbImage) -> (Vec<u8>, Vec<u8>) {
    // Simple palette: white and red
    let palette = vec![
        255, 255, 255, // White
        255, 0, 0,     // Red
    ];
    
    let mut indexed_data = Vec::new();
    
    for pixel in img.pixels() {
        let [r, g, b] = pixel.0;
        if r > 200 && g < 50 && b < 50 {
            indexed_data.push(1); // Red
        } else {
            indexed_data.push(0); // White
        }
    }
    
    (palette, indexed_data)
}
