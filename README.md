# Atmo Rust - Atmospheric Simulation Visualizer

A Rust-based atmospheric simulation with a React Three Fiber web visualizer for tectonic plates on a 3D globe.

## Project Structure

```
atmo-rust/
├── src/                           # Rust simulation code
│   ├── sim_manager.rs            # Main simulation manager with JSON export
│   ├── planet.rs                 # Planet data structures
│   ├── plate.rs                  # Tectonic plate structures
│   └── ...
├── examples/
│   ├── export_to_json.rs         # Script to generate JSON data
│   └── plate_generation/         # React Three Fiber visualizer
│       ├── src/
│       │   ├── components/
│       │   │   ├── GlobeVisualizer.jsx
│       │   │   ├── PlateRenderer.jsx
│       │   │   └── InfoPanel.jsx
│       │   ├── App.jsx
│       │   └── main.jsx
│       ├── public/
│       │   └── simulation_data.json
│       ├── package.json
│       └── vite.config.js
└── Cargo.toml
```

## Features

### Rust Simulation
- **Planet Generation**: Create planets with configurable radius and density
- **Tectonic Plate Simulation**: Generate realistic tectonic plates with varying sizes, densities, and thicknesses
- **JSON Export**: Export simulation data for visualization
- **Database Storage**: Persistent storage using RocksDB

### React Three Fiber Visualizer (Vite)
- **3D Globe**: Interactive sphere representing the planet
- **Tectonic Plates**: Color-coded plates rendered on the globe surface
- **Interactive Controls**:
  - Mouse drag to rotate
  - Scroll to zoom
  - Click plates for detailed information
- **Information Panel**: Real-time display of simulation statistics
- **Responsive Design**: Works on desktop and mobile
- **Fast Development**: Powered by Vite for instant hot reload

## Getting Started

### Prerequisites
- Rust (latest stable)
- Node.js (v16 or later)
- npm or yarn

### Installation

1. **Clone the repository**
   ```bash
   git clone <repository-url>
   cd atmo-rust
   ```

2. **Build the Rust simulation**
   ```bash
   cargo build
   ```

### Running the Visualizer

1. **Navigate to the visualizer**
   ```bash
   cd examples/plate_generation
   ```

2. **Install dependencies**
   ```bash
   npm install
   ```

3. **Generate simulation data**
   ```bash
   npm run generate-data
   ```

4. **Start the development server**
   ```bash
   npm run dev
   ```
   Open [http://localhost:3000](http://localhost:3000) to view the visualizer

### Manual Data Generation

You can also generate data manually with custom parameters:

```bash
# Generate with default parameters
cargo run --example export_to_json

# Generate to specific file
cargo run --example export_to_json -- custom_output.json
```

## Usage

### Web Interface
- **Rotate**: Click and drag to rotate the globe
- **Zoom**: Use mouse wheel to zoom in/out
- **Select Plates**: Click on colored regions to see plate details
- **View Statistics**: Check the info panel for simulation data

### Customizing Simulations

Edit `examples/export_to_json.rs` to modify simulation parameters:

```rust
let planet_config = SimPlanetParams {
    radius: 6372,           // Planet radius in km
    mantle_density_gcm3: Some(4.5), // Density in g/cm³
};

// Adjust plate coverage (0.0 to 1.0)
manager.make_plates(0.6);
```

## Data Format

The exported JSON contains:
```json
{
  "sim": {
    "id": "uuid",
    "planet_id": "uuid",
    "sim_id": "uuid"
  },
  "planet": {
    "id": "uuid",
    "sim_id": "uuid", 
    "radius_km": 6372,
    "mantle_density_gcm3": 4.5,
    "plate_ids": ["uuid1", "uuid2", ...]
  },
  "plates": [
    {
      "id": "uuid",
      "center": {"x": 0.5, "y": 0.3, "z": 0.8},
      "radius_km": 1500,
      "thickness_km": 45,
      "density": 3.2,
      "planet_id": "uuid",
      "platelet_ids": []
    }
  ]
}
```

## Development

### Adding New Features

1. **Rust Side**: Modify simulation logic in `src/`
2. **Web Side**: Add React components in `examples/plate_generation/src/components/`
3. **Data Export**: Update `SimExportData` structure if needed

### Building the Visualizer for Production

```bash
cd examples/plate_generation
npm run build
```

This creates an optimized build in the `dist/` folder.

## License

[Add your license here]
