# Plate Generation Visualizer

A React Three Fiber visualizer for exploring tectonic plate generation from the atmospheric simulation.

## 🌍 Overview

This example demonstrates the plate generation capabilities of the atmo-rust simulation through an interactive 3D visualization. View and interact with tectonic plates on a realistic globe representation.

## ✨ Features

- **3D Globe Visualization**: Interactive planet with realistic rendering
- **Tectonic Plate Display**: Color-coded plates based on density and thickness
- **Interactive Controls**: Mouse rotation, zoom, and click-to-select functionality
- **Real-time Statistics**: Planet properties and plate information
- **Responsive Design**: Modern UI that works on desktop and mobile

## 🚀 Quick Start

### Prerequisites
- Node.js (v16 or later)
- Rust (for data generation)

### Installation

1. **Navigate to the example directory**
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

5. **Open your browser**
   Visit [http://localhost:3000](http://localhost:3000)

## 🎮 Usage

### Controls
- **🖱️ Drag**: Rotate the globe
- **🔍 Scroll**: Zoom in/out
- **🎯 Click**: Select plates for detailed information

### Interface
- **Info Panel**: Shows planet statistics and selected plate details
- **Controls Panel**: Quick reference for interaction methods
- **3D Globe**: Auto-rotating planet with color-coded tectonic plates

## 🛠️ Development

### Available Scripts

- `npm run dev` - Start development server with hot reload
- `npm run build` - Build for production
- `npm run preview` - Preview production build
- `npm run generate-data` - Generate new simulation data

### Customizing the Simulation

To modify the simulation parameters, edit the data generation script:

```bash
# Edit the export script in the parent directory
../../examples/export_to_json.rs
```

You can adjust:
- Planet radius and density
- Plate coverage percentage
- Generation algorithms

### Customizing the Visualization

The visualizer is built with modular React components:

- `src/components/GlobeVisualizer.jsx` - Main 3D globe
- `src/components/PlateRenderer.jsx` - Plate visualization logic
- `src/components/InfoPanel.jsx` - Information display
- `src/index.css` - Styling and theme

## 📊 Data Format

The visualizer expects JSON data in this format:

```json
{
  "sim": { "id": "...", "planet_id": "...", "sim_id": "..." },
  "planet": {
    "id": "...",
    "radius_km": 6372,
    "mantle_density_gcm3": 4.5,
    "plate_ids": ["..."]
  },
  "plates": [
    {
      "id": "...",
      "center": {"x": 0.5, "y": 0.3, "z": 0.8},
      "radius_km": 1500,
      "thickness_km": 45,
      "density": 3.2,
      "planet_id": "...",
      "platelet_ids": []
    }
  ]
}
```

## 🎨 Visual Features

- **Color Coding**: Plates are colored based on their physical properties
- **Interactive Selection**: Click any plate to see detailed information
- **Smooth Animations**: Auto-rotation and smooth transitions
- **Space Theme**: Dark background with starfield for immersion
- **Responsive Layout**: Adapts to different screen sizes

## 🔧 Technical Details

- **Built with Vite**: Fast development and optimized builds
- **React Three Fiber**: Declarative 3D graphics with React
- **Three.js**: Powerful 3D rendering engine
- **Modern JavaScript**: ES modules and latest React features

## 🚀 Deployment

To build for production:

```bash
npm run build
```

The built files will be in the `dist/` directory, ready for deployment to any static hosting service.

## 🤝 Contributing

This example is part of the larger atmo-rust project. To contribute:

1. Make changes to the visualizer
2. Test with `npm run dev`
3. Build with `npm run build`
4. Submit a pull request to the main repository

## 📝 License

Same as the parent atmo-rust project.
