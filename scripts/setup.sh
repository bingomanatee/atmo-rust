#\!/bin/bash

echo "🌍 Atmo Rust Visualizer Setup"
echo "=============================="

# Check if Rust is installed
if \! command -v cargo &> /dev/null; then
    echo "❌ Rust/Cargo not found"
    echo "📥 Please install Rust from: https://rustup.rs/"
    echo ""
    echo "After installing Rust, run this script again."
    exit 1
fi

echo "✅ Rust found: $(cargo --version)"

# Check if Node.js is installed
if \! command -v node &> /dev/null; then
    echo "❌ Node.js not found"
    echo "📥 Please install Node.js from: https://nodejs.org/"
    exit 1
fi

echo "✅ Node.js found: $(node --version)"

# Build Rust project
echo ""
echo "🔨 Building Rust simulation..."
cargo build --release

if [ $? -ne 0 ]; then
    echo "❌ Rust build failed"
    exit 1
fi

echo "✅ Rust build successful"

# Install Node.js dependencies
echo ""
echo "📦 Installing Node.js dependencies..."
npm install

if [ $? -ne 0 ]; then
    echo "❌ npm install failed"
    exit 1
fi

echo "✅ Node.js dependencies installed"

# Generate initial simulation data
echo ""
echo "🌍 Generating simulation data..."
cargo run --example export_to_json public/simulation_data.json

if [ $? -ne 0 ]; then
    echo "❌ Failed to generate simulation data"
    echo "ℹ️  Using sample data instead"
else
    echo "✅ Simulation data generated"
fi

echo ""
echo "🎉 Setup complete\!"
echo ""
echo "To start the visualizer:"
echo "  npm start"
echo ""
echo "To regenerate data:"
echo "  npm run generate-data"
echo ""
echo "Open http://localhost:3000 to view the visualizer"
