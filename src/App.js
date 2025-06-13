import React, { useState, useEffect } from 'react';
import { Canvas } from '@react-three/fiber';
import { OrbitControls, Stars } from '@react-three/drei';
import GlobeVisualizer from './components/GlobeVisualizer';
import InfoPanel from './components/InfoPanel';
import './App.css';

function App() {
  const [simulationData, setSimulationData] = useState(null);
  const [selectedPlate, setSelectedPlate] = useState(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);

  useEffect(() => {
    // Load simulation data
    fetch('/simulation_data.json')
      .then(response => {
        if (!response.ok) {
          throw new Error(`HTTP error! status: ${response.status}`);
        }
        return response.json();
      })
      .then(data => {
        console.log('Loaded simulation data:', data);
        setSimulationData(data);
        setLoading(false);
      })
      .catch(err => {
        console.error('Failed to load simulation data:', err);
        setError(err.message);
        setLoading(false);
      });
  }, []);

  if (loading) {
    return (
      <div className="loading">
        <h2>Loading simulation data...</h2>
        <p>Make sure to run: <code>npm run generate-data</code></p>
      </div>
    );
  }

  if (error) {
    return (
      <div className="error">
        <h2>Error loading simulation data</h2>
        <p>{error}</p>
        <p>Make sure to run: <code>npm run generate-data</code></p>
      </div>
    );
  }

  return (
    <div className="App">
      <Canvas
        camera={{ position: [0, 0, 8], fov: 60 }}
        style={{ background: '#000011' }}
      >
        <ambientLight intensity={0.3} />
        <pointLight position={[10, 10, 10]} intensity={1} />
        <Stars 
          radius={300} 
          depth={60} 
          count={20000} 
          factor={7} 
          saturation={0} 
          fade={true}
        />
        
        <GlobeVisualizer 
          simulationData={simulationData}
          onPlateSelect={setSelectedPlate}
          selectedPlate={selectedPlate}
        />
        
        <OrbitControls 
          enablePan={true}
          enableZoom={true}
          enableRotate={true}
          minDistance={3}
          maxDistance={15}
        />
      </Canvas>
      
      <InfoPanel
        simulationData={simulationData}
        selectedPlate={selectedPlate}
        onPlateSelect={setSelectedPlate}
      />

      <div className="controls">
        <h4>Controls</h4>
        <p>🖱️ Drag to rotate</p>
        <p>🔍 Scroll to zoom</p>
        <p>🎯 Click plates for details</p>
      </div>
    </div>
  );
}

export default App;
