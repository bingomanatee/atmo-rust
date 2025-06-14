import React from 'react';

function InfoPanel({ simulationData, selectedPlate, onPlateSelect }) {
  if (!simulationData) {
    return null;
  }

  const { sim, planet, plates } = simulationData;

  // Calculate total plate coverage
  const totalPlateArea = plates.reduce((sum, plate) => {
    return sum + Math.PI * Math.pow(plate.radius_km, 2);
  }, 0);
  
  const planetSurfaceArea = 4 * Math.PI * Math.pow(planet.radius_km, 2);
  const coverage = (totalPlateArea / planetSurfaceArea) * 100;

  // Find largest and smallest plates
  const sortedPlates = [...plates].sort((a, b) => b.radius_km - a.radius_km);
  const largestPlate = sortedPlates[0];
  const smallestPlate = sortedPlates[sortedPlates.length - 1];

  return (
    <div className="info-panel">
      <h3>Plate Generation</h3>
      
      <div className="stat">
        <span className="stat-label">Planet Radius:</span>
        <span className="stat-value">{planet.radius_km.toLocaleString()} km</span>
      </div>
      
      <div className="stat">
        <span className="stat-label">Mantle Density:</span>
        <span className="stat-value">{planet.mantle_density_gcm3.toFixed(2)} g/cm³</span>
      </div>
      
      <div className="stat">
        <span className="stat-label">Total Plates:</span>
        <span className="stat-value">{plates.length}</span>
      </div>
      
      <div className="stat">
        <span className="stat-label">Plate Coverage:</span>
        <span className="stat-value">{coverage.toFixed(1)}%</span>
      </div>

      {selectedPlate ? (
        <div style={{ marginTop: '20px', borderTop: '1px solid #333', paddingTop: '15px' }}>
          <h3>Selected Plate</h3>
          <div className="stat">
            <span className="stat-label">Radius:</span>
            <span className="stat-value">{selectedPlate.radius_km.toLocaleString()} km</span>
          </div>
          <div className="stat">
            <span className="stat-label">Thickness:</span>
            <span className="stat-value">{selectedPlate.thickness_km} km</span>
          </div>
          <div className="stat">
            <span className="stat-label">Density:</span>
            <span className="stat-value">{selectedPlate.density.toFixed(2)} g/cm³</span>
          </div>
          <div className="stat">
            <span className="stat-label">Center:</span>
            <span className="stat-value">
              ({selectedPlate.center.x.toFixed(2)}, {selectedPlate.center.y.toFixed(2)}, {selectedPlate.center.z.toFixed(2)})
            </span>
          </div>
          <button 
            onClick={() => onPlateSelect(null)}
            style={{
              marginTop: '10px',
              padding: '5px 10px',
              background: '#4a9eff',
              border: 'none',
              borderRadius: '4px',
              color: 'white',
              cursor: 'pointer'
            }}
          >
            Deselect
          </button>
        </div>
      ) : (
        <div style={{ marginTop: '20px', borderTop: '1px solid #333', paddingTop: '15px' }}>
          <h4>Plate Statistics</h4>
          <div className="stat">
            <span className="stat-label">Largest Plate:</span>
            <span className="stat-value">{largestPlate.radius_km.toLocaleString()} km</span>
          </div>
          <div className="stat">
            <span className="stat-label">Smallest Plate:</span>
            <span className="stat-value">{smallestPlate.radius_km.toLocaleString()} km</span>
          </div>
          <p style={{ fontSize: '0.9rem', color: '#ccc', marginTop: '10px' }}>
            Click on a plate to see details
          </p>
        </div>
      )}
    </div>
  );
}

export default InfoPanel;
