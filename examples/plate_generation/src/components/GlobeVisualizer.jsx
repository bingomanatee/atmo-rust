import React, { useRef, useMemo } from 'react';
import { useFrame } from '@react-three/fiber';
import * as THREE from 'three';
import PlateRenderer from './PlateRenderer';

function GlobeVisualizer({ simulationData, onPlateSelect, selectedPlate }) {
  const globeRef = useRef();
  
  // Create globe geometry and material
  const globeGeometry = useMemo(() => new THREE.SphereGeometry(1, 64, 32), []);
  
  const globeMaterial = useMemo(() => {
    return new THREE.MeshPhongMaterial({
      color: 0x2233aa,
      transparent: true,
      opacity: 0.6,
      wireframe: false,
    });
  }, []);

  // Convert planet radius to normalized scale (planet radius in km to Three.js units)
  const planetRadius = simulationData?.planet?.radius_km || 6372;
  const scale = 1.0; // Keep globe at unit scale, plates will be positioned relative to this

  // Auto-rotate the globe slowly
  useFrame((state, delta) => {
    if (globeRef.current) {
      globeRef.current.rotation.y += delta * 0.1;
    }
  });

  if (!simulationData) {
    return null;
  }

  return (
    <group>
      {/* Main planet sphere */}
      <mesh ref={globeRef} geometry={globeGeometry} material={globeMaterial} scale={scale} />
      
      {/* Render tectonic plates */}
      <PlateRenderer 
        plates={simulationData.plates}
        planetRadius={planetRadius}
        onPlateSelect={onPlateSelect}
        selectedPlate={selectedPlate}
        globeScale={scale}
      />
    </group>
  );
}

export default GlobeVisualizer;
