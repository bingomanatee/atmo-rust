import React, { useMemo } from 'react';
import * as THREE from 'three';

function PlateRenderer({ plates, planetRadius, onPlateSelect, selectedPlate, globeScale }) {
  
  // Generate colors for plates based on their properties
  const plateColors = useMemo(() => {
    if (!plates) return [];
    
    return plates.map((plate, index) => {
      // Color based on density and thickness
      const densityNorm = (plate.density - 2.0) / 3.0; // Normalize density roughly 2-5 g/cm³
      const thicknessNorm = plate.thickness_km / 100.0; // Normalize thickness
      
      // Create color variation
      const hue = (index * 137.5) % 360; // Golden angle for good distribution
      const saturation = 0.6 + densityNorm * 0.4;
      const lightness = 0.4 + thicknessNorm * 0.3;
      
      return new THREE.Color().setHSL(hue / 360, saturation, lightness);
    });
  }, [plates]);

  // Create plate geometries
  const plateGeometries = useMemo(() => {
    if (!plates) return [];
    
    return plates.map((plate) => {
      // Convert plate center from Vec3 to normalized position on sphere
      const center = new THREE.Vector3(plate.center.x, plate.center.y, plate.center.z);
      center.normalize();
      
      // Calculate plate radius as a fraction of planet radius
      const plateRadiusRatio = plate.radius_km / planetRadius;
      const angularRadius = Math.asin(plateRadiusRatio); // Angular radius in radians
      
      // Create a circular geometry for the plate
      const segments = Math.max(8, Math.floor(angularRadius * 32)); // More segments for larger plates
      const geometry = new THREE.CircleGeometry(angularRadius, segments);
      
      // Position the circle on the sphere surface
      // Create a quaternion to rotate the circle to face the center position
      const up = new THREE.Vector3(0, 0, 1);
      const quaternion = new THREE.Quaternion().setFromUnitVectors(up, center);
      
      return {
        geometry,
        position: center.clone().multiplyScalar(globeScale + 0.001), // Slightly above surface
        quaternion,
        plate
      };
    });
  }, [plates, planetRadius, globeScale]);

  const handlePlateClick = (event, plate) => {
    event.stopPropagation();
    onPlateSelect(plate);
  };

  if (!plates || plates.length === 0) {
    return null;
  }

  return (
    <group>
      {plateGeometries.map((plateGeo, index) => {
        const isSelected = selectedPlate && selectedPlate.id === plateGeo.plate.id;
        const color = plateColors[index];
        
        return (
          <mesh
            key={plateGeo.plate.id}
            geometry={plateGeo.geometry}
            position={plateGeo.position}
            quaternion={plateGeo.quaternion}
            onClick={(event) => handlePlateClick(event, plateGeo.plate)}
          >
            <meshPhongMaterial
              color={color}
              transparent={true}
              opacity={isSelected ? 0.9 : 0.7}
              side={THREE.DoubleSide}
              emissive={isSelected ? color.clone().multiplyScalar(0.3) : new THREE.Color(0x000000)}
            />
          </mesh>
        );
      })}
    </group>
  );
}

export default PlateRenderer;
