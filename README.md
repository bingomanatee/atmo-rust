# Plate Simulation System: Design Overview

This document outlines the core simulation architecture for a tectonic and geodynamic model based on a hexagonal grid system (H3) and pressure-driven mechanics.

## Key Grid Assumptions

* **Planetary resolution base**: H3 resolution **L2**
* **Platelets**: Spawn from **L3** cells

    * Initially snapped to L3 grid but become **free-floating 3D agents**
    * Store: momentum (Vec3), thickness, density, compression, cohesion
* **Pressure cells**: Use **L3** for pressure simulation (matches platelets)

## Simulation Layers

### 1. Platelet Model

* Each L3 cell starts with an **oceanic platelet**
* Platelets have a stack of 1–3 materials:

    * Oceanic
    * Mid-crust (thicker intermediate zones)
    * Continental
* Compression/tension modifies:

    * Stack thickness and density
    * Lateral and vertical momentum

### 2. Cohesion and Neighbor System

* Platelets track neighbor links (initially via H3 L3 neighbors)
* Cohesion between neighbors evolves over time:

    * High compression → stronger bonds
    * Tension or shear → weakening, potential fracture

### 3. Pressure Model

* Each L3 cell stores a scalar **pressure value**
* **Uniform base pressure** across planet
* Pressure flows **laterally** to lower-pressure neighbors
* If pressure **accumulates beyond a threshold**, it triggers an **upwell** (volcano event)

### 4. Volcano / Upwell System (Edge-Free)

* Upwells occur at pressure cells, not along explicit edges
* Volcano event triggers if:

    * Pressure > threshold (e.g. ambient + crust thickness factor)
    * Crust above is not too thick
* Upwell creates or thickens a platelet stack

    * Adds mass to local platelet
    * Slightly boosts vertical/lateral momentum

### 5. Downwells / Absorption

* High downward pressure zones with:

    * **Dense crust above**
    * **Low upward pressure**
* May cause platelet removal or density compression
* Simulates subduction without explicit edge modeling

### 6. Plate Motion Mechanics

* Each platelet's momentum adjusted by:

    * Pressure gradients (sideways push)
    * Crust collisions (repulsion, compression)
* Platelets with aligned vectors may bond into **plates**
* Cohesion governs bond stability and momentum transfer

## Simulation Cycle

1. **Partition platelets** into H3 L2 bins
2. **Compute neighbor-based interactions** (tension, compression)
3. **Move platelets** by momentum
4. **Update cohesion** and potential bonds
5. **Adjust global pressure field**
6. **Check for volcanoes (upwells) and downwells**
7. **Form or fracture plates**
8. Repeat

## Parallelization Plan

* Partition planet into L2 zones
* Each worker processes platelets within zone and 1-ring neighbors
* Store data in RocksDB column families for parallel-safe read/write
* Pressure model operates on L3 grid and is globally synchronized each cycle

## Long-Term Dynamics

* Continental plates emerge through collision/volcano events
* Platelets can dissolve due to downwells under continental mass
* Plate boundaries are emergent, not explicitly modeled
* Volcano distribution is a product of crust/pressure dynamics

## Advantages

* Simplifies traditional ridge/trench modeling
* Unifies volcanoes, crust creation, and plate growth
* Enables efficient parallel execution
* Scales naturally with H3 grid resolution

---

This model allows realism and procedural emergence with reasonable compute cost and no hard-coded fault lines. Plates are not primary; **platelets and pressure are the true actors**.
