use glam::Vec3;
use h3o::{CellIndex, DirectedEdgeIndex, LatLng, Resolution};
use rand::Rng;
use rand::seq::SliceRandom;

pub struct CellAndBase(CellIndex, CellIndex);

pub struct PointSampler {
    primary_cells: Vec<CellIndex>,
    current_index: usize,
    planet_radius_km: i32,
}

impl PointSampler {
    /// Create a new sampler for one plate
    pub fn new(planet_radius_km: i32) -> Self {
        PointSampler {
            primary_cells: CellIndex::base_cells().collect(),
            current_index: 0,
            planet_radius_km,
        }
    }

    fn ll_tuple(cell: &CellIndex) -> (f64, f64) {
        let ll = LatLng::from(*cell);
        (ll.lat_radians(), ll.lng_radians())
    }

    fn random_child(parent: &CellIndex) -> CellIndex {
        let mut rng = rand::rng();

        let children: Vec<CellIndex> = parent.children(Resolution::Three).collect();
        let idx = rng.random_range(0..children.len());
        children[idx]
    }

    fn cell_to_unit(cell: &CellIndex) -> Vec3 {
        let (lat, lng) = PointSampler::ll_tuple(&cell);
        let (cos_lat, sin_lat) = (lat.cos(), lat.sin());
        let (cos_lng, sin_lng) = (lng.cos(), lng.sin());
        // Convert spherical coordinates to Cartesian coordinates on unit sphere
        Vec3::new(
            (cos_lat * cos_lng) as f32,
            (cos_lat * sin_lng) as f32,
            sin_lat as f32,
        )
    }

    /*
    returns a "semi-random" point on the planet.
    each time it is called it returns a random h3 location in a different
    primary cell region.

    The points are "spaced out" and not perfectly overalpping
    which is _not_ a trait of true randomness but prevents plates from
    being perfectly on top of each other.

    In theory they will be as far apart from each other as a primary radii
    but in truth they will be from (level 3 radii ... (2 x primary radii - level 3 radii))
    apart from each other.
     */
    pub fn random_point_on_planet(&mut self) -> Vec3 {
        if self.current_index == 0 {
            let mut rng = rand::rng();
            self.primary_cells.shuffle(&mut rng);
        }

        let primary_cell = self.primary_cells[self.current_index];
        let cell = PointSampler::random_child(&primary_cell);

        let unit = PointSampler::cell_to_unit(&cell);
        let point = unit * self.planet_radius_km as f32;

        // 4) Advance the index, wrapping around
        self.current_index = (self.current_index + 1) % self.primary_cells.len();

        point
    }
}

pub struct H3Utils;

impl H3Utils {
    // note - given the opaccity of rust closures this is not a useful pattern -
    // use the below method in the future
    /// Iterate over all H3 cells at the given resolution, calling `callback` for each cell.
    pub fn iter_at<F>(resolution: Resolution, mut callback: F)
    where
        F: FnMut(CellIndex),
    {
        for base_cell in CellIndex::base_cells() {
            for cell in base_cell.children(resolution) {
                callback(cell);
            }
        }
    }

    pub fn iter_cells_with_base(
        resolution: Resolution,
    ) -> impl Iterator<Item = (CellIndex, CellIndex)> {
        CellIndex::base_cells()
            .into_iter()
            .flat_map(move |base_cell| {
                base_cell
                    .children(resolution)
                    .into_iter()
                    .map(move |child| (child, base_cell))
            })
    }

    /// Get all neighboring cells for a given cell index
    pub fn neighbors_for(cell_index: CellIndex) -> Vec<CellIndex> {
        // Get all neighboring cells
        let mut neighbors: Vec<CellIndex> = cell_index
            .grid_disk::<Vec<_>>(1)
            .into_iter()
            .filter(|&neighbor| neighbor != cell_index)
            .collect();
        neighbors
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::planet::EARTH_RADIUS_KM;
    use glam::Vec3;

    #[test]
    fn test_point_sampler_creation() {
        let radius_km = EARTH_RADIUS_KM;
        let sampler = PointSampler::new(radius_km);
        assert_eq!(sampler.planet_radius_km, radius_km);
        assert_eq!(sampler.current_index, 0);
        assert_eq!(sampler.primary_cells.len(), 122); // h3 base cells count
    }

    #[test]
    fn test_random_point_on_planet_radius() {
        let radius_km = 6371;
        let mut sampler = PointSampler::new(radius_km);

        let point = sampler.random_point_on_planet();
        let length = point.length();

        // The length should be close to radius_km (allow small floating error)
        let radius_f64 = radius_km as f32;
        let epsilon = 1.0; // 1 km tolerance
        assert!(
            (length - radius_f64).abs() < epsilon,
            "Point length {} not close to radius {}",
            length,
            radius_f64
        );
    }

    #[test]
    fn test_random_points_are_different() {
        let radius_km = 6371;
        let mut sampler = PointSampler::new(radius_km);

        let point1 = sampler.random_point_on_planet();
        let point2 = sampler.random_point_on_planet();

        // Points should not be exactly equal
        assert_ne!(point1, point2, "Two consecutive points should not be equal");
    }

    #[test]
    fn test_points_are_on_sphere_surface() {
        let radius_km = 6371;
        let mut sampler = PointSampler::new(radius_km);

        for _ in 0..10 {
            let point = sampler.random_point_on_planet();
            let length = point.length();
            let radius_f64 = radius_km as f32;
            let epsilon = 1.0;
            assert!(
                (length - radius_f64).abs() < epsilon,
                "Point length {} not close to radius {}",
                length,
                radius_f64
            );
        }
    }

    #[test]
    fn test_closest_neighbor_distances_statistics() {
        let radius_km = 6371;
        let mut sampler = PointSampler::new(radius_km);

        let num_points = 30;
        let mut points: Vec<Vec3> = Vec::with_capacity(num_points);

        for _ in 0..num_points {
            points.push(sampler.random_point_on_planet());
        }

        // For each point, find the closest neighbor distance
        let mut closest_distances = Vec::with_capacity(num_points);

        for (i, p) in points.iter().enumerate() {
            let mut min_dist = f32::MAX;
            for (j, q) in points.iter().enumerate() {
                if i == j {
                    continue;
                }
                let dist = (*p - *q).length();
                if dist < min_dist {
                    min_dist = dist;
                }
            }
            closest_distances.push(min_dist);
        }

        // Compute statistics: min, average, std deviation
        let min_distance = closest_distances.iter().cloned().fold(f32::MAX, f32::min);
        let sum: f32 = closest_distances.iter().sum();
        let avg_distance = sum / closest_distances.len() as f32;

        let variance = closest_distances
            .iter()
            .map(|d| {
                let diff = d - avg_distance as f32;
                diff * diff
            })
            .sum::<f32>()
            / closest_distances.len() as f32;
        let std_dev = variance.sqrt();

        println!(
            "Closest neighbor distances (km): min = {:.2}, avg = {:.2}, std dev = {:.2}",
            min_distance, avg_distance, std_dev
        );

        // Basic sanity checks
        assert!(
            min_distance > 100.0,
            "Minimum closest neighbor distance should be > 100"
        );
        assert!(
            avg_distance > 1000.0,
            "Average closest neighbor distance should be > 1000"
        );
        assert!(std_dev >= 300.0, "Standard deviation should be >= 300");
    }
}
