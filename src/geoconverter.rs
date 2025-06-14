use glam::Vec3;
use h3o::{CellIndex, LatLng, Resolution};
use rand::Rng;
use rand::seq::SliceRandom;

/// Utility to convert between Vec3 and CellIndex for a spherical planet
pub struct GeoCellConverter {
    pub radius_km: f64,
    pub resolution: Resolution,
}

impl GeoCellConverter {
    pub fn new(radius_km: f64, resolution: Resolution) -> GeoCellConverter {
        GeoCellConverter {
            radius_km,
            resolution,
        }
    }

    fn cell_radius_km(&self) -> f32 {
        let edge_len = self.resolution.edge_length_km() as f32;
        edge_len * 1.1547
    }

    pub fn randomize_point_on_sphere(&self, base: Vec3) -> Vec3 {
        let mut rng = rand::rng();
        let max_jitter = self.cell_radius_km();
        let jitter = Vec3::new(
            rng.random_range(-max_jitter..max_jitter),
            rng.random_range(-max_jitter..max_jitter),
            rng.random_range(-max_jitter..max_jitter),
        );

        (base + jitter).normalize() * (self.radius_km as f32)
    }

    fn cell_to_radians(cell: &CellIndex) -> (f64, f64) {
        let ll = LatLng::from(*cell);
        (ll.lat_radians(), ll.lng_radians())
    }
    fn vec3_to_lat_lng(v: Vec3) -> Option<LatLng> {
        let lat = v.z.asin().to_degrees() as f64;
        let lon = v.y.atan2(v.x).to_degrees() as f64;
        LatLng::new(lat, lon).ok()
    }
    fn cell_to_unit_vec3(cell: &CellIndex) -> Vec3 {
        let (lat, lng) = GeoCellConverter::cell_to_radians(&cell);
        let (cos_lat, sin_lat) = (lat.cos(), lat.sin());
        let (cos_lng, sin_lng) = (lng.cos(), lng.sin());
        // Convert spherical coordinates to Cartesian coordinates on unit sphere
        Vec3::new(
            (cos_lat * cos_lng) as f32,
            (cos_lat * sin_lng) as f32,
            sin_lat as f32,
        )
    }
    pub fn vec3_to_cell(&self, point: Vec3) -> Option<CellIndex> {
        let unit_vec = point.normalize();
        let lat_lng = GeoCellConverter::vec3_to_lat_lng(unit_vec)?;
        Some(lat_lng.to_cell(self.resolution))
    }

    pub fn cell_to_vec3(&self, cell: CellIndex) -> Vec3 {
        let unit = GeoCellConverter::cell_to_unit_vec3(&cell);
        unit * self.radius_km as f32
    }
}
