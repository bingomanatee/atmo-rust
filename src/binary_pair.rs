use h3o::CellIndex;

/// Binary pair for levelling operations
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct BinaryPair {
    pub cell_a: CellIndex,
    pub cell_b: CellIndex,
}

impl BinaryPair {
    /// Create a new binary pair with consistent ordering
    pub fn new(cell_a: CellIndex, cell_b: CellIndex) -> Self {
        // Ensure consistent ordering: higher ID first
        if cell_a > cell_b {
            Self { cell_a, cell_b }
        } else {
            Self { cell_a: cell_b, cell_b: cell_a }
        }
    }

    /// Generate a unique string ID for this pair
    pub fn to_string_id(&self) -> String {
        format!("{}:{}", self.cell_a, self.cell_b)
    }

    /// Parse a binary pair from a string ID
    pub fn from_string_id(id: &str) -> Result<Self, String> {
        let parts: Vec<&str> = id.split(':').collect();
        if parts.len() != 2 {
            return Err(format!("Invalid binary pair ID format: {}", id));
        }

        let cell_a = parts[0].parse::<u64>()
            .map_err(|_| format!("Invalid cell A ID: {}", parts[0]))?;
        let cell_b = parts[1].parse::<u64>()
            .map_err(|_| format!("Invalid cell B ID: {}", parts[1]))?;

        let cell_a = CellIndex::try_from(cell_a)
            .map_err(|_| format!("Invalid H3 cell A: {}", cell_a))?;
        let cell_b = CellIndex::try_from(cell_b)
            .map_err(|_| format!("Invalid H3 cell B: {}", cell_b))?;

        Ok(Self::new(cell_a, cell_b))
    }

    /// Get the other cell in the pair
    pub fn other_cell(&self, cell: CellIndex) -> Option<CellIndex> {
        if self.cell_a == cell {
            Some(self.cell_b)
        } else if self.cell_b == cell {
            Some(self.cell_a)
        } else {
            None
        }
    }

    /// Check if this pair contains the given cell
    pub fn contains(&self, cell: CellIndex) -> bool {
        self.cell_a == cell || self.cell_b == cell
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_binary_pair_creation() {
        let cell_a = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let cell_b = CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap();

        let pair1 = BinaryPair::new(cell_a, cell_b);
        let pair2 = BinaryPair::new(cell_b, cell_a);

        // Should be the same regardless of order
        assert_eq!(pair1, pair2);
        assert_eq!(pair1.to_string_id(), pair2.to_string_id());
    }

    #[test]
    fn test_string_id_roundtrip() {
        let cell_a = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let cell_b = CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap();

        let pair = BinaryPair::new(cell_a, cell_b);
        let id = pair.to_string_id();
        let parsed_pair = BinaryPair::from_string_id(&id).unwrap();

        assert_eq!(pair, parsed_pair);
    }

    #[test]
    fn test_other_cell() {
        let cell_a = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let cell_b = CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap();

        let pair = BinaryPair::new(cell_a, cell_b);

        assert_eq!(pair.other_cell(cell_a), Some(cell_b));
        assert_eq!(pair.other_cell(cell_b), Some(cell_a));

        let cell_c = CellIndex::try_from(0x8a1fb46622f7fff_u64).unwrap();
        assert_eq!(pair.other_cell(cell_c), None);
    }

    #[test]
    fn test_contains() {
        let cell_a = CellIndex::try_from(0x8a1fb46622dffff_u64).unwrap();
        let cell_b = CellIndex::try_from(0x8a1fb46622e7fff_u64).unwrap();
        let cell_c = CellIndex::try_from(0x8a1fb46622f7fff_u64).unwrap();

        let pair = BinaryPair::new(cell_a, cell_b);

        assert!(pair.contains(cell_a));
        assert!(pair.contains(cell_b));
        assert!(!pair.contains(cell_c));
    }
}
