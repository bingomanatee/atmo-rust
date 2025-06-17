use crate::asthenosphere::AsthenosphereCell;
use crate::rock_store::RockStore;
use h3o::{CellIndex, Resolution};
use rocksdb::Error;

struct CellUpdate {
    cell: AsthenosphereCell,
    energy_k: f64,
    volume: f64,
}

impl CellUpdate {
    fn new(cell: AsthenosphereCell) -> CellUpdate {
        CellUpdate {
            cell: cell.clone(),
            energy_k: cell.energy_k,
            volume: cell.volume,
        }
    }
}

pub fn level_cells(base_cell: CellIndex, store: &RockStore, step: u32, res: Resolution) {
    let child_cells = base_cell.children(res);

    let cells = child_cells
        .filter_map(|cell| match store.get_asth(cell, step) {
            Ok(Some(old_cell)) => Ok(CellUpdate::new(old_cell)),
            Ok(None) => Err(format!("Cell not found: {}, {}", cell, step)),
            Err(e) => Err(format!("Error fetching cell: {:?}", e)),
        })
        .collect();
}
