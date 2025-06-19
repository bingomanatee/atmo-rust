use crate::asthenosphere::AsthenosphereCell;
use crate::rock_store::RockStore;
use h3o::{CellIndex, Resolution};
use std::collections::HashMap;
use crate::deprecated::asth_cell_update::AsthCellUpdate;

pub fn erode_cells(
    base_cell: CellIndex,
    store: &RockStore,
    step: u32,
    res: Resolution,
) -> Result<(), String> {
    let child_cells = base_cell.children(res);

    // Collect cells into a HashMap indexed by cell.cell (CellIndex)
    let cells: Result<HashMap<CellIndex, AsthCellUpdate>, String> = child_cells
        .map(|cell| match store.get_asth(cell, step) {
            Ok(Some(old_cell)) => Ok((cell, AsthCellUpdate::new(old_cell))),
            Ok(None) => Err(format!("Cell not found: {}, {}", cell, step)),
            Err(e) => Err(format!("Error fetching cell: {:?}", e)),
        })
        .collect();

    match cells {
        Ok(mut cell_map) => {
            // Collect transfer operations first to avoid borrowing conflicts
            let mut transfers_to_make = Vec::new();

            // Iterate over references to find transfers needed
            for (&cell_index, cell_update) in &cell_map {
                let mut lowest_index = cell_index;
                let mut lowest_volume = cell_update.volume;
                for &neighbor_cell in &cell_update.cell.neighbors {
                    if cell_map.contains_key(&neighbor_cell) {
                        let neighbor = cell_map.get(&neighbor_cell).unwrap();
                        if neighbor.cell.volume < lowest_volume {
                            lowest_index = neighbor_cell;
                            lowest_volume = neighbor.cell.volume;
                        }
                    }
                }
                if lowest_volume < cell_update.cell.volume {
                    transfers_to_make.push((cell_index, lowest_index));
                }
            }

            // Now execute the transfers
            for (from_index, to_index) in transfers_to_make {
                transfer(&mut cell_map, from_index, to_index);
            }

            let mut updated_cells: Vec<AsthenosphereCell>= Vec::new();

            for (cell, update) in cell_map {
                if (update.energy_k != update.cell.energy_j) || update.volume != update.cell.volume {
                    updated_cells.push(AsthenosphereCell {
                        energy_j: update.energy_k,
                        volume: update.volume,
                        ..update.cell
                    })
                }
            }

            store.put_asth_batch(&updated_cells);

            Ok(())
        }
        Err(e) => Err(e),
    }
}

fn transfer(cell_map: &mut HashMap<CellIndex, AsthCellUpdate>, from_index: CellIndex, to_index: CellIndex) {
    if !cell_map.contains_key(&from_index) || !cell_map.contains_key(&to_index) {
        return;
    }

    let source_volume = cell_map.get(&from_index).unwrap().volume;
    let dest_volume = cell_map.get(&to_index).unwrap().volume;

    let volume_difference = source_volume - dest_volume;
    if volume_difference <= 0.0 {
        return;
    }

    let amount_to_move = TRANSFER_RATIO * volume_difference;

    // Update source cell (remove volume and proportional energy)
    if let Some(source) = cell_map.get_mut(&from_index) {
        let energy_ratio = source.energy_k / source.volume;
        source.volume -= amount_to_move;
        source.energy_k -= amount_to_move * energy_ratio;
        source.cell.volume = source.volume;
        source.cell.energy_j = source.energy_k;
    }

    // Calculate source energy ratio before borrowing the destination cell
    let source_energy_ratio = cell_map.get(&from_index).unwrap().energy_k / cell_map.get(&from_index).unwrap().volume;

    // Update destination cell (add volume and proportional energy)
    if let Some(dest) = cell_map.get_mut(&to_index) {
        dest.volume += amount_to_move;
        dest.energy_k += amount_to_move * source_energy_ratio;
        dest.cell.volume = dest.volume;
        dest.cell.energy_j = dest.energy_k;
    }
}

const TRANSFER_RATIO: f64 = 0.125;
