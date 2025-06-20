pub mod sim;
mod planet;
mod plate_generator;
mod vary;
pub mod rock_store;
mod helpers;
pub mod sim_manager;
mod h3_utils;
mod geoconverter;
mod thermal_plate_generator;
mod asthenosphere;
mod steps;
pub mod constants;
pub mod gif_exporter;
pub mod png_exporter;
mod erode_cell;
// New simulation system
pub mod asth_cell_next;
pub mod binary_pair;
pub mod sim_next;

#[cfg(test)]
mod sim_next_test;

// Deprecated linked system (preserved for reference)
mod deprecated;
mod asthenosphere_linked;
mod asth_mutate;
// pub mod asth_sim_linked; // Temporarily disabled due to compilation issues
// mod asl_leveller; // Depends on asth_sim_linked
// mod asl_direct_transfer_leveller; // Depends on asth_sim_linked
mod asl_binary_pair_leveller;