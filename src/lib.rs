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
pub mod convection;

// New simulation system
pub mod asth_cell_next;
pub mod binary_pair;
pub mod sim_next;

#[cfg(test)]
mod sim_next_test;

#[cfg(test)]
mod volcano_bias_test;

// Deprecated linked system (preserved for reference)
mod asth_mutate;
mod plate;
