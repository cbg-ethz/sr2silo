mod config;
mod process;
#[cfg(test)]
mod test;
mod types;
mod util;

use crate::config::input_transformer_config::Config;
use anyhow::Result;
use process::merge_batches::merge_batches;
use process::merge_final_result::merge_final_result;
use process::sort::sort_batches;
use process::split_data_into_batches::split_data_into_batches;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use std::fs;

#[pyfunction]
fn run_with_config(config_path: &str) -> PyResult<()> {
    run(config_path).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
}

fn run(config_path: &str) -> Result<()> {
    println!("Reading config from file: {}", config_path);
    let config = Config::read_from_file(config_path.to_string())?;
    println!("Config read successfully: {:?}", config);

    if fs::metadata(&config.output_dir).is_ok() {
        println!("Output directory exists, removing: {}", config.output_dir);
        fs::remove_dir_all(&config.output_dir)?;
    }

    println!("Splitting data into batches...");
    let batch_id_to_batch = split_data_into_batches(&config)?;
    println!("Data split into batches successfully.");

    println!("Sorting batches...");
    let batch_id_to_sorted_batch = sort_batches(batch_id_to_batch, &config)?;
    println!("Batches sorted successfully.");

    println!("Merging batches...");
    let merged_batch_paths = merge_batches(batch_id_to_sorted_batch, &config)?;
    println!("Batches merged successfully.");

    println!("Merging final result...");
    merge_final_result(merged_batch_paths, &config)?;
    println!("Final result merged successfully.");

    Ok(())
}

#[pymodule]
fn silo_input_transformer(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(run_with_config, m)?)?;
    Ok(())
}
