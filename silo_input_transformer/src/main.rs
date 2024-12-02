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
use std::fs;

fn main() -> Result<()> {
    let config = Config::read_from_file("config.yaml".to_string())?;

    if fs::metadata(&config.output_dir).is_ok() {
        fs::remove_dir_all(&config.output_dir)?;
    }

    let batch_id_to_batch = split_data_into_batches(&config)?;
    let batch_id_to_sorted_batch = sort_batches(batch_id_to_batch, &config)?;
    let merged_batch_paths = merge_batches(batch_id_to_sorted_batch, &config)?;

    merge_final_result(merged_batch_paths, &config)?;

    Ok(())
}
