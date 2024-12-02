use crate::config::input_transformer_config::Config;
use crate::util::create_new_file;
use anyhow::Result;
use std::fs::File;
use std::path::{Path, PathBuf};

pub fn merge_final_result(merge_batch_paths: Vec<PathBuf>, config: &Config) -> Result<()> {
    let mut merged_batches_file =
        create_new_file(&Path::new(&config.output_dir).join("silo_input.ndjson"))?;

    for merge_batch_path in merge_batch_paths {
        let mut merge_batch_file = File::open(&merge_batch_path)?;
        std::io::copy(&mut merge_batch_file, &mut merged_batches_file)?;
    }

    Ok(())
}
