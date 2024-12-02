use serde::Deserialize;
use std::fs::File;
use std::io::read_to_string;
use std::path::PathBuf;

pub fn read_ndjson<T: for<'a> Deserialize<'a>>(merged_batch_paths: &PathBuf) -> Vec<T> {
    read_to_string(File::open(&merged_batch_paths).unwrap())
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str::<T>(&line).unwrap())
        .collect()
}
