use anyhow::Context;
use anyhow::Result;
use std::fs::File;
use std::path::Path;

pub fn create_new_file(path: &Path) -> Result<File> {
    File::create(path).with_context(|| format!("Failed to create file: {:?}", path))
}
