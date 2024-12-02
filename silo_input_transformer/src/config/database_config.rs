use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct DatabaseConfig {
    pub schema: DatabaseSchema,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct DatabaseSchema {
    pub primary_key: String,
}

impl DatabaseConfig {
    pub fn read_from_file(file: &str) -> Result<DatabaseConfig> {
        let f = std::fs::File::open(file)
            .context(format!("Tried to read database config from {file}"))?;
        serde_yaml::from_reader(f).map_err(Into::into)
    }
}
