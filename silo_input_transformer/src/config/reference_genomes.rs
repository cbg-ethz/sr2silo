use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ReferenceGenomes {
    pub nucleotide_sequences: Vec<ReferenceSequence>,
    pub genes: Vec<ReferenceSequence>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ReferenceSequence {
    pub name: String,
    pub sequence: String,
}

impl ReferenceGenomes {
    pub fn read_from_file(file: &str) -> anyhow::Result<ReferenceGenomes> {
        let f = std::fs::File::open(file)?;
        serde_json::from_reader(f).map_err(Into::into)
    }
}
