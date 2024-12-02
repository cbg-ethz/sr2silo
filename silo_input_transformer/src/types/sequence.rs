use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize)]
pub struct Sequence {
    pub primary_key: String,
    pub sequence: String,
}
