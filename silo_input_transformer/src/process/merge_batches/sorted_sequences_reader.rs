use crate::types::sequence::Sequence;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::mem;
use std::path::PathBuf;

pub struct SortedSequencesReader {
    sequence_file: BufReader<File>,
    sequence_buffer: String,
    next: Option<Result<Sequence>>,
}

impl SortedSequencesReader {
    pub fn new(sequence_file_path: &PathBuf) -> Result<Self> {
        let sequence_file = BufReader::new(
            File::open(sequence_file_path)
                .context(format!("Reading sequences from {sequence_file_path:?}"))?,
        );
        let mut reader = SortedSequencesReader {
            sequence_file,
            sequence_buffer: String::new(),
            next: None,
        };
        reader.next = reader.get_next_sequence();
        Ok(reader)
    }

    pub fn get_next_entry_for_primary_key(
        &mut self,
        primary_key: &str,
    ) -> Option<Result<Sequence>> {
        match &self.next {
            Some(Ok(sequence)) if sequence.primary_key != primary_key => return None,
            _ => {}
        };
        let new_next = self.get_next_sequence();
        mem::replace(&mut self.next, new_next)
    }

    fn get_next_sequence(&mut self) -> Option<Result<Sequence>> {
        self.sequence_buffer.clear();
        match self.sequence_file.read_line(&mut self.sequence_buffer) {
            Ok(0) => None,
            Ok(_) => Some(serde_json::from_str(&self.sequence_buffer).context(format!(
                "Failed to parse sequence: {}",
                self.sequence_buffer
            ))),
            Err(err) => Some(Err(err.into())),
        }
    }
}
