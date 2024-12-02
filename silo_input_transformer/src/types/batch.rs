use crate::config::input_transformer_config::SequenceType;
use crate::config::reference_genomes::{ReferenceGenomes, ReferenceSequence};
use serde::Serialize;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

#[derive(Debug, Serialize, Eq, PartialEq, Hash, Clone)]
pub struct BatchId(pub String);

impl From<usize> for BatchId {
    fn from(id: usize) -> Self {
        BatchId(format!("batch_{id}"))
    }
}

#[derive(Debug)]
pub struct Batch {
    pub metadata_file: PathBuf,
    pub sequence_files: SequenceFiles,
}

impl Batch {
    pub fn new(path: &Path, reference_genomes: &ReferenceGenomes) -> anyhow::Result<Self> {
        std::fs::create_dir_all(path)?;

        Ok(Batch {
            metadata_file: path.join("metadata.ndjson"),
            sequence_files: SequenceFiles::new(path, reference_genomes),
        })
    }
}

#[derive(Debug)]
pub struct SequenceFiles {
    pub aligned_nucleotide_sequence: HashMap<String, PathBuf>,
    pub amino_acid_sequence: HashMap<String, PathBuf>,
    pub unaligned_nucleotide_sequence: HashMap<String, PathBuf>,
}

impl SequenceFiles {
    pub fn new(path: &Path, reference_genomes: &ReferenceGenomes) -> Self {
        SequenceFiles {
            aligned_nucleotide_sequence: get_sequence_files_from_reference_sequences(
                path,
                &reference_genomes.nucleotide_sequences,
                "aligned_",
            ),
            amino_acid_sequence: get_sequence_files_from_reference_sequences(
                path,
                &reference_genomes.genes,
                "gene_",
            ),
            unaligned_nucleotide_sequence: get_sequence_files_from_reference_sequences(
                path,
                &reference_genomes.nucleotide_sequences,
                "unaligned_",
            ),
        }
    }

    pub fn get_sequence_files(&self, sequence_type: SequenceType) -> &HashMap<String, PathBuf> {
        match sequence_type {
            SequenceType::AminoAcid => &self.amino_acid_sequence,
            SequenceType::Nucleotide => &self.aligned_nucleotide_sequence,
            SequenceType::UnalignedNucleotide => &self.unaligned_nucleotide_sequence,
        }
    }
}

fn get_sequence_files_from_reference_sequences(
    path: &Path,
    reference_sequences: &[ReferenceSequence],
    prefix: &str,
) -> HashMap<String, PathBuf> {
    reference_sequences
        .iter()
        .map(|sequence| {
            (
                sequence.name.clone(),
                path.join(format!("{prefix}{}.ndjson", sequence.name)),
            )
        })
        .collect()
}
