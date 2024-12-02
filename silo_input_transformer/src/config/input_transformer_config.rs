use anyhow::Result;
use serde::Deserialize;
use std::path::PathBuf;

#[derive(Debug, Deserialize)]
pub struct Config {
    pub file_inputs: FileInputConfig,
    pub file_prefixes: FilePrefixesConfig,
    pub output_dir: String,
    pub batch_size: usize,
}

#[derive(Copy, Clone, Debug)]
pub enum SequenceType {
    AminoAcid,
    Nucleotide,
    UnalignedNucleotide,
}

impl Config {
    pub fn read_from_file(file: String) -> Result<Config> {
        let f = std::fs::File::open(file)?;
        serde_yaml::from_reader(f).map_err(Into::into)
    }
    pub fn get_sequence_file(&self, sequence_type: SequenceType, sequence_name: &str) -> PathBuf {
        let sequence_prefix = match sequence_type {
            SequenceType::AminoAcid => &self.file_prefixes.amino_acid_sequence,
            SequenceType::Nucleotide => &self.file_prefixes.nucleotide_sequence,
            SequenceType::UnalignedNucleotide => &self.file_prefixes.unaligned_nucleotide_sequence,
        };
        PathBuf::from(&self.file_inputs.sequence_file_directory)
            .join(format!("{sequence_prefix}{sequence_name}.fasta"))
    }
}

#[derive(Debug, Deserialize)]
pub struct FileInputConfig {
    pub metadata: String,
    pub database_config: String,
    pub reference_genomes: String,
    pub sequence_file_directory: String,
}

#[derive(Debug, Deserialize)]
pub struct FilePrefixesConfig {
    pub amino_acid_sequence: String,
    pub nucleotide_sequence: String,
    pub unaligned_nucleotide_sequence: String,
}
