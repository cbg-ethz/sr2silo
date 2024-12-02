use crate::config::database_config::{DatabaseConfig, DatabaseSchema};
use crate::config::input_transformer_config::{Config, FileInputConfig, FilePrefixesConfig};
use crate::config::reference_genomes::{ReferenceGenomes, ReferenceSequence};
use crate::types::batch::{Batch, BatchId};
use crate::util::create_new_file;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use tempdir::TempDir;

#[derive(Debug)]
pub struct PreparedBatch {
    pub batch_id_to_batch: HashMap<BatchId, Batch>,
    pub config: Config,
    /** This directory will be deleted when the PreparedBatch goes out of scope */
    _input_tempdir: TempDir,
    /** This directory will be deleted when the PreparedBatch goes out of scope */
    pub _output_tempdir: TempDir,
}

#[derive(Debug, Clone)]
pub struct BatchContent {
    metadata: serde_json::Value,
    unaligned_nucleotide_sequence_1: Option<serde_json::Value>,
    unaligned_nucleotide_sequence_2: Option<serde_json::Value>,
    aligned_nucleotide_sequences_1: Option<serde_json::Value>,
    aligned_nucleotide_sequences_2: Option<serde_json::Value>,
    gene_1: Option<serde_json::Value>,
    gene_2: Option<serde_json::Value>,
}

impl BatchContent {
    pub fn new(metadata: serde_json::Value) -> Self {
        BatchContent {
            metadata,
            unaligned_nucleotide_sequence_1: None,
            unaligned_nucleotide_sequence_2: None,
            aligned_nucleotide_sequences_1: None,
            aligned_nucleotide_sequences_2: None,
            gene_1: None,
            gene_2: None,
        }
    }

    pub fn with_unaligned_nucleotide_sequence_1(mut self, sequence: serde_json::Value) -> Self {
        self.unaligned_nucleotide_sequence_1 = Some(sequence);
        self
    }

    pub fn with_unaligned_nucleotide_sequence_2(mut self, sequence: serde_json::Value) -> Self {
        self.unaligned_nucleotide_sequence_2 = Some(sequence);
        self
    }

    pub fn with_aligned_nucleotide_sequence_1(mut self, sequence: serde_json::Value) -> Self {
        self.aligned_nucleotide_sequences_1 = Some(sequence);
        self
    }

    pub fn with_aligned_nucleotide_sequence_2(mut self, sequence: serde_json::Value) -> Self {
        self.aligned_nucleotide_sequences_2 = Some(sequence);
        self
    }

    pub fn with_gene_1(mut self, sequence: serde_json::Value) -> Self {
        self.gene_1 = Some(sequence);
        self
    }

    pub fn with_gene_2(mut self, sequence: serde_json::Value) -> Self {
        self.gene_2 = Some(sequence);
        self
    }
}

pub const PRIMARY_KEY: &'static str = "my_primary_key";

pub fn prepare_batch(batch_contents: Vec<BatchContent>) -> PreparedBatch {
    let input_dir = TempDir::new("batch").unwrap();
    let output_dir = TempDir::new("output").unwrap();

    let config = Config {
        output_dir: path_to_string(output_dir.path()),
        batch_size: 3,
        file_inputs: FileInputConfig {
            metadata: "".to_string(),
            database_config: path_to_string(input_dir.path().join("database_config.yaml")),
            reference_genomes: path_to_string(input_dir.path().join("reference_genomes.json")),
            sequence_file_directory: "".to_string(),
        },
        file_prefixes: FilePrefixesConfig {
            amino_acid_sequence: "".to_string(),
            nucleotide_sequence: "".to_string(),
            unaligned_nucleotide_sequence: "".to_string(),
        },
    };

    write_database_config(&config);
    write_reference_genomes(&config);

    let reference_genomes = get_test_reference_genomes();

    let batch_id = BatchId::from(0);
    let batch = Batch::new(&output_dir.path().join(&batch_id.0), &reference_genomes).unwrap();

    write_batch_files(batch_contents, &batch);

    PreparedBatch {
        batch_id_to_batch: HashMap::from([(batch_id, batch)]),
        config,
        _input_tempdir: input_dir,
        _output_tempdir: output_dir,
    }
}

fn path_to_string<P: AsRef<Path>>(path: P) -> String {
    path.as_ref().to_str().unwrap().to_string()
}

fn write_reference_genomes(config: &Config) {
    serde_json::to_writer(
        File::create(&config.file_inputs.reference_genomes).unwrap(),
        &get_test_reference_genomes(),
    )
    .unwrap()
}

fn write_database_config(config: &Config) {
    serde_yaml::to_writer(
        File::create(&config.file_inputs.database_config).unwrap(),
        &DatabaseConfig {
            schema: DatabaseSchema {
                primary_key: PRIMARY_KEY.to_string(),
            },
        },
    )
    .unwrap()
}

fn write_batch_files(batch_contents: Vec<BatchContent>, batch: &Batch) {
    let mut metadata_file = create_new_file(&batch.metadata_file).unwrap();
    let mut nucleotide_sequence_1_file = create_new_file(
        &batch.sequence_files.unaligned_nucleotide_sequence["nucleotide_sequence_1"],
    )
    .unwrap();
    let mut nucleotide_sequence_2_file = create_new_file(
        &batch.sequence_files.unaligned_nucleotide_sequence["nucleotide_sequence_2"],
    )
    .unwrap();
    let mut aligned_nucleotide_sequence_1_file =
        create_new_file(&batch.sequence_files.aligned_nucleotide_sequence["nucleotide_sequence_1"])
            .unwrap();
    let mut aligned_nucleotide_sequence_2_file =
        create_new_file(&batch.sequence_files.aligned_nucleotide_sequence["nucleotide_sequence_2"])
            .unwrap();
    let mut gene_1_file =
        create_new_file(&batch.sequence_files.amino_acid_sequence["gene_1"]).unwrap();
    let mut gene_2_file =
        create_new_file(&batch.sequence_files.amino_acid_sequence["gene_2"]).unwrap();

    for batch_content in batch_contents {
        serde_json::to_writer(&metadata_file, &batch_content.metadata).unwrap();
        metadata_file.write_all(b"\n").unwrap();

        write_sequence(
            &mut nucleotide_sequence_1_file,
            batch_content.unaligned_nucleotide_sequence_1,
        );
        write_sequence(
            &mut nucleotide_sequence_2_file,
            batch_content.unaligned_nucleotide_sequence_2,
        );
        write_sequence(
            &mut aligned_nucleotide_sequence_1_file,
            batch_content.aligned_nucleotide_sequences_1,
        );
        write_sequence(
            &mut aligned_nucleotide_sequence_2_file,
            batch_content.aligned_nucleotide_sequences_2,
        );
        write_sequence(&mut gene_1_file, batch_content.gene_1);
        write_sequence(&mut gene_2_file, batch_content.gene_2);
    }
}

fn write_sequence(sequence_file: &mut File, sequence: Option<serde_json::Value>) {
    if let Some(sequence) = sequence {
        serde_json::to_writer(&*sequence_file, &sequence).unwrap();
        sequence_file.write_all(b"\n").unwrap();
    }
}

fn get_test_reference_genomes() -> ReferenceGenomes {
    ReferenceGenomes {
        nucleotide_sequences: vec![
            ReferenceSequence {
                name: "nucleotide_sequence_1".to_string(),
                sequence: "ATGC".to_string(),
            },
            ReferenceSequence {
                name: "nucleotide_sequence_2".to_string(),
                sequence: "GGAATTCC".to_string(),
            },
        ],
        genes: vec![
            ReferenceSequence {
                name: "gene_1".to_string(),
                sequence: "MAAA*".to_string(),
            },
            ReferenceSequence {
                name: "gene_2".to_string(),
                sequence: "MTTTTT*".to_string(),
            },
        ],
    }
}
