use crate::config::database_config::DatabaseConfig;
use crate::config::input_transformer_config::{Config, SequenceType};
use crate::config::reference_genomes::{ReferenceGenomes, ReferenceSequence};
use crate::types::batch::{Batch, BatchId};
use crate::types::sequence::Sequence;
use crate::util::create_new_file;
use anyhow::bail;
use anyhow::Result;
use csv::ReaderBuilder;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

pub fn split_data_into_batches(config: &Config) -> Result<HashMap<BatchId, Batch>> {
    let batches_path = Path::new(&config.output_dir).join("01_batches");
    std::fs::create_dir_all(&batches_path)?;

    let (primary_key_to_batch, batch_id_to_batch) =
        split_metadata_into_batches(&config, &batches_path)?;
    split_sequence_files_into_batches(&primary_key_to_batch, &config, &batch_id_to_batch)?;

    let primary_key_to_batch_file =
        create_new_file(&batches_path.join("primary_key_to_batch_id.json"))?;

    serde_json::to_writer(&primary_key_to_batch_file, &primary_key_to_batch)?;

    Ok(batch_id_to_batch)
}

fn split_sequence_files_into_batches(
    primary_key_to_batch_id: &HashMap<String, BatchId>,
    config: &Config,
    batch_id_to_batch: &HashMap<BatchId, Batch>,
) -> Result<()> {
    let reference_genomes =
        ReferenceGenomes::read_from_file(&config.file_inputs.reference_genomes)?;

    let sequences_and_type = vec![
        (SequenceType::AminoAcid, &reference_genomes.genes),
        (
            SequenceType::UnalignedNucleotide,
            &reference_genomes.nucleotide_sequences,
        ),
        (
            SequenceType::Nucleotide,
            &reference_genomes.nucleotide_sequences,
        ),
    ];

    for (sequence_type, sequences) in sequences_and_type {
        for sequence in sequences {
            let sequence_file_path = config.get_sequence_file(sequence_type, &sequence.name);
            let mut fasta_reader = needletail::parse_fastx_file(&sequence_file_path)?;

            let batch_id_to_sequence_file =
                get_batch_id_to_sequence_file_map(batch_id_to_batch, sequence_type, &sequence)?;

            while let Some(record) = fasta_reader.next() {
                let record = record?;

                let sequence_id = String::from_utf8_lossy(record.id());

                let Some(batch_id) = primary_key_to_batch_id.get(&sequence_id.to_string()) else {
                    bail!("Primary key not found in metadata file for sequence {sequence_id:?} in file {sequence_file_path:?}")
                };

                let sequence_entry = Sequence {
                    primary_key: sequence_id.to_string(),
                    sequence: String::from_utf8_lossy(&record.seq()).to_string(),
                };

                let mut sequence_file = &batch_id_to_sequence_file[batch_id];
                serde_json::to_writer(sequence_file, &sequence_entry)?;
                sequence_file.write_all(b"\n")?;
            }
        }
    }
    Ok(())
}

fn get_batch_id_to_sequence_file_map(
    batch_id_to_batch: &HashMap<BatchId, Batch>,
    sequence_type: SequenceType,
    sequence: &&ReferenceSequence,
) -> Result<HashMap<BatchId, File>> {
    batch_id_to_batch
        .iter()
        .map(|(batch_id, batch)| {
            let sequence_file_path =
                batch.sequence_files.get_sequence_files(sequence_type)[&sequence.name].clone();
            (batch_id, sequence_file_path)
        })
        .map(|(batch_id, sequence_file_path)| {
            let file = create_new_file(&sequence_file_path)?;
            Ok((batch_id.clone(), file))
        })
        .collect::<Result<HashMap<BatchId, File>>>()
}

fn split_metadata_into_batches(
    config: &Config,
    batches_path: &PathBuf,
) -> Result<(HashMap<String, BatchId>, HashMap<BatchId, Batch>)> {
    let database_config = DatabaseConfig::read_from_file(&config.file_inputs.database_config)?;
    let reference_genomes =
        ReferenceGenomes::read_from_file(&config.file_inputs.reference_genomes)?;

    let mut primary_key_to_batch = HashMap::<String, BatchId>::new();
    let mut batch_id_to_batch = HashMap::<BatchId, Batch>::new();

    let mut tsv_reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(&config.file_inputs.metadata)?;

    let mut current_batch = 0;
    let mut entries_in_current_batch = 0;

    let mut current_batch_file = insert_new_batch(
        &batches_path,
        &reference_genomes,
        &mut batch_id_to_batch,
        current_batch,
    )?;

    for tsv_entry in tsv_reader.deserialize() {
        if entries_in_current_batch >= config.batch_size {
            current_batch += 1;
            entries_in_current_batch = 0;
            current_batch_file = insert_new_batch(
                &batches_path,
                &reference_genomes,
                &mut batch_id_to_batch,
                current_batch,
            )?;
        }

        let record: HashMap<String, String> = tsv_entry?;

        let primary_key = match record.get(&database_config.schema.primary_key) {
            None => bail!("Primary key not found in metadata file for entry {record:?}"),
            Some(primary_key) => primary_key,
        };

        primary_key_to_batch.insert(primary_key.to_string(), current_batch.into());
        entries_in_current_batch += 1;

        serde_json::to_writer(&current_batch_file, &record)?;
        current_batch_file.write_all(b"\n")?;
    }

    Ok((primary_key_to_batch, batch_id_to_batch))
}

fn insert_new_batch(
    batches_path: &&PathBuf,
    reference_genomes: &ReferenceGenomes,
    batch_id_to_batch: &mut HashMap<BatchId, Batch>,
    current_batch: usize,
) -> Result<File> {
    batch_id_to_batch.insert(
        current_batch.into(),
        Batch::new(
            &batches_path.join(current_batch.to_string()),
            &reference_genomes,
        )?,
    );
    let metadata_file_path = &batch_id_to_batch[&current_batch.into()].metadata_file;
    Ok(create_new_file(&metadata_file_path)?)
}
