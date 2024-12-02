mod sorted_sequences_reader;

use crate::config::database_config::DatabaseConfig;
use crate::config::input_transformer_config::Config;
use crate::config::reference_genomes::ReferenceGenomes;
use crate::process::merge_batches::sorted_sequences_reader::SortedSequencesReader;
use crate::types::batch::{Batch, BatchId};
use crate::util::create_new_file;
use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

pub fn merge_batches(
    batch_id_to_sorted_batch: HashMap<BatchId, Batch>,
    config: &Config,
) -> Result<Vec<PathBuf>> {
    let database_config = DatabaseConfig::read_from_file(&config.file_inputs.database_config)?;
    let reference_genomes =
        ReferenceGenomes::read_from_file(&config.file_inputs.reference_genomes)?;
    let merged_output_path = Path::new(&config.output_dir).join("03_merged");
    std::fs::create_dir_all(&merged_output_path)?;

    let mut merged_batch_paths = vec![];

    for (batch_id, batch) in batch_id_to_sorted_batch {
        let mut metadata_file = BufReader::new(File::open(batch.metadata_file)?);
        let mut amino_acid_sequence_files =
            create_sequence_readers(batch.sequence_files.amino_acid_sequence)?;
        let mut aligned_nucleotide_sequence_files =
            create_sequence_readers(batch.sequence_files.aligned_nucleotide_sequence)?;
        let mut unaligned_nucleotide_sequence_files =
            create_sequence_readers(batch.sequence_files.unaligned_nucleotide_sequence)?;

        let target_file_path = merged_output_path.join(format!("{}.ndjson", batch_id.0));
        let mut target_file = create_new_file(&target_file_path)?;

        let mut buffer = String::new();
        while metadata_file.read_line(&mut buffer)? != 0 {
            let metadata: HashMap<String, String> = serde_json::from_str(&buffer)?;
            let metadata_primary_key = &metadata[&database_config.schema.primary_key].clone();

            let mut sequence_entry = SequenceEntry {
                metadata,
                nucleotide_insertions: HashMap::new(),
                amino_acid_insertions: HashMap::new(),
                aligned_nucleotide_sequences: HashMap::new(),
                unaligned_nucleotide_sequences: HashMap::new(),
                aligned_amino_acid_sequences: HashMap::new(),
            };

            for nucleotide_sequence in &reference_genomes.nucleotide_sequences {
                sequence_entry
                    .nucleotide_insertions
                    .insert(nucleotide_sequence.name.clone(), Vec::new());
            }
            for gene in &reference_genomes.genes {
                sequence_entry
                    .amino_acid_insertions
                    .insert(gene.name.clone(), Vec::new());
            }

            process_sequence_files(
                &mut amino_acid_sequence_files,
                metadata_primary_key,
                &mut sequence_entry.aligned_amino_acid_sequences,
            )?;
            process_sequence_files(
                &mut aligned_nucleotide_sequence_files,
                metadata_primary_key,
                &mut sequence_entry.aligned_nucleotide_sequences,
            )?;
            process_sequence_files(
                &mut unaligned_nucleotide_sequence_files,
                metadata_primary_key,
                &mut sequence_entry.unaligned_nucleotide_sequences,
            )?;

            serde_json::to_writer(&target_file, &sequence_entry)?;
            target_file.write_all(b"\n")?;

            buffer.clear();
        }

        merged_batch_paths.push(target_file_path);
    }

    Ok(merged_batch_paths)
}

fn process_sequence_files(
    amino_acid_sequence_files: &mut HashMap<String, SortedSequencesReader>,
    metadata_primary_key: &String,
    target_map: &mut HashMap<String, Option<String>>,
) -> Result<()> {
    for (sequence_name, sequence_file_reader) in amino_acid_sequence_files {
        let next_sequence =
            sequence_file_reader.get_next_entry_for_primary_key(metadata_primary_key);
        let sequence = match next_sequence {
            None => None,
            Some(sequence_result) => Some(sequence_result?.sequence),
        };
        target_map.insert((*sequence_name).clone(), sequence);
    }
    Ok(())
}

fn create_sequence_readers(
    sequence_file_paths: HashMap<String, PathBuf>,
) -> Result<HashMap<String, SortedSequencesReader>> {
    sequence_file_paths
        .iter()
        .map(|(sequence_name, sequence_file_path)| {
            Ok((
                sequence_name.clone(),
                SortedSequencesReader::new(sequence_file_path)?,
            ))
        })
        .collect()
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
struct SequenceEntry {
    metadata: HashMap<String, String>,
    nucleotide_insertions: HashMap<String, Vec<String>>,
    amino_acid_insertions: HashMap<String, Vec<String>>,
    aligned_nucleotide_sequences: HashMap<String, Option<String>>,
    unaligned_nucleotide_sequences: HashMap<String, Option<String>>,
    aligned_amino_acid_sequences: HashMap<String, Option<String>>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test::prepare_batch::{prepare_batch, BatchContent, PRIMARY_KEY};
    use crate::test::read_ndjson::read_ndjson;
    use serde_json::json;

    #[test]
    fn given_missing_sequences_then_merged_batch_contains_null() {
        let prepared_batch = prepare_batch(vec![
            BatchContent::new(json!({PRIMARY_KEY: "1", "someValue": "value1"}))
                .with_aligned_nucleotide_sequence_1(json!({"primary_key": "1", "sequence": "ATGC"}))
                .with_unaligned_nucleotide_sequence_1(json!({"primary_key": "1", "sequence": "A"})),
            BatchContent::new(json!({PRIMARY_KEY: "2", "someValue": "value2"})),
            BatchContent::new(json!({PRIMARY_KEY: "3", "someValue": "value3"}))
                .with_gene_1(json!({"primary_key": "3", "sequence": "MAAA*"})),
        ]);

        let merged_batch_paths =
            merge_batches(prepared_batch.batch_id_to_batch, &prepared_batch.config).unwrap();

        let deserialized_result: HashMap<_, _> =
            read_ndjson::<SequenceEntry>(&merged_batch_paths[0])
                .into_iter()
                .map(|entry| (entry.metadata[PRIMARY_KEY].clone(), entry))
                .collect();

        let empty_nucleotide_sequences = HashMap::from([
            ("nucleotide_sequence_1".to_string(), None),
            ("nucleotide_sequence_2".to_string(), None),
        ]);
        let empty_genes =
            HashMap::from([("gene_1".to_string(), None), ("gene_2".to_string(), None)]);

        assert_eq!(
            deserialized_result["1"].metadata,
            HashMap::from([
                (PRIMARY_KEY.to_string(), "1".to_string()),
                ("someValue".to_string(), "value1".to_string())
            ])
        );
        assert_eq!(
            deserialized_result["1"].aligned_nucleotide_sequences,
            HashMap::from([
                (
                    "nucleotide_sequence_1".to_string(),
                    Some("ATGC".to_string())
                ),
                ("nucleotide_sequence_2".to_string(), None)
            ])
        );
        assert_eq!(
            deserialized_result["1"].unaligned_nucleotide_sequences,
            HashMap::from([
                ("nucleotide_sequence_1".to_string(), Some("A".to_string())),
                ("nucleotide_sequence_2".to_string(), None)
            ])
        );
        assert_eq!(
            deserialized_result["1"].aligned_amino_acid_sequences,
            empty_genes
        );

        assert_eq!(
            deserialized_result["2"].metadata,
            HashMap::from([
                (PRIMARY_KEY.to_string(), "2".to_string()),
                ("someValue".to_string(), "value2".to_string())
            ])
        );
        assert_eq!(
            deserialized_result["2"].aligned_nucleotide_sequences,
            empty_nucleotide_sequences
        );
        assert_eq!(
            deserialized_result["2"].unaligned_nucleotide_sequences,
            empty_nucleotide_sequences
        );
        assert_eq!(
            deserialized_result["2"].aligned_amino_acid_sequences,
            empty_genes
        );

        assert_eq!(
            deserialized_result["3"].metadata,
            HashMap::from([
                (PRIMARY_KEY.to_string(), "3".to_string()),
                ("someValue".to_string(), "value3".to_string())
            ])
        );
        assert_eq!(
            deserialized_result["3"].aligned_nucleotide_sequences,
            empty_nucleotide_sequences
        );
        assert_eq!(
            deserialized_result["3"].unaligned_nucleotide_sequences,
            empty_nucleotide_sequences
        );
        assert_eq!(
            deserialized_result["3"].aligned_amino_acid_sequences,
            HashMap::from([
                ("gene_1".to_string(), Some("MAAA*".to_string())),
                ("gene_2".to_string(), None)
            ])
        );
    }
}
