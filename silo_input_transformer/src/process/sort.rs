use crate::config::database_config::DatabaseConfig;
use crate::config::input_transformer_config::Config;
use crate::config::reference_genomes::ReferenceGenomes;
use crate::types::batch::{Batch, BatchId};
use crate::util::create_new_file;
use anyhow::Result;
use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

pub fn sort_batches(
    batch_id_to_batch: HashMap<BatchId, Batch>,
    config: &Config,
) -> Result<HashMap<BatchId, Batch>> {
    let sorted_batches_path = Path::new(&config.output_dir).join("02_sorted");
    fs::create_dir_all(&sorted_batches_path)?;

    let database_config = DatabaseConfig::read_from_file(&config.file_inputs.database_config)?;
    let reference_genomes =
        ReferenceGenomes::read_from_file(&config.file_inputs.reference_genomes)?;

    let batch_id_to_sorted_batch = batch_id_to_batch
        .into_iter()
        .map(|(batch_id, batch)| {
            let buf = sorted_batches_path.join(batch_id.0.clone());
            let sorted_batch = Batch::new(&buf, &reference_genomes)?;

            let sequence_files = &batch.sequence_files;
            let sorted_sequence_files = &sorted_batch.sequence_files;
            let sequence_file_paths_to_sort = vec![
                zip_sequence_files(
                    &sequence_files.amino_acid_sequence,
                    &sorted_sequence_files.amino_acid_sequence,
                ),
                zip_sequence_files(
                    &sequence_files.unaligned_nucleotide_sequence,
                    &sorted_sequence_files.unaligned_nucleotide_sequence,
                ),
                zip_sequence_files(
                    &sequence_files.aligned_nucleotide_sequence,
                    &sorted_sequence_files.aligned_nucleotide_sequence,
                ),
            ]
            .into_iter()
            .flatten()
            .map(|(path, target_path)| (path, target_path, "primary_key"))
            .chain(vec![(
                batch.metadata_file.clone(),
                sorted_batch.metadata_file.clone(),
                database_config.schema.primary_key.as_str(),
            )]);

            for (unsorted_file_path, target_file_path, key) in sequence_file_paths_to_sort {
                sort_file_by_key(&unsorted_file_path, key, &target_file_path)?;
            }

            Ok((batch_id, sorted_batch))
        })
        .collect::<Result<HashMap<BatchId, Batch>>>()?;

    Ok(batch_id_to_sorted_batch)
}

fn zip_sequence_files(
    sequence_files: &HashMap<String, PathBuf>,
    sorted_sequence_files: &HashMap<String, PathBuf>,
) -> Vec<(PathBuf, PathBuf)> {
    sequence_files
        .iter()
        .map(|(sequence_name, path)| (path.clone(), sorted_sequence_files[sequence_name].clone()))
        .collect()
}

fn sort_file_by_key(
    unsorted_file_path: &PathBuf,
    key: &str,
    target_file_path: &PathBuf,
) -> Result<()> {
    let mut deserialized_lines = fs::read_to_string(unsorted_file_path)?
        .lines()
        .map(|line| Ok(serde_json::from_str::<HashMap<String, String>>(line)?))
        .collect::<Result<Vec<HashMap<String, String>>>>()?;

    deserialized_lines.sort_by(|a, b| a[key].cmp(&b[key]));

    let mut target_file = create_new_file(&target_file_path)?;
    for deserialized_line in deserialized_lines {
        serde_json::to_writer(&target_file, &deserialized_line)?;
        target_file.write_all(b"\n")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::config::input_transformer_config::SequenceType;
    use crate::test::prepare_batch::{prepare_batch, BatchContent, PRIMARY_KEY};
    use crate::test::read_ndjson::read_ndjson;
    use crate::types::batch::{Batch, BatchId};
    use crate::types::sequence::Sequence;
    use serde_json::json;
    use std::collections::HashMap;

    #[test]
    fn test_sorting() {
        let prepared_batch = prepare_batch(
            [9, 1, 5, 8, 0, 3, 6, 2, 4, 7, 10]
                .map(generate_batch_content)
                .to_vec(),
        );

        let sorted_batches =
            super::sort_batches(prepared_batch.batch_id_to_batch, &prepared_batch.config).unwrap();

        let batch = &sorted_batches[&BatchId::from(0)];
        assert_metadata_looks_good(batch);
        assert_sequence_files_look_good(
            batch,
            SequenceType::UnalignedNucleotide,
            "nucleotide_sequence_1",
            "unaligned_nucleotide_sequence_1",
        );
        assert_sequence_files_look_good(
            batch,
            SequenceType::UnalignedNucleotide,
            "nucleotide_sequence_2",
            "unaligned_nucleotide_sequence_2",
        );
        assert_sequence_files_look_good(
            batch,
            SequenceType::Nucleotide,
            "nucleotide_sequence_1",
            "aligned_nucleotide_sequence_1",
        );
        assert_sequence_files_look_good(
            batch,
            SequenceType::Nucleotide,
            "nucleotide_sequence_2",
            "aligned_nucleotide_sequence_2",
        );
        assert_sequence_files_look_good(batch, SequenceType::AminoAcid, "gene_1", "gene_1");
        assert_sequence_files_look_good(batch, SequenceType::AminoAcid, "gene_2", "gene_2");
    }

    fn generate_batch_content(index: usize) -> BatchContent {
        BatchContent::new(
            json!({PRIMARY_KEY: index.to_string(), "someValue": format!("value{}", index)}),
        )
        .with_unaligned_nucleotide_sequence_1(json!({
            "primary_key": index.to_string(),
            "sequence": format!("unaligned_nucleotide_sequence_1_{index}")
        }))
        .with_unaligned_nucleotide_sequence_2(json!({
            "primary_key": index.to_string(),
            "sequence": format!("unaligned_nucleotide_sequence_2_{index}")
        }))
        .with_aligned_nucleotide_sequence_1(json!({
            "primary_key": index.to_string(),
            "sequence": format!("aligned_nucleotide_sequence_1_{index}")
        }))
        .with_aligned_nucleotide_sequence_2(json!({
            "primary_key": index.to_string(),
            "sequence": format!("aligned_nucleotide_sequence_2_{index}")
        }))
        .with_gene_1(json!({
            "primary_key": index.to_string(),
            "sequence": format!("gene_1_{index}")}
        ))
        .with_gene_2(json!({
            "primary_key": index.to_string(),
            "sequence": format!("gene_2_{index}")
        }))
    }

    fn assert_metadata_looks_good(batch: &Batch) {
        let metadata = read_ndjson::<HashMap<String, String>>(&batch.metadata_file);
        assert_eq!(metadata.len(), 11);
        assert_eq!(
            metadata[0],
            HashMap::from([
                (PRIMARY_KEY.to_string(), "0".to_string()),
                ("someValue".to_string(), "value0".to_string())
            ])
        );
        assert_eq!(
            metadata[2],
            HashMap::from([
                (PRIMARY_KEY.to_string(), "10".to_string()),
                ("someValue".to_string(), "value10".to_string())
            ])
        );
        assert_eq!(
            metadata[10],
            HashMap::from([
                (PRIMARY_KEY.to_string(), "9".to_string()),
                ("someValue".to_string(), "value9".to_string())
            ])
        );
    }

    fn assert_sequence_files_look_good(
        batch: &Batch,
        sequence_type: SequenceType,
        sequence_name: &str,
        sequence_content_prefix: &str,
    ) {
        let sequences = read_ndjson::<Sequence>(
            &batch.sequence_files.get_sequence_files(sequence_type)[sequence_name],
        );

        assert_eq!(sequences.len(), 11);
        assert_eq!(sequences[0].primary_key, "0");
        assert_eq!(
            sequences[0].sequence,
            format!("{sequence_content_prefix}_0")
        );
        assert_eq!(sequences[2].primary_key, "10");
        assert_eq!(
            sequences[2].sequence,
            format!("{sequence_content_prefix}_10")
        );
        assert_eq!(sequences[10].primary_key, "9");
        assert_eq!(
            sequences[10].sequence,
            format!("{sequence_content_prefix}_9")
        );
    }
}
