# Usage of `process-from-vpipe` Command

The `process-from-vpipe` command is used to convert BAM alignments to SILO format. It reads metadata from a timeline file.

## Important Note

Before using the `process-from-vpipe` command, ensure that primers have been trimmed off from the input BAM file. This means that either soft-clippings or hard-clippings must be applied to the input file to remove primer sequences.

## Command Overview

The basic usage of the command is as follows:

```bash
sr2silo process-from-vpipe \
    --input-file INPUT.bam \
    --sample-id SAMPLE_ID \
    --batch-id BATCH_ID \
    --timeline-file TIMELINE.tsv \
    --output-fp OUTPUT.ndjson \
    --reference sars-cov-2
```

### Required Arguments

- `--input-file, -i`: Path to the input BAM alignment file.
- `--sample-id, -s`: Sample ID to use for metadata.
- `--batch-id, -b`: Batch ID to use for metadata.
- `--timeline-file, -t`: Path to the timeline metadata file.
- `--output-fp, -o`: Path for the output file (will be auto-suffixed with `.ndjson.zst`).

### Optional Arguments

- `--reference, -r`: The nucleotide/amino acid reference from the resources folder. Default is `sars-cov-2`.
- `--skip-merge/--no-skip-merge`: Whether to skip merging of paired-end reads. Default is `--no-skip-merge`.

## Example Usage

```bash
sr2silo process-from-vpipe \
    --input-file ./data/sample/alignments/REF_aln_trim.bam \
    --sample-id "A1_05_2024_10_08" \
    --batch-id "20241024_2411515907" \
    --timeline-file ./data/timeline.tsv \
    --output-fp ./results/output.ndjson \
    --reference sars-cov-2
```

To upload the results to SILO, use the separate `submit-to-loculus` command:

```bash
sr2silo submit-to-loculus \
    --processed-file ./results/output.ndjson.zst \
    --sample-id "A1_05_2024_10_08"
```
