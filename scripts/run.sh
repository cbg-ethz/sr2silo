#!/bin/bash

# Define the arguments
INPUT_FILE="./tests/data/samples/A1_05_2024_10_08/20241024_2411515907/alignments/REF_aln_trim.bam"
#INPUT_FILE="../../../data/sr2silo/demo_real/A1_10_2024_09_30/20241018_AAG55WNM5/alignments/REF_aln_trim.bam"
SAMPLE_ID="A1_05_2024_10_08"
BATCH_ID="20241024_2411515907"
OUTPUT_FILE="./results_neo/silo_input.ndjson"
TIMELINE_FILE="./tests/data/samples/timeline_A1_05_2024_10_08.tsv"
PRIMERS_FILE="./tests/data/samples_large/primers.yaml"
NUC_REFERENCE="sars-cov-2"

# Run using sr2silo CLI
# Add --upload flag if you want to upload and submit to SILO
sr2silo import-to-loculus \
    --input-file "$INPUT_FILE" \
    --sample-id "$SAMPLE_ID" \
    --batch-id "$BATCH_ID" \
    --timeline-file "$TIMELINE_FILE" \
    --primer-file "$PRIMERS_FILE" \
    --output-fp "$OUTPUT_FILE" \
    --reference "$NUC_REFERENCE" \
    # --upload # Uncomment to enable upload to S3 and submission to SILO
