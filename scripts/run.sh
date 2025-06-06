#!/bin/bash

# Define the arguments
INPUT_FILE="tests/data/samples_large/A1_05_2024_10_08/20241024_2411515907/alignments/REF_aln_trim.bam"
#INPUT_FILE="../../../data/sr2silo/demo_real/A1_10_2024_09_30/20241018_AAG55WNM5/alignments/REF_aln_trim.bam"
SAMPLE_ID="A1_05_2024_10_08"
BATCH_ID="20241024_2411515907"
OUTPUT_FILE="./results_neo/silo_input.ndjson"
TIMELINE_FILE="tests/data/samples_large/timeline_A1_05_2024_10_08.tsv"
PRIMERS_FILE="./resources/sars-cov-2/primers/primers.yaml"
NUC_REFERENCE="sars-cov-2"

# Run using sr2silo CLI
# Use submit-to-loculus command to upload and submit processed files to SILO
sr2silo process-from-vpipe \
    --input-file "$INPUT_FILE" \
    --sample-id "$SAMPLE_ID" \
    --batch-id "$BATCH_ID" \
    --timeline-file "$TIMELINE_FILE" \
    --primer-file "$PRIMERS_FILE" \
    --output-fp "$OUTPUT_FILE" \
    --reference "$NUC_REFERENCE"

# Uncomment the following lines to upload and submit the processed file to SILO
sr2silo submit-to-loculus \
    --processed-file "$OUTPUT_FILE.zst" \
    --sample-id "$SAMPLE_ID"
