#!/bin/bash

# Define the arguments
INPUT_FILE="tests/data/samples_large/A1_05_2024_10_08/20241024_2411515907/alignments/REF_aln_trim.bam"
#INPUT_FILE="../../../data/sr2silo/demo_real/A1_10_2024_09_30/20241018_AAG55WNM5/alignments/REF_aln_trim.bam"
SAMPLE_ID="A1_05_2024_10_08"
OUTPUT_FILE="./results_neo/silo_input.ndjson"
TIMELINE_FILE="tests/data/samples_large/timeline_A1_05_2024_10_08.tsv"
# Optional: specify LAPIS URL for custom reference genomes
# If not provided, uses default SARS-CoV-2 references (NCBI Reference Sequence: NC_045512.2)
LAPIS_URL="https://lapis.cov-spectrum.org/open/v2"

# Run using sr2silo CLI
# Use submit-to-loculus command to upload and submit processed files to SILO

# Example 1: Use custom LAPIS server for references
sr2silo process-from-vpipe \
    --input-file "$INPUT_FILE" \
    --sample-id "$SAMPLE_ID" \
    --timeline-file "$TIMELINE_FILE" \
    --lapis-url "$LAPIS_URL" \
    --output-fp "$OUTPUT_FILE"

# Example 2: Use default SARS-CoV-2 references (recommended for most use cases)
# sr2silo process-from-vpipe \
#     --input-file "$INPUT_FILE" \
#     --sample-id "$SAMPLE_ID" \
#     --timeline-file "$TIMELINE_FILE" \
#     --output-fp "$OUTPUT_FILE"

# Uncomment the following lines to upload and submit the processed file to SILO
sr2silo submit-to-loculus \
    --processed-file "$OUTPUT_FILE.zst" \
    --sample-id "$SAMPLE_ID"
