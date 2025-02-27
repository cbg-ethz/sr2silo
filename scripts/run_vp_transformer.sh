#!/bin/bash

# Define the arguments
INPUT_FILE="./tests/data/samples/A1_05_2024_10_08/20241024_2411515907/alignments/REF_aln_trim.bam"
SAMPLE_ID="A1_05_2024_10_08"
BATCH_ID="20241024_2411515907"
OUTPUT_FILE="./results_neo/silo_input.ndjson"
TIMELINE_FILE="./tests/data/samples/timeline_A1_05_2024_10_08.tsv"
PRIMERS_FILE="./tests/data/samples_large/primers.yaml"
NUC_REFERENCE="sars-cov-2"

# Run the Python script with the arguments
python scripts/vp_transformer.py \
    --input_file "$INPUT_FILE" \
    --sample_id "$SAMPLE_ID" \
    --batch_id "$BATCH_ID" \
    --timeline_file "$TIMELINE_FILE" \
    --primer_file "$PRIMERS_FILE" \
    --output_fp "$OUTPUT_FILE" \
    --reference "$NUC_REFERENCE"
