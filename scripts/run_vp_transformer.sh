#!/bin/bash

# Define the arguments
INPUT_FILE="./tests/data/samples/A1_05_2024_10_08/20241024_2411515907/alignments/REF_aln_trim.bam"
SAMPLE_ID="A1_05_2024_10_08"
BATCH_ID="20241024_2411515907"
RESULTS_DIR="./results_neo"
TIMELINE_FILE="./tests/data/samples/timeline_A1_05_2024_10_08.tsv"
PRIMERS_FILE="./tests/data/samples_large/primers.yaml"
NUC_REFERENCE="sars-cov-2"
DATABASE_CONFIG="./scripts/database_config.yaml"

# Run the Python script with the arguments
python scripts/vp_transformer.py \
    --input_file "$INPUT_FILE" \
    --sample_id "$SAMPLE_ID" \
    --batch_id "$BATCH_ID" \
    --result_dir "$RESULTS_DIR" \
    --timeline_file "$TIMELINE_FILE" \
    --primer_file "$PRIMERS_FILE" \
    --reference "$NUC_REFERENCE" \
    --database_config "$DATABASE_CONFIG"
