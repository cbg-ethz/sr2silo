#!/bin/bash

# Define the arguments
SAMPLE_DIR="./tests/data/samples/A1_05_2024_10_08/20241024_2411515907/alignments"
SAMPLE_ID="A1_05_2024_10_08"
BATCH_ID="20241024_2411515907"
RESULTS_DIR="./results_neo"
TIMELINE_FILE="./tests/data/samples/timeline_A1_05_2024_10_08.tsv"
PRIMERS_FILE="./tests/data/samples_large/primers.yaml"
NUC_REFERENCE="sars-cov-2"
DATABASE_CONFIG="./scripts/database_config.yaml"

# Run the Python script with the arguments
python scripts/vp_transformer.py \
    --sample_dir "$SAMPLE_DIR" \
    --sample_id "$SAMPLE_ID" \
    --batch_id "$BATCH_ID" \
    --result_dir "$RESULTS_DIR" \
    --timeline_file "$TIMELINE_FILE" \
    --primer_file "$PRIMERS_FILE" \
    --nuc_reference "$NUC_REFERENCE" \
