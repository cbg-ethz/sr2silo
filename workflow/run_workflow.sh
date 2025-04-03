#!/bin/bash
# Run Snakemake workflow with conda environment support

set -e

# Change to the workflow directory
cd "$(dirname "$0")"

# Create or update conda environments
echo "Setting up conda environments..."
snakemake --use-conda --conda-create-envs-only

# Run the workflow
echo "Running workflow..."
snakemake --use-conda -c 1 "$@"

echo "Workflow completed."
