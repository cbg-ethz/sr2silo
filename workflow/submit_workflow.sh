#!/bin/bash
# =============================================================================
# Sr2Silo Workflow Submission Script
# =============================================================================
#
# HOW TO USE:
# -----------
# This script submits the sr2silo workflow as a SLURM job with configurable options.
#
# Basic usage:
#   ./submit_workflow.sh
#       - Submits job with 4 cores using workflow-defined conda environments
#
# With options:
#   ./submit_workflow.sh -c 8
#       - Uses 8 cores instead of the default 4
#
#   ./submit_workflow.sh -e
#       - Uses your currently active conda environment instead of creating new ones
#
#   ./submit_workflow.sh -e -c 12
#       - Uses your current conda environment with 12 cores
#
# The script will submit a SLURM job that:
#   - Allocates the specified number of cores
#   - Assigns 8GB memory per core
#   - Sets 160GB temporary disk space
#   - Sets a time limit of 12 hours
#   - Names the job "sr2silo" for easy identification
#   - Outputs logs to /cluster/home/koehng/logs/
#
# =============================================================================

# Submit Snakemake workflow to SLURM with conda environment options

# Parse command-line arguments
USE_EXISTING_ENV=false
CORES=4  # Default number of cores

function show_usage {
  echo "Usage: $0 [-e] [-c CORES]"
  echo "  -e            Use existing conda environment (default: use workflow-defined conda envs)"
  echo "  -c CORES      Number of cores to use (default: 4)"
  exit 1
}

while getopts "ec:" opt; do
  case $opt in
    e) USE_EXISTING_ENV=true ;;
    c) CORES=$OPTARG ;;
    *) show_usage ;;
  esac
done

# Prepare Snakemake command
if [ "$USE_EXISTING_ENV" = true ]; then
  # Use the active conda environment
  echo "Using existing conda environment"
  SNAKEMAKE_CMD="snakemake -c $CORES"
else
  # Use the conda environments defined in the workflow
  echo "Using workflow-defined conda environments"
  SNAKEMAKE_CMD="snakemake --use-conda -c $CORES"
fi

echo "Running with $CORES cores"

# Submit job to SLURM
# -J sr2silo : Sets the job name to "sr2silo" for easy identification in the SLURM queue
sbatch \
    --mail-type=END \
    --ntasks=1 \
    --cpus-per-task=$CORES \
    --mem-per-cpu=8G \
    --tmp=160G \
    --time=12:00:00 \
    -o /cluster/home/koehng/logs/sr2silo.out \
    -e /cluster/home/koehng/logs/sr2silo.err \
    -J sr2silo \
    --wrap "$SNAKEMAKE_CMD"
