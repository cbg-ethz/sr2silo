#!/bin/bash
# =============================================================================
# Sr2Silo Workflow SLURM Batch Script - Custom Configuration Template
# =============================================================================
#
# This version allows easy customization of parameters
# Copy this file and modify the parameters as needed
#
# HOW TO USE:
# -----------
# 1. Copy this file: cp sr2silo_workflow_template.sbatch my_workflow.sbatch
# 2. Edit the SBATCH directives and variables below as needed
# 3. Submit with: sbatch my_workflow.sbatch
#
# =============================================================================

# SLURM Job Configuration - Modify these as needed
#SBATCH --job-name=sr2silo
#SBATCH --mail-type=END
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=CORES_PLACEHOLDER
#SBATCH --mem-per-cpu=8G
#SBATCH --tmp=160G
#SBATCH --time=12:00:00
#SBATCH --output=LOG_DIR_PLACEHOLDER/sr2silo.out
#SBATCH --error=LOG_DIR_PLACEHOLDER/sr2silo.err

# =============================================================================
# Configuration Variables - Edit these as needed
# =============================================================================

# Number of cores to use (should match --cpus-per-task above)
CORES=CORES_PLACEHOLDER

# Use existing conda environment (true) or workflow-defined conda envs (false)
USE_EXISTING_ENV=USE_ENV_PLACEHOLDER

# Log directory (should match --output and --error paths above)
LOG_DIR="LOG_DIR_PLACEHOLDER"

# =============================================================================
# Script Execution - Usually no need to modify below this line
# =============================================================================

echo "=============================================="
echo "Sr2Silo Workflow Starting"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Cores: $CORES"
echo "Use existing env: $USE_EXISTING_ENV"
echo "Log directory: $LOG_DIR"
echo "=============================================="

# Load eth_proxy module for internet access on ETH cluster
echo "Loading eth_proxy module..."
module load eth_proxy
echo "Proxy loaded: ${https_proxy}"

# Create log directory if it doesn't exist
mkdir -p "$LOG_DIR"

# Prepare and execute Snakemake command
if [ "$USE_EXISTING_ENV" = true ]; then
    # Use the active conda environment
    echo "Using existing conda environment"
    echo "Running: snakemake -c $CORES"
    snakemake -c $CORES
else
    # Use the conda environments defined in the workflow
    echo "Using workflow-defined conda environments"
    echo "Running: snakemake --use-conda --conda-frontend conda -c $CORES"
    snakemake --use-conda --conda-frontend conda -c $CORES
fi

SNAKEMAKE_EXIT_CODE=$?

echo "=============================================="
echo "Sr2Silo Workflow Completed"
echo "Exit code: $SNAKEMAKE_EXIT_CODE"
echo "=============================================="

exit $SNAKEMAKE_EXIT_CODE
