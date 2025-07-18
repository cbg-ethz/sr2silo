#!/bin/bash
# =============================================================================
# Sr2Silo Workflow SLURM Batch Script
# =============================================================================
#
# This is a standalone .sbatch file for running the sr2silo workflow
#
# PREREQUISITES:
# --------------
# - The sr2silo-workflow conda environment must be created first:
#   make create-workflow-env
#
# HOW TO USE:
# -----------
# 1. Edit the CORES variable below to set the number of cores (default: 10)
# 2. Edit the USE_EXISTING_ENV variable to choose conda environment mode:
#    - false: Use workflow-defined conda environments (default)
#    - true:  Use your currently active conda environment
# 3. Edit the LOG_DIR variable to set custom log directory if needed
# 4. Edit the CONDA_ENV variable if you want to use a different environment
# 5. Submit with: sbatch sr2silo_workflow.sbatch
#
# =============================================================================

#SBATCH --job-name=sr2silo
#SBATCH --mail-type=END
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=8G
#SBATCH --tmp=160G
#SBATCH --time=12:00:00
#SBATCH --output=/cluster/work/bewi/members/koehng/logs/sr2silo.out
#SBATCH --error=/cluster/work/bewi/members/koehng/logs/sr2silo.err

# =============================================================================
# Configuration Variables - Edit these as needed
# =============================================================================

# Number of cores to use (should match --cpus-per-task above)
CORES=20

# Use existing conda environment (true) or workflow-defined conda envs (false)
USE_EXISTING_ENV=false

# Log directory (should match --output and --error paths above)
LOG_DIR="/cluster/work/bewi/members/koehng/logs"

# Conda environment name (change this if you want to use a different environment)
CONDA_ENV="sr2silo-workflow"

# =============================================================================
# Script Execution
# =============================================================================

echo "=============================================="
echo "Sr2Silo Workflow Starting"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Cores: $CORES"
echo "Use existing env: $USE_EXISTING_ENV"
echo "Log directory: $LOG_DIR"
echo "Conda environment: $CONDA_ENV"
echo "=============================================="

# Load eth_proxy module for internet access on ETH cluster
echo "Loading eth_proxy module..."
module load eth_proxy
echo "Proxy loaded: ${https_proxy}"

# Initialize conda for bash shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Check if the conda environment exists
if ! conda env list | grep -q "^${CONDA_ENV} "; then
    echo "ERROR: Conda environment '${CONDA_ENV}' not found!"
    echo "Please create it first with: make create-workflow-env"
    exit 1
fi

# Activate the conda environment that contains snakemake
echo "Activating conda environment: $CONDA_ENV"
conda activate "$CONDA_ENV"

# Verify snakemake is available
echo "Snakemake version: $(snakemake --version)"

# Create log directory if it doesn't exist
mkdir -p "$LOG_DIR"

# Prepare and execute Snakemake command
if [ "$USE_EXISTING_ENV" = true ]; then
    # Use the active conda environment
    echo "Using existing conda environment"
    echo "Running: snakemake -c $CORES" --rerun-incomplete
    snakemake -c $CORES
else
    # Use the conda environments defined in the workflow
    echo "Using workflow-defined conda environments"
    echo "Running: snakemake --use-conda --conda-frontend conda -c $CORES"
    snakemake --use-conda --conda-frontend conda -c $CORES --rerun-incomplete
fi

SNAKEMAKE_EXIT_CODE=$?

echo "=============================================="
echo "Sr2Silo Workflow Completed"
echo "Exit code: $SNAKEMAKE_EXIT_CODE"
echo "=============================================="

exit $SNAKEMAKE_EXIT_CODE
