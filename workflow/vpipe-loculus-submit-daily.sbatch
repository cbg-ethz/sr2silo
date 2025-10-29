#!/bin/bash
# =============================================================================
# VPipe → Loculus Submit Daily Workflow SLURM Batch Script (Self-Submitting)
# =============================================================================
#
# This script runs the sr2silo workflow daily at 23:00 (11 PM) and automatically
# resubmits itself for the next day.
#
# PREREQUISITES:
# --------------
# - The sr2silo-workflow conda environment must be created first:
#   make create-workflow-env
#
# HOW TO USE:
# -----------
# 1. Edit the configuration variables below as needed
# 2. Submit the job once: sbatch vpipe-loculus-submit.sbatch
# 3. The job will run at 23:00 and automatically resubmit itself
# 4. To stop the cycle, use: scancel <job_id>
#
# =============================================================================

#SBATCH --job-name=vpipe-loculus-submit-daily
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=8G
#SBATCH --tmp=800G
#SBATCH --time=12:00:00
#SBATCH --begin=23:00
#SBATCH --output=/cluster/project/pangolin/research/W-ASAP/logs/vpipe-loculus-submit-daily-%j.out
#SBATCH --error=/cluster/project/pangolin/research/W-ASAP/logs/vpipe-loculus-submit-daily-%j.err

# =============================================================================
# Configuration Variables - Edit these as needed
# =============================================================================

# Number of cores to use (should match --cpus-per-task above)
CORES=20

# Use existing conda environment (true) or workflow-defined conda envs (false)
USE_EXISTING_ENV=true

# Log directory (should match --output and --error paths above)
LOG_DIR="/cluster/project/pangolin/research/W-ASAP/logs"

# Conda environment name
CONDA_ENV="sr2silo-workflow"

# Path to the workflow sbatch script
WORKFLOW_SCRIPT="./sr2silo_workflow.sbatch"

# Path to this script (for resubmission)
DAILY_SCRIPT_NAME="vpipe-loculus-submit-daily.sbatch"

# =============================================================================
# Script Execution
# =============================================================================

echo "=============================================="
echo "VPipe → Loculus Submit Daily Workflow Starting"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Date: $(date)"
echo "Cores: $CORES"
echo "Use existing env: $USE_EXISTING_ENV"
echo "Log directory: $LOG_DIR"
echo "Conda environment: $CONDA_ENV"
echo "=============================================="

# Create log directory if it doesn't exist
mkdir -p "$LOG_DIR"

# Determine the workflow directory.
# Prefer SLURM_SUBMIT_DIR because sbatch copies the script to /var/spool, which lacks the Snakefile.
if [ -n "$SLURM_SUBMIT_DIR" ] && [ -d "$SLURM_SUBMIT_DIR" ]; then
    CANDIDATE_DIR="$SLURM_SUBMIT_DIR"
else
    CANDIDATE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

# If the Snakefile is not directly in the submit directory, fall back to a workflow/ subdirectory.
if [ ! -f "$CANDIDATE_DIR/Snakefile" ] && [ -f "$CANDIDATE_DIR/workflow/Snakefile" ]; then
    CANDIDATE_DIR="$CANDIDATE_DIR/workflow"
fi

WORKFLOW_DIR="$CANDIDATE_DIR"
DAILY_SCRIPT="$WORKFLOW_DIR/$DAILY_SCRIPT_NAME"

echo "Workflow directory: $WORKFLOW_DIR"
cd "$WORKFLOW_DIR" || { echo "ERROR: Failed to cd to $WORKFLOW_DIR"; exit 1; }

# Quick sanity checks
if [ ! -f "Snakefile" ]; then
    echo "ERROR: No Snakefile found in $WORKFLOW_DIR"
    echo "Contents of directory:" && ls -la
    exit 1
fi

# Prefer explicit snakefile path to avoid CWD ambiguity
SNAKEFILE="$WORKFLOW_DIR/Snakefile"

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

# Run the workflow
echo ""
echo "=============================================="
echo "Running sr2silo workflow..."
echo "=============================================="

if [ "$USE_EXISTING_ENV" = true ]; then
    # Use the active conda environment
    echo "Using existing conda environment"
    echo "Running: snakemake -s $SNAKEFILE -c $CORES --rerun-incomplete --keep-going"
    snakemake -s "$SNAKEFILE" -c $CORES --rerun-incomplete --keep-going
else
    # Use the conda environments defined in the workflow
    echo "Using workflow-defined conda environments"
    echo "Running: snakemake -s $SNAKEFILE --use-conda --conda-frontend conda -c $CORES --rerun-incomplete --keep-going"
    snakemake -s "$SNAKEFILE" --use-conda --conda-frontend conda -c $CORES --rerun-incomplete --keep-going
fi

WORKFLOW_EXIT_CODE=$?

echo ""
echo "=============================================="
echo "VPipe → Loculus Submit Workflow Completed"
echo "Exit code: $WORKFLOW_EXIT_CODE"
echo "Completed at: $(date)"
echo "=============================================="

# Resubmit this script for tomorrow at 23:00
echo ""
echo "=============================================="
echo "Resubmitting for tomorrow at 23:00..."
echo "=============================================="

if [ ! -f "$DAILY_SCRIPT" ]; then
    echo "ERROR: Could not locate daily script at $DAILY_SCRIPT"
    echo "The daily workflow cycle has been interrupted."
    exit 1
fi

NEW_JOB_ID=$(sbatch "$DAILY_SCRIPT" | awk '{print $4}')

if [ $? -eq 0 ]; then
    echo "Successfully resubmitted as job ID: $NEW_JOB_ID"
    echo "Next run scheduled for tomorrow at 23:00"
else
    echo "ERROR: Failed to resubmit job!"
    echo "The daily workflow cycle has been interrupted."
    exit 1
fi

echo "=============================================="
echo "Daily workflow cycle continuing..."
echo "=============================================="

exit $WORKFLOW_EXIT_CODE
