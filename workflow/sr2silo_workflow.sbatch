#!/bin/bash
# =============================================================================
# Sr2Silo Workflow SLURM Batch Script
# =============================================================================
#
# This is a standalone .sbatch file for running the sr2silo workflow
# Each Snakemake rule is submitted as its own Slurm job with resources
# mapped from rule resources (threads, mem_mb, tmp_gb).
# =============================================================================

#SBATCH --job-name=sr2silo
#SBATCH --mail-type=END
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=12:00:00
#SBATCH --output=/cluster/project/pangolin/W-ASAP/logs/sr2silo.out
#SBATCH --error=/cluster/project/pangolin/W-ASAP/logs/sr2silo.err

# =============================================================================
# Configuration Variables - Edit these as needed
# =============================================================================

# Number of parallel jobs to submit via Slurm
CORES=20

# Use existing conda environment (true) or workflow-defined conda envs (false)
USE_EXISTING_ENV=false

# Log directory for Snakemake job stdout/stderr
LOG_DIR="/cluster/project/pangolin/W-ASAP/logs"

# Conda environment name (contains snakemake)
CONDA_ENV="sr2silo-workflow"

# =============================================================================
# Script Execution
# =============================================================================

echo "=============================================="
echo "Sr2Silo Workflow Starting"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Jobs (parallel submissions): $CORES"
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

# Cluster submission command mapping rule resources to Slurm
CLUSTER_CMD="sbatch --parsable \
  --job-name=sr2silo_{rule} \
  --cpus-per-task={threads} \
  --mem={resources.mem_mb}M \
  --tmp={resources.tmp_gb}G \
  --output=$LOG_DIR/sr2silo_{rule}.%j.out \
  --error=$LOG_DIR/sr2silo_{rule}.%j.err"

# Prepare and execute Snakemake command
if [ "$USE_EXISTING_ENV" = true ]; then
    echo "Using existing conda environment"
    echo "Running: snakemake -j $CORES --rerun-incomplete --keep-going --restart-times 2 --latency-wait 120 \
      --default-resources mem_mb=12000 tmp_gb=40 \
      --cluster \"$CLUSTER_CMD\""
    snakemake -j $CORES --rerun-incomplete --keep-going --restart-times 2 --latency-wait 120 \
      --default-resources mem_mb=12000 tmp_gb=40 \
      --cluster "$CLUSTER_CMD"
else
    echo "Using workflow-defined conda environments"
    echo "Running: snakemake --use-conda --conda-frontend conda -j $CORES --rerun-incomplete --keep-going --restart-times 2 --latency-wait 120 \
      --default-resources mem_mb=12000 tmp_gb=40 \
      --cluster \"$CLUSTER_CMD\""
    snakemake --use-conda --conda-frontend conda -j $CORES --rerun-incomplete --keep-going --restart-times 2 --latency-wait 120 \
      --default-resources mem_mb=12000 tmp_gb=40 \
      --cluster "$CLUSTER_CMD"
fi

SNAKEMAKE_EXIT_CODE=$?

echo "=============================================="
echo "Sr2Silo Workflow Completed"
echo "Exit code: $SNAKEMAKE_EXIT_CODE"
echo "=============================================="

exit $SNAKEMAKE_EXIT_CODE
