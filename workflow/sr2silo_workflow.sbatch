#!/bin/bash
# =============================================================================
# Sr2Silo Workflow SLURM Batch Script
# =============================================================================
#
# This is a standalone .sbatch file for running the sr2silo workflow
# Each Snakemake rule is submitted as its own Slurm job with resources
# mapped from rule resources (threads, mem_per_cpu_mb, tmp_gb).
# =============================================================================

#SBATCH --job-name=sr2silo
#SBATCH --mail-type=END
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
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

# Create Snakemake profile directory
mkdir -p .snakemake/profiles/slurm

# Create a simple profile config that uses local execution with resource limits
cat > .snakemake/profiles/slurm/config.yaml << 'EOF'
jobs: 20
rerun-incomplete: true
keep-going: true
latency-wait: 120
restart-times: 2
use-conda: true
conda-frontend: conda
default-resources:
  mem_mb: 20000
  mem_per_cpu_mb: 4000
  tmp_gb: 40
  runtime: 120
EOF

# Function to calculate available tmp space and limit concurrent heavy jobs
calculate_tmp_limit() {
    # Get available tmp space in GB (rough estimate)
    if [ -n "$SLURM_TMPDIR" ]; then
        AVAIL_TMP=$(df -BG "$SLURM_TMPDIR" 2>/dev/null | awk 'NR==2{gsub("G","",$4); print $4}' || echo "200")
    else
        AVAIL_TMP=$(df -BG /tmp 2>/dev/null | awk 'NR==2{gsub("G","",$4); print $4}' || echo "200")
    fi

    # Reserve 20GB buffer and ensure we don't go below 40GB per heavy job
    USABLE_TMP=$((AVAIL_TMP - 20))
    MAX_HEAVY_JOBS=$((USABLE_TMP / 40))

    # Ensure at least 1 heavy job can run, max 10 to avoid over-submission
    if [ $MAX_HEAVY_JOBS -lt 1 ]; then MAX_HEAVY_JOBS=1; fi
    if [ $MAX_HEAVY_JOBS -gt 10 ]; then MAX_HEAVY_JOBS=10; fi

    echo $MAX_HEAVY_JOBS
}

MAX_TMP_JOBS=$(calculate_tmp_limit)
echo "Available tmp space allows max $MAX_TMP_JOBS concurrent heavy jobs"

# Prepare and execute Snakemake command with resource constraints
if [ "$USE_EXISTING_ENV" = true ]; then
    echo "Using existing conda environment"
    echo "Running: snakemake --cores $CORES --rerun-incomplete --keep-going \
      --resources tmp_gb=$((MAX_TMP_JOBS * 40)) \
      --default-resources mem_mb=10000 mem_per_cpu_mb=2500 tmp_gb=40 runtime=120"
    snakemake --cores $CORES --rerun-incomplete --keep-going \
      --resources tmp_gb=$((MAX_TMP_JOBS * 40)) \
      --default-resources mem_mb=10000 mem_per_cpu_mb=2500 tmp_gb=40 runtime=120
else
    echo "Using workflow-defined conda environments"
    echo "Running: snakemake --use-conda --conda-frontend conda --cores $CORES --rerun-incomplete --keep-going \
      --resources tmp_gb=$((MAX_TMP_JOBS * 40)) \
      --default-resources mem_mb=10000 mem_per_cpu_mb=2500 tmp_gb=40 runtime=120"
    snakemake --use-conda --conda-frontend conda --cores $CORES --rerun-incomplete --keep-going \
      --resources tmp_gb=$((MAX_TMP_JOBS * 40)) \
      --default-resources mem_mb=10000 mem_per_cpu_mb=2500 tmp_gb=40 runtime=120
fi

SNAKEMAKE_EXIT_CODE=$?

echo "=============================================="
echo "Sr2Silo Workflow Completed"
echo "Exit code: $SNAKEMAKE_EXIT_CODE"
echo "=============================================="

exit $SNAKEMAKE_EXIT_CODE
