#!/bin/bash
# =============================================================================
# Daily sr2silo submission - pass virus name as argument
# Usage: sbatch --job-name=sr2silo-covid --export=VIRUS=covid deployments/submit-daily.sbatch
# =============================================================================

#SBATCH --job-name=sr2silo-daily
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=8G
#SBATCH --tmp=800G
#SBATCH --time=12:00:00
#SBATCH --begin=23:00
#SBATCH --output=/cluster/project/pangolin/research/W-ASAP/logs/%x-%j.out
#SBATCH --error=/cluster/project/pangolin/research/W-ASAP/logs/%x-%j.err

set -euo pipefail

# Require VIRUS to be set
: "${VIRUS:?Set VIRUS via --export=VIRUS=covid}"

PROJECT_ROOT="/cluster/project/pangolin/research/W-ASAP"
CONDA_ENV="base"
CORES="${SLURM_CPUS_PER_TASK:-20}"

echo "=== sr2silo daily: $VIRUS (job $SLURM_JOB_ID) ==="
echo "Node: $SLURM_NODELIST | Cores: $CORES | $(date)"

# Setup environment
module load eth_proxy 2>/dev/null || true

# Initialize conda (using hook to preserve system PATH including sbatch)
CONDA_EXE="/cluster/project/pangolin/resources/miniconda3/bin/conda"
eval "$("$CONDA_EXE" shell.bash hook)"
conda activate "$CONDA_ENV"

# Load credentials
eval "$(python "$PROJECT_ROOT/sr2silo/deployments/secrets.py")"

# Extract resource requirements from config for next run
CONFIG_FILE="$PROJECT_ROOT/sr2silo/deployments/$VIRUS/config.yaml"
NEXT_CPUS=$(grep "^SLURM_CPUS:" "$CONFIG_FILE" | awk '{print $2}' | tr -d '"' || echo "20")
NEXT_MEM=$(grep "^SLURM_MEM:" "$CONFIG_FILE" | awk '{print $2}' | tr -d '"' || echo "160G")

# Calculate mem-per-cpu for sbatch (total mem / cpus)
# Handle G suffix
if [[ "$NEXT_MEM" == *"G" ]]; then
    TOTAL_MEM_GB=${NEXT_MEM%G}
    MEM_PER_CPU_GB=$(( TOTAL_MEM_GB / NEXT_CPUS ))
    MEM_PER_CPU="${MEM_PER_CPU_GB}G"
else
    # Assume MB if no suffix, or handle appropriately. For now assume G if simple number
    # Fallback to safe default if parsing fails
    MEM_PER_CPU="8G"
fi

echo "Next run resources: CPUS=$NEXT_CPUS, MEM=$NEXT_MEM (MEM_PER_CPU=$MEM_PER_CPU)"

# Run workflow and capture exit code
cd "$PROJECT_ROOT/sr2silo/workflow"
set +e  # Temporarily disable exit on error
snakemake --configfile "../deployments/$VIRUS/config.yaml" -j"$CORES" --rerun-incomplete --keep-going --rerun-trigger mtime --conda-frontend conda --conda-prefix "/cluster/project/pangolin/resources/snake-envs" --use-conda
SNAKEMAKE_EXIT=$?
set -e  # Re-enable exit on error

echo "=== Completed $(date) with exit code $SNAKEMAKE_EXIT ==="

# Decide whether to resubmit based on exit code:
# - 0: Success, resubmit
# - 1: Some jobs failed but workflow ran (with --keep-going), resubmit
# - Other: Critical failure (config error, missing files, etc.), do NOT resubmit
if [[ $SNAKEMAKE_EXIT -gt 1 ]]; then
    echo "ERROR: Snakemake exited with code $SNAKEMAKE_EXIT (critical failure). NOT resubmitting."
    exit $SNAKEMAKE_EXIT
fi

# Resubmit for tomorrow
/cluster/apps/slurm/bin/sbatch \
    --job-name="sr2silo-$VIRUS" \
    --cpus-per-task="$NEXT_CPUS" \
    --mem-per-cpu="$MEM_PER_CPU" \
    --export=VIRUS="$VIRUS",GPG_PASSPHRASE="$GPG_PASSPHRASE" \
    "$PROJECT_ROOT/sr2silo/deployments/submit-daily.sbatch"
