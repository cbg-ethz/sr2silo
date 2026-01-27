#!/bin/bash
# =============================================================================
# Historical sr2silo submission - one-off run for specific date range
# Usage: sbatch --export=GPG_PASSPHRASE="your_passphrase" deployments/submit-historical.sbatch
# =============================================================================

#SBATCH --job-name=sr2silo-covid-historical
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=8G
#SBATCH --time=4:00:00
#SBATCH --output=/cluster/project/pangolin/research/W-ASAP/logs/%x-%j.out
#SBATCH --error=/cluster/project/pangolin/research/W-ASAP/logs/%x-%j.err

set -euo pipefail

# Require GPG_PASSPHRASE to be set for credential decryption
: "${GPG_PASSPHRASE:?Set GPG_PASSPHRASE via --export=GPG_PASSPHRASE=your_passphrase}"

PROJECT_ROOT="/cluster/project/pangolin/research/W-ASAP"
CONDA_ENV="sr2silo-workflow"
CORES="${SLURM_CPUS_PER_TASK:-20}"

# Historical date range - edit these as needed
START_DATE="2022-02-05"
END_DATE="2022-02-09"
VIRUS="covid"

echo "=== sr2silo historical run: $VIRUS ($START_DATE to $END_DATE) ==="
echo "Job: $SLURM_JOB_ID | Node: $SLURM_NODELIST | Cores: $CORES | $(date)"

# Setup environment
module load eth_proxy 2>/dev/null || true

# Initialize conda
CONDA_EXE="/cluster/work/bewi/members/koehng/miniconda3/bin/conda"
eval "$("$CONDA_EXE" shell.bash hook)"
conda activate "$CONDA_ENV"

# Load credentials
eval "$(python "$PROJECT_ROOT/sr2silo/deployments/secrets.py")"

# Run workflow with historical date range
cd "$PROJECT_ROOT/sr2silo/workflow"
snakemake --configfile "../deployments/$VIRUS/config.yaml" \
    --config START_DATE="$START_DATE" END_DATE="$END_DATE" \
    -j"$CORES" --rerun-incomplete --keep-going

echo "=== Completed $(date) ==="
