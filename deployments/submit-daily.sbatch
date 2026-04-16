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
# preset default (to keep unbound detection happy)
declare -x ftp_proxy="http://proxy.service.consul:3128"
declare -x http_proxy="http://proxy.service.consul:3128"
declare -x https_proxy="http://proxy.service.consul:3128"
declare -x no_proxy="api.wandb.ai,app.neptune.ai,jupyter.eu-dev.hpc.ethz.ch,jupyter-staging.euler.hpc.ethz.ch,jupyter.euler.hpc.ethz.ch,.consul,localhost,127.0.0.1,127.0.0.0/8,169.254.0.0/16"
# overwrite the above with actual module's setting
source /etc/profile.d/lmod.sh || true
source /etc/profile.d/module_path.sh || true
module load eth_proxy || true
#module load eth_proxy 2>/dev/null || true
echo "https_proxy: ${https_proxy}"

# Initialize conda (using hook to preserve system PATH including sbatch)
CONDA_EXE="/cluster/project/pangolin/resources/miniconda3/bin/conda"
eval "$("$CONDA_EXE" shell.bash hook)"
conda activate "$CONDA_ENV"

# Load credentials
eval "$(python "$PROJECT_ROOT/sr2silo/deployments/secrets.py")"

# Extract resource requirements from config for next run
CONFIG_FILE="$PROJECT_ROOT/sr2silo/deployments/$VIRUS/config.yaml"
# HACK parse configs with Gawk. Should switch to 'yq'
NEXT_CPUS=$(gawk '/^SLURM_CPUS:/ {print gensub(/"/, "", "g", $2)}' "$CONFIG_FILE");	NEXT_CPUS=${NEXT_CPUS:-20}
NEXT_MEM=$(gawk '/^SLURM_MEM:/ {print gensub(/"/, "", "g", $2)}' "$CONFIG_FILE");	NEXT_MEM=${NEXT_MEM:-160G}

# Calculate mem-per-cpu for sbatch (total mem / cpus)
# Handle G suffix
if [[ "$NEXT_MEM" =~ ^([[:digit:]]+)([MG]?)$ ]]; then
    # If no suffix, handle appropriately: skip the suffix too (i.e.: "?")
    TOTAL_MEM_GB=${BASH_REMATCH[1]}
    MEM_PER_CPU_GB=$(( TOTAL_MEM_GB / NEXT_CPUS ))
    MEM_PER_CPU="${MEM_PER_CPU_GB}${BASH_REMATCH[2]}"
else
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
