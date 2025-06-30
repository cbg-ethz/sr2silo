#!/bin/bash
# =============================================================================
# Sr2Silo Workflow .sbatch File Generator
# =============================================================================
#
# This script generates customized .sbatch files for the sr2silo workflow
#
# HOW TO USE:
# -----------
# Basic usage:
#   ./generate_sbatch.sh
#       - Creates sr2silo_10cores.sbatch with 10 cores using workflow-defined conda environments
#
# With options:
#   ./generate_sbatch.sh -c 8
#       - Creates sr2silo_8cores.sbatch with 8 cores
#
#   ./generate_sbatch.sh -e
#       - Creates sr2silo_existing_env.sbatch using your currently active conda environment
#
#   ./generate_sbatch.sh -c 12 -e -o my_custom_job.sbatch
#       - Creates my_custom_job.sbatch with 12 cores using existing environment
#
# =============================================================================

# Default values
CORES=10
USE_EXISTING_ENV=false
LOG_DIR="/cluster/work/bewi/members/koehng/logs"
OUTPUT_FILE=""

function show_usage {
  echo "Usage: $0 [-e] [-c CORES] [-l LOG_DIR] [-o OUTPUT_FILE]"
  echo "  -e            Use existing conda environment (default: use workflow-defined conda envs)"
  echo "  -c CORES      Number of cores to use (default: 10)"
  echo "  -l LOG_DIR    Directory for log files (default: /cluster/work/bewi/members/koehng/logs)"
  echo "  -o OUTPUT     Output .sbatch filename (default: auto-generated based on options)"
  exit 1
}

while getopts "ec:l:o:" opt; do
  case $opt in
    e) USE_EXISTING_ENV=true ;;
    c) CORES=$OPTARG ;;
    l) LOG_DIR=$OPTARG ;;
    o) OUTPUT_FILE=$OPTARG ;;
    *) show_usage ;;
  esac
done

# Generate default output filename if not specified
if [ -z "$OUTPUT_FILE" ]; then
  if [ "$USE_EXISTING_ENV" = true ]; then
    OUTPUT_FILE="sr2silo_${CORES}cores_existing_env.sbatch"
  else
    OUTPUT_FILE="sr2silo_${CORES}cores.sbatch"
  fi
fi

# Check if template exists
TEMPLATE_FILE="sr2silo_workflow_template.sbatch"
if [ ! -f "$TEMPLATE_FILE" ]; then
  echo "Error: Template file '$TEMPLATE_FILE' not found!"
  echo "Make sure you're running this script from the workflow directory."
  exit 1
fi

echo "Generating .sbatch file: $OUTPUT_FILE"
echo "  Cores: $CORES"
echo "  Use existing env: $USE_EXISTING_ENV"
echo "  Log directory: $LOG_DIR"

# Create the customized .sbatch file
sed -e "s|CORES_PLACEHOLDER|$CORES|g" \
    -e "s|USE_ENV_PLACEHOLDER|$USE_EXISTING_ENV|g" \
    -e "s|LOG_DIR_PLACEHOLDER|$LOG_DIR|g" \
    "$TEMPLATE_FILE" > "$OUTPUT_FILE"

# Make the generated file executable
chmod +x "$OUTPUT_FILE"

echo ""
echo "✓ Generated: $OUTPUT_FILE"
echo ""
echo "To submit the job, run:"
echo "  sbatch $OUTPUT_FILE"
echo ""
echo "To check job status:"
echo "  squeue -u \$USER"
echo ""
echo "To check job output:"
echo "  tail -f $LOG_DIR/sr2silo.out"
echo "  tail -f $LOG_DIR/sr2silo.err"
