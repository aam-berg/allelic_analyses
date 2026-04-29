#!/bin/bash
#SBATCH --job-name=proseq_pause_%j
#SBATCH --output=logs/proseq_pause_%A_%a.out
#SBATCH --error=logs/proseq_pause_%A_%a.err
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

# =============================================================================
# run_motif.sh — SLURM job script for a single motif (v2)
#
# CHANGED from v1:
#   - Default: merged bigwigs only (no --skip_replicates needed)
#   - Use --include_replicates to also process individual replicate bigwigs
#
# Usage (direct):
#   sbatch run_motif.sh AC0001
#   sbatch --export=INCLUDE_REPS=1 run_motif.sh AC0001  # with replicates
#
# Usage (array, called by submit_all.sh):
#   Reads motif ID from motif_list.txt based on SLURM_ARRAY_TASK_ID.
# =============================================================================

set -euo pipefail

# --- Load modules (adjust to your cluster) ---
source activate test_m

# --- Determine script directory ---
SCRIPT_DIR="."

# --- Get motif ID ---
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    # Array job mode: read motif ID from the motif list file
    MOTIF_LIST="${SCRIPT_DIR}/motif_list.txt"
    MOTIF_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MOTIF_LIST")
else
    # Direct submission mode
    MOTIF_ID="${1:?Usage: sbatch run_motif.sh <motif_id>}"
fi

echo "============================================"
echo "Processing motif: ${MOTIF_ID}"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Node: $(hostname)"
echo "Time: $(date)"
echo "============================================"

# --- Build extra arguments ---
# Default: merged only (no flag needed; run_motif.R defaults to merged only).
# Set INCLUDE_REPS=1 to also run individual replicates.
EXTRA_ARGS=""
if [[ "${INCLUDE_REPS:-0}" == "1" ]]; then
    EXTRA_ARGS="--include_replicates"
    echo "Mode: including individual replicates"
else
    echo "Mode: merged bigwigs only (default)"
fi

# --- Run ---
Rscript "${SCRIPT_DIR}/run_motif.R" \
    --motif_id "${MOTIF_ID}" \
    --config "${SCRIPT_DIR}/config.R" \
    ${EXTRA_ARGS}

echo ""
echo "Done: $(date)"