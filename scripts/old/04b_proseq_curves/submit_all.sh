#!/bin/bash
# =============================================================================
# submit_all.sh — Discover all motif archetypes and submit as a SLURM array job
#
# CHANGED from v1:
#   - Removed --merged-only (merged is now the default in run_motif.R)
#   - Added --with-replicates to opt INTO replicate processing
#
# Usage:
#   bash submit_all.sh [--dry-run] [--with-replicates] [--motifs AC0001,AC0002,...] [--max-array N]
#
# Options:
#   --dry-run          Print what would be submitted without actually submitting
#   --with-replicates  Also process individual replicate bigwigs (default: merged only)
#   --motifs LIST      Comma-separated list of specific motif IDs to process
#   --max-array N      Maximum concurrent array tasks [default: 1000]
# =============================================================================

set -euo pipefail

SCRIPT_DIR="."
ANNOT_DIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/annot_motifs"
OUTPUT_BASE="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/proseq_pausing_analysis"

# --- Parse arguments ---
DRY_RUN=false
WITH_REPLICATES=false
SPECIFIC_MOTIFS=""
MAX_ARRAY=1000

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)           DRY_RUN=true;          shift ;;
        --with-replicates)   WITH_REPLICATES=true;   shift ;;
        --motifs)            SPECIFIC_MOTIFS="$2";   shift 2 ;;
        --max-array)         MAX_ARRAY="$2";         shift 2 ;;
        *)                   echo "Unknown option: $1"; exit 1 ;;
    esac
done

# --- Build motif list ---
MOTIF_LIST="${SCRIPT_DIR}/motif_list.txt"

if [[ -n "$SPECIFIC_MOTIFS" ]]; then
    echo "$SPECIFIC_MOTIFS" | tr ',' '\n' > "$MOTIF_LIST"
else
    # Discover all motif IDs from annotated files
    ls "${ANNOT_DIR}"/*_annotated.tsv.gz 2>/dev/null \
        | xargs -n1 basename \
        | sed 's/_annotated\.tsv\.gz$//' \
        | sort -V \
        > "$MOTIF_LIST"
fi

N_MOTIFS=$(wc -l < "$MOTIF_LIST")

if [[ "$N_MOTIFS" -eq 0 ]]; then
    echo "ERROR: No motif files found in ${ANNOT_DIR}"
    exit 1
fi

echo "Found ${N_MOTIFS} motif archetypes to process."
echo "First few: $(head -5 "$MOTIF_LIST" | tr '\n' ' ')"
echo ""

# --- Create directories ---
mkdir -p "${SCRIPT_DIR}/logs"
mkdir -p "${OUTPUT_BASE}"

# --- Set export variable for replicates ---
EXPORT_VARS=""
if $WITH_REPLICATES; then
    echo "Note: Including individual replicates (--with-replicates)"
    EXPORT_VARS="--export=ALL,INCLUDE_REPS=1"
else
    echo "Note: Merged bigwigs only (default)"
fi

# --- Submit ---
SLURM_SCRIPT="${SCRIPT_DIR}/run_motif.sh"

# SBATCH_CMD="sbatch \
#     --array=1-${N_MOTIFS}%${MAX_ARRAY} \
#     --job-name=proseq_pause \
#     --output=${SCRIPT_DIR}/logs/proseq_pause_%A_%a.out \
#     --error=${SCRIPT_DIR}/logs/proseq_pause_%A_%a.err \
#     ${EXPORT_VARS} \
#     ${SLURM_SCRIPT}"

SBATCH_CMD="sbatch \
    --array=1-65 \
    --job-name=proseq_pause \
    --output=${SCRIPT_DIR}/logs/proseq_pause_%A_%a.out \
    --error=${SCRIPT_DIR}/logs/proseq_pause_%A_%a.err \
    ${EXPORT_VARS} \
    ${SLURM_SCRIPT}"

echo "Submission command:"
echo "  ${SBATCH_CMD}"
echo ""

if $DRY_RUN; then
    echo "[DRY RUN] Would submit ${N_MOTIFS} jobs. Motif list saved to: ${MOTIF_LIST}"
else
    eval "$SBATCH_CMD"
    echo ""
    echo "Submitted ${N_MOTIFS} array jobs."
    echo "Monitor with: squeue -u \$USER -n proseq_pause"
    echo "Output will be in: ${OUTPUT_BASE}/<motif_id>/"
fi