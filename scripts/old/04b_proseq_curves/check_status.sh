#!/bin/bash
# =============================================================================
# check_status.sh — Check run status across all motif archetypes
#
# Categorizes each motif as:
#   COMPLETED     — Has COMPLETED.txt (plots were generated)
#   ALL_SKIPPED   — Has ALL_SKIPPED.txt (all presets had too few instances)
#   PARTIAL_SKIP  — Has SKIPPED_PRESETS.txt but also COMPLETED.txt
#   NO_OUTPUT     — Output directory exists but no status files (likely failed)
#   NOT_RUN       — No output directory at all
#
# Usage:
#   bash check_status.sh [output_base_dir] [motif_list.txt]
# =============================================================================

set -euo pipefail

OUTPUT_BASE="${1:-/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/proseq_pausing_analysis}"
MOTIF_LIST="${2:-motif_list.txt}"

if [[ ! -f "$MOTIF_LIST" ]]; then
    echo "ERROR: Motif list not found: $MOTIF_LIST"
    exit 1
fi

echo "============================================"
echo "Motif PRO-seq run status check"
echo "Output base: $OUTPUT_BASE"
echo "Motif list:  $MOTIF_LIST"
echo "============================================"
echo ""

n_completed=0
n_all_skipped=0
n_partial=0
n_no_output=0
n_not_run=0

while IFS= read -r motif_id; do
    dir="${OUTPUT_BASE}/${motif_id}"

    if [[ ! -d "$dir" ]]; then
        echo "NOT_RUN       ${motif_id}"
        ((n_not_run++)) || true
    elif [[ -f "${dir}/ALL_SKIPPED.txt" ]]; then
        echo "ALL_SKIPPED   ${motif_id}  — $(grep 'MIN_INSTANCES' "${dir}/ALL_SKIPPED.txt" | head -1)"
        ((n_all_skipped++)) || true
    elif [[ -f "${dir}/COMPLETED.txt" && -f "${dir}/SKIPPED_PRESETS.txt" ]]; then
        echo "PARTIAL_SKIP  ${motif_id}  — $(grep 'Plots generated' "${dir}/COMPLETED.txt")"
        ((n_partial++)) || true
    elif [[ -f "${dir}/COMPLETED.txt" ]]; then
        echo "COMPLETED     ${motif_id}  — $(grep 'Plots generated' "${dir}/COMPLETED.txt")"
        ((n_completed++)) || true
    else
        echo "NO_OUTPUT     ${motif_id}  — directory exists but no status file (check logs)"
        ((n_no_output++)) || true
    fi
done < "$MOTIF_LIST"

echo ""
echo "============================================"
echo "Summary"
echo "============================================"
echo "  COMPLETED:     ${n_completed}"
echo "  ALL_SKIPPED:   ${n_all_skipped}"
echo "  PARTIAL_SKIP:  ${n_partial}"
echo "  NO_OUTPUT:     ${n_no_output}  (likely script failures — check SLURM logs)"
echo "  NOT_RUN:       ${n_not_run}"
echo "  Total:         $(wc -l < "$MOTIF_LIST")"
