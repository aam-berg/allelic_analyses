#!/bin/bash
# =============================================================================
# find_failed_motifs.sh — Identify motifs that failed due to technical issues
#
# Distinguishes between:
#   OK           — COMPLETED.txt exists (plots generated successfully)
#   SKIPPED      — ALL_SKIPPED.txt exists (too few instances, not a bug)
#   PARTIAL      — COMPLETED.txt + SKIPPED_PRESETS.txt (some presets OK, some skipped)
#   FAILED       — Output dir exists but no COMPLETED/ALL_SKIPPED (technical failure)
#   FAILED_EMPTY — Output dir exists but data/ is empty (crashed during extraction)
#   NOT_RUN      — No output dir at all (never submitted or SLURM didn't start)
#
# Output:
#   failed_motifs.txt          — motif IDs to re-run (FAILED + FAILED_EMPTY)
#   not_run_motifs.txt         — motif IDs never attempted
#   run_status_summary.tsv     — full status table for all motifs
#
# Usage:
#   bash find_failed_motifs.sh [motif_list.txt] [output_base_dir]
# =============================================================================

set -euo pipefail

MOTIF_LIST="${1:-motif_list.txt}"
OUTPUT_BASE="${2:-/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/proseq_pausing_analysis}"
LOG_DIR="${3:-./logs}"

if [[ ! -f "$MOTIF_LIST" ]]; then
    echo "ERROR: Motif list not found: $MOTIF_LIST"
    echo "Usage: bash find_failed_motifs.sh [motif_list.txt] [output_base_dir] [log_dir]"
    exit 1
fi

FAILED_FILE="failed_motifs.txt"
NOT_RUN_FILE="not_run_motifs.txt"
SUMMARY_FILE="run_status_summary.tsv"

> "$FAILED_FILE"
> "$NOT_RUN_FILE"

# Header for summary
echo -e "motif_id\tstatus\tn_profiles\tn_plots\terror_hint" > "$SUMMARY_FILE"

n_ok=0; n_skipped=0; n_partial=0; n_failed=0; n_failed_empty=0; n_not_run=0

while IFS= read -r motif_id; do
    [[ -z "$motif_id" || "$motif_id" == \#* ]] && continue

    dir="${OUTPUT_BASE}/${motif_id}"

    if [[ ! -d "$dir" ]]; then
        # --- NOT_RUN ---
        echo "$motif_id" >> "$NOT_RUN_FILE"
        echo -e "${motif_id}\tNOT_RUN\t0\t0\tno output directory" >> "$SUMMARY_FILE"
        ((n_not_run++)) || true
        continue
    fi

    has_completed=false
    has_all_skipped=false
    has_skipped_presets=false
    [[ -f "${dir}/COMPLETED.txt" ]]        && has_completed=true
    [[ -f "${dir}/ALL_SKIPPED.txt" ]]      && has_all_skipped=true
    [[ -f "${dir}/SKIPPED_PRESETS.txt" ]]  && has_skipped_presets=true

    # Count actual outputs
    n_profiles=$(find "${dir}/data" -name "*_profiles.tsv.gz" 2>/dev/null | wc -l)
    n_plots=$(find "${dir}/plots" -name "*.pdf" 2>/dev/null | wc -l)

    if $has_completed && ! $has_skipped_presets; then
        # --- OK ---
        echo -e "${motif_id}\tOK\t${n_profiles}\t${n_plots}\t" >> "$SUMMARY_FILE"
        ((n_ok++)) || true

    elif $has_completed && $has_skipped_presets; then
        # --- PARTIAL ---
        echo -e "${motif_id}\tPARTIAL\t${n_profiles}\t${n_plots}\tsome presets skipped (low n)" >> "$SUMMARY_FILE"
        ((n_partial++)) || true

    elif $has_all_skipped; then
        # --- SKIPPED (not a failure) ---
        echo -e "${motif_id}\tSKIPPED\t${n_profiles}\t${n_plots}\tall presets below MIN_INSTANCES" >> "$SUMMARY_FILE"
        ((n_skipped++)) || true

    else
        # --- FAILED (technical issue) ---
        # Try to extract error hint from SLURM logs
        error_hint=""

        # Look for the most recent log file mentioning this motif
        log_file=""
        if [[ -d "$LOG_DIR" ]]; then
            log_file=$(grep -rl "Processing motif: ${motif_id}" "${LOG_DIR}"/*.err 2>/dev/null | tail -1 || true)
            if [[ -z "$log_file" ]]; then
                log_file=$(grep -rl "Processing motif: ${motif_id}" "${LOG_DIR}"/*.out 2>/dev/null | tail -1 || true)
            fi
        fi

        if [[ -n "$log_file" ]]; then
            # Extract last error-like lines
            error_hint=$(grep -i "error\|fatal\|cannot\|failed\|traceback\|killed\|segfault\|oom\|memory" \
                         "$log_file" 2>/dev/null | tail -3 | tr '\n' ' | ' | head -c 200 || true)
        fi

        # Distinguish between empty dir (crashed early) and partial data (crashed during plotting)
        if [[ "$n_profiles" -eq 0 ]]; then
            echo "$motif_id" >> "$FAILED_FILE"
            echo -e "${motif_id}\tFAILED_EMPTY\t${n_profiles}\t${n_plots}\t${error_hint:-no profiles generated}" >> "$SUMMARY_FILE"
            ((n_failed_empty++)) || true
        else
            echo "$motif_id" >> "$FAILED_FILE"
            echo -e "${motif_id}\tFAILED\t${n_profiles}\t${n_plots}\t${error_hint:-profiles exist but no status file}" >> "$SUMMARY_FILE"
            ((n_failed++)) || true
        fi
    fi

done < "$MOTIF_LIST"

# --- Summary ---
echo ""
echo "============================================"
echo "Run status summary"
echo "============================================"
echo "  OK:            ${n_ok}"
echo "  PARTIAL:       ${n_partial}  (some presets skipped for low n — not failures)"
echo "  SKIPPED:       ${n_skipped}  (all presets below MIN_INSTANCES — not failures)"
echo "  FAILED:        ${n_failed}  (data extracted but crashed before completion)"
echo "  FAILED_EMPTY:  ${n_failed_empty}  (crashed before/during data extraction)"
echo "  NOT_RUN:       ${n_not_run}  (never started)"
echo "  ---"
echo "  Total motifs:  $(wc -l < "$MOTIF_LIST")"
echo ""

n_rerun=$(wc -l < "$FAILED_FILE")
if [[ "$n_rerun" -gt 0 ]]; then
    echo "Motifs to re-run: ${n_rerun}"
    echo "  List:  ${FAILED_FILE}"
    echo ""
    echo "To re-submit failed motifs:"
    echo "  bash submit_all.sh --motifs \$(paste -sd, ${FAILED_FILE})"
    echo ""
    echo "Or manually:"
    echo "  while read m; do ./run_motif.sh \$m; done < ${FAILED_FILE}"
else
    echo "No technical failures detected."
fi

n_notrun=$(wc -l < "$NOT_RUN_FILE")
if [[ "$n_notrun" -gt 0 ]]; then
    echo ""
    echo "Motifs never attempted: ${n_notrun}"
    echo "  List: ${NOT_RUN_FILE}"
fi

echo ""
echo "Full status table: ${SUMMARY_FILE}"