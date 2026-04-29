#!/bin/bash
# =============================================================================
# 0x_find_missing_motifs.sh — Find motifs that didn't produce annotated output
#
# Compares the motif manifest (built by 02_annotate_motifs.sh) against the
# files actually present in OUTDIR. Writes a list of missing motif IDs that
# can be passed to 0x_resubmit_missing_motifs.sh.
#
# A motif is "missing" if:
#   - The manifest lists it AND
#   - ${OUTDIR}/${motif_id}_annotated.tsv.gz does NOT exist
#
# This is the recovery path for tasks that failed due to timeout, OOM, or
# transient errors. Re-run 02_annotate_motifs.sh first if the manifest itself
# is stale.
#
# USAGE:
#   bash 0x_find_missing_motifs.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="."
source "${SCRIPT_DIR}/config.sh"

MANIFEST="${OUTDIR}/_motif_manifest.txt"
MISSING_LIST="${OUTDIR}/_missing_motifs.txt"

if [[ ! -f "${MANIFEST}" ]]; then
    echo "[ERROR] Manifest not found: ${MANIFEST}" >&2
    echo "        Run 02_annotate_motifs.sh first to create it." >&2
    exit 1
fi

> "${MISSING_LIST}"

TOTAL=0
MISSING_COUNT=0
DONE_COUNT=0

while IFS=$'\t' read -r task_id motif_id bed_path; do
    TOTAL=$((TOTAL + 1))
    OUT="${OUTDIR}/${motif_id}_annotated.tsv.gz"
    if [[ -f "${OUT}" ]]; then
        DONE_COUNT=$((DONE_COUNT + 1))
    else
        MISSING_COUNT=$((MISSING_COUNT + 1))
        echo -e "${task_id}\t${motif_id}\t${bed_path}" >> "${MISSING_LIST}"
    fi
done < "${MANIFEST}"

echo "================================================================"
echo "Missing motif scan"
echo "================================================================"
echo "  Total in manifest:   ${TOTAL}"
echo "  Successfully done:   ${DONE_COUNT}"
echo "  Missing:             ${MISSING_COUNT}"
echo "  Missing list:        ${MISSING_LIST}"
echo ""

if (( MISSING_COUNT == 0 )); then
    echo "[OK] All motifs are annotated. Nothing to retry."
    rm -f "${MISSING_LIST}"
    exit 0
fi

echo "First few missing:"
head -10 "${MISSING_LIST}" | cut -f2 | sed 's/^/    /'
if (( MISSING_COUNT > 10 )); then
    echo "    ... and $((MISSING_COUNT - 10)) more"
fi

echo ""
echo "Next: bash 0x_resubmit_missing_motifs.sh --apply"
