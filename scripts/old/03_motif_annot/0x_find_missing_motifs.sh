#!/bin/bash
# =============================================================================
# 0x_find_missing_motifs.sh — Identify motifs with missing/empty annotated output
#
# Compares input BED files in MOTIF_DIR (config.sh) against
# <MOTIF>_annotated.tsv.gz files in OUTDIR. Missing or empty outputs are
# written to ${OUTDIR}/missing_motifs.txt for later resubmission via
# 0x_resubmit_missing_motifs.sh.
#
# Idempotent and side-effect-free except for that one file.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

MISSING_LIST="${OUTDIR}/missing_motifs.txt"

total=0
missing=0
> "${MISSING_LIST}"  # truncate

for bed in "${MOTIF_DIR}"/*.bed; do
    motif_id=$(basename "${bed}" .bed)
    out_file="${OUTDIR}/${motif_id}_annotated.tsv.gz"
    total=$((total + 1))
    if [[ ! -f "${out_file}" || ! -s "${out_file}" ]]; then
        echo "${bed}" >> "${MISSING_LIST}"
        missing=$((missing + 1))
    fi
done

echo "=============================================="
echo "Missing motif check"
echo "  Total input BEDs:   ${total}"
echo "  Successful outputs: $((total - missing))"
echo "  Missing/empty:      ${missing}"
echo "=============================================="

if [[ ${missing} -gt 0 ]]; then
    echo ""
    echo "List written to: ${MISSING_LIST}"
    echo ""
    echo "First 10 missing:"
    head -10 "${MISSING_LIST}" | while read -r f; do
        echo "  $(basename "$f")"
    done
    [[ ${missing} -gt 10 ]] && echo "  ... and $((missing - 10)) more"
    echo ""
    echo "Resubmit with: bash 0x_resubmit_missing_motifs.sh"
else
    echo ""
    echo "All motifs completed successfully!"
fi
