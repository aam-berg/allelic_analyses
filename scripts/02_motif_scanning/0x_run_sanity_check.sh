#!/bin/bash
# =============================================================================
# 0x_run_sanity_check.sh — Pre-flight check on the motif PWM set.
#
# Reads MEME_FILE + METADATA_FILE (from config.sh), produces:
#   ${SANITYDIR}/motif_width_distribution.png
#   ${SANITYDIR}/motif_orientation_analysis.png
#   ${SANITYDIR}/orientation_analysis.tsv     ← used by 02_postprocess.py
#
# Cheap (~1 min). Login-node OK; no SLURM needed.
# Run before run_pipeline.sh to enable palindrome dedup at step 2.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

mkdir -p "${SANITYDIR}"

set +u
source activate "${CONDA_ENV}"
set -u

python "${SCRIPT_DIR}/0x_sanity_check_motifs.py" \
    --meme "${MEME_FILE}" \
    --metadata "${METADATA_FILE}" \
    --output-dir "${SANITYDIR}" \
    --width-cutoff 8

echo "Sanity check done. Output: ${SANITYDIR}"
echo "  - motif_width_distribution.png"
echo "  - motif_orientation_analysis.png"
echo "  - orientation_analysis.tsv  (palindrome table for step 2)"
