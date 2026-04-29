#!/bin/bash
#SBATCH --job-name=moods_postprocess
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/postprocess_%j.out
#SBATCH --error=logs/postprocess_%j.err

# =============================================================================
# 02_postprocess.sbatch — Merge per-chrom BEDs, palindrome dedup, summary
#
# If ${PALINDROME_TABLE} (defined in config.sh) exists, it is passed to
# 02_postprocess.py for palindrome-strand deduplication. To force no dedup,
# delete or rename the file before submitting.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="."
source "${SCRIPT_DIR}/config.sh"

mkdir -p logs "${MERGEDIR}"

if [[ -f "${PALINDROME_TABLE}" ]]; then
    PALINDROME_ARG="--palindrome-table ${PALINDROME_TABLE}"
    echo "Palindrome table found; will dedupe palindromic motifs."
else
    PALINDROME_ARG=""
    echo "No palindrome table; no dedup will be applied."
fi

set +u
source activate "${CONDA_ENV}"
set -u

python "${SCRIPT_DIR}/02_postprocess.py" \
    --input-dir "${CHROMDIR}" \
    --output-dir "${MERGEDIR}" \
    --background "${PREPDIR}/background.txt" \
    --meme "${MEME_FILE}" \
    --pvalue "${PVALUE}" \
    --chrom-sizes "${CHROM_SIZES}" \
    ${PALINDROME_ARG}
