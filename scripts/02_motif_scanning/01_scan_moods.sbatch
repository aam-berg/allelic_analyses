#!/bin/bash
#SBATCH --job-name=moods_scan
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/moods_scan_%A_%a.out
#SBATCH --error=logs/moods_scan_%A_%a.err

# =============================================================================
# 01_scan_moods.sbatch — One MOODS scan per chromosome (array task)
#
# Submitted by run_pipeline.sh, but can also be submitted directly:
#   N=$(wc -l < ${PREPDIR}/chromosomes.txt)
#   sbatch --array=1-$N 01_scan_moods.sbatch
# =============================================================================

set -euo pipefail

SCRIPT_DIR="."
source "${SCRIPT_DIR}/config.sh"

mkdir -p logs "${CHROMDIR}"

CHROM_FILE="${PREPDIR}/chromosomes.txt"
if [[ ! -f "${CHROM_FILE}" ]]; then
    echo "ERROR: Chromosome file not found: ${CHROM_FILE}" >&2
    echo "       Run 00_prepare.py first." >&2
    exit 1
fi

CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CHROM_FILE}")
if [[ -z "${CHROM}" ]]; then
    echo "ERROR: No chromosome for task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

echo "MOODS scan: ${CHROM} (task ${SLURM_ARRAY_TASK_ID})"
echo "  Genome:   ${GENOME}"
echo "  MEME:     ${MEME_FILE}"
echo "  P-value:  ${PVALUE}"
echo "  Output:   ${CHROMDIR}/${CHROM}.bed"

set +u
source activate "${CONDA_ENV}"
set -u

python "${SCRIPT_DIR}/01_scan_moods.py" \
    --genome "${GENOME}" \
    --chrom "${CHROM}" \
    --meme "${MEME_FILE}" \
    --background "${PREPDIR}/background.txt" \
    --pvalue "${PVALUE}" \
    --output "${CHROMDIR}/${CHROM}.bed"

echo "Done: ${CHROM}"
