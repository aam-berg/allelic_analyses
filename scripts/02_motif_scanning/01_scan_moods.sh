#!/bin/bash
#SBATCH --job-name=moods_scan
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/moods_scan_%A_%a.out
#SBATCH --error=logs/moods_scan_%A_%a.err

# ============================================================================
# 01_scan_moods.sbatch — SLURM array job: one task per chromosome
#
# Submit with:
#   N=$(wc -l < ${PREPDIR}/chromosomes.txt)
#   sbatch --array=1-${N} 01_scan_moods.sbatch
#
# Environment variables (set in run_pipeline.sh or export before submitting):
#   GENOME       — path to mm39.fa
#   MEME_FILE    — path to consensus_pwms.meme
#   PREPDIR      — directory containing background.txt and chromosomes.txt
#   OUTDIR       — output directory for BED files
#   PVALUE       — p-value threshold (default: 1e-5)
#   SCAN_SCRIPT  — path to 01_scan_moods.py
# ============================================================================

#set -euo pipefail

# Defaults (override by exporting before sbatch)
GENOME="${GENOME:-/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/genome/mm39.fa}"
MEME_FILE="${MEME_FILE:-/home/alb1273/pausing_phase_project/resources/consensus_pwms.meme}"
PREPDIR="${PREPDIR:-/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/moods_scan/prep}"
OUTDIR="${OUTDIR:-/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/moods_scan/per_chrom}"
PVALUE="${PVALUE:-1e-5}"
SCAN_SCRIPT="01_scan_moods.py"

# Get chromosome for this array task
CHROM_FILE="${PREPDIR}/chromosomes.txt"
if [[ ! -f "${CHROM_FILE}" ]]; then
    echo "ERROR: Chromosome file not found: ${CHROM_FILE}" >&2
    exit 1
fi

CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CHROM_FILE}")
if [[ -z "${CHROM}" ]]; then
    echo "ERROR: No chromosome for task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

echo "=========================================="
echo "MOODS scan: ${CHROM}"
echo "Task ID:    ${SLURM_ARRAY_TASK_ID}"
echo "Genome:     ${GENOME}"
echo "MEME file:  ${MEME_FILE}"
echo "P-value:    ${PVALUE}"
echo "Output dir: ${OUTDIR}"
echo "=========================================="

mkdir -p "${OUTDIR}"

# Activate your conda environment (adjust as needed)
# module load conda2/4.2.13   # or whatever module system O2 uses
source activate motif_scanning    # your environment with MOODS-python installed

# If using a module-based Python with MOODS installed:
# module load python/3.10.11

# Run the scan
python "${SCAN_SCRIPT}" \
    --genome "${GENOME}" \
    --chrom "${CHROM}" \
    --meme "${MEME_FILE}" \
    --background "${PREPDIR}/background.txt" \
    --pvalue "${PVALUE}" \
    --output "${OUTDIR}/${CHROM}.bed"

echo "Done: ${CHROM}"
