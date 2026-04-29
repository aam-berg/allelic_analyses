#!/bin/bash
#SBATCH --job-name=moods_qc
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/moods_qc_%j.out
#SBATCH --error=logs/moods_qc_%j.err

# ============================================================================
# 03_qc_report.sh — Comprehensive QC for MOODS genome-wide scan
#
# This runs the full QC analysis including:
#   - Genome coverage (all motifs + per-motif)
#   - Genomic context classification (promoter/intragenic/intergenic)
#   - TSS distance distributions
#   - Information content analysis
#   - Palindrome detection
#   - Publication-quality plots
#
# Requirements:
#   - bedtools in PATH
#   - Python 3.8+ with numpy, matplotlib, pandas
#   - MOODS-python (for IC/p-value estimation — optional)
#
# Submit with:
#   sbatch 03_qc_report.sbatch
# ============================================================================

set -euo pipefail

# ---- Paths (edit these) ----
BASE="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/moods_scan"
MOTIF_DIR="${BASE}/merged/per_motif"
MERGED_BED="${BASE}/merged/all_motifs_merged.bed.gz"
SUMMARY_TSV="${BASE}/merged/motif_hit_summary.tsv"
GTF="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/annotation/gencode.vM38.chr_patch_hapl_scaff.annotation.gtf.gz"
CHROM_SIZES="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/genome/mm39.chrom.sizes"
MEME="/home/alb1273/pausing_phase_project/resources/consensus_pwms.meme"
BACKGROUND="${BASE}/prep/background.txt"
PVALUE="1e-5"
OUTPUT_DIR="${BASE}/qc_report"

QC_SCRIPT="03_qc_report.py"

# ---- Setup ----
mkdir -p logs "${OUTPUT_DIR}"

echo "=========================================="
echo "MOODS QC Report"
echo "Started: $(date)"
echo "=========================================="

# Activate environment (adjust as needed)
source activate motif_scanning

# Ensure bedtools is available
if ! command -v bedtools &> /dev/null; then
    echo "ERROR: bedtools not found in PATH" >&2
    # Try loading module
    module load bedtools/2.31.0 2>/dev/null || true
    if ! command -v bedtools &> /dev/null; then
        echo "ERROR: Still can't find bedtools. Please install or load module." >&2
        exit 1
    fi
fi
echo "bedtools: $(which bedtools)"
echo "Python:   $(which python)"

# ---- Run QC ----
python "${QC_SCRIPT}" \
    --motif-dir "${MOTIF_DIR}" \
    --merged-bed "${MERGED_BED}" \
    --summary-tsv "${SUMMARY_TSV}" \
    --gtf "${GTF}" \
    --chrom-sizes "${CHROM_SIZES}" \
    --meme "${MEME}" \
    --background "${BACKGROUND}" \
    --pvalue "${PVALUE}" \
    --output-dir "${OUTPUT_DIR}" \
    --tss-sample-n 50000 \
    --tss-max-dist 100000

echo ""
echo "=========================================="
echo "Done: $(date)"
echo "Output: ${OUTPUT_DIR}"
echo "=========================================="
