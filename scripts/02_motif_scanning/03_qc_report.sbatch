#!/bin/bash
#SBATCH --job-name=moods_qc
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/moods_qc_%j.out
#SBATCH --error=logs/moods_qc_%j.err

# =============================================================================
# 03_qc_report.sbatch — Comprehensive QC report
# =============================================================================

set -euo pipefail

SCRIPT_DIR="."
source "${SCRIPT_DIR}/config.sh"

mkdir -p logs "${QCDIR}"

# bedtools sanity check
if ! command -v bedtools &>/dev/null; then
    module load bedtools/2.31.0 2>/dev/null || true
    if ! command -v bedtools &>/dev/null; then
        echo "ERROR: bedtools not found; install or load module." >&2
        exit 1
    fi
fi

set +u
source activate "${CONDA_ENV}"
set -u

python "${SCRIPT_DIR}/03_qc_report.py" \
    --motif-dir "${MERGEDIR}/per_motif" \
    --merged-bed "${MERGEDIR}/all_motifs_merged.bed.gz" \
    --summary-tsv "${MERGEDIR}/motif_hit_summary.tsv" \
    --gtf "${GTF_FILE}" \
    --chrom-sizes "${CHROM_SIZES}" \
    --meme "${MEME_FILE}" \
    --background "${PREPDIR}/background.txt" \
    --pvalue "${PVALUE}" \
    --promoter-upstream "${PROMOTER_UPSTREAM}" \
    --promoter-downstream "${PROMOTER_DOWNSTREAM}" \
    --output-dir "${QCDIR}" \
    --tss-sample-n 50000 \
    --tss-max-dist 100000
