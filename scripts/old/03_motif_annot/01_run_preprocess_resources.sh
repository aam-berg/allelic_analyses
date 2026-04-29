#!/bin/bash
#SBATCH --job-name=preprocess_resources
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/preprocess_resources_%j.out
#SBATCH --error=logs/preprocess_resources_%j.err

# =============================================================================
# 01_run_preprocess_resources.sh — SLURM wrapper for 01_preprocess_resources.R
#
# Reads paths from config.sh.
#
# NOTE: The chain file for mm10→mm39 liftOver will be auto-downloaded by the
# R script if missing. If your compute nodes lack internet, pre-download on
# a login node:
#   wget -P "$(dirname ${CHAIN_FILE})" \
#     https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

mkdir -p logs "${PREPROCESSED_DIR}"

echo "Job ID: ${SLURM_JOB_ID:-N/A}"
echo "Node:   $(hostname)"
echo "Date:   $(date)"
echo ""
echo "  GTF:                ${GTF_FILE}"
echo "  ATAC dir:           ${ATAC_DIR}"
echo "  ATAC files:         ${ATAC_FILES}"
echo "  Chain file:         ${CHAIN_FILE}"
echo "  Output:             ${PREPROCESSED_DIR}"
echo ""

set +u
source activate "${CONDA_ENV}"
set -u

Rscript "${SCRIPT_DIR}/01_preprocess_resources.R" \
    --gtf "${GTF_FILE}" \
    --atac_dir "${ATAC_DIR}" \
    --atac_files "${ATAC_FILES}" \
    --chain_file "${CHAIN_FILE}" \
    --outdir "${PREPROCESSED_DIR}"

echo ""
echo "Preprocessing finished: $(date)"
