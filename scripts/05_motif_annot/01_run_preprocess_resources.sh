#!/bin/bash
#SBATCH --job-name=annot_v2_preprocess
#SBATCH --output=logs/preprocess_%j.out
#SBATCH --error=logs/preprocess_%j.err

# =============================================================================
# 01_run_preprocess_resources.sh — SLURM wrapper for 01_preprocess_resources.R
#
# This must be run ONCE before submitting per-motif annotation jobs.
# Reads paths and parameters from config.sh.
#
# USAGE:
#   sbatch 01_run_preprocess_resources.sh
#
# OR locally (no SLURM):
#   bash 01_run_preprocess_resources.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="."
source "${SCRIPT_DIR}/config.sh"

# Apply SLURM resource defaults from config if running interactively
# (when submitted via sbatch, these go in the #SBATCH header, but we use
# CLI sbatch flags from run_pipeline.sh instead — see SLURM directives in
# this file are minimal as a result).

mkdir -p logs "${PREPROCESSED_DIR}"

echo "================================================================"
echo "05_motif_annot/ preprocessing"
echo "Job ID: ${SLURM_JOB_ID:-N/A}"
echo "Node:   $(hostname)"
echo "Date:   $(date)"
echo "================================================================"
echo ""
echo "Inputs:"
echo "  GTF:                    ${GTF_FILE}"
echo "  ATAC consensus BED:     ${ATAC_CONSENSUS_BED}"
echo "  RNA-seq expression:     ${RNASEQ_EXPRESSION_TSV}"
echo "  RNA-seq SJ glob:        ${RNASEQ_SJ_GLOB}"
echo "  SJ min unique reads:    ${SJ_MIN_UNIQUE_READS}"
echo "  Output:                 ${PREPROCESSED_DIR}"
echo ""

set +u
source activate "${CONDA_ENV}"
set -u

Rscript "${SCRIPT_DIR}/01_preprocess_resources.R" \
    --gtf "${GTF_FILE}" \
    --atac_consensus_bed "${ATAC_CONSENSUS_BED}" \
    --rnaseq_expression_tsv "${RNASEQ_EXPRESSION_TSV}" \
    --rnaseq_sj_glob "${RNASEQ_SJ_GLOB}" \
    --sj_min_unique_reads "${SJ_MIN_UNIQUE_READS}" \
    --outdir "${PREPROCESSED_DIR}"

echo ""
echo "Preprocessing finished: $(date)"
