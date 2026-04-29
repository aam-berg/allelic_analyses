#!/bin/bash
#SBATCH --job-name=moods_scan
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-5 ####700%50
#SBATCH --output=logs/moods_%A_%a.out
#SBATCH --error=logs/moods_%A_%a.err

# =============================================================================
# MOODS Array Job: Scan Genome for Motif Occurrences
# =============================================================================
# This script runs MOODS for each motif cluster in parallel using SLURM arrays.
#
# Usage:
#   sbatch --array=1-N%50 02_run_moods_array.sh <p_threshold> <manifest> <genome> <outdir>
#
# Example:
#   sbatch 02_run_moods_array.sh 1e-4 \
#       /path/to/cluster_manifest.txt \
#       /path/to/mm39.fa \
#       /path/to/output
#
# Parameters:
#   p_threshold: MOODS p-value threshold (e.g., 1e-4 or 1e-5)
#   manifest: Path to cluster_manifest.txt (one cluster ID per line)
#   genome: Path to mm39.fa reference genome
#   outdir: Base output directory
# =============================================================================

set -euo pipefail

# Parse arguments
P_THRESH="${1:-1e-5}"
MANIFEST="${2:-/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/motif_splitting/cluster_manifest.txt}"
GENOME="${3:-/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/genome/mm39.fa}"
BASE_OUTDIR="${4:-/n/scratch/users/a/alb1273/pausing_phase_project_intermediates}"
MEME_FILE="${5:-/home/alb1273/pausing_phase_project/resources/consensus_pwms.meme}"
METADATA_FILE="${6:-/home/alb1273/pausing_phase_project/resources/metadata.tsv}"

# Get cluster ID for this array task
CLUSTER=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${MANIFEST})

if [[ -z "${CLUSTER}" ]]; then
    echo "Error: No cluster found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "=============================================="
echo "MOODS Motif Scanning Array Job"
echo "=============================================="
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Cluster: ${CLUSTER}"
echo "P-value threshold: ${P_THRESH}"
echo "Output directory: ${BASE_OUTDIR}"
echo "Start time: $(date)"
echo "=============================================="

# Paths
SCRIPTS_DIR="/home/alb1273/pausing_phase_project/scripts"

# Load conda environment
source activate motif_scanning

# =============================================================================
# Run MOODS
# =============================================================================
echo ""
echo "Running MOODS for cluster ${CLUSTER}..."

python ${SCRIPTS_DIR}/scan_motifs_moods_improved.py \
    --meme "${MEME_FILE}" \
    --genome "${GENOME}" \
    --output "${BASE_OUTDIR}" \
    --pvalue ${P_THRESH} \
    --metadata "${METADATA_FILE}" \
    --cluster "${CLUSTER}"

echo ""
echo "=============================================="
echo "MOODS scan complete for ${CLUSTER}"
echo "End time: $(date)"
echo "=============================================="
