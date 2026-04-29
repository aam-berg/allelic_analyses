#!/bin/bash
# =============================================================================
# Master Script: Run Complete FIMO Motif Finding Workflow
# =============================================================================
# This script orchestrates:
# 1. Splitting the MEME file by cluster
# 2. Running FIMO at p < 1e-4
# 3. Running FIMO at p < 1e-5
# 4. Aggregating results
#
# Usage:
#   ./run_fimo_workflow.sh
#
# Prerequisites:
#   - conda environment 'pausing_phase' created
#   - Resource files downloaded (genome, etc.)
# =============================================================================

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================
PROJECT_DIR="/home/alb1273/pausing_phase_project"
SCRATCH_DIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates"
SCRIPTS_DIR="${PROJECT_DIR}/scripts"

# Input files
MEME_FILE="${PROJECT_DIR}/resources/consensus_pwms.meme"
METADATA_FILE="${PROJECT_DIR}/resources/metadata.tsv"
GENOME="${SCRATCH_DIR}/resources/genome/mm39.fa"    # Note change from mm10 to mm39

# Output directories
SPLIT_DIR="${SCRATCH_DIR}/motif_splitting"
LOG_DIR="${SCRATCH_DIR}/logs"

# Create directories
mkdir -p "${SPLIT_DIR}" "${LOG_DIR}"

echo "=============================================="
echo "FIMO Motif Finding Workflow"
echo "=============================================="
echo "Project: ${PROJECT_DIR}"
echo "Scratch: ${SCRATCH_DIR}"
echo "Start time: $(date)"
echo "=============================================="

# =============================================================================
# Step 0: Check prerequisites
# =============================================================================
echo ""
echo "[0/4] Checking prerequisites..."

# Check conda environment exists
if ! conda env list | grep -q "pausing_phase"; then
    echo "  Creating conda environment..."
    conda env create -f ${PROJECT_DIR}/envs/pausing_phase.yaml
fi

# Check genome exists
if [[ ! -f "${GENOME}" ]]; then
    echo "Error: Genome file not found: ${GENOME}"
    echo "Please run 01_download_resources.sh first"
    exit 1
fi

# Check MEME file exists
if [[ ! -f "${MEME_FILE}" ]]; then
    echo "Error: MEME file not found: ${MEME_FILE}"
    exit 1
fi

echo "  All prerequisites satisfied"

# =============================================================================
# Step 1: Split MEME file by cluster
# =============================================================================
echo ""
echo "[1/4] Splitting MEME file by cluster..."

# Activate environment
source activate pausing_phase 2>/dev/null || conda activate pausing_phase

python ${SCRIPTS_DIR}/split_meme_by_cluster.py \
    --meme "${MEME_FILE}" \
    --metadata "${METADATA_FILE}" \
    --outdir "${SPLIT_DIR}"

# Get number of clusters
N_CLUSTERS=$(wc -l < "${SPLIT_DIR}/cluster_manifest.txt")
echo "  Found ${N_CLUSTERS} clusters"

# =============================================================================
# Step 2: Submit FIMO array jobs
# =============================================================================
echo ""
echo "[2/4] Submitting FIMO array jobs..."

# Create log directories
mkdir -p ${LOG_DIR}/fimo_em4
mkdir -p ${LOG_DIR}/fimo_em5

cd ${LOG_DIR}/fimo_em4

# Submit for p < 1e-4
echo "  Submitting FIMO jobs for p < 1e-4..."
JOB_ID_EM4=$(sbatch \
    --array=1-${N_CLUSTERS}%50 \
    --output=fimo_%A_%a.out \
    --error=fimo_%A_%a.err \
    ${SCRIPTS_DIR}/02_run_fimo_array.sh \
    1e-4 \
    "${SPLIT_DIR}/cluster_manifest.txt" \
    "${GENOME}" \
    "${SCRATCH_DIR}" \
    | awk '{print $4}')

echo "    Submitted job ${JOB_ID_EM4} (p < 1e-4)"

cd ${LOG_DIR}/fimo_em5

# Submit for p < 1e-5
echo "  Submitting FIMO jobs for p < 1e-5..."
JOB_ID_EM5=$(sbatch \
    --array=1-${N_CLUSTERS}%50 \
    --output=fimo_%A_%a.out \
    --error=fimo_%A_%a.err \
    ${SCRIPTS_DIR}/02_run_fimo_array.sh \
    1e-5 \
    "${SPLIT_DIR}/cluster_manifest.txt" \
    "${GENOME}" \
    "${SCRATCH_DIR}" \
    | awk '{print $4}')

echo "    Submitted job ${JOB_ID_EM5} (p < 1e-5)"

# =============================================================================
# Step 3: Submit aggregation job (runs after FIMO completes)
# =============================================================================
echo ""
echo "[3/4] Submitting aggregation job..."

cd ${LOG_DIR}

# Submit aggregation job with dependency on both FIMO jobs
AGG_JOB=$(sbatch \
    --dependency=afterok:${JOB_ID_EM4}:${JOB_ID_EM5} \
    --output=aggregate_%j.out \
    --error=aggregate_%j.err \
    ${SCRIPTS_DIR}/03_aggregate_fimo_results.sh \
    "${SCRATCH_DIR}" \
    | awk '{print $4}')

echo "  Submitted aggregation job ${AGG_JOB} (depends on ${JOB_ID_EM4} and ${JOB_ID_EM5})"

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "=============================================="
echo "Workflow submitted successfully!"
echo "=============================================="
echo ""
echo "Jobs submitted:"
echo "  FIMO p<1e-4: ${JOB_ID_EM4} (${N_CLUSTERS} tasks)"
echo "  FIMO p<1e-5: ${JOB_ID_EM5} (${N_CLUSTERS} tasks)"
echo "  Aggregation: ${AGG_JOB} (depends on FIMO jobs)"
echo ""
echo "Monitor with:"
echo "  squeue -u \$USER"
echo "  sacct -j ${JOB_ID_EM4},${JOB_ID_EM5},${AGG_JOB}"
echo ""
echo "Logs will be in:"
echo "  ${LOG_DIR}/fimo_em4/"
echo "  ${LOG_DIR}/fimo_em5/"
echo ""
echo "Results will be in:"
echo "  ${SCRATCH_DIR}/motif_archetypes_em4/motif_locations/"
echo "  ${SCRATCH_DIR}/motif_archetypes_em5/motif_locations/"