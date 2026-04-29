#!/usr/bin/env bash
# =============================================================================
# submit_01.sh — Submit 01_download_and_process.sh as a SLURM job array
# =============================================================================
#
# Each sample runs as a separate SLURM task, so all samples process in parallel.
#
# USAGE:
#   sbatch submit_01.sh            # Submit all samples
#
# To rerun just one failed task (e.g. task 1):
#   sbatch --array=1 submit_01.sh
#
# =============================================================================

#SBATCH --job-name=proseq_proc
#SBATCH --array=0-3
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --output=logs/proseq_proc_%A_%a.out
#SBATCH --error=logs/proseq_proc_%A_%a.err

# --- Map SLURM_ARRAY_TASK_ID → SRR accession ---
# UPDATE THIS ARRAY if you add/remove samples.
# Order doesn't need to match config.sh; the processing script finds the
# matching sample name automatically.
SRRS=(
    "SRR17640044"   # index 0 → WT_mESC_PROseq_Rep1
    "SRR17640045"   # index 1 → WT_mESC_PROseq_Rep2
    "SRR17640046"   # index 0 → WT_mESC_PROseq_Rep1
    "SRR17640047"   # index 1 → WT_mESC_PROseq_Rep2
)

SRR="${SRRS[$SLURM_ARRAY_TASK_ID]}"
echo "[submit] Array task ${SLURM_ARRAY_TASK_ID} → ${SRR}"
echo "[submit] Node: $(hostname), CPUs: ${SLURM_CPUS_PER_TASK}, Time limit: $(squeue -j ${SLURM_JOB_ID} -h -o %l)"

mkdir -p logs

# Run the processing script for this single sample
bash 01_download_and_process.sh "${SRR}"