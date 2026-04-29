#!/bin/bash
#SBATCH --job-name=motif_snp_cov
#SBATCH --partition=short
#SBATCH --time=0-02:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/motif_snp_cov_%a.out
#SBATCH --error=logs/motif_snp_cov_%a.err
# ==============================================================================
# run_motif_snp_coverage.sh
# ==============================================================================
# Usage:
#   1. First, generate the cluster list:
#        mkdir -p logs
#        tail -n +2 /home/alb1273/pausing_phase_project/resources/metadata.tsv \
#          | cut -f2 | sort -u > cluster_list.txt
#        N=$(wc -l < cluster_list.txt)
#        echo "Total clusters: $N"
#
#   2. Submit the array:
#        sbatch --array=0-$((N-1)) run_motif_snp_coverage.sh
#
#   Or to limit concurrency (e.g. max 100 at a time):
#        sbatch --array=0-$((N-1))%100 run_motif_snp_coverage.sh
#
#   3. After all jobs finish, merge results:
#        Rscript merge_and_plot.R
# ==============================================================================

# ---- Paths ----
CLUSTER_LIST="cluster_list.txt"
RSCRIPT="process_single_cluster.R"

# ---- Get cluster ID for this array task ----
# Arrays are 0-indexed; sed is 1-indexed
LINE=$((SLURM_ARRAY_TASK_ID + 1))
CLUSTER_ID=$(sed -n "${LINE}p" "$CLUSTER_LIST")

if [ -z "$CLUSTER_ID" ]; then
    echo "[ERROR] No cluster found at line $LINE of $CLUSTER_LIST"
    exit 1
fi

echo "[INFO] Task $SLURM_ARRAY_TASK_ID → cluster $CLUSTER_ID"

# ---- Load R ----
module load R/4.3.1 2>/dev/null || module load R 2>/dev/null || true

# ---- Run ----
Rscript "$RSCRIPT" "$CLUSTER_ID"