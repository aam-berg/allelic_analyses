#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# run_qc_plots.sh — Submit the two QC/prep plotting scripts
# =============================================================================
#
# USAGE:
#   bash run_qc_plots.sh
#
# DEPENDENCIES:
#   R packages: data.table, ggplot2, patchwork, scales, viridis
#   Install if missing:
#     Rscript -e 'install.packages(c("patchwork", "viridis"), repos="https://cloud.r-project.org")'
#   (data.table, ggplot2, scales are likely already installed)
#
# SBATCH RESOURCES:
#   The motif archetype script reads ~hundreds of gzipped TSVs, so give it
#   enough memory and time:
#     QC plots:      ~5 min, 4G RAM
#     Archetype plots: ~10-30 min (depending on number of archetypes), 16G RAM
# =============================================================================

SCRIPT_DIR="."

echo "================================================================"
echo "Submitting QC plotting jobs"
echo "================================================================"

# --- Check R packages ---
echo "[INFO] Checking R package availability..."
Rscript -e '
pkgs <- c("data.table", "ggplot2", "patchwork", "scales", "viridis")
missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
    cat("Installing missing packages:", paste(missing, collapse = ", "), "\n")
    install.packages(missing, repos = "https://cloud.r-project.org", quiet = TRUE)
} else {
    cat("All required packages available.\n")
}
'

# --- Plot 1: PRO-seq + WASP QC ---
echo ""
echo "[1] PRO-seq + WASP QC plots..."
sbatch --partition=short --time=1:00:00 --mem=4G --cpus-per-task=1 \
    -o "${SCRIPT_DIR}/logs/plot_qc_%j.out" \
    -e "${SCRIPT_DIR}/logs/plot_qc_%j.err" \
    --wrap="Rscript ${SCRIPT_DIR}/plot_proseq_wasp_qc.R"

# --- Plot 2: Motif archetype prep ---
echo "[2] Motif archetype preparatory plots..."
sbatch --partition=short --time=2:00:00 --mem=16G --cpus-per-task=1 \
    -o "${SCRIPT_DIR}/logs/plot_motif_%j.out" \
    -e "${SCRIPT_DIR}/logs/plot_motif_%j.err" \
    --wrap="Rscript ${SCRIPT_DIR}/plot_motif_archetype_prep.R"

echo ""
echo "[DONE] Jobs submitted. Check logs/ for output."
echo ""
echo "Output locations:"
echo "  QC plots:      /n/scratch/.../proseq/qc/plots/"
echo "  Archetype plots: /n/scratch/.../annot_motifs/prep_plots/"