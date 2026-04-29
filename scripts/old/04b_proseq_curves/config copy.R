#!/usr/bin/env Rscript
# =============================================================================
# config.R — Central configuration for PRO-seq × TFBS pausing analysis pipeline
#
# v3: Two-encounter model (co-directional / anti-directional)
# =============================================================================

# ----- Paths -----------------------------------------------------------------
PROSEQ_BW_DIR   <- "/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/proseq/bigwig"
ANNOT_MOTIF_DIR <- "/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/annot_motifs"
METADATA_FILE   <- "/home/alb1273/pausing_phase_project/resources/metadata.tsv"
PWM_MEME_FILE   <- "/home/alb1273/pausing_phase_project/resources/consensus_pwms.meme"

# Output base directory
OUTPUT_BASE     <- "/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/proseq_pausing_analysis"

# ----- Analysis parameters ---------------------------------------------------
FLANK_BP        <- 500L    # bp flanking each side of motif center
BIN_SIZE        <- 1L      # resolution in bp (1 = single-nt)

# Normalization: "minmax" (0-1 per instance), "zscore", or "none"
NORM_METHOD     <- "minmax"

# Minimum number of motif instances required to produce a plot
MIN_INSTANCES   <- 100L

# ----- BigWig file sets ------------------------------------------------------
get_bigwig_sets <- function() {
  bw_dir <- PROSEQ_BW_DIR
  sets <- list()

  # Merged bigwigs
  for (end_type in c("3prime", "5prime")) {
    sets[[paste0("merged_", end_type)]] <- list(
      plus  = file.path(bw_dir, "merged", sprintf("WT_PROseq_merged_%s_plus.bw",  end_type)),
      minus = file.path(bw_dir, "merged", sprintf("WT_PROseq_merged_%s_minus.bw", end_type)),
      label = paste0("merged_", end_type)
    )
  }

  # Individual replicates
  for (rep in 1:4) {
    for (end_type in c("3prime", "5prime")) {
      tag <- sprintf("rep%d_%s", rep, end_type)
      sets[[tag]] <- list(
        plus  = file.path(bw_dir, "individual", sprintf("WT_PROseq_rep%d_%s_plus.bw",  rep, end_type)),
        minus = file.path(bw_dir, "individual", sprintf("WT_PROseq_rep%d_%s_minus.bw", rep, end_type)),
        label = tag
      )
    }
  }
  return(sets)
}

# ----- TFBS filter presets ---------------------------------------------------
FILTER_PRESETS <- list(
  promoter_accessible_pc_lnc = list(
    is_promoter       = TRUE,
    atac_accessible   = TRUE,
    gene_type_include = c("protein_coding", "lncRNA"),
    chroms_exclude    = c("chrM", "chrY")
  )
)

ACTIVE_PRESETS <- names(FILTER_PRESETS)

# ----- Overlay plot pairs (presets to overlay on same axes) -------------------
OVERLAY_PAIRS <- list()   # no overlay pairs when only one active preset

# ----- Plot aesthetics -------------------------------------------------------
PLOT_WIDTH       <- 16    # inches
PLOT_HEIGHT      <- 14    # 2 panels + 2 logos (smaller than 4-panel version)
PLOT_DPI         <- 200
MOTIF_LINE_COLOR <- "red"
MOTIF_LINE_TYPE  <- "dashed"

OVERLAY_COLORS <- c(
  "#D55E00",  # vermillion (first preset in each pair)
  "#0072B2"   # blue       (second preset in each pair)
)

# ----- SLURM parameters -----------------------------------------------------
SLURM_PARTITION <- "short"
SLURM_TIME      <- "01:00:00"
SLURM_MEM       <- "16G"
SLURM_CPUS      <- 1
SLURM_MODULES   <- c("R/4.3.1", "gcc/9.2.0")