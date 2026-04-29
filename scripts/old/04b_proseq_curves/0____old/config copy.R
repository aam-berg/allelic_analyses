#!/usr/bin/env Rscript
# =============================================================================
# config.R — Central configuration for PRO-seq × TFBS pausing analysis pipeline
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

# Normalization: "minmax" (0-1 per instance), "zscore", "rpm", or "none"
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
# Available filter keys:
#   is_intragenic        : logical
#   is_promoter          : logical
#   atac_accessible      : logical (TRUE = overlap with EITHER ATAC dataset)
#   gene_type_include    : character vector of patterns to grep in intragenic_gene_types
#   gene_type_exclude    : character vector of patterns to exclude
#   genic_region_include : character vector of patterns to grep in genic_regions
#   genic_region_exclude : character vector of patterns to exclude
#   no_snp               : logical (TRUE = exclude SNP-overlapping sites)
#   min_score            : numeric, minimum motif match score
#   chroms_exclude       : character vector of chromosomes to exclude

FILTER_PRESETS <- list(
  # intragenic_accessible_pc_lnc = list(
  #   is_intragenic     = TRUE,
  #   is_promoter       = FALSE,
  #   atac_accessible   = TRUE,
  #   gene_type_include = c("protein_coding", "lncRNA"),
  #   chroms_exclude    = c("chrM", "chrY")
  # )

  # #### IGNORE for now; might be relevant soon
  promoter_accessible_pc_lnc = list(
    is_promoter       = TRUE,
    atac_accessible   = TRUE,
    gene_type_include = c("protein_coding", "lncRNA"),
    chroms_exclude    = c("chrM", "chrY")
  )
  #### IGNORE the below (probably entirely)
  # intragenic_pc_lnc = list(
  #   is_intragenic     = TRUE,
  #   is_promoter       = FALSE,
  #   gene_type_include = c("protein_coding", "lncRNA"),
  #   chroms_exclude    = c("chrM", "chrY")
  # ),
  # promoter_pc_lnc = list(
  #   is_promoter       = TRUE,
  #   gene_type_include = c("protein_coding", "lncRNA"),
  #   chroms_exclude    = c("chrM", "chrY")
  # ),
  # intragenic_accessible_pc_intron = list(
  #   is_intragenic        = TRUE,
  #   is_promoter          = FALSE,
  #   atac_accessible      = TRUE,
  #   gene_type_include    = c("protein_coding"),
  #   genic_region_include = c("intron"),
  #   genic_region_exclude = c("exon", "CDS", "UTR"),
  #   chroms_exclude       = c("chrM", "chrY")
  # ),
  # intragenic_accessible_pc_exon = list(
  #   is_intragenic        = TRUE,
  #   is_promoter          = FALSE,
  #   atac_accessible      = TRUE,
  #   gene_type_include    = c("protein_coding"),
  #   genic_region_include = c("exon"),
  #   chroms_exclude       = c("chrM", "chrY")
  # )
)

ACTIVE_PRESETS <- names(FILTER_PRESETS)

# ----- Overlay plot pairs (presets to overlay on same axes) -------------------
# OVERLAY_PAIRS <- list(
#   c("intragenic_accessible_pc_lnc")
#   c("promoter_accessible_pc_lnc")
#   #c("intragenic_accessible_pc_intron", "intragenic_accessible_pc_exon")
# )
OVERLAY_PAIRS <- list()   # no overlay pairs (use when only one active preset)

# ----- Plot aesthetics -------------------------------------------------------
MOTIF_LINE_COLOR <- "red"
MOTIF_LINE_TYPE  <- "twodash"


# ----- SLURM parameters -----------------------------------------------------
SLURM_PARTITION <- "short"
SLURM_TIME      <- "01:00:00"
SLURM_MEM       <- "16G"
SLURM_CPUS      <- 1
SLURM_MODULES   <- c("R/4.3.1", "gcc/9.2.0")


# =============================================================================
# config_v2_additions.R — New/changed settings for v2
#
# Add these to your existing config.R (or adjust existing values).
# =============================================================================

# --- CHANGED: Minimum instances lowered from 1000 to 100 ---
MIN_INSTANCES <- 1   # was 1000 in v1

# --- NEW: Path to TF metadata file ---
# Contains motif_id <-> tf_name mapping (the "cluster" column = motif archetype ID)
METADATA_FILE <- "/home/alb1273/pausing_phase_project/resources/metadata.tsv"

# --- NEW/UPDATED: Plot dimensions ---
# Increased to accommodate larger fonts and sequence logo
PLOT_WIDTH  <- 18   # inches (was ~12)
PLOT_HEIGHT <- 16   # inches (was ~10)
PLOT_DPI    <- 200

# --- NEW: Overlay colors and labels ---
# Named vector: preset_name -> color
# These are used for multi-curve overlay plots.
# The names also serve as legend labels.
# To customize legend labels, use descriptive names as keys:
#
# Or with custom display labels, rename your presets in OVERLAY_PAIRS accordingly.
OVERLAY_COLORS <- c(
  "#D55E00",  # vermillion (first preset in each pair)
  "#0072B2"   # blue       (second preset in each pair)
)
# NOTE: If OVERLAY_COLORS is unnamed, it will be auto-mapped to preset names.
# For explicit control, use a named vector matching your preset names.