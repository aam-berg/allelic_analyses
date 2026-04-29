#!/usr/bin/env Rscript
# =============================================================================
# replot.R — Re-generate plots from saved intermediate data
#
# This avoids re-extracting signal from bigwig files. Useful for:
# - Changing plot aesthetics (colors, titles, etc.)
# - Generating overlay plots for new preset combinations
# - Adjusting normalization (if raw matrices were saved)
#
# Usage:
#   Rscript replot.R --motif_id AC0001 [--bw_set merged_3prime] [--preset all]
#   Rscript replot.R --motif_id AC0001 --renormalize zscore
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--motif_id",     type = "character", help = "Motif archetype ID"),
  make_option("--bw_set",       type = "character", default = "all"),
  make_option("--preset",       type = "character", default = "all"),
  make_option("--renormalize",  type = "character", default = "",
              help = "Re-normalize from raw matrices: minmax, zscore, or none"),
  make_option("--config",       type = "character", default = "config.R"),
  make_option("--overlay_pairs", type = "character", default = "",
              help = "Comma-separated pairs for overlay, e.g. 'p1:p2,p3:p4'")
)

opt <- parse_args(OptionParser(option_list = option_list))

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
if (length(script_dir) == 0 || script_dir == "") script_dir <- "."

source(file.path(script_dir, opt$config))
source(file.path(script_dir, "functions.R"))

motif_id <- opt$motif_id
out_base <- file.path(OUTPUT_BASE, motif_id)
data_dir <- file.path(out_base, "data")

if (!dir.exists(data_dir)) stop("No data directory found for ", motif_id)

# Load PWM
all_pwms <- parse_meme_pwms(PWM_MEME_FILE)
pwm <- all_pwms[[motif_id]]
motif_width <- nrow(pwm)
motif_half_width <- motif_width / 2

# Find all saved profile files
profile_files <- list.files(data_dir, pattern = "_profiles\\.tsv\\.gz$", full.names = TRUE)

if (length(profile_files) == 0) stop("No profile data found in ", data_dir)

message(sprintf("Found %d profile files for %s", length(profile_files), motif_id))

# --- Re-normalize if requested -----------------------------------------------
if (opt$renormalize != "") {
  message(sprintf("Re-normalizing with method: %s", opt$renormalize))

  matrix_files <- list.files(data_dir, pattern = "_matrices\\.rds$", full.names = TRUE)

  for (mf in matrix_files) {
    # Parse preset and bw_label from filename
    bn <- basename(mf)
    parts <- sub(paste0("^", motif_id, "_"), "", sub("_matrices\\.rds$", "", bn))
    # Find the bw_label (last part after the last underscore pair matching a bw set name)
    # This is a bit fragile; we'll try to reconstruct
    bw_sets <- get_bigwig_sets()

    found <- FALSE
    for (bw_name in names(bw_sets)) {
      if (grepl(paste0("_", bw_name, "$"), parts)) {
        preset_name <- sub(paste0("_", bw_name, "$"), "", parts)
        bw_label <- bw_name

        matrices <- readRDS(mf)
        new_profiles_list <- list()

        for (ori_name in names(matrices)) {
          mat <- matrices[[ori_name]]
          if (is.null(mat)) next

          mat <- abs(mat)
          mat_norm <- normalize_matrix(mat, method = opt$renormalize)

          valid <- rowSums(is.na(mat_norm)) < ncol(mat_norm)
          mat_valid <- mat_norm[valid, , drop = FALSE]

          mean_sig <- colMeans(mat_valid, na.rm = TRUE)
          se_sig <- apply(mat_valid, 2, function(x) sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x))))
          positions <- seq(-FLANK_BP, FLANK_BP, by = 1L)

          ori_label <- ORIENTATION_DEFS[[ori_name]]$label

          new_profiles_list[[ori_name]] <- data.table(
            orientation = ori_label,
            position    = positions,
            mean_signal = mean_sig,
            se_signal   = se_sig,
            n_instances = nrow(mat_valid)
          )
        }

        new_profiles <- rbindlist(new_profiles_list, fill = TRUE)

        # Overwrite profile file
        prof_file <- file.path(data_dir,
                               sprintf("%s_%s_%s_profiles.tsv.gz",
                                       motif_id, preset_name, bw_label))
        fwrite(new_profiles, prof_file, sep = "\t")
        message(sprintf("  Re-normalized: %s", basename(prof_file)))

        found <- TRUE
        break
      }
    }
    if (!found) message(sprintf("  Could not parse filename: %s", bn))
  }
}

# --- Re-plot from profile files ----------------------------------------------
dir.create(file.path(out_base, "plots_replot"), recursive = TRUE, showWarnings = FALSE)

for (pf in profile_files) {
  bn <- basename(pf)
  profiles <- fread(pf)

  # Parse filename to get preset_name and bw_label
  parts <- sub(paste0("^", motif_id, "_"), "", sub("_profiles\\.tsv\\.gz$", "", bn))

  # Try to match bw set names
  bw_sets <- get_bigwig_sets()
  matched <- FALSE

  for (bw_name in names(bw_sets)) {
    if (grepl(paste0("_", bw_name, "$"), parts)) {
      preset_name <- sub(paste0("_", bw_name, "$"), "", parts)
      bw_label <- bw_name

      # Load counts
      counts_file <- file.path(data_dir,
                               sprintf("%s_%s_%s_counts.tsv",
                                       motif_id, preset_name, bw_label))
      counts <- c()
      if (file.exists(counts_file)) {
        counts_dt <- fread(counts_file)
        counts <- setNames(counts_dt$n_instances, counts_dt$orientation)
      }

      p <- make_quadrant_plot(
        profiles_dt     = profiles,
        counts          = counts,
        pwm_matrix      = pwm,
        motif_id        = motif_id,
        bw_label        = bw_label,
        filter_label    = preset_name,
        motif_half_width = motif_half_width
      )

      plot_file <- file.path(out_base, "plots_replot",
                             sprintf("%s_%s_%s.pdf", motif_id, preset_name, bw_label))
      ggsave(plot_file, p, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
      ggsave(sub("\\.pdf$", ".png", plot_file), p,
             width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)

      message(sprintf("  Replotted: %s", basename(plot_file)))
      matched <- TRUE
      break
    }
  }

  if (!matched) message(sprintf("  Could not match: %s", bn))
}

message("\nReplotting complete.")