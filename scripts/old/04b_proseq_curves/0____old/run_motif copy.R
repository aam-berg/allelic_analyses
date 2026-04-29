#!/usr/bin/env Rscript
# =============================================================================
# run_motif.R — Process a single motif archetype (v2)
#
# Changes from v1:
#   - Default: only merged bigwigs (use --include_replicates to add individual)
#   - TF name subtitle from metadata.tsv
#   - Sequence logos (ggseqlogo) instead of color bars
#   - Greatly increased font sizes (documented with FONT_SIZE comments)
#   - Shared y-lim across all subplots for a motif archetype
#   - Cumulative SNP overlap stats displayed on plots
#   - MIN_INSTANCES default = 100
#   - Writes skip note file when too few instances
#   - Color-coded overlay curves with legend labels
#
# Usage:
#   Rscript run_motif.R --motif_id AC0001
#   Rscript run_motif.R --motif_id AC0001 --include_replicates
#   Rscript run_motif.R --motif_id AC0001 --bw_set merged_3prime --preset intragenic_accessible_pc
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(ggplot2)
  library(ggseqlogo)    # for sequence logo rendering
  library(grid)
  library(gridExtra)
  library(patchwork)    # for combining plots with shared axes
})

# --- Parse arguments ---------------------------------------------------------
option_list <- list(
  make_option("--motif_id", type = "character", help = "Motif archetype ID (e.g. AC0001)"),
  make_option("--bw_set",   type = "character", default = "all",
              help = "BigWig set label or 'all' [default: all]"),
  make_option("--preset",   type = "character", default = "all",
              help = "Filter preset name or 'all' [default: all]"),
  make_option("--config",   type = "character", default = "config.R",
              help = "Path to config.R [default: config.R]"),
  # CHANGED: default is now TRUE (skip replicates by default; merged only)
  make_option("--include_replicates", type = "logical", default = FALSE,
              action = "store_true",
              help = "Include individual replicate bigwigs (default: merged only)"),
  make_option("--metadata", type = "character", default = NULL,
              help = "Path to metadata.tsv with TF names [default: from config]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --- Source config and functions ---------------------------------------------
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
if (length(script_dir) == 0 || script_dir == "") script_dir <- "."

source(file.path(script_dir, opt$config))
source(file.path(script_dir, "functions.R"))
source(file.path(script_dir, "plotting_functions.R"))

motif_id <- opt$motif_id
if (is.null(motif_id)) stop("Must specify --motif_id")

message(sprintf("\n========== Processing motif: %s ==========\n", motif_id))

# --- Setup output directories ------------------------------------------------
out_base <- file.path(OUTPUT_BASE, motif_id)
dir.create(file.path(out_base, "plots"),         recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_base, "data"),           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_base, "rederived_pwms"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_base, "plots_overlay"),  recursive = TRUE, showWarnings = FALSE)

# --- Load TF metadata --------------------------------------------------------
metadata_path <- if (!is.null(opt$metadata)) opt$metadata else METADATA_FILE
if (file.exists(metadata_path)) {
  message("Loading TF metadata...")
  tf_meta <- fread(metadata_path)
  # Get unique TF names for this motif archetype (cluster column = motif_id)
  # NEW (use a temporary variable with a name that doesn't collide):
  target_id <- motif_id
  tf_names_for_motif <- sort(unique(
      tf_meta[cluster == target_id, tf_name]
  ))
    
  if (length(tf_names_for_motif) == 0) {
    tf_subtitle <- "(no associated TFs in metadata)"
  } else {
    tf_subtitle <- paste(tf_names_for_motif, collapse = ", ")
  }
  message(sprintf("  TFs for %s: %s", motif_id, tf_subtitle))
} else {
  warning(sprintf("Metadata file not found: %s", metadata_path))
  tf_names_for_motif <- character(0)
  tf_subtitle <- "(metadata not available)"
}

# --- Load annotated TFBS data ------------------------------------------------
annot_file <- file.path(ANNOT_MOTIF_DIR, paste0(motif_id, "_annotated.tsv.gz"))
if (!file.exists(annot_file)) {
  stop(sprintf("Annotated TFBS file not found: %s", annot_file))
}

message("Loading annotated TFBS data...")
dt <- fread(annot_file)
message(sprintf("  Loaded %d rows (%d unique positions)",
                nrow(dt), length(unique(paste(dt$chrom, dt$start, dt$end)))))

# --- Load PWMs ---------------------------------------------------------------
message("Parsing PWMs...")
all_pwms <- parse_meme_pwms(PWM_MEME_FILE)
if (!motif_id %in% names(all_pwms)) {
  stop(sprintf("Motif %s not found in PWM file", motif_id))
}
pwm <- all_pwms[[motif_id]]
motif_width <- nrow(pwm)
motif_half_width <- motif_width / 2
consensus <- pwm_consensus(pwm)
message(sprintf("  Motif %s: width=%d, consensus=%s", motif_id, motif_width, consensus))

# --- Determine which BigWig sets and presets to run --------------------------
bw_sets <- get_bigwig_sets()
if (opt$bw_set != "all") {
  bw_sets <- bw_sets[opt$bw_set]
}
# CHANGED: Default is merged-only. Use --include_replicates to include individual.
if (!opt$include_replicates) {
  bw_sets <- bw_sets[grep("^merged_", names(bw_sets))]
  message("  Running merged bigwigs only (use --include_replicates for all)")
}

presets <- FILTER_PRESETS[ACTIVE_PRESETS]
if (opt$preset != "all") {
  presets <- presets[opt$preset]
}

# --- Verify bigwig files exist -----------------------------------------------
for (bw_name in names(bw_sets)) {
  bw <- bw_sets[[bw_name]]
  if (!file.exists(bw$plus))  warning(sprintf("BigWig not found: %s", bw$plus))
  if (!file.exists(bw$minus)) warning(sprintf("BigWig not found: %s", bw$minus))
}

# --- Storage for overlay plots and shared y-lim computation ------------------
# Key: bw_label, Value: list of (preset_name -> profiles, counts, snp_stats)
overlay_cache <- list()

# Collect all y-values across ALL preset × bw combos for shared y-lim
global_y_values <- numeric(0)

# Track which presets were skipped
skipped_presets <- list()

# =============================================================================
# HELPER: Compute cumulative SNP stats from filtered data
# =============================================================================
# The annotation stores SNP overlaps in ANNULAR (non-overlapping) bins:
#   snp_F121-9_overlap       = direct motif overlap
#   snp_F121-9_flank10bp_overlap  = 1–10bp from edge  (annular)
#   snp_F121-9_flank25bp_overlap  = 11–25bp from edge (annular)
#   snp_F121-9_flank50bp_overlap  = 26–50bp from edge (annular)
#   snp_F121-9_flank100bp_overlap = 51–100bp from edge (annular)
#
# We compute CUMULATIVE counts: "within Xbp" = has SNP in motif OR any
# annular bin up to X bp.

compute_cumulative_snp_stats <- function(filtered_dt, hybrid = "F121-9") {
  prefix <- paste0("snp_", hybrid)

  # Column names for direct and annular overlaps
  col_direct <- paste0(prefix, "_overlap")
  flank_dists <- c(10, 25, 50, 100)
  col_flanks <- paste0(prefix, "_flank", flank_dists, "bp_overlap")

  # Check which columns exist
  available <- c(col_direct, col_flanks)
  available <- available[available %in% names(filtered_dt)]

  if (length(available) == 0) {
    return(data.table(
      label = c("direct overlap", "within 10bp", "within 25bp",
                "within 50bp", "within 100bp"),
      count = rep(NA_integer_, 5),
      total = nrow(filtered_dt)
    ))
  }

  n_total <- nrow(filtered_dt)

  # Direct overlap
  n_direct <- if (col_direct %in% names(filtered_dt)) {
    sum(filtered_dt[[col_direct]], na.rm = TRUE)
  } else { NA_integer_ }

  # Cumulative: within Xbp = direct OR any annular bin with d_outer <= X
  # Build cumulative logical vectors
  cum_cols <- list(
    "within 10bp"  = c(col_direct, paste0(prefix, "_flank10bp_overlap")),
    "within 25bp"  = c(col_direct, paste0(prefix, "_flank10bp_overlap"),
                        paste0(prefix, "_flank25bp_overlap")),
    "within 50bp"  = c(col_direct, paste0(prefix, "_flank10bp_overlap"),
                        paste0(prefix, "_flank25bp_overlap"),
                        paste0(prefix, "_flank50bp_overlap")),
    "within 100bp" = c(col_direct, paste0(prefix, "_flank10bp_overlap"),
                        paste0(prefix, "_flank25bp_overlap"),
                        paste0(prefix, "_flank50bp_overlap"),
                        paste0(prefix, "_flank100bp_overlap"))
  )

  results <- data.table(
    label = c("direct overlap", names(cum_cols)),
    count = NA_integer_,
    total = n_total
  )
  results[label == "direct overlap", count := n_direct]

  for (lbl in names(cum_cols)) {
    cols_needed <- cum_cols[[lbl]]
    cols_present <- cols_needed[cols_needed %in% names(filtered_dt)]
    if (length(cols_present) > 0) {
      # A motif has SNP within X bp if ANY of the relevant columns is TRUE
      has_snp <- Reduce(`|`, lapply(cols_present, function(cc) {
        filtered_dt[[cc]]
      }))
      results[label == lbl, count := sum(has_snp, na.rm = TRUE)]
    }
  }

  return(results)
}

# =============================================================================
# MAIN LOOP: iterate over filter presets × bigwig sets
# =============================================================================

# --- PASS 1: Compute all profiles, collect y-ranges -------------------------
message("\n========== PASS 1: Computing profiles ==========\n")

all_results <- list()  # store everything for pass 2

for (preset_name in names(presets)) {
  filters <- presets[[preset_name]]

  message(sprintf("\n--- Filter preset: %s ---", preset_name))

  # Filter TFBSs
  filtered <- filter_tfbs(dt, filters)
  n_plus  <- sum(filtered$strand == "+")
  n_minus <- sum(filtered$strand == "-")
  message(sprintf("  After filtering: %d rows (%d on + strand, %d on - strand)",
                  nrow(filtered), n_plus, n_minus))

  if (nrow(filtered) < MIN_INSTANCES) {
    message(sprintf("  Too few instances (%d < %d), skipping this preset.",
                    nrow(filtered), MIN_INSTANCES))
    skipped_presets[[preset_name]] <- list(
      reason = "too_few_instances",
      n_found = nrow(filtered),
      min_required = MIN_INSTANCES
    )
    next
  }

  # Compute cumulative SNP stats for this filtered subset
  snp_stats <- compute_cumulative_snp_stats(filtered)

  # Re-derive PWMs from filtered subsets (sanity check)
  for (s in c("+", "-")) {
    sub <- filtered[strand == s]
    if (nrow(sub) > 0) {
      rederived <- rederive_pwm(sub, motif_width)
      if (!is.null(rederived)) {
        rederived_file <- file.path(out_base, "rederived_pwms",
                                     sprintf("%s_%s_%s_strand.tsv",
                                             motif_id, preset_name,
                                             ifelse(s == "+", "plus", "minus")))
        fwrite(as.data.table(rederived), rederived_file, sep = "\t")
        message(sprintf("  Re-derived PWM (%s strand, n=%d): %s",
                        s, nrow(sub), pwm_consensus(rederived)))
      }
    }
  }

  # Save filtered TFBS info
  filtered_meta <- filtered[, .(chrom, start, end, strand, score, motif_sequence,
                                 is_intragenic, is_promoter, genic_regions,
                                 intragenic_gene_types)]
  fwrite(filtered_meta,
         file.path(out_base, "data",
                   sprintf("%s_%s_filtered_sites.tsv.gz", motif_id, preset_name)),
         sep = "\t")

  # --- Process each BigWig set -----------------------------------------------
  for (bw_name in names(bw_sets)) {
    bw <- bw_sets[[bw_name]]
    bw_label <- bw$label

    message(sprintf("\n  BigWig set: %s", bw_label))

    if (!file.exists(bw$plus) || !file.exists(bw$minus)) {
      message("    BigWig files not found, skipping.")
      next
    }

    # Compute metaprofiles for all 4 orientations
    result <- compute_metaprofiles(
      filtered_dt   = filtered,
      bw_plus       = bw$plus,
      bw_minus      = bw$minus,
      flank_bp      = FLANK_BP,
      norm_method   = NORM_METHOD,
      min_instances = MIN_INSTANCES
    )

    profiles <- result$profiles
    counts   <- result$counts

    if (nrow(profiles) == 0) {
      message("    No profiles computed (insufficient data).")
      next
    }

    # Collect y-values for global y-lim
    if ("mean_signal" %in% names(profiles)) {
      global_y_values <- c(global_y_values, profiles$mean_signal)
    }

    # --- Save data -----------------------------------------------------------
    data_file <- file.path(out_base, "data",
                           sprintf("%s_%s_%s_profiles.tsv.gz",
                                   motif_id, preset_name, bw_label))
    fwrite(profiles, data_file, sep = "\t")
    message(sprintf("    Saved profile data: %s", basename(data_file)))

    counts_dt <- data.table(
      orientation = names(counts),
      n_instances = counts
    )
    counts_file <- file.path(out_base, "data",
                             sprintf("%s_%s_%s_counts.tsv",
                                     motif_id, preset_name, bw_label))
    fwrite(counts_dt, counts_file, sep = "\t")

    matrices_file <- file.path(out_base, "data",
                               sprintf("%s_%s_%s_matrices.rds",
                                       motif_id, preset_name, bw_label))
    saveRDS(result$matrices, matrices_file)

    # Store for pass 2
    key <- paste(preset_name, bw_label, sep = "||")
    all_results[[key]] <- list(
      profiles    = profiles,
      counts      = counts,
      preset_name = preset_name,
      bw_label    = bw_label,
      snp_stats   = snp_stats
    )

    # Cache for overlay plots
    cache_key <- bw_label
    if (is.null(overlay_cache[[cache_key]])) {
      overlay_cache[[cache_key]] <- list()
    }
    overlay_cache[[cache_key]][[preset_name]] <- list(
      profiles  = profiles,
      counts    = counts,
      snp_stats = snp_stats
    )
  }  # end bw_sets loop
}  # end presets loop


# =============================================================================
# Compute shared y-lim across ALL subplots for this motif archetype
# =============================================================================

if (length(global_y_values) > 0) {
  shared_ylim <- c(
    min(global_y_values, na.rm = TRUE),
    max(global_y_values, na.rm = TRUE) * 1.05  # 5% headroom
  )
  # If min is positive, start at 0 for cleaner plots
  if (shared_ylim[1] >= 0) shared_ylim[1] <- 0
  message(sprintf("\nShared y-lim: [%.4f, %.4f]", shared_ylim[1], shared_ylim[2]))
} else {
  shared_ylim <- NULL
}


# =============================================================================
# PASS 2: Generate plots with shared y-lim
# =============================================================================

message("\n========== PASS 2: Generating plots ==========\n")

for (key in names(all_results)) {
  res <- all_results[[key]]

  plot_file <- file.path(out_base, "plots",
                         sprintf("%s_%s_%s.pdf",
                                 motif_id, res$preset_name, res$bw_label))

  p <- make_quadrant_plot(
    profiles_dt      = res$profiles,
    counts           = res$counts,
    pwm_matrix       = pwm,
    motif_id         = motif_id,
    bw_label         = res$bw_label,
    filter_label     = res$preset_name,
    motif_half_width = motif_half_width,
    tf_subtitle      = tf_subtitle,
    snp_stats        = res$snp_stats,
    shared_ylim      = shared_ylim
  )

  ggsave(plot_file, p, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
  message(sprintf("  Saved plot: %s", basename(plot_file)))

  png_file <- sub("\\.pdf$", ".png", plot_file)
  ggsave(png_file, p, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
}


# =============================================================================
# OVERLAY PLOTS (pairs of presets on same axes)
# =============================================================================

message("\n\n========== Generating overlay plots ==========\n")

for (pair in OVERLAY_PAIRS) {
  p1 <- pair[1]
  p2 <- pair[2]

  for (cache_key in names(overlay_cache)) {
    cache <- overlay_cache[[cache_key]]

    if (!p1 %in% names(cache) || !p2 %in% names(cache)) {
      message(sprintf("  Skipping overlay %s vs %s for %s (data not available)",
                      p1, p2, cache_key))
      next
    }

    profiles_list <- list()
    profiles_list[[p1]] <- cache[[p1]]$profiles
    profiles_list[[p2]] <- cache[[p2]]$profiles

    counts_list <- list()
    counts_list[[p1]] <- cache[[p1]]$counts
    counts_list[[p2]] <- cache[[p2]]$counts

    snp_stats_list <- list()
    snp_stats_list[[p1]] <- cache[[p1]]$snp_stats
    snp_stats_list[[p2]] <- cache[[p2]]$snp_stats

    p <- make_overlay_quadrant(
      profiles_list    = profiles_list,
      counts_list      = counts_list,
      pwm_matrix       = pwm,
      motif_id         = motif_id,
      bw_label         = cache_key,
      preset_names     = c(p1, p2),
      motif_half_width = motif_half_width,
      tf_subtitle      = tf_subtitle,
      snp_stats_list   = snp_stats_list,
      shared_ylim      = shared_ylim,
      # Colors and labels: users can customize via OVERLAY_COLORS in config.R
      # Format: named vector c("preset1" = "#color1", "preset2" = "#color2")
      colors           = OVERLAY_COLORS
    )

    overlay_file <- file.path(out_base, "plots_overlay",
                              sprintf("%s_%s_vs_%s_%s.pdf",
                                      motif_id, p1, p2, cache_key))
    ggsave(overlay_file, p, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)

    png_file <- sub("\\.pdf$", ".png", overlay_file)
    ggsave(png_file, p, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)

    message(sprintf("  Saved overlay: %s", basename(overlay_file)))
  }
}


# =============================================================================
# WRITE SKIP NOTES (for motifs with too few instances)
# =============================================================================
# This file lets you quickly check which motifs were skipped due to
# insufficient data vs. actual script failures.

if (length(skipped_presets) > 0) {
  skip_note_file <- file.path(out_base, "SKIPPED_PRESETS.txt")
  skip_lines <- c(
    sprintf("Motif: %s", motif_id),
    sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("MIN_INSTANCES threshold: %d", MIN_INSTANCES),
    ""
  )
  for (pn in names(skipped_presets)) {
    info <- skipped_presets[[pn]]
    skip_lines <- c(skip_lines,
      sprintf("Preset '%s': SKIPPED — %s (found %d, need >= %d)",
              pn, info$reason, info$n_found, info$min_required)
    )
  }
  writeLines(skip_lines, skip_note_file)
  message(sprintf("\nWrote skip note: %s", skip_note_file))
}

# Also write a SUCCESS marker if at least some plots were generated
if (length(all_results) > 0) {
  writeLines(
    c(sprintf("Motif: %s", motif_id),
      sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      sprintf("Plots generated: %d", length(all_results)),
      sprintf("Presets skipped: %d", length(skipped_presets))),
    file.path(out_base, "COMPLETED.txt")
  )
} else if (length(skipped_presets) == length(presets)) {
  # ALL presets were skipped — write a clear note
  writeLines(
    c(sprintf("Motif: %s", motif_id),
      sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      "STATUS: ALL presets skipped (too few instances for every preset)",
      sprintf("MIN_INSTANCES threshold: %d", MIN_INSTANCES)),
    file.path(out_base, "ALL_SKIPPED.txt")
  )
  message(sprintf("\nAll presets skipped for %s — wrote ALL_SKIPPED.txt", motif_id))
}


# =============================================================================
# SUMMARY
# =============================================================================

message(sprintf("\n========== Done processing motif %s ==========", motif_id))
message(sprintf("Output directory: %s", out_base))
message("Contents:")
message(sprintf("  plots/           — individual preset × bigwig plots"))
message(sprintf("  plots_overlay/   — overlay comparison plots"))
message(sprintf("  data/            — profiles (.tsv.gz), counts (.tsv), matrices (.rds)"))
message(sprintf("  rederived_pwms/  — re-derived PWMs from filtered subsets"))
if (length(skipped_presets) > 0) {
  message(sprintf("  SKIPPED_PRESETS.txt — presets with too few instances"))
}