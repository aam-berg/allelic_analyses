#!/usr/bin/env Rscript
# =============================================================================
# 06_plot_pwm_score_dists.R — Plot (iv): PWM score distributions
# =============================================================================
#
# Two variants:
#   (a) Per-archetype faceted violins for top-N archetypes:
#       MOODS PWM score distribution, comparing all-TFBSs vs qualifying-TFBSs.
#       Tests whether SNP-overlapping TFBSs are enriched/depleted at strong
#       motif positions (would indicate purifying selection on motif strength).
#   (b) Pooled global panel:
#       Within-archetype z-normalized scores; pooled across archetypes; same
#       all vs qualifying comparison. Cross-archetype-comparable.
#
# USAGE:
#   Rscript 06_plot_pwm_score_dists.R \
#     --summary_dir .../per_archetype_summaries/ \
#     --qualifying_table .../tables/qualifying_per_archetype.tsv \
#     --plots_dir .../plots/ \
#     --top_n 25 \
#     --min_motif_width 8 \
#     --qualifying_col qualifying_overlap
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(ggplot2)
    library(scales)
})

option_list <- list(
    make_option("--summary_dir",       type = "character"),
    make_option("--qualifying_table",  type = "character"),
    make_option("--plots_dir",         type = "character"),
    make_option("--top_n",             type = "integer", default = 25),
    make_option("--min_motif_width",   type = "integer", default = 8),
    make_option("--qualifying_col",    type = "character", default = "qualifying_overlap"),
    make_option("--lib_path",          type = "character", default = "")
)
opt <- parse_args(OptionParser(option_list = option_list))

lib_path <- if (nchar(opt$lib_path) > 0) opt$lib_path else
            file.path(dirname(sub("--file=", "",
                grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))),
                "lib_diagnostics.R")
source(lib_path)

cat("============================================================\n")
cat("06_plot_pwm_score_dists.R\n")
cat("============================================================\n\n")

dir.create(opt$plots_dir, recursive = TRUE, showWarnings = FALSE)

# Determine top archetypes by n_qualifying_<thresh> matching opt$qualifying_col
qual_dt <- fread(opt$qualifying_table)
qual_dt <- qual_dt[!is.na(motif_width) & motif_width > opt$min_motif_width]
n_col <- paste0("n_", opt$qualifying_col)
if (!n_col %in% names(qual_dt)) {
    stop("Column not found in qualifying table: ", n_col)
}
setorderv(qual_dt, n_col, order = -1L)
top_motifs <- head(qual_dt$motif_id[!is.na(qual_dt[[n_col]]) & qual_dt[[n_col]] > 0],
                    opt$top_n)
cat("Top", length(top_motifs), "archetypes selected\n")

# -----------------------------------------------------------------------------
# Read per-TFBS scores for top motifs (and a sample for global plot)
# -----------------------------------------------------------------------------
read_scores <- function(motif_id, summary_dir) {
    p <- file.path(summary_dir, paste0(motif_id, "_scores.tsv.gz"))
    if (!file.exists(p)) return(data.table())
    fread(p)
}

cat("Reading scores for top motifs...\n")
top_scores_list <- lapply(top_motifs, read_scores, summary_dir = opt$summary_dir)
top_scores <- rbindlist(top_scores_list, fill = TRUE)

if (nrow(top_scores) == 0) {
    cat("No score files found in ", opt$summary_dir, "\n")
    quit(save = "no", status = 0)
}

q_flag <- opt$qualifying_col
if (!q_flag %in% names(top_scores)) {
    stop("Qualifying flag column not in scores: ", q_flag)
}
top_scores[, group := ifelse(get(q_flag), "qualifying", "all_other")]
top_scores[, motif_id := factor(motif_id, levels = top_motifs)]

# -----------------------------------------------------------------------------
# Plot (iv.a) — per-archetype faceted violins for top motifs
# -----------------------------------------------------------------------------
# Build a "all" vs "qualifying" duplication: every TFBS appears once in "all",
# and qualifying ones also appear in "qualifying" so we get true "all" baseline
plot_dt <- rbind(
    copy(top_scores)[, group := "all"],
    top_scores[get(q_flag) == TRUE][, group := "qualifying"]
)
plot_dt[, group := factor(group, levels = c("all", "qualifying"))]

p_facets <- ggplot(plot_dt, aes(x = group, y = score, fill = group)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.3, alpha = 0.85) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white",
                  linewidth = 0.3) +
    scale_fill_manual(values = c(all = "#bdbdbd", qualifying = "#08519c"),
                       guide = "none") +
    facet_wrap(~ motif_id, scales = "free_y", ncol = 5) +
    labs(title = sprintf("PWM score distributions: top %d archetypes",
                            length(top_motifs)),
         subtitle = sprintf("Qualifying = %s; archetypes width > %d bp",
                              opt$qualifying_col, opt$min_motif_width),
         x = NULL, y = "MOODS score") +
    DIAG_THEME +
    theme(strip.text = element_text(size = 8, face = "bold"),
          axis.text.x = element_text(size = 8))

n_rows <- ceiling(length(top_motifs) / 5)
save_diag_plot(p_facets,
    file.path(opt$plots_dir,
        sprintf("pwm_score_dists_per_archetype_top%d.pdf", length(top_motifs))),
    width = 14, height = max(4, 1.8 * n_rows))

# -----------------------------------------------------------------------------
# Plot (iv.b) — global pooled with within-archetype z-normalization
# -----------------------------------------------------------------------------
# For ALL archetypes (not just top), z-normalize within archetype,
# then pool. Compare "all" vs "qualifying".
cat("\nReading scores for all archetypes...\n")
all_files <- list.files(opt$summary_dir, pattern = "_scores\\.tsv\\.gz$",
                         full.names = TRUE)
# Filter to width > min_motif_width archetypes
keep_motifs <- qual_dt$motif_id  # already filtered by width
keep_files <- all_files[basename(all_files) %in% paste0(keep_motifs, "_scores.tsv.gz")]
cat("Reading", length(keep_files), "score files...\n")

if (length(keep_files) > 0) {
    all_scores <- rbindlist(lapply(keep_files, fread), fill = TRUE)
    cat("Total TFBSs:", nrow(all_scores), "\n")

    if (q_flag %in% names(all_scores)) {
        all_scores[, score_z := (score - mean(score, na.rm = TRUE)) /
                                  sd(score, na.rm = TRUE), by = motif_id]
        # Build all + qualifying duplication
        all_z <- rbind(
            copy(all_scores)[, group := "all"],
            all_scores[get(q_flag) == TRUE][, group := "qualifying"]
        )
        all_z[, group := factor(group, levels = c("all", "qualifying"))]
        all_z <- all_z[!is.na(score_z) & is.finite(score_z)]

        # Plot: density of z-scores
        p_global <- ggplot(all_z, aes(x = score_z, fill = group, color = group)) +
            geom_density(alpha = 0.4, linewidth = 0.7) +
            scale_fill_manual(values = c(all = "#bdbdbd", qualifying = "#08519c"),
                               name = NULL) +
            scale_color_manual(values = c(all = "grey40", qualifying = "#08306b"),
                                guide = "none") +
            labs(title = "PWM score distributions (within-archetype z-normalized)",
                 subtitle = sprintf("Pooled across %d archetypes (width > %d bp). Qualifying = %s.",
                                      length(keep_files),
                                      opt$min_motif_width, opt$qualifying_col),
                 x = "Within-archetype z(score)",
                 y = "Density") +
            geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
            DIAG_THEME

        save_diag_plot(p_global,
            file.path(opt$plots_dir, "pwm_score_dists_global_zscore.pdf"),
            width = 8, height = 5)

        # Print summary stats
        cat("\nGlobal z-score summary:\n")
        cat("  All TFBSs:        mean z = 0 by construction; n =",
            sum(all_z$group == "all"), "\n")
        med_q <- median(all_z$score_z[all_z$group == "qualifying"], na.rm = TRUE)
        cat("  Qualifying TFBSs: median z =", round(med_q, 3),
            "; n =", sum(all_z$group == "qualifying"), "\n")
        if (med_q < -0.05) {
            cat("  -> Qualifying TFBSs are ENRICHED at LOWER scores (purifying selection sign).\n")
        } else if (med_q > 0.05) {
            cat("  -> Qualifying TFBSs are ENRICHED at HIGHER scores (positive selection sign).\n")
        } else {
            cat("  -> No clear enrichment.\n")
        }
    }
}

cat("\nDone.\n")
