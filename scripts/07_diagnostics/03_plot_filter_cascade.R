#!/usr/bin/env Rscript
# =============================================================================
# 03_plot_filter_cascade.R — Plot (i): the filtering cascade
# =============================================================================
#
# Reads the aggregated archetype_summary.tsv and produces:
#
#   plots/cascade_n_tfbs_all.pdf            n TFBSs at each cascade stage, all archetypes
#   plots/cascade_n_tfbs_widthgt8.pdf        n TFBSs, archetypes width > 8 bp only
#   plots/cascade_unique_bp_all.pdf          unique bp covered, all archetypes
#   plots/cascade_unique_bp_widthgt8.pdf     unique bp covered, width > 8 only
#   plots/cascade_per_archetype_heatmap.pdf  survival fraction per archetype × stage
#
# USAGE:
#   Rscript 03_plot_filter_cascade.R \
#     --archetype_summary .../tables/archetype_summary.tsv \
#     --plots_dir .../plots/ \
#     --min_motif_width 8
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(ggplot2)
    library(scales)
})

option_list <- list(
    make_option("--archetype_summary", type = "character"),
    make_option("--plots_dir",          type = "character"),
    make_option("--min_motif_width",    type = "integer", default = 8),
    make_option("--lib_path",           type = "character", default = "")
)
opt <- parse_args(OptionParser(option_list = option_list))

lib_path <- if (nchar(opt$lib_path) > 0) opt$lib_path else
            file.path(dirname(sub("--file=", "",
                grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))),
                "lib_diagnostics.R")
source(lib_path)

cat("============================================================\n")
cat("03_plot_filter_cascade.R\n")
cat("============================================================\n\n")

dir.create(opt$plots_dir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(opt$archetype_summary)
cat("Loaded archetype summary:", nrow(dt), "rows\n")

# Identify cascade stages from column names (n_<stage>)
stage_cols <- grep("^n_", names(dt), value = TRUE)
# Exclude n_qualifying_* and n_mapped_tfs / n_with_tpm
stage_cols <- setdiff(stage_cols, c("n_mapped_tfs", "n_with_tpm"))
stage_cols <- stage_cols[!grepl("^n_qualifying_", stage_cols)]
stage_labels <- sub("^n_", "", stage_cols)
cat("Cascade stages:", paste(stage_labels, collapse = " -> "), "\n\n")

# Pretty labels for plotting
pretty_labels <- c(
    all              = "All scanned",
    standard_chroms  = "Standard chrom",
    intragenic_sense = "Intragenic",
    pc_or_lncrna     = "PC / lncRNA",
    expressed        = "Expressed (TPM≥1)",
    not_promoter     = "Not promoter",
    atac_accessible  = "ATAC-accessible"
)

# -----------------------------------------------------------------------------
# Pivot to long format
# -----------------------------------------------------------------------------
build_long <- function(d) {
    n_long <- melt(d, id.vars = c("motif_id", "motif_width"),
                    measure.vars = stage_cols,
                    variable.name = "stage_col", value.name = "n_tfbs")
    n_long[, stage := sub("^n_", "", as.character(stage_col))]

    ubp_cols <- paste0("ubp_", stage_labels)
    ubp_cols_present <- intersect(ubp_cols, names(d))
    if (length(ubp_cols_present) > 0) {
        u_long <- melt(d, id.vars = c("motif_id", "motif_width"),
                        measure.vars = ubp_cols_present,
                        variable.name = "stage_col", value.name = "ubp")
        u_long[, stage := sub("^ubp_", "", as.character(stage_col))]
        n_long[u_long[, .(motif_id, stage, ubp)],
                ubp := i.ubp, on = c("motif_id", "stage")]
    }
    n_long[, stage := factor(stage, levels = stage_labels)]
    n_long[, stage_pretty := factor(pretty_labels[as.character(stage)],
                                      levels = pretty_labels[stage_labels])]
    n_long
}

# -----------------------------------------------------------------------------
# Plot helper: stacked / aggregated bar at each cascade stage
# -----------------------------------------------------------------------------
plot_cascade <- function(long_dt, metric, title, subtitle) {
    agg <- long_dt[, .(total = sum(get(metric), na.rm = TRUE)), by = stage_pretty]
    if (sum(agg$total) == 0) {
        return(ggplot() +
            annotate("text", x = 0.5, y = 0.5, label = "No data") +
            theme_void())
    }

    # Drop ratio for context
    agg[, frac_of_first := total / total[1]]

    ggplot(agg, aes(x = stage_pretty, y = total)) +
        geom_col(fill = "#4a90d9", color = "white", linewidth = 0.3) +
        geom_text(aes(label = paste0(comma(total), "\n(",
                                       sprintf("%.1f%%", 100 * frac_of_first), ")")),
                  vjust = -0.3, size = 3, lineheight = 0.85) +
        scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.18))) +
        labs(title = title, subtitle = subtitle,
             x = NULL,
             y = if (metric == "n_tfbs") "Number of TFBSs"
                  else "Unique bp covered") +
        DIAG_THEME +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# -----------------------------------------------------------------------------
# Generate the four cascade plots: 2 archetype-sets × 2 metrics
# -----------------------------------------------------------------------------
archetype_subsets <- list(
    list(label = "all",
         dt = dt,
         desc = sprintf("All archetypes (n=%d)", nrow(dt))),
    list(label = "widthgt8",
         dt = dt[!is.na(motif_width) & motif_width > opt$min_motif_width],
         desc = sprintf("Archetypes width > %d bp (n=%d)",
                          opt$min_motif_width,
                          nrow(dt[!is.na(motif_width) & motif_width > opt$min_motif_width])))
)

for (subset in archetype_subsets) {
    if (nrow(subset$dt) == 0) {
        cat("Skipping subset", subset$label, "(empty)\n")
        next
    }
    long <- build_long(subset$dt)

    # n_tfbs
    p1 <- plot_cascade(long, "n_tfbs",
        title = sprintf("Filter cascade — TFBS counts (%s)", subset$label),
        subtitle = subset$desc)
    save_diag_plot(p1, file.path(opt$plots_dir,
        sprintf("cascade_n_tfbs_%s.pdf", subset$label)),
        width = 9, height = 5.5)

    # unique bp
    if (sum(long$ubp, na.rm = TRUE) > 0) {
        p2 <- plot_cascade(long, "ubp",
            title = sprintf("Filter cascade — Unique bp covered (%s)", subset$label),
            subtitle = subset$desc)
        save_diag_plot(p2, file.path(opt$plots_dir,
            sprintf("cascade_unique_bp_%s.pdf", subset$label)),
            width = 9, height = 5.5)
    }
}

# -----------------------------------------------------------------------------
# Per-archetype heatmap: survival fraction at each stage
# -----------------------------------------------------------------------------
sub_dt <- dt[!is.na(motif_width) & motif_width > opt$min_motif_width]
if (nrow(sub_dt) > 0) {
    long <- build_long(sub_dt)
    # Compute fraction of TFBSs surviving at each stage (vs n_all)
    long[, n_all := n_tfbs[stage == "all"], by = motif_id]
    long[, frac := ifelse(n_all > 0, n_tfbs / n_all, NA_real_)]

    # Order motifs by final-stage retention
    final_stage <- tail(stage_labels, 1)
    motif_order <- long[stage == final_stage][order(-frac), motif_id]
    long[, motif_id := factor(motif_id, levels = motif_order)]

    # Cap at top 100 most-retained for readability
    keep_motifs <- head(motif_order, 100)
    long_sub <- long[motif_id %in% keep_motifs]
    long_sub[, motif_id := factor(motif_id, levels = keep_motifs)]

    p3 <- ggplot(long_sub, aes(x = stage_pretty, y = motif_id, fill = frac)) +
        geom_tile(color = "white", linewidth = 0.1) +
        scale_fill_gradient(low = "#f7fbff", high = "#08519c",
                            limits = c(0, 1), na.value = "grey90",
                            name = "Survival\nfraction") +
        labs(title = "Per-archetype TFBS retention through cascade",
             subtitle = sprintf("Top 100 archetypes by final retention (width > %d bp)",
                                  opt$min_motif_width),
             x = NULL, y = NULL) +
        DIAG_THEME +
        theme(axis.text.y = element_text(size = 6),
              axis.text.x = element_text(angle = 30, hjust = 1))
    save_diag_plot(p3,
        file.path(opt$plots_dir, "cascade_per_archetype_heatmap.pdf"),
        width = 9, height = 14)
}

cat("\nDone. Wrote plots to ", opt$plots_dir, "\n")
