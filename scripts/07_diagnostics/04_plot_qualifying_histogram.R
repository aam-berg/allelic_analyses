#!/usr/bin/env Rscript
# =============================================================================
# 04_plot_qualifying_histogram.R — Plot (ii): qualifying TFBSs per archetype
# =============================================================================
#
# For each "qualifying" definition (overlap, flank10, ..., flank800), produces
# a histogram of n_qualifying_<thresh> across motif archetypes. Adds:
#   - Vertical red dashed line at the median
#   - Vertical green-gradient dashed lines at cumulative thresholds
#     (≥250, ≥500, ≥1000, ≥2000), with annotated counts
#   - ECDF panel below the histogram for direct read-off
#
# Default: only plot for archetypes with motif_width > 8 (toggleable).
# Aesthetic mirrors the user's notebook Cell 19 with minor refinements.
#
# USAGE:
#   Rscript 04_plot_qualifying_histogram.R \
#     --qualifying_table .../tables/qualifying_per_archetype.tsv \
#     --plots_dir .../plots/ \
#     --min_motif_width 8
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(ggplot2)
    library(scales)
    library(patchwork)
})

option_list <- list(
    make_option("--qualifying_table",  type = "character"),
    make_option("--plots_dir",         type = "character"),
    make_option("--min_motif_width",   type = "integer", default = 8),
    make_option("--lib_path",          type = "character", default = "")
)
opt <- parse_args(OptionParser(option_list = option_list))

lib_path <- if (nchar(opt$lib_path) > 0) opt$lib_path else
            file.path(dirname(sub("--file=", "",
                grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))),
                "lib_diagnostics.R")
source(lib_path)

cat("============================================================\n")
cat("04_plot_qualifying_histogram.R\n")
cat("============================================================\n\n")

dir.create(opt$plots_dir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(opt$qualifying_table)
cat("Loaded qualifying table:", nrow(dt), "rows\n")

q_cols <- grep("^n_qualifying_", names(dt), value = TRUE)

# Threshold labels for plotting
threshold_label_map <- c(
    "n_qualifying_overlap"    = "Direct body overlap",
    "n_qualifying_flank10"    = "Within 10 bp of motif edge",
    "n_qualifying_flank25"    = "Within 25 bp",
    "n_qualifying_flank50"    = "Within 50 bp",
    "n_qualifying_flank100"   = "Within 100 bp",
    "n_qualifying_flank200"   = "Within 200 bp",
    "n_qualifying_flank400"   = "Within 400 bp",
    "n_qualifying_flank800"   = "Within 800 bp"
)

# -----------------------------------------------------------------------------
# Plot helper for one threshold
# -----------------------------------------------------------------------------
make_qual_hist <- function(d, q_col, label, archetype_subset_label) {
    counts <- d[[q_col]]
    counts <- counts[!is.na(counts) & counts > 0]
    n_arch <- length(counts)
    if (n_arch == 0) {
        return(ggplot() +
            annotate("text", x = 0.5, y = 0.5,
                     label = paste0(label, "\nNo archetypes with ≥1 qualifying TFBS"),
                     size = 4) +
            theme_void())
    }
    max_val <- max(counts)
    med_val <- median(counts)
    bw <- if (max_val <= 50) 5 else if (max_val <= 200) 10 else
          if (max_val <= 500) 25 else if (max_val <= 2000) 50 else 100

    cumul_thresholds <- c(250, 500, 1000, 2000)
    cumul_dt <- data.table(
        threshold = cumul_thresholds,
        n_archetypes = sapply(cumul_thresholds, function(t) sum(counts >= t))
    )
    cumul_dt <- cumul_dt[threshold <= max_val & n_archetypes > 0]

    if (nrow(cumul_dt) > 0) {
        green_pal <- colorRampPalette(c("#a1d99b", "#006d2c"))(nrow(cumul_dt))
        cumul_dt[, line_color := green_pal]
    }

    p_hist <- ggplot(data.table(count = counts), aes(x = count)) +
        geom_histogram(binwidth = bw, fill = "#4a90d9", color = "white",
                       linewidth = 0.3) +
        geom_vline(xintercept = med_val, linetype = "dashed",
                   color = "#d62728", linewidth = 0.5)

    if (nrow(cumul_dt) > 0) {
        p_hist <- p_hist + geom_vline(data = cumul_dt,
            aes(xintercept = threshold, color = line_color),
            linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
            scale_color_identity()
    }

    # Build legend text
    legend_lines <- c(sprintf("Median = %d", round(med_val)))
    legend_colors <- c("#d62728")
    if (nrow(cumul_dt) > 0) {
        for (i in seq_len(nrow(cumul_dt))) {
            legend_lines <- c(legend_lines,
                sprintf("≥%s: %d archetypes",
                        comma(cumul_dt$threshold[i]),
                        cumul_dt$n_archetypes[i]))
            legend_colors <- c(legend_colors, cumul_dt$line_color[i])
        }
    }

    p_hist <- p_hist +
        scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
        labs(title = sprintf("Archetype yield: %s", label),
             subtitle = sprintf("%d archetypes ≥1 qualifying TFBS (%s); total = %s",
                                  n_arch, archetype_subset_label, comma(sum(counts))),
             x = "Number of qualifying TFBSs per archetype",
             y = "Number of archetypes") +
        DIAG_THEME

    # Annotate legend in upper-right
    p_hist_built <- ggplot_build(p_hist)
    y_max <- max(p_hist_built$data[[1]]$count)
    line_spacing <- y_max * 0.07
    y_start <- y_max * 0.95
    n_lines <- length(legend_lines)

    # Background rectangle
    p_hist <- p_hist + annotate("rect",
        xmin = max_val * 0.55, xmax = max_val * 1.0,
        ymin = y_start - n_lines * line_spacing - line_spacing * 0.4,
        ymax = y_start + line_spacing * 0.5,
        fill = "white", alpha = 0.85, color = "grey70", linewidth = 0.4)
    for (i in seq_len(n_lines)) {
        p_hist <- p_hist + annotate("text",
            x = max_val * 0.97,
            y = y_start - (i - 1) * line_spacing,
            label = legend_lines[i], hjust = 1, vjust = 0.5,
            size = 3.2, fontface = "bold", color = legend_colors[i])
    }

    # ECDF panel
    p_ecdf <- ggplot(data.table(count = counts), aes(x = count)) +
        stat_ecdf(geom = "step", color = "#08306b", linewidth = 0.7) +
        geom_vline(xintercept = med_val, linetype = "dashed",
                   color = "#d62728", linewidth = 0.4) +
        labs(x = "Number of qualifying TFBSs per archetype",
             y = "ECDF (fraction of archetypes ≤ x)") +
        DIAG_THEME +
        theme(plot.margin = margin(0, 5, 5, 5))

    # Combine via patchwork
    combined <- p_hist / p_ecdf + plot_layout(heights = c(3, 1))
    combined
}

# -----------------------------------------------------------------------------
# Apply width filter (default: motif_width > 8)
# -----------------------------------------------------------------------------
archetype_subsets <- list(
    all       = list(d = dt, label = "all archetypes"),
    widthgt8  = list(d = dt[!is.na(motif_width) & motif_width > opt$min_motif_width],
                      label = sprintf("width > %d bp", opt$min_motif_width))
)

for (subset_name in names(archetype_subsets)) {
    subset <- archetype_subsets[[subset_name]]
    cat("\nSubset:", subset_name, "—", nrow(subset$d), "archetypes\n")
    for (q_col in q_cols) {
        if (!q_col %in% names(subset$d)) next
        label <- threshold_label_map[q_col]
        if (is.na(label)) label <- q_col

        p <- make_qual_hist(subset$d, q_col, label, subset$label)
        out_path <- file.path(opt$plots_dir,
            sprintf("qualifying_hist_%s_%s.pdf",
                    sub("^n_qualifying_", "", q_col),
                    subset_name))
        save_diag_plot(p, out_path, width = 7.5, height = 6)
        cat("  Wrote:", basename(out_path), "\n")
    }
}

cat("\nDone.\n")
