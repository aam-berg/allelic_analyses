#!/usr/bin/env Rscript
# =============================================================================
# 05_plot_top_archetypes.R — Plot (iii): top archetypes panel
# =============================================================================
#
# For top-N archetypes (by qualifying-TFBS count), produces a panel:
#   - Horizontal bar = qualifying-TFBS count
#   - Bar color = log2(max TPM across mapped TFs + 1)
#   - To the right of each bar: ggseqlogo of the archetype's PWM
#   - Right-most label: the TF name with the max TPM
#
# Uses qualifying definition = "overlap" (direct body overlap) by default,
# toggleable.
#
# USAGE:
#   Rscript 05_plot_top_archetypes.R \
#     --qualifying_table .../tables/qualifying_per_archetype.tsv \
#     --pwm_meme .../consensus_pwms.meme \
#     --plots_dir .../plots/ \
#     --top_n_list 10,25,50 \
#     --min_motif_width 8 \
#     --qualifying_col n_qualifying_overlap
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(ggplot2)
    library(scales)
    library(patchwork)
    library(ggseqlogo)
})

option_list <- list(
    make_option("--qualifying_table",  type = "character"),
    make_option("--pwm_meme",          type = "character"),
    make_option("--plots_dir",         type = "character"),
    make_option("--top_n_list",        type = "character", default = "10,25,50"),
    make_option("--min_motif_width",   type = "integer", default = 8),
    make_option("--qualifying_col",    type = "character", default = "n_qualifying_overlap"),
    make_option("--lib_path",          type = "character", default = "")
)
opt <- parse_args(OptionParser(option_list = option_list))

lib_path <- if (nchar(opt$lib_path) > 0) opt$lib_path else
            file.path(dirname(sub("--file=", "",
                grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))),
                "lib_diagnostics.R")
source(lib_path)

cat("============================================================\n")
cat("05_plot_top_archetypes.R\n")
cat("============================================================\n\n")

dir.create(opt$plots_dir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(opt$qualifying_table)
pwms <- parse_meme_pwms(opt$pwm_meme)
cat("Loaded archetypes:", nrow(dt), " | PWMs:", length(pwms), "\n")

# Apply width filter
dt <- dt[!is.na(motif_width) & motif_width > opt$min_motif_width]
cat("After width filter (>", opt$min_motif_width, "):", nrow(dt), "\n")

q_col <- opt$qualifying_col
if (!q_col %in% names(dt)) {
    stop("Qualifying column not in table: ", q_col)
}

# Sort by qualifying count (descending)
dt <- dt[!is.na(get(q_col)) & get(q_col) > 0]
setorderv(dt, q_col, order = -1L)
cat("Archetypes with ≥1 qualifying TFBS:", nrow(dt), "\n")

if (nrow(dt) == 0) {
    cat("No archetypes to plot. Exiting.\n")
    quit(save = "no", status = 0)
}

top_ns <- as.integer(strsplit(opt$top_n_list, ",")[[1]])

# -----------------------------------------------------------------------------
# Make a single plot for top-N
# -----------------------------------------------------------------------------
make_top_panel <- function(d, n, q_col_local) {
    sub <- head(d, n)
    sub[, motif_id_factor := factor(motif_id, levels = rev(motif_id))]
    sub[, log2_tpm := log2(pmax(max_tpm_across_mapped, 0.01) + 1)]

    # Bar plot
    p_bars <- ggplot(sub, aes(y = motif_id_factor, x = get(q_col_local),
                                fill = log2_tpm)) +
        geom_col(width = 0.8) +
        geom_text(aes(label = comma(get(q_col_local))),
                  hjust = -0.15, size = 2.8) +
        geom_text(aes(label = ifelse(is.na(max_tf_name), "(no TF)",
                                       sprintf("%s (TPM=%.1f)",
                                               max_tf_name, max_tpm_across_mapped))),
                  x = 0, hjust = 1.05, size = 2.6, color = "grey25") +
        scale_x_continuous(labels = comma, expand = expansion(mult = c(0.50, 0.20))) +
        scale_fill_gradient(low = "#fee5d9", high = "#a50f15",
                            na.value = "grey70",
                            name = "log2(max\nTPM + 1)") +
        labs(title = sprintf("Top %d archetypes by qualifying TFBS count", n),
             subtitle = sprintf("Qualifying = %s; archetypes width > %d bp",
                                  sub("^n_qualifying_", "", q_col_local),
                                  opt$min_motif_width),
             x = "Qualifying TFBSs", y = NULL) +
        DIAG_THEME +
        theme(panel.grid.major.y = element_blank())

    # Logo for each archetype, right of bar — use ggseqlogo facet
    # Build a named list of PWM matrices in the order of motif_id (top to bottom)
    motif_order <- as.character(sub$motif_id)
    logos_data <- list()
    for (mid in motif_order) {
        if (!is.null(pwms[[mid]])) {
            mat <- t(pwms[[mid]]$ppm)   # ggseqlogo wants rows = bases
            rownames(mat) <- c("A", "C", "G", "T")
            logos_data[[mid]] <- mat
        }
    }
    if (length(logos_data) > 0) {
        p_logos <- ggseqlogo(logos_data, ncol = 1, method = "bits") +
            facet_wrap(~seq_group, ncol = 1, strip.position = "left") +
            theme(strip.text.y.left = element_text(angle = 0, size = 7),
                  axis.text.y = element_blank(),
                  axis.text.x = element_blank(),
                  axis.title  = element_blank(),
                  panel.spacing = unit(0.05, "lines"),
                  strip.background = element_blank())
        # Reorder facets to match bars (top-down)
        # ggseqlogo's facet order is fixed by list order — we constructed it
        # in motif_order (descending qualifying count, so top motif first)
        # which matches the bar order top-to-bottom (since y = factor with
        # levels = rev(motif_id)).
        combined <- p_bars + p_logos + plot_layout(widths = c(2.3, 1))
    } else {
        combined <- p_bars
    }
    combined
}

for (n in top_ns) {
    if (n > nrow(dt)) {
        cat("Top-", n, ": only ", nrow(dt), " archetypes available.\n", sep = "")
        n <- nrow(dt)
    }
    p <- make_top_panel(dt, n, q_col)
    h <- max(4, 0.25 * n + 1.5)
    save_diag_plot(p,
        file.path(opt$plots_dir,
            sprintf("top_archetypes_n%d_%s.pdf", n,
                    sub("^n_qualifying_", "", q_col))),
        width = 11, height = h)
    cat("  Wrote: top_archetypes_n", n, "_*.pdf\n", sep = "")
}

cat("\nDone.\n")
