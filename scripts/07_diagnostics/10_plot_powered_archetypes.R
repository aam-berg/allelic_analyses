#!/usr/bin/env Rscript
# =============================================================================
# 10_plot_powered_archetypes.R — Joint power filter (critical diagnostic)
# =============================================================================
#
# This is plot (ii) re-applied with an additional joint coverage filter:
# count, per archetype, the number of TFBSs that pass the full cascade AND
# have at least one paired het-SNP with sufficient allele-specific reads in
# each assay (ATAC, PRO-seq, RNA-seq). This tells us how many archetypes
# we can actually power for downstream analyses.
#
# Reads:
#   - 06_allele_pairs/pair_tables/{motif_id}_pair_table.tsv.gz
#     (uses pause_min_reads_per_allele, atac_ref/alt_total, etc.)
#   - tables/qualifying_per_archetype.tsv (for cascade context)
#
# Three filters applied (toggleable via thresholds in config.sh):
#   ATAC powered:    min(atac_ref_total, atac_alt_total) >= POWER_ATAC_MIN_PER_ALLELE
#   PRO-seq powered: pause_min_reads_per_allele >= POWER_PROSEQ_PAUSE_MIN_PER_ALLELE
#   RNA-seq powered: at least one junction cell with both ref_total + alt_total
#                    >= POWER_RNASEQ_AT_JUNCTION_MIN_PER_ALLELE
#
# For each combination of {ATAC-only, PRO-seq-only, ALL three}, plot a
# histogram in the same aesthetic as plot (ii).
#
# USAGE:
#   Rscript 10_plot_powered_archetypes.R \
#     --pair_table_dir .../pair_tables/ \
#     --qualifying_table .../tables/qualifying_per_archetype.tsv \
#     --plots_dir .../plots/ \
#     --tables_dir .../tables/ \
#     --power_atac_min 10 \
#     --power_proseq_pause_min 5 \
#     --power_rnaseq_min 4 \
#     --min_motif_width 8 \
#     --pwm_meme .../consensus_pwms.meme
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(ggplot2)
    library(scales)
    library(patchwork)
})

option_list <- list(
    make_option("--pair_table_dir",         type = "character"),
    make_option("--qualifying_table",        type = "character"),
    make_option("--plots_dir",               type = "character"),
    make_option("--tables_dir",              type = "character"),
    make_option("--power_atac_min",          type = "integer", default = 10),
    make_option("--power_proseq_pause_min",  type = "integer", default = 5),
    make_option("--power_rnaseq_min",        type = "integer", default = 4),
    make_option("--min_motif_width",         type = "integer", default = 8),
    make_option("--pwm_meme",                type = "character"),
    make_option("--lib_path",                type = "character", default = "")
)
opt <- parse_args(OptionParser(option_list = option_list))

lib_path <- if (nchar(opt$lib_path) > 0) opt$lib_path else
            file.path(dirname(sub("--file=", "",
                grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))),
                "lib_diagnostics.R")
source(lib_path)

cat("============================================================\n")
cat("10_plot_powered_archetypes.R\n")
cat("============================================================\n\n")

dir.create(opt$plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$tables_dir, recursive = TRUE, showWarnings = FALSE)

# Width filter
widths <- get_archetype_pwm_widths(opt$pwm_meme)
keep_motifs <- widths[motif_width > opt$min_motif_width, motif_id]
cat("Archetypes width >", opt$min_motif_width, ":", length(keep_motifs), "\n")

# Pair tables
files <- list.files(opt$pair_table_dir, pattern = "_pair_table\\.tsv\\.gz$",
                     full.names = TRUE)
keep_files <- files[basename(files) %in% paste0(keep_motifs, "_pair_table.tsv.gz")]
cat("Pair tables:", length(keep_files), "\n")
if (length(keep_files) == 0) {
    cat("[SKIP] No pair tables found.\n")
    quit(save = "no", status = 0)
}

# Columns to read; not all may exist depending on pipeline status
needed_cols <- c("motif_id", "motif_hit_id", "pair_id",
                  "atac_ref_total", "atac_alt_total",
                  "pause_min_reads_per_allele")
# Plus per-junction-cell ref/alt totals (we don't know exactly which cells
# are present, so read everything matching the pattern)

read_pair_for_power <- function(f) {
    # Read header first, then select what we need
    h <- names(fread(f, nrows = 0))
    junction_cols <- grep("^(donor|acceptor)_(sense|antisense)_(upstream|downstream)_(ref|alt)_total$",
                           h, value = TRUE)
    sel <- intersect(c(needed_cols, junction_cols), h)
    fread(f, select = sel, fill = TRUE)
}

cat("Reading pair tables (only needed columns)...\n")
all_pairs <- rbindlist(lapply(keep_files, function(f) {
    tryCatch(read_pair_for_power(f), error = function(e) {
        warning("Skipping ", f, ": ", conditionMessage(e))
        data.table()
    })
}), fill = TRUE)
cat("Total pairs:", nrow(all_pairs), "\n")
if (nrow(all_pairs) == 0) {
    cat("[SKIP] No pair data.\n")
    quit(save = "no", status = 0)
}

# -----------------------------------------------------------------------------
# Per-pair powered flags
# -----------------------------------------------------------------------------
all_pairs[, atac_min_per_allele :=
             pmin(as.numeric(atac_ref_total), as.numeric(atac_alt_total),
                   na.rm = FALSE)]
all_pairs[, atac_powered := !is.na(atac_min_per_allele) &
                              atac_min_per_allele >= opt$power_atac_min]

all_pairs[, proseq_powered :=
             !is.na(pause_min_reads_per_allele) &
             pause_min_reads_per_allele >= opt$power_proseq_pause_min]

# RNA-seq: any junction cell with both ref_total and alt_total >= threshold
junction_cells <- unique(sub("_(ref|alt)_total$", "",
    grep("_(ref|alt)_total$", names(all_pairs), value = TRUE)))
junction_cells <- setdiff(junction_cells, "atac")  # exclude atac

if (length(junction_cells) > 0) {
    rnaseq_powered <- rep(FALSE, nrow(all_pairs))
    for (cell in junction_cells) {
        ref_col <- paste0(cell, "_ref_total")
        alt_col <- paste0(cell, "_alt_total")
        if (ref_col %in% names(all_pairs) && alt_col %in% names(all_pairs)) {
            ref_v <- as.numeric(all_pairs[[ref_col]])
            alt_v <- as.numeric(all_pairs[[alt_col]])
            cell_pow <- !is.na(ref_v) & !is.na(alt_v) &
                         ref_v >= opt$power_rnaseq_min &
                         alt_v >= opt$power_rnaseq_min
            rnaseq_powered <- rnaseq_powered | cell_pow
        }
    }
    all_pairs[, rnaseq_powered := rnaseq_powered]
} else {
    all_pairs[, rnaseq_powered := NA]
}

cat("Powered fractions:\n")
cat("  ATAC:    ",   sprintf("%.1f%%", 100 * mean(all_pairs$atac_powered, na.rm = TRUE)), "\n")
cat("  PRO-seq: ",   sprintf("%.1f%%", 100 * mean(all_pairs$proseq_powered, na.rm = TRUE)), "\n")
cat("  RNA-seq: ",   sprintf("%.1f%%", 100 * mean(all_pairs$rnaseq_powered, na.rm = TRUE)), "\n")

# -----------------------------------------------------------------------------
# Per-archetype counts: how many UNIQUE motif_hit_ids have at least one
# powered pair under each definition?
# -----------------------------------------------------------------------------
power_defs <- list(
    atac_only      = quote(atac_powered),
    proseq_only    = quote(proseq_powered),
    rnaseq_only    = quote(rnaseq_powered),
    atac_and_proseq = quote(atac_powered & proseq_powered),
    all_three      = quote(atac_powered & proseq_powered & rnaseq_powered)
)

per_arch_list <- list()
for (def_name in names(power_defs)) {
    expr <- power_defs[[def_name]]
    powered <- all_pairs[eval(expr) & !is.na(eval(expr))]
    if (nrow(powered) == 0) {
        per_arch_list[[def_name]] <- data.table(motif_id = character(0),
                                                  n_powered_motif_hits = integer(0),
                                                  power_def = def_name)
        next
    }
    counts <- powered[, .(n_powered_motif_hits = uniqueN(motif_hit_id)),
                       by = motif_id]
    counts[, power_def := def_name]
    per_arch_list[[def_name]] <- counts
}
per_arch <- rbindlist(per_arch_list, fill = TRUE)

# Save table
out_tab <- file.path(opt$tables_dir, "powered_per_archetype.tsv")
fwrite(per_arch, out_tab, sep = "\t", na = "NA")
cat("Wrote:", out_tab, "\n")

# -----------------------------------------------------------------------------
# Histograms (one per power_def), same aesthetic as plot (ii)
# -----------------------------------------------------------------------------
make_powered_hist <- function(d, def_name) {
    counts <- d[power_def == def_name]$n_powered_motif_hits
    counts <- counts[!is.na(counts) & counts > 0]
    n_arch <- length(counts)
    if (n_arch == 0) {
        return(ggplot() +
            annotate("text", x = 0.5, y = 0.5,
                     label = paste0(def_name, "\nNo archetypes powered"),
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

    p <- ggplot(data.table(count = counts), aes(x = count)) +
        geom_histogram(binwidth = bw, fill = "#54278f", color = "white",
                       linewidth = 0.3) +
        geom_vline(xintercept = med_val, linetype = "dashed",
                   color = "#d62728", linewidth = 0.5)
    if (nrow(cumul_dt) > 0) {
        p <- p + geom_vline(data = cumul_dt,
            aes(xintercept = threshold, color = line_color),
            linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
            scale_color_identity()
    }
    p <- p +
        scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
        labs(title = sprintf("Powered archetypes: %s", def_name),
             subtitle = sprintf("%d archetypes ≥1 powered TFBS; total = %s pairs",
                                  n_arch, comma(sum(counts))),
             x = "Number of powered motif hits per archetype",
             y = "Number of archetypes") +
        DIAG_THEME

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
    pb <- ggplot_build(p)
    y_max <- max(pb$data[[1]]$count)
    line_spacing <- y_max * 0.07
    y_start <- y_max * 0.95
    n_lines <- length(legend_lines)
    p <- p + annotate("rect",
        xmin = max_val * 0.55, xmax = max_val * 1.0,
        ymin = y_start - n_lines * line_spacing - line_spacing * 0.4,
        ymax = y_start + line_spacing * 0.5,
        fill = "white", alpha = 0.85, color = "grey70", linewidth = 0.4)
    for (i in seq_len(n_lines)) {
        p <- p + annotate("text",
            x = max_val * 0.97,
            y = y_start - (i - 1) * line_spacing,
            label = legend_lines[i], hjust = 1, vjust = 0.5,
            size = 3.2, fontface = "bold", color = legend_colors[i])
    }
    p
}

for (def_name in names(power_defs)) {
    p <- make_powered_hist(per_arch, def_name)
    save_diag_plot(p,
        file.path(opt$plots_dir,
            sprintf("powered_archetypes_%s.pdf", def_name)),
        width = 7.5, height = 5)
    cat("  Wrote: powered_archetypes_", def_name, ".pdf\n", sep = "")
}

cat("\nDone.\n")
