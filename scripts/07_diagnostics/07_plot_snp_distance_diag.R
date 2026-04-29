#!/usr/bin/env Rscript
# =============================================================================
# 07_plot_snp_distance_diag.R — ATAC coupling distance diagnostic
# =============================================================================
#
# For each (motif, SNP) pair from 06_allele_pairs/ pair tables, compute:
#   abs_log2_atac_ratio = |log2((atac_ref_total + 1) / (atac_alt_total + 1))|
#   abs_snp_dist = |snp_distance_from_motif_edge_bp|  (0 if SNP in body)
#
# Bin pairs by abs_snp_dist into [0, 1-25, 25-50, 50-100, 100-200, 200-400, 400-800].
# Plot:
#   (a) Distribution of |log2 ratio| per bin (boxplot/violin), with a "background"
#       distribution from random SNP↔motif shuffling for comparison.
#   (b) Median/IQR per bin overlaid on background — empirical effect distance is
#       where bin curve merges with background.
#
# REQUIRES: 06_allele_pairs/pair_tables/*_pair_table.tsv.gz with extended
# MAX_PAIR_FLANK_BP=800 (apply patch first).
#
# USAGE:
#   Rscript 07_plot_snp_distance_diag.R \
#     --pair_table_dir .../pair_tables/ \
#     --plots_dir .../plots/ \
#     --tables_dir .../tables/ \
#     --bin_edges 0,1,25,50,100,200,400,800 \
#     --n_permutations 20 \
#     --min_motif_width 8 \
#     --pwm_meme .../consensus_pwms.meme
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(ggplot2)
    library(scales)
})

option_list <- list(
    make_option("--pair_table_dir",   type = "character"),
    make_option("--plots_dir",         type = "character"),
    make_option("--tables_dir",        type = "character"),
    make_option("--bin_edges",         type = "character",
                default = "0,1,25,50,100,200,400,800"),
    make_option("--n_permutations",    type = "integer", default = 20),
    make_option("--min_motif_width",   type = "integer", default = 8),
    make_option("--pwm_meme",          type = "character"),
    make_option("--lib_path",          type = "character", default = "")
)
opt <- parse_args(OptionParser(option_list = option_list))

lib_path <- if (nchar(opt$lib_path) > 0) opt$lib_path else
            file.path(dirname(sub("--file=", "",
                grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))),
                "lib_diagnostics.R")
source(lib_path)

cat("============================================================\n")
cat("07_plot_snp_distance_diag.R\n")
cat("============================================================\n\n")

dir.create(opt$plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$tables_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Filter motif IDs by width
# -----------------------------------------------------------------------------
widths <- get_archetype_pwm_widths(opt$pwm_meme)
keep_motifs <- widths[motif_width > opt$min_motif_width, motif_id]
cat("Archetypes with width >", opt$min_motif_width, ":", length(keep_motifs), "\n")

# -----------------------------------------------------------------------------
# Read pair tables — only the columns needed
# -----------------------------------------------------------------------------
files <- list.files(opt$pair_table_dir, pattern = "_pair_table\\.tsv\\.gz$",
                     full.names = TRUE)
keep_files <- files[basename(files) %in% paste0(keep_motifs, "_pair_table.tsv.gz")]
cat("Pair tables to read:", length(keep_files), "\n")
if (length(keep_files) == 0) {
    cat("No pair tables found. Skipping diagnostic.\n")
    quit(save = "no", status = 0)
}

needed_cols <- c("pair_id", "motif_id", "motif_hit_id",
                  "snp_distance_from_motif_edge_bp",
                  "snp_in_motif_body",
                  "atac_ref_total", "atac_alt_total")

read_pair_table <- function(f) {
    dt <- fread(f, select = needed_cols, fill = TRUE)
    dt
}
all_pairs <- rbindlist(lapply(keep_files, function(f) {
    tryCatch(read_pair_table(f), error = function(e) {
        warning("Skipping ", f, ": ", conditionMessage(e))
        data.table()
    })
}), fill = TRUE)
cat("Total pairs read:", nrow(all_pairs), "\n")

if (nrow(all_pairs) == 0 ||
    !all(c("atac_ref_total", "atac_alt_total") %in% names(all_pairs))) {
    cat("No ATAC count data available. Skipping.\n")
    quit(save = "no", status = 0)
}

# Compute |log2 ratio| with pseudocount
all_pairs[, atac_ref_total := as.numeric(atac_ref_total)]
all_pairs[, atac_alt_total := as.numeric(atac_alt_total)]
all_pairs[, abs_log2_atac_ratio := abs(log2((atac_ref_total + 1) /
                                              (atac_alt_total + 1)))]
all_pairs[, abs_snp_dist := abs(as.integer(snp_distance_from_motif_edge_bp))]

# Filter: must have ATAC reads
all_pairs <- all_pairs[(atac_ref_total + atac_alt_total) >= 5]
cat("After ATAC coverage filter (≥5 reads total):", nrow(all_pairs), "\n")

if (nrow(all_pairs) == 0) {
    cat("Nothing to plot.\n")
    quit(save = "no", status = 0)
}

# -----------------------------------------------------------------------------
# Bin distances
# -----------------------------------------------------------------------------
bin_edges <- as.integer(strsplit(opt$bin_edges, ",")[[1]])
# Cut: half-open [a, b). Use right=FALSE for [a,b).
bin_labels <- character(length(bin_edges))
for (i in seq_along(bin_edges)) {
    if (i == length(bin_edges)) {
        bin_labels[i] <- paste0("≥", bin_edges[i])
    } else if (bin_edges[i] == bin_edges[i + 1] - 1) {
        bin_labels[i] <- as.character(bin_edges[i])
    } else {
        bin_labels[i] <- sprintf("%d-%d", bin_edges[i], bin_edges[i + 1] - 1)
    }
}
# Add a final +Inf so cut covers everything
all_pairs[, dist_bin := cut(abs_snp_dist,
                              breaks = c(bin_edges, Inf),
                              labels = bin_labels[seq_along(bin_edges)],
                              right = FALSE,
                              include.lowest = TRUE)]
all_pairs <- all_pairs[!is.na(dist_bin)]
cat("Pairs after binning:", nrow(all_pairs), "\n")
print(table(all_pairs$dist_bin))

# -----------------------------------------------------------------------------
# Build "background" distribution by permuting SNP-motif assignments.
# Per permutation: shuffle the abs_snp_dist column independently of
# abs_log2_atac_ratio. This breaks the (motif × SNP) coupling but preserves
# both marginal distributions, giving a per-bin background.
# -----------------------------------------------------------------------------
set.seed(42)
n_perm <- opt$n_permutations
cat("Building background with", n_perm, "permutations...\n")

bg_list <- vector("list", n_perm)
for (i in seq_len(n_perm)) {
    shuffled <- copy(all_pairs)
    shuffled[, abs_snp_dist_perm := sample(abs_snp_dist)]
    shuffled[, dist_bin_perm := cut(abs_snp_dist_perm,
                                     breaks = c(bin_edges, Inf),
                                     labels = bin_labels[seq_along(bin_edges)],
                                     right = FALSE, include.lowest = TRUE)]
    bg_list[[i]] <- shuffled[, .(perm = i, dist_bin = dist_bin_perm,
                                    abs_log2_atac_ratio)]
}
bg_dt <- rbindlist(bg_list)

# -----------------------------------------------------------------------------
# Compute observed and background medians + IQRs per bin
# -----------------------------------------------------------------------------
obs_summary <- all_pairs[, .(
    n_pairs = .N,
    median_abs_log2 = median(abs_log2_atac_ratio, na.rm = TRUE),
    p25 = as.numeric(quantile(abs_log2_atac_ratio, 0.25, na.rm = TRUE)),
    p75 = as.numeric(quantile(abs_log2_atac_ratio, 0.75, na.rm = TRUE))
), by = dist_bin]
obs_summary[, source := "observed"]

bg_summary <- bg_dt[, .(
    median_abs_log2 = median(abs_log2_atac_ratio, na.rm = TRUE),
    p25 = as.numeric(quantile(abs_log2_atac_ratio, 0.25, na.rm = TRUE)),
    p75 = as.numeric(quantile(abs_log2_atac_ratio, 0.75, na.rm = TRUE))
), by = dist_bin]
bg_summary[, source := "background"]
bg_summary[, n_pairs := obs_summary$n_pairs[match(dist_bin, obs_summary$dist_bin)]]

combined <- rbind(obs_summary, bg_summary, fill = TRUE)
combined[, source := factor(source, levels = c("observed", "background"))]

# -----------------------------------------------------------------------------
# Plot 1: per-bin distributions (violins) with overlaid background median
# -----------------------------------------------------------------------------
p1 <- ggplot(all_pairs, aes(x = dist_bin, y = abs_log2_atac_ratio)) +
    geom_violin(fill = "#4a90d9", alpha = 0.5, linewidth = 0.3, scale = "width") +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white",
                 linewidth = 0.3) +
    geom_point(data = bg_summary, aes(x = dist_bin, y = median_abs_log2),
               color = "#d62728", size = 2.5, shape = 18) +
    geom_text(data = obs_summary,
              aes(x = dist_bin, label = paste0("n=", comma(n_pairs))),
              y = max(all_pairs$abs_log2_atac_ratio, na.rm = TRUE) * 0.95,
              size = 2.8, color = "grey40") +
    coord_cartesian(ylim = c(0,
        quantile(all_pairs$abs_log2_atac_ratio, 0.99, na.rm = TRUE))) +
    labs(title = "ATAC allelic asymmetry vs SNP distance from motif edge",
         subtitle = sprintf("Red diamonds = background median (%d permutations)",
                              n_perm),
         x = "|SNP–motif edge distance| (bp)",
         y = "|log2((ref+1)/(alt+1))|  ATAC reads at motif body") +
    DIAG_THEME +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_diag_plot(p1,
    file.path(opt$plots_dir, "snp_distance_atac_asymmetry.pdf"),
    width = 9, height = 5.5)

# -----------------------------------------------------------------------------
# Plot 2: signal vs distance line plot (cleaner)
# -----------------------------------------------------------------------------
combined_clean <- combined[!is.na(median_abs_log2)]
p2 <- ggplot(combined_clean, aes(x = dist_bin, y = median_abs_log2,
                                    color = source, group = source)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = p25, ymax = p75), width = 0.2, linewidth = 0.5) +
    scale_color_manual(values = c(observed = "#08519c", background = "#d62728"),
                        name = NULL) +
    labs(title = "ATAC asymmetry vs SNP distance",
         subtitle = sprintf("Median ± IQR per bin. Where observed merges with background = effective coupling distance. (%d perms)",
                              n_perm),
         x = "|SNP–motif edge distance| (bp)",
         y = "Median |log2 ATAC ref/alt|") +
    DIAG_THEME +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "top")

save_diag_plot(p2,
    file.path(opt$plots_dir, "snp_distance_atac_summary.pdf"),
    width = 8, height = 5)

# -----------------------------------------------------------------------------
# Save the bin summary table
# -----------------------------------------------------------------------------
fwrite(combined,
       file.path(opt$tables_dir, "snp_distance_atac_diagnostic.tsv"),
       sep = "\t", na = "NA")
cat("Wrote:", file.path(opt$tables_dir, "snp_distance_atac_diagnostic.tsv"), "\n")

cat("\nDone.\n")
