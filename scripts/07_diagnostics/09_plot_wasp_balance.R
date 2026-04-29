#!/usr/bin/env Rscript
# =============================================================================
# 09_plot_wasp_balance.R — WASP cascade + ref:alt balance + per-SNP coverage
# =============================================================================
#
# Three plots from 08_compute_wasp_balance.sh outputs:
#
#   (A) plots/wasp_cascade.pdf
#       Stacked bars per assay × sample showing reads at each WASP stage.
#
#   (B) plots/allele_balance.pdf
#       Per assay × sample: ref vs alt totals as side-by-side bars,
#       with reference-fraction overlaid (should be ~0.5).
#
#   (C) plots/per_snp_coverage_hist.pdf
#       Histogram of per-SNP read depth (ref+alt summed across reps within
#       each assay), with vertical lines at 10/25/50 read thresholds and
#       counts of SNPs above each.
#
# USAGE:
#   Rscript 09_plot_wasp_balance.R \
#     --wasp_dir .../wasp/ \
#     --plots_dir .../plots/
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(ggplot2)
    library(scales)
    library(patchwork)
})

option_list <- list(
    make_option("--wasp_dir",   type = "character"),
    make_option("--plots_dir",  type = "character"),
    make_option("--lib_path",   type = "character", default = "")
)
opt <- parse_args(OptionParser(option_list = option_list))

lib_path <- if (nchar(opt$lib_path) > 0) opt$lib_path else
            file.path(dirname(sub("--file=", "",
                grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))),
                "lib_diagnostics.R")
source(lib_path)

cat("============================================================\n")
cat("09_plot_wasp_balance.R\n")
cat("============================================================\n\n")

dir.create(opt$plots_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# (A) WASP cascade
# -----------------------------------------------------------------------------
cascade_path <- file.path(opt$wasp_dir, "wasp_cascade.tsv")
if (!file.exists(cascade_path)) {
    cat("[SKIP] No wasp_cascade.tsv found; run 08_compute_wasp_balance.sh first.\n")
} else {
    cas <- fread(cascade_path)
    cas[, n_reads_num := suppressWarnings(as.numeric(n_reads))]
    cas <- cas[!is.na(n_reads_num)]
    if (nrow(cas) == 0) {
        cat("[WARN] cascade table empty (all values NA).\n")
    } else {
        # Keep only the cumulative cascade stages (not the ref/alt rows)
        cas_stages <- cas[stage %in% c("total_input", "wasp_passed", "allele_split_total")]
        cas_stages[, stage := factor(stage,
            levels = c("total_input", "wasp_passed", "allele_split_total"),
            labels = c("Total input", "WASP-passed (vW=1)", "Allele-split total"))]

        p_a <- ggplot(cas_stages, aes(x = sample, y = n_reads_num, fill = stage)) +
            geom_col(position = position_dodge(width = 0.8), width = 0.7) +
            scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.10))) +
            scale_fill_manual(values = c("Total input" = "#969696",
                                          "WASP-passed (vW=1)" = "#4292c6",
                                          "Allele-split total" = "#08519c"),
                               name = NULL) +
            facet_wrap(~ assay, scales = "free", ncol = 3) +
            labs(title = "WASP filter cascade per assay × replicate",
                 subtitle = "Reads at each filtering stage. WASP-passed reads should be the\nmajority of input; allele-split should ≈ WASP-passed.",
                 x = NULL, y = "Number of reads") +
            DIAG_THEME +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "top")

        save_diag_plot(p_a,
            file.path(opt$plots_dir, "wasp_cascade.pdf"),
            width = 12, height = 5)
        cat("Wrote: wasp_cascade.pdf\n")
    }
}

# -----------------------------------------------------------------------------
# (B) Allele balance
# -----------------------------------------------------------------------------
balance_path <- file.path(opt$wasp_dir, "allele_balance.tsv")
if (file.exists(balance_path)) {
    bal <- fread(balance_path)
    bal[, n_ref := suppressWarnings(as.numeric(n_ref))]
    bal[, n_alt := suppressWarnings(as.numeric(n_alt))]
    bal[, ref_frac := suppressWarnings(as.numeric(ref_frac))]
    bal <- bal[!is.na(n_ref) & !is.na(n_alt)]

    if (nrow(bal) > 0) {
        bal_long <- melt(bal, id.vars = c("assay", "sample"),
                          measure.vars = c("n_ref", "n_alt"),
                          variable.name = "allele", value.name = "n_reads")
        bal_long[, allele := factor(allele, levels = c("n_ref", "n_alt"),
                                       labels = c("Ref", "Alt"))]

        p_bars <- ggplot(bal_long, aes(x = sample, y = n_reads, fill = allele)) +
            geom_col(position = position_dodge(width = 0.8), width = 0.7) +
            scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.10))) +
            scale_fill_manual(values = c(Ref = "#3182bd", Alt = "#e6550d"),
                               name = NULL) +
            facet_wrap(~ assay, scales = "free", ncol = 3) +
            labs(title = "Global ref vs alt read counts (post-WASP)",
                 x = NULL, y = "Number of reads") +
            DIAG_THEME +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "top")

        # Reference fraction panel
        p_frac <- ggplot(bal, aes(x = sample, y = ref_frac)) +
            geom_col(fill = "#54278f", width = 0.6) +
            geom_hline(yintercept = 0.5, linetype = "dashed", color = "red",
                        linewidth = 0.5) +
            geom_text(aes(label = sprintf("%.3f", ref_frac)),
                      vjust = -0.3, size = 2.8) +
            scale_y_continuous(limits = c(0, 0.65),
                                expand = expansion(mult = c(0, 0.10))) +
            facet_wrap(~ assay, scales = "free_x", ncol = 3) +
            labs(subtitle = "Reference fraction = ref / (ref + alt). Ideal = 0.5 (red).",
                 x = NULL, y = "Ref / (ref+alt)") +
            DIAG_THEME +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))

        combined <- p_bars / p_frac + plot_layout(heights = c(2, 1))
        save_diag_plot(combined,
            file.path(opt$plots_dir, "allele_balance.pdf"),
            width = 12, height = 8)
        cat("Wrote: allele_balance.pdf\n")
    }
}

# -----------------------------------------------------------------------------
# (C) Per-SNP coverage histograms
# -----------------------------------------------------------------------------
per_snp_path <- file.path(opt$wasp_dir, "per_snp_coverage.tsv.gz")
if (!file.exists(per_snp_path)) {
    cat("[SKIP] No per_snp_coverage.tsv.gz; run 08 first (or it skipped if pair indices missing).\n")
} else {
    pile <- fread(per_snp_path)
    if (nrow(pile) == 0) {
        cat("[WARN] per-SNP coverage table empty.\n")
    } else {
        # Sum reads per SNP per assay (combining ref+alt and across reps)
        pile[, n_reads := as.integer(n_reads)]
        per_snp <- pile[, .(total_reads = sum(n_reads, na.rm = TRUE)),
                          by = .(chrom, pos, assay)]

        thresholds <- c(10, 25, 50)
        # For annotation per assay
        thresh_dt <- per_snp[, {
            tots <- total_reads
            .(threshold = thresholds,
               n_snps_above = sapply(thresholds, function(t) sum(tots >= t)))
        }, by = assay]

        # Plot histograms per assay
        p_c <- ggplot(per_snp[total_reads > 0], aes(x = total_reads)) +
            geom_histogram(fill = "#4a90d9", color = "white", linewidth = 0.2,
                           bins = 60) +
            geom_vline(data = data.table(threshold = thresholds),
                       aes(xintercept = threshold), linetype = "dashed",
                       color = "#d62728", linewidth = 0.4) +
            facet_wrap(~ assay, scales = "free", ncol = 3) +
            scale_x_log10(labels = comma,
                          breaks = c(1, 10, 100, 1000, 10000)) +
            scale_y_continuous(labels = comma) +
            labs(title = "Per-SNP read depth (summed across reps + alleles per assay)",
                 subtitle = "Red dashed lines = read-depth thresholds (10 / 25 / 50). x-axis log10 scale.",
                 x = "Total reads at SNP position",
                 y = "Number of SNPs") +
            DIAG_THEME

        # Build per-panel annotation text
        annot_dt <- thresh_dt[, .(label = paste(sprintf("\u2265%d: %s SNPs",
                                                          threshold,
                                                          comma(n_snps_above)),
                                                  collapse = "\n")),
                               by = assay]

        p_c <- p_c + geom_text(data = annot_dt,
            aes(x = Inf, y = Inf, label = label),
            hjust = 1.05, vjust = 1.2, size = 3, color = "grey20")

        save_diag_plot(p_c,
            file.path(opt$plots_dir, "per_snp_coverage_hist.pdf"),
            width = 12, height = 5)
        cat("Wrote: per_snp_coverage_hist.pdf\n")
    }
}

cat("\nDone.\n")
