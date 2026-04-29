#!/usr/bin/env Rscript
# =============================================================================
# 07_pair_qc.R — Cross-motif QC report
# =============================================================================
#
# Walks all per-motif pair tables and produces:
#   1. pair_qc_summary.tsv   — per-motif summary stats (n_pairs, n with PWM
#                              scoring, n with PRO-seq coverage, etc.)
#   2. pair_qc_global.txt    — text report with histograms and global stats
#   3. coverage_histograms.tsv — flat tables for matplotlib later
#
# Doesn't make plots (text-only). Plots are an analysis-level concern.
#
# USAGE:
#   Rscript 07_pair_qc.R \
#     --pair_table_dir /path/to/pair_tables \
#     --outdir /path/to/qc/
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option("--pair_table_dir", type = "character"),
    make_option("--outdir",         type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

cat("============================================================\n")
cat("07_pair_qc.R\n")
cat("============================================================\n\n")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(opt$pair_table_dir, pattern = "_pair_table\\.tsv\\.gz$",
                    full.names = TRUE)
cat("Per-motif pair tables found:", length(files), "\n")
if (length(files) == 0) {
    stop("No per-motif pair tables in ", opt$pair_table_dir)
}

# -----------------------------------------------------------------------------
# 1. Per-motif summary
# -----------------------------------------------------------------------------
summary_rows <- list()
all_pairs <- list()
for (f in files) {
    motif_id <- sub("_pair_table\\.tsv\\.gz$", "", basename(f))
    dt <- fread(f)
    if (nrow(dt) == 0) {
        summary_rows[[motif_id]] <- data.table(
            motif_id = motif_id,
            n_pairs = 0L
        )
        next
    }

    sr <- list(motif_id = motif_id, n_pairs = nrow(dt))
    if ("snp_relative_position" %in% names(dt)) {
        tab <- table(dt$snp_relative_position)
        sr$n_body              <- as.integer(tab["body"] %||% 0)
        sr$n_upstream_flank    <- as.integer(tab["upstream_flank"] %||% 0)
        sr$n_downstream_flank  <- as.integer(tab["downstream_flank"] %||% 0)
    }
    if ("pwm_scoring_applicable" %in% names(dt)) {
        sr$n_pwm_scored <- sum(dt$pwm_scoring_applicable, na.rm = TRUE)
    }
    if ("delta_pwm_score" %in% names(dt)) {
        sr$median_abs_delta_pwm <-
            round(median(abs(dt$delta_pwm_score), na.rm = TRUE), 3)
        sr$n_strong_pwm_change <-
            sum(abs(dt$delta_pwm_score) > 1.0, na.rm = TRUE)
    }
    if ("atac_ref_total" %in% names(dt)) {
        n_pairs_with_atac <- sum((dt$atac_ref_total + dt$atac_alt_total) > 0,
                                  na.rm = TRUE)
        sr$n_with_atac_coverage <- n_pairs_with_atac
    }
    if ("pause_total_reads" %in% names(dt)) {
        sr$n_with_pause_reads <- sum(dt$pause_total_reads > 0, na.rm = TRUE)
        sr$n_with_pause_min10 <-
            sum(dt$pause_min_reads_per_allele >= 10, na.rm = TRUE)
    }
    if ("log2_pause_index_ratio_alt_over_ref" %in% names(dt)) {
        v <- dt$log2_pause_index_ratio_alt_over_ref
        sr$n_with_log2_ratio <- sum(!is.na(v))
        sr$median_log2_ratio <- round(median(v, na.rm = TRUE), 3)
    }

    summary_rows[[motif_id]] <- as.data.table(sr)
    all_pairs[[motif_id]] <- dt
}
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

summary_dt <- rbindlist(summary_rows, fill = TRUE)
summary_path <- file.path(opt$outdir, "pair_qc_summary.tsv")
fwrite(summary_dt, summary_path, sep = "\t", na = "NA")
cat("  Wrote per-motif summary:", summary_path, "\n")

# -----------------------------------------------------------------------------
# 2. Coverage histograms
# -----------------------------------------------------------------------------
all_dt <- rbindlist(all_pairs, fill = TRUE)
cat("\nGlobal pair count:", nrow(all_dt), "\n")

bin_edges_atac    <- c(-1, 0, 5, 10, 25, 50, 100, 1e9)
bin_edges_pause   <- c(-1, 0, 5, 10, 25, 50, 100, 1e9)
bin_edges_psi_min <- c(-1, 0, 4, 10, 25, 50, 1e9)
hist_dt <- list()

if ("atac_ref_total" %in% names(all_dt)) {
    h <- table(cut(all_dt$atac_ref_total + all_dt$atac_alt_total,
                    bin_edges_atac))
    hist_dt$atac_total <- data.table(metric = "atac_total",
                                     bin = names(h), n_pairs = as.integer(h))
}

if ("pause_total_reads" %in% names(all_dt)) {
    h <- table(cut(all_dt$pause_total_reads, bin_edges_pause))
    hist_dt$pause_total <- data.table(metric = "pause_total",
                                       bin = names(h), n_pairs = as.integer(h))
}

if ("pause_min_reads_per_allele" %in% names(all_dt)) {
    h <- table(cut(all_dt$pause_min_reads_per_allele, bin_edges_psi_min))
    hist_dt$pause_min <- data.table(metric = "pause_min_per_allele",
                                     bin = names(h), n_pairs = as.integer(h))
}

if (length(hist_dt) > 0) {
    hist_combined <- rbindlist(hist_dt)
    hist_path <- file.path(opt$outdir, "coverage_histograms.tsv")
    fwrite(hist_combined, hist_path, sep = "\t")
    cat("  Wrote coverage histograms:", hist_path, "\n")
}

# -----------------------------------------------------------------------------
# 3. Text report
# -----------------------------------------------------------------------------
report_path <- file.path(opt$outdir, "pair_qc_global.txt")
sink(report_path)
cat("================================================================\n")
cat("F121-9 Allele Pair QC Report\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

cat("Motifs processed:                ", length(files), "\n")
cat("Total pairs across motifs:       ", nrow(all_dt), "\n\n")

if ("snp_relative_position" %in% names(all_dt)) {
    cat("Pair distribution by relative position:\n")
    print(all_dt[, .N, by = snp_relative_position])
    cat("\n")
}

if ("pwm_scoring_applicable" %in% names(all_dt)) {
    cat("PWM scoring applicability:\n")
    print(all_dt[, .N, by = pwm_scoring_applicable])
    cat("\n")
}

if ("snp_offset_from_motif_center_bp" %in% names(all_dt)) {
    cat("Distribution of |SNP offset from motif center|:\n")
    abs_off <- abs(all_dt$snp_offset_from_motif_center_bp)
    cat("  median: ", median(abs_off), "bp\n")
    cat("  IQR:    ", quantile(abs_off, 0.25), "-", quantile(abs_off, 0.75), "bp\n")
    cat("  90%-le: ", quantile(abs_off, 0.90), "bp\n")
    cat("\n")
}

if ("delta_pwm_score" %in% names(all_dt)) {
    cat("PWM ref->alt impact (body-SNPs only):\n")
    v <- all_dt$delta_pwm_score
    cat("  Pairs scored:", sum(!is.na(v)), "\n")
    cat("  Median |delta|:", round(median(abs(v), na.rm = TRUE), 3), "\n")
    cat("  Pairs with |delta| > 1.0:", sum(abs(v) > 1.0, na.rm = TRUE), "\n")
    cat("  Pairs with |delta| > 2.0:", sum(abs(v) > 2.0, na.rm = TRUE), "\n")
    cat("\n")
}

if ("atac_ref_total" %in% names(all_dt)) {
    cat("ATAC coverage at motifs:\n")
    cat("  Pairs with ANY allelic ATAC reads:",
        sum((all_dt$atac_ref_total + all_dt$atac_alt_total) > 0), "\n")
    cat("  Pairs with >=10 reads/allele (ref AND alt):",
        sum(pmin(all_dt$atac_ref_total, all_dt$atac_alt_total) >= 10), "\n")
    cat("\n")
}

if ("pause_total_reads" %in% names(all_dt)) {
    cat("PRO-seq pause-window coverage:\n")
    cat("  Pairs with ANY pause-window reads:",
        sum(all_dt$pause_total_reads > 0, na.rm = TRUE), "\n")
    cat("  Pairs with >=10 reads/allele in pause window:",
        sum(all_dt$pause_min_reads_per_allele >= 10, na.rm = TRUE), "\n")
    cat("  Pairs with valid log2(pause_index alt/ref):",
        sum(!is.na(all_dt$log2_pause_index_ratio_alt_over_ref)), "\n")
    cat("\n")
}

# PSI cell-by-cell summary
psi_cells <- c("donor_sense_upstream", "donor_sense_downstream",
                "donor_antisense_upstream", "donor_antisense_downstream",
                "acceptor_sense_upstream", "acceptor_sense_downstream",
                "acceptor_antisense_upstream", "acceptor_antisense_downstream")
for (cell in psi_cells) {
    cn_psi_ref <- paste0(cell, "_psi_ref")
    cn_psi_alt <- paste0(cell, "_psi_alt")
    if (cn_psi_ref %in% names(all_dt)) {
        n_with_both <- sum(!is.na(all_dt[[cn_psi_ref]]) &
                           !is.na(all_dt[[cn_psi_alt]]))
        cat("Cell ", cell, ":\n", sep = "")
        cat("  Pairs with PSI ref AND alt:", n_with_both, "\n")
    }
}
sink()
cat("\n  Wrote text report:", report_path, "\n\n")

cat("==========================================\n")
cat("QC complete. Files in:", opt$outdir, "\n")
cat("==========================================\n")
