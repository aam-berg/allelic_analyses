#!/usr/bin/env Rscript
# =============================================================================
# 04c_compute_pause_indices.R — Compute local pause indices per allele
# =============================================================================
#
# Reads per-motif PRO-seq window counts and computes:
#
#   pause_index_ref / pause_index_alt
#     = (proseq_pause_<allele> / pause_window_width_bp) /
#       max(proseq_genebody_<allele> / gene_body_width_bp, eps)
#   log2_pause_index_ratio_alt_over_ref
#     = log2( (pause_index_alt + eps) / (pause_index_ref + eps) )
#
# Plus power columns:
#   pause_total_reads = proseq_pause_ref + proseq_pause_alt
#   genebody_total_reads = proseq_genebody_ref + proseq_genebody_alt
#   pause_min_reads_per_allele = min(pause_ref, pause_alt)
#
# These let downstream code filter by power before plotting.
#
# IMPORTANT — RAW COUNT COLUMNS PRESERVED:
#   The raw count columns from 04 are passed through unchanged. They're
#   the substrate for proper binomial / negative-binomial tests. The
#   log-ratios are the convenient summary; raw counts are the truth.
#
# USAGE:
#   Rscript 04c_compute_pause_indices.R \
#     --window_counts /path/to/proseq_counts/AC0001_proseq_window_counts.tsv.gz \
#     --pause_window_half_bp 10 \
#     --gene_body_flank_bp 2000 \
#     --motif_id AC0001 \
#     --outdir /path/to/pause_indices/
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option("--window_counts",       type = "character"),
    make_option("--pause_window_half_bp", type = "integer", default = 10),
    make_option("--gene_body_flank_bp",   type = "integer", default = 2000),
    make_option("--motif_id",            type = "character"),
    make_option("--outdir",              type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

cat("============================================================\n")
cat("04c_compute_pause_indices.R — motif:", opt$motif_id, "\n")
cat("============================================================\n\n")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(opt$window_counts)
if (nrow(dt) == 0 || !"proseq_pause_ref" %in% names(dt)) {
    cat("[INFO] Empty window counts; writing empty pause index table.\n")
    fwrite(data.table(),
        file.path(opt$outdir, paste0(opt$motif_id, "_pause_indices.tsv.gz")),
        sep = "\t", compress = "gzip")
    quit(save = "no", status = 0)
}

cat("Motifs:", nrow(dt), "\n")

eps <- 1e-3

pause_w_bp    <- 2 * opt$pause_window_half_bp + 1   # inclusive width
genebody_w_bp <- opt$gene_body_flank_bp

# Per-bp normalized counts
dt[, pause_density_ref    := proseq_pause_ref    / pause_w_bp]
dt[, pause_density_alt    := proseq_pause_alt    / pause_w_bp]
dt[, genebody_density_ref := proseq_genebody_ref / genebody_w_bp]
dt[, genebody_density_alt := proseq_genebody_alt / genebody_w_bp]

# Pause index = pause_density / gene_body_density
dt[, pause_index_ref := pause_density_ref / pmax(genebody_density_ref, eps)]
dt[, pause_index_alt := pause_density_alt / pmax(genebody_density_alt, eps)]

# Log-ratio (alt over ref): positive = more paused on alt allele
dt[, log2_pause_index_ratio_alt_over_ref :=
       log2((pause_index_alt + eps) / (pause_index_ref + eps))]

# Power columns
dt[, pause_total_reads        := proseq_pause_ref + proseq_pause_alt]
dt[, genebody_total_reads     := proseq_genebody_ref + proseq_genebody_alt]
dt[, pause_min_reads_per_allele :=
       pmin(proseq_pause_ref, proseq_pause_alt)]
dt[, genebody_min_reads_per_allele :=
       pmin(proseq_genebody_ref, proseq_genebody_alt)]

# Sanity logs
cat("\nPause index summary (across motifs):\n")
cat("  pause_index_ref:    median=",
    round(median(dt$pause_index_ref, na.rm = TRUE), 2),
    "  (>1 = paused, <1 = not paused)\n")
cat("  pause_index_alt:    median=",
    round(median(dt$pause_index_alt, na.rm = TRUE), 2), "\n")
cat("  log2_ratio_alt/ref: median=",
    round(median(dt$log2_pause_index_ratio_alt_over_ref, na.rm = TRUE), 3), "\n")
cat("  Motifs with >=10 reads in pause window (both alleles): ",
    sum(dt$pause_min_reads_per_allele >= 10), "/", nrow(dt), "\n", sep = "")
cat("  Motifs with >=20 reads in gene body (both alleles):    ",
    sum(dt$genebody_min_reads_per_allele >= 20), "/", nrow(dt), "\n", sep = "")

# Output: just the new columns + motif_hit_id (the wide_counts file already
# has the raw count columns; we don't duplicate them)
out_cols <- c("motif_hit_id",
              "pause_density_ref", "pause_density_alt",
              "genebody_density_ref", "genebody_density_alt",
              "pause_index_ref", "pause_index_alt",
              "log2_pause_index_ratio_alt_over_ref",
              "pause_total_reads", "genebody_total_reads",
              "pause_min_reads_per_allele", "genebody_min_reads_per_allele")

out_path <- file.path(opt$outdir,
    paste0(opt$motif_id, "_pause_indices.tsv.gz"))
fwrite(dt[, ..out_cols], out_path, sep = "\t", na = "NA", compress = "gzip")

cat("\nWrote:", out_path, "\n")
cat("Dimensions:", nrow(dt), "rows x", length(out_cols), "cols\n")
