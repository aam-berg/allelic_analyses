#!/usr/bin/env Rscript
# =============================================================================
# 05b_compute_psi_per_allele.R — PSI per allele for each motif's 8 junction cells
# =============================================================================
#
# PSI definition (per allele):
#   PSI = use_reads / (use_reads + skip_reads)
#
# Where "use" = reads spliced through the site (an N skip with one boundary
# at the site), "skip" = continuous reads spanning the site (no N intron
# overlaps it).
#
# This is the "junction-centric" PSI variant: we focus on a specific splice
# site (donor or acceptor position) and ask, of the reads that COULD inform
# splicing at this site, what fraction actually splices through it. This is
# more conservative than the standard "exon inclusion PSI" but is the right
# measure for our question — does pausing at a TFBS change usage of the
# nearest splice site?
#
# OUTPUT (wide, per motif × cell):
#   motif_hit_id, junction_cell, site_chrom, site_pos, site_strand,
#   site_type, source, motif_type, distance,
#   use_ref_total, use_alt_total, skip_ref_total, skip_alt_total,
#   psi_ref, psi_alt, delta_psi (= psi_alt - psi_ref),
#   total_reads (= sum of use+skip across alleles),
#   psi_min_reads_per_allele
#
# A wide motif × cells table (8 rows per motif) and a wide motif × all_cells
# table (1 row per motif, 8×N columns) are both produced — the second is for
# direct join into the master pair table.
#
# USAGE:
#   Rscript 05b_compute_psi_per_allele.R \
#     --junction_counts /path/to/AC0001_rnaseq_junctions.tsv.gz \
#     --psi_min_total_reads 4 \
#     --motif_id AC0001 \
#     --outdir /path/to/psi/
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option("--junction_counts",     type = "character"),
    make_option("--psi_min_total_reads", type = "integer", default = 4),
    make_option("--motif_id",            type = "character"),
    make_option("--outdir",              type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

cat("============================================================\n")
cat("05b_compute_psi_per_allele.R — motif:", opt$motif_id, "\n")
cat("============================================================\n\n")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(opt$junction_counts)
if (nrow(dt) == 0) {
    cat("[INFO] Empty junction counts; writing empty PSI tables.\n")
    fwrite(data.table(),
        file.path(opt$outdir, paste0(opt$motif_id, "_psi_long.tsv.gz")),
        sep = "\t", compress = "gzip")
    fwrite(data.table(),
        file.path(opt$outdir, paste0(opt$motif_id, "_psi_wide.tsv.gz")),
        sep = "\t", compress = "gzip")
    quit(save = "no", status = 0)
}

cat("Junction-cell rows:", nrow(dt), "\n")
cat("Unique motifs:", length(unique(dt$motif_hit_id)), "\n")

# ---- PSI per allele ----
dt[, total_use_ref  := use_ref_total]
dt[, total_use_alt  := use_alt_total]
dt[, total_skip_ref := skip_ref_total]
dt[, total_skip_alt := skip_alt_total]

dt[, ref_total := total_use_ref + total_skip_ref]
dt[, alt_total := total_use_alt + total_skip_alt]
dt[, total_reads_both_alleles := ref_total + alt_total]
dt[, psi_min_reads_per_allele := pmin(ref_total, alt_total)]

# PSI requires at least psi_min_total_reads in (use+skip) per allele to be reported
dt[, psi_ref := ifelse(ref_total >= opt$psi_min_total_reads,
                        total_use_ref / ref_total, NA_real_)]
dt[, psi_alt := ifelse(alt_total >= opt$psi_min_total_reads,
                        total_use_alt / alt_total, NA_real_)]
dt[, delta_psi := psi_alt - psi_ref]

# ---- LONG OUTPUT (one row per motif × cell) ----
long_cols <- c("motif_hit_id", "junction_cell",
               "site_chrom", "site_pos", "site_strand",
               "site_type", "source", "motif_type", "distance",
               "total_use_ref", "total_use_alt",
               "total_skip_ref", "total_skip_alt",
               "ref_total", "alt_total",
               "total_reads_both_alleles", "psi_min_reads_per_allele",
               "psi_ref", "psi_alt", "delta_psi")

long_path <- file.path(opt$outdir, paste0(opt$motif_id, "_psi_long.tsv.gz"))
fwrite(dt[, ..long_cols], long_path, sep = "\t", na = "NA", compress = "gzip")
cat("\nWrote LONG (motif × cell):", long_path, "\n")

# ---- WIDE OUTPUT (one row per motif; 8 cells × ~6 cols = ~48 cols) ----
# Pivot so each cell's columns get prefixed
wide_dt <- dcast(
    dt,
    motif_hit_id ~ junction_cell,
    value.var = c("site_chrom", "site_pos", "site_strand", "site_type",
                  "source", "motif_type", "distance",
                  "total_use_ref", "total_use_alt",
                  "total_skip_ref", "total_skip_alt",
                  "ref_total", "alt_total",
                  "psi_ref", "psi_alt", "delta_psi"),
    fun.aggregate = function(x) if (length(x) == 0) NA else x[1]
)

# Reorder columns: keep value-type prefixes grouped per cell for readability.
# data.table dcast names columns like "site_pos_donor_sense_upstream".
# We want them like "donor_sense_upstream_site_pos".
# This is cosmetic but helps for inspection.
old_names <- names(wide_dt)[-1]
new_names <- vapply(old_names, function(n) {
    # Find cell suffix among the 8 known cells
    cells <- c("donor_sense_upstream", "donor_sense_downstream",
                "donor_antisense_upstream", "donor_antisense_downstream",
                "acceptor_sense_upstream", "acceptor_sense_downstream",
                "acceptor_antisense_upstream", "acceptor_antisense_downstream")
    for (c in cells) {
        suffix <- paste0("_", c)
        if (endsWith(n, suffix)) {
            value_part <- substr(n, 1, nchar(n) - nchar(suffix))
            return(paste(c, value_part, sep = "_"))
        }
    }
    n
}, character(1))

setnames(wide_dt, c("motif_hit_id", new_names))
setcolorder(wide_dt, c("motif_hit_id", sort(new_names)))

wide_path <- file.path(opt$outdir, paste0(opt$motif_id, "_psi_wide.tsv.gz"))
fwrite(wide_dt, wide_path, sep = "\t", na = "NA", compress = "gzip")
cat("Wrote WIDE (motif × all-cells flat):", wide_path, "\n")
cat("  Wide dims:", nrow(wide_dt), "x", ncol(wide_dt), "\n")

# ---- Sanity report ----
cat("\nPSI summary across cells:\n")
for (cell in unique(dt$junction_cell)) {
    sub <- dt[junction_cell == cell]
    n_with_psi <- sum(!is.na(sub$psi_ref) & !is.na(sub$psi_alt))
    if (n_with_psi == 0) next
    median_dpsi <- median(abs(sub$delta_psi), na.rm = TRUE)
    n_strong <- sum(abs(sub$delta_psi) >= 0.1, na.rm = TRUE)
    cat(sprintf("  %s: n_with_psi=%d, median |delta_psi|=%.3f, n |dpsi|>=0.1: %d\n",
                cell, n_with_psi, median_dpsi, n_strong))
}
