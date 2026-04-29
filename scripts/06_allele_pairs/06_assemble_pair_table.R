#!/usr/bin/env Rscript
# =============================================================================
# 06_assemble_pair_table.R — Master join of all per-pair / per-motif features
# =============================================================================
#
# For ONE motif, joins together:
#   - pair_index (geometry, per pair)
#   - pwm_scores (per pair, PWM ref/alt)
#   - atac_counts (per motif → broadcast to all pairs of that motif)
#   - proseq_window_counts (per motif → broadcast)
#   - pause_indices (per motif → broadcast)
#   - psi_wide (per motif → broadcast; 8 cells × ~16 cols)
#
# Produces ONE per-motif table:
#   {motif_id}_pair_table.tsv.gz
#
# IMPORTANT: this table does NOT include the static-genomic features from
# 05_motif_annot/. Those are in {motif_id}_annotated.tsv.gz and joined by
# motif_hit_id at analysis time. See README for the recommended join.
#
# USAGE:
#   Rscript 06_assemble_pair_table.R \
#     --pair_index .../AC0001_pair_index.tsv.gz \
#     --pwm_scores .../AC0001_pwm_scores.tsv.gz \
#     --atac_counts .../AC0001_atac_counts.tsv.gz \
#     --proseq_window .../AC0001_proseq_window_counts.tsv.gz \
#     --pause_indices .../AC0001_pause_indices.tsv.gz \
#     --psi_wide .../AC0001_psi_wide.tsv.gz \
#     --motif_id AC0001 \
#     --outdir /path/to/pair_tables/
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option("--pair_index",     type = "character"),
    make_option("--pwm_scores",     type = "character"),
    make_option("--atac_counts",    type = "character"),
    make_option("--proseq_window",  type = "character"),
    make_option("--pause_indices",  type = "character"),
    make_option("--psi_wide",       type = "character"),
    make_option("--motif_id",       type = "character"),
    make_option("--outdir",         type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

cat("============================================================\n")
cat("06_assemble_pair_table.R — motif:", opt$motif_id, "\n")
cat("============================================================\n\n")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

read_or_empty <- function(p, label) {
    if (!file.exists(p)) {
        cat("  [WARN]", label, "missing:", p, "\n")
        return(data.table())
    }
    dt <- fread(p)
    cat("  [OK]   ", label, ":", nrow(dt), "rows x", ncol(dt), "cols\n")
    dt
}

cat("Loading inputs:\n")
pair_dt    <- read_or_empty(opt$pair_index,    "pair_index    ")
pwm_dt     <- read_or_empty(opt$pwm_scores,    "pwm_scores    ")
atac_dt    <- read_or_empty(opt$atac_counts,   "atac_counts   ")
window_dt  <- read_or_empty(opt$proseq_window, "proseq_window ")
pause_dt   <- read_or_empty(opt$pause_indices, "pause_indices ")
psi_dt     <- read_or_empty(opt$psi_wide,      "psi_wide      ")

if (nrow(pair_dt) == 0) {
    cat("\n[INFO] Empty pair index; writing empty assembled table.\n")
    fwrite(data.table(),
        file.path(opt$outdir, paste0(opt$motif_id, "_pair_table.tsv.gz")),
        sep = "\t", compress = "gzip")
    quit(save = "no", status = 0)
}

# --- Join PWM scores (per pair_id) ---
if (nrow(pwm_dt) > 0) {
    pwm_keep <- setdiff(names(pwm_dt), c("motif_hit_id"))
    pair_dt <- pwm_dt[, ..pwm_keep][pair_dt, on = "pair_id"]
    cat("\nAfter PWM join:", ncol(pair_dt), "cols\n")
}

# --- Join ATAC (per motif_hit_id; broadcast) ---
if (nrow(atac_dt) > 0) {
    pair_dt <- atac_dt[pair_dt, on = "motif_hit_id"]
    cat("After ATAC join:", ncol(pair_dt), "cols\n")
}

# --- Join PRO-seq window counts (per motif_hit_id; broadcast) ---
if (nrow(window_dt) > 0) {
    pair_dt <- window_dt[pair_dt, on = "motif_hit_id"]
    cat("After PRO-seq window join:", ncol(pair_dt), "cols\n")
}

# --- Join pause indices (per motif_hit_id; broadcast) ---
if (nrow(pause_dt) > 0) {
    pair_dt <- pause_dt[pair_dt, on = "motif_hit_id"]
    cat("After pause-index join:", ncol(pair_dt), "cols\n")
}

# --- Join PSI wide (per motif_hit_id; broadcast) ---
if (nrow(psi_dt) > 0) {
    pair_dt <- psi_dt[pair_dt, on = "motif_hit_id"]
    cat("After PSI join:", ncol(pair_dt), "cols\n")
}

# --- Reorder canonical columns first ---
canonical_first <- c(
    "pair_id", "motif_hit_id", "motif_id",
    "motif_chrom", "motif_start", "motif_end", "motif_strand", "motif_score",
    "motif_width",
    "snp_chrom", "snp_pos", "snp_ref", "snp_alt",
    "snp_offset_from_motif_center_bp",
    "snp_offset_from_motif_center_abs_bp",
    "snp_distance_from_motif_edge_bp",
    "snp_in_motif_body", "snp_upstream_of_motif",
    "snp_relative_position", "snp_position_in_motif",
    "pwm_scoring_applicable",
    "motif_sequence_ref", "motif_sequence_alt",
    "pwm_score_ref", "pwm_score_alt", "delta_pwm_score"
)
canonical_present <- intersect(canonical_first, names(pair_dt))
other_cols <- setdiff(names(pair_dt), canonical_present)
setcolorder(pair_dt, c(canonical_present, other_cols))

out_path <- file.path(opt$outdir,
    paste0(opt$motif_id, "_pair_table.tsv.gz"))
fwrite(pair_dt, out_path, sep = "\t", na = "NA", compress = "gzip")

cat("\n=========================================\n")
cat("Wrote master pair table:\n")
cat("  ", out_path, "\n")
cat("  Dimensions:", nrow(pair_dt), "rows x", ncol(pair_dt), "cols\n")
cat("=========================================\n")

# --- Sanity print ---
cat("\nPair categories (snp_relative_position × pwm applicable):\n")
print(pair_dt[, .N, by = .(snp_relative_position, pwm_scoring_applicable)])

if ("delta_pwm_score" %in% names(pair_dt)) {
    n_strong <- sum(abs(pair_dt$delta_pwm_score) > 1.0, na.rm = TRUE)
    cat("\nPairs with |delta_pwm_score| > 1.0:", n_strong,
        "(", round(100 * n_strong / nrow(pair_dt), 1), "% of all)\n")
}

if ("log2_pause_index_ratio_alt_over_ref" %in% names(pair_dt)) {
    cat("Pairs with non-NA pause-index log-ratio:",
        sum(!is.na(pair_dt$log2_pause_index_ratio_alt_over_ref)), "\n")
}
