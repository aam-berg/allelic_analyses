#!/usr/bin/env Rscript
# =============================================================================
# 01_build_summary.R — Per-motif filter cascade + qualifying stats
# =============================================================================
#
# For ONE motif's annotated TSV, applies the filter cascade and writes:
#
#   {SUMMARY_DIR}/{motif_id}_summary.tsv     one row, ~50 cols of stats
#   {SUMMARY_DIR}/{motif_id}_scores.tsv.gz   per-TFBS scores at each stage
#                                            (only kept TFBSs at the final stage,
#                                             plus a 'qualifying' column)
#
# The summary row contains:
#   motif_id, motif_width
#   For each cascade stage:
#     n_<stage>            number of TFBSs surviving
#     ubp_<stage>          unique bp covered by surviving TFBSs (after reduce)
#     score_mean_<stage>   mean MOODS score
#     score_median_<stage> median MOODS score
#     score_p90_<stage>    90th percentile MOODS score
#   For each qualifying-distance threshold (overlap, flank10, ..., flank800):
#     n_qualifying_<thresh>
#     score_mean_qualifying_<thresh>
#     score_median_qualifying_<thresh>
#
# USAGE:
#   Rscript 01_build_summary.R \
#     --annot_tsv .../{motif_id}_annotated.tsv.gz \
#     --motif_id AC0001 \
#     --motif_width 6 \
#     --apply_expression_filter TRUE \
#     --expression_tpm_min 1.0 \
#     --hybrid F121-9 \
#     --outdir .../per_archetype_summaries/
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option("--annot_tsv",                type = "character"),
    make_option("--motif_id",                  type = "character"),
    make_option("--motif_width",               type = "integer"),
    make_option("--apply_expression_filter",   type = "logical", default = TRUE),
    make_option("--expression_tpm_min",        type = "double",  default = 1.0),
    make_option("--standard_chroms",           type = "character",
                default = paste0("chr", c(1:19, "X"), collapse = ",")),
    make_option("--gene_types_include",        type = "character",
                default = "protein_coding,lncRNA"),
    make_option("--hybrid",                    type = "character", default = "F121-9"),
    make_option("--outdir",                    type = "character"),
    make_option("--lib_path",                  type = "character", default = "")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Source library (relative to this script unless overridden)
lib_path <- if (nchar(opt$lib_path) > 0) opt$lib_path else
            file.path(dirname(sub("--file=", "",
                grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))),
                "lib_diagnostics.R")
if (!file.exists(lib_path)) {
    # Fall back to env var or hard path
    lib_path <- Sys.getenv("LIB_DIAGNOSTICS_R", lib_path)
}
source(lib_path)

cat("01_build_summary.R — motif:", opt$motif_id,
    " | width:", opt$motif_width, "\n")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Read annotated TSV
# -----------------------------------------------------------------------------
if (!file.exists(opt$annot_tsv)) {
    cat("  [WARN] Missing annotated TSV; writing empty summary.\n")
    fwrite(data.table(motif_id = opt$motif_id,
                      motif_width = opt$motif_width,
                      missing_input = TRUE),
           file.path(opt$outdir, paste0(opt$motif_id, "_summary.tsv")),
           sep = "\t")
    quit(save = "no", status = 0)
}

dt <- fread(opt$annot_tsv)
n_total <- nrow(dt)
cat("  Read", n_total, "TFBS rows.\n")

if (n_total == 0) {
    cat("  Empty annotation; writing empty summary.\n")
    fwrite(data.table(motif_id = opt$motif_id,
                      motif_width = opt$motif_width,
                      empty_input = TRUE),
           file.path(opt$outdir, paste0(opt$motif_id, "_summary.tsv")),
           sep = "\t")
    quit(save = "no", status = 0)
}

# Annotated TSV has BED 0-based half-open: convert start/end to 1-based for our use
dt[, start_1b := start + 1L]
dt[, end_1b := end]

# -----------------------------------------------------------------------------
# Build cascade
# -----------------------------------------------------------------------------
std_chroms <- strsplit(opt$standard_chroms, ",")[[1]]
gene_types <- strsplit(opt$gene_types_include, ",")[[1]]

cascade <- build_cascade_steps(
    standard_chroms        = std_chroms,
    gene_types_include     = gene_types,
    expression_tpm_min     = opt$expression_tpm_min,
    apply_expression_filter = opt$apply_expression_filter
)

# -----------------------------------------------------------------------------
# Apply cascade cumulatively, recording stats per stage
# -----------------------------------------------------------------------------
summary_row <- list(motif_id = opt$motif_id, motif_width = opt$motif_width,
                     missing_input = FALSE, empty_input = FALSE)

current_idx <- seq_len(n_total)

for (s in cascade) {
    keep_in_full <- s$fn(dt)
    current_idx <- intersect(current_idx, keep_in_full)
    sub <- dt[current_idx]

    n_s <- length(current_idx)
    summary_row[[paste0("n_", s$label)]] <- n_s

    if (n_s == 0) {
        summary_row[[paste0("ubp_", s$label)]]          <- 0L
        summary_row[[paste0("score_mean_", s$label)]]    <- NA_real_
        summary_row[[paste0("score_median_", s$label)]]  <- NA_real_
        summary_row[[paste0("score_p90_", s$label)]]     <- NA_real_
    } else {
        summary_row[[paste0("ubp_", s$label)]] <-
            as.integer(count_unique_bp(sub$chrom, sub$start_1b, sub$end_1b))
        summary_row[[paste0("score_mean_", s$label)]] <-
            mean(sub$score, na.rm = TRUE)
        summary_row[[paste0("score_median_", s$label)]] <-
            as.numeric(median(sub$score, na.rm = TRUE))
        summary_row[[paste0("score_p90_", s$label)]] <-
            as.numeric(quantile(sub$score, 0.9, na.rm = TRUE))
    }
}

# Cache the post-cascade subset for qualifying calculations
final_idx <- current_idx
post <- dt[final_idx]
n_post <- nrow(post)
cat("  After full cascade:", n_post, "TFBSs.\n")

# -----------------------------------------------------------------------------
# Qualifying TFBSs at each distance threshold
# -----------------------------------------------------------------------------
# The 05_motif_annot schema gives:
#   snp_<hybrid>_overlap                         direct body overlap
#   snp_<hybrid>_flank10bp_overlap               1-10 bp from edge
#   snp_<hybrid>_flank25bp_overlap               11-25
#   snp_<hybrid>_flank50bp_overlap               26-50
#   snp_<hybrid>_flank100bp_overlap              51-100
#   snp_<hybrid>_flank200bp_overlap              101-200
#   snp_<hybrid>_flank400bp_overlap              201-400
#   snp_<hybrid>_flank800bp_overlap              401-800
#
# A "qualifying TFBS at threshold T" means: the TFBS has a SNP within T bp of
# its edge (or in its body for T=0). This is a CUMULATIVE test, so we OR the
# overlap columns up to and including the threshold.

build_qualifying_mask <- function(post_dt, hybrid, threshold) {
    if (nrow(post_dt) == 0) return(logical(0))

    cols_needed <- character(0)
    direct_col <- paste0("snp_", hybrid, "_overlap")
    cols_needed <- c(cols_needed, direct_col)

    flank_dists <- c(10, 25, 50, 100, 200, 400, 800)
    for (fd in flank_dists) {
        if (fd <= threshold || threshold == 0 && fd == 0) {
            cols_needed <- c(cols_needed,
                paste0("snp_", hybrid, "_flank", fd, "bp_overlap"))
        }
    }
    if (threshold == 0) {
        cols_needed <- direct_col
    }

    cols_present <- intersect(cols_needed, names(post_dt))
    if (length(cols_present) == 0) return(rep(FALSE, nrow(post_dt)))

    mask <- rep(FALSE, nrow(post_dt))
    for (cn in cols_present) {
        mask <- mask | isTRUE_safe(post_dt[[cn]])
    }
    mask
}

q_thresholds <- c(0, 10, 25, 50, 100, 200, 400, 800)
q_labels <- ifelse(q_thresholds == 0, "overlap",
                    paste0("flank", q_thresholds))

for (i in seq_along(q_thresholds)) {
    th <- q_thresholds[i]
    lab <- q_labels[i]
    mask <- build_qualifying_mask(post, opt$hybrid, th)
    n_q <- sum(mask)
    summary_row[[paste0("n_qualifying_", lab)]] <- n_q
    if (n_q == 0) {
        summary_row[[paste0("score_mean_qualifying_", lab)]]   <- NA_real_
        summary_row[[paste0("score_median_qualifying_", lab)]] <- NA_real_
    } else {
        scores <- post$score[mask]
        summary_row[[paste0("score_mean_qualifying_", lab)]]   <- mean(scores, na.rm = TRUE)
        summary_row[[paste0("score_median_qualifying_", lab)]] <-
            as.numeric(median(scores, na.rm = TRUE))
    }
}

cat("  Qualifying TFBSs (overlap):", summary_row$n_qualifying_overlap,
    " (flank200):", summary_row$n_qualifying_flank200, "\n")

# -----------------------------------------------------------------------------
# Write summary row
# -----------------------------------------------------------------------------
summary_dt <- as.data.table(summary_row)
fwrite(summary_dt,
       file.path(opt$outdir, paste0(opt$motif_id, "_summary.tsv")),
       sep = "\t", na = "NA")
cat("  Wrote summary:", file.path(opt$outdir, paste0(opt$motif_id, "_summary.tsv")), "\n")

# -----------------------------------------------------------------------------
# Per-TFBS scores at the post-cascade stage (for PWM dist plots)
# Save lightweight: chrom, start, end, strand, score, motif_hit_id (if present),
# qualifying flags at each threshold.
# -----------------------------------------------------------------------------
score_dt <- post[, .(
    chrom, start = start_1b, end = end_1b, strand, score
)]
if ("motif_hit_id" %in% names(post)) {
    score_dt[, motif_hit_id := post$motif_hit_id]
}
for (i in seq_along(q_thresholds)) {
    th <- q_thresholds[i]
    lab <- q_labels[i]
    score_dt[[paste0("qualifying_", lab)]] <- build_qualifying_mask(post, opt$hybrid, th)
}
score_dt[, motif_id := opt$motif_id]
fwrite(score_dt,
       file.path(opt$outdir, paste0(opt$motif_id, "_scores.tsv.gz")),
       sep = "\t", compress = "gzip", na = "NA")
cat("  Wrote per-TFBS scores: ", nrow(score_dt), "rows\n")
