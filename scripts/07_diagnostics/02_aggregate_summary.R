#!/usr/bin/env Rscript
# =============================================================================
# 02_aggregate_summary.R — Merge per-motif summary rows into the master table
# =============================================================================
#
# Produces:
#   {TABLES_DIR}/archetype_summary.tsv
#       One row per motif archetype. Columns:
#         motif_id, motif_width, missing_input, empty_input,
#         n_<stage>, ubp_<stage>, score_mean_<stage>, ... (per cascade stage)
#         n_qualifying_<thresh>, score_mean_qualifying_<thresh>, ...
#         max_tpm_across_mapped, max_tf_name, n_mapped_tfs   (joined from metadata)
#
# Also writes:
#   {TABLES_DIR}/qualifying_per_archetype.tsv     compact: motif_id, motif_width,
#                                                  + n_qualifying_<thresh> for each
#                                                  threshold + max_tpm.
#
# USAGE:
#   Rscript 02_aggregate_summary.R \
#     --summary_dir .../per_archetype_summaries/ \
#     --pwm_meme .../consensus_pwms.meme \
#     --metadata .../metadata.tsv \
#     --expression .../gene_expression_summary.tsv \
#     --outdir .../tables/
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option("--summary_dir",   type = "character"),
    make_option("--pwm_meme",      type = "character"),
    make_option("--metadata",      type = "character"),
    make_option("--expression",    type = "character"),
    make_option("--outdir",        type = "character"),
    make_option("--lib_path",      type = "character", default = "")
)
opt <- parse_args(OptionParser(option_list = option_list))

lib_path <- if (nchar(opt$lib_path) > 0) opt$lib_path else
            file.path(dirname(sub("--file=", "",
                grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))),
                "lib_diagnostics.R")
source(lib_path)

cat("============================================================\n")
cat("02_aggregate_summary.R\n")
cat("============================================================\n\n")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Read all per-motif summaries
# -----------------------------------------------------------------------------
files <- list.files(opt$summary_dir, pattern = "_summary\\.tsv$",
                     full.names = TRUE)
cat("Found", length(files), "per-motif summary files\n")
if (length(files) == 0) {
    stop("No summary files in ", opt$summary_dir)
}
all_summaries <- rbindlist(lapply(files, fread), fill = TRUE)
cat("Aggregated:", nrow(all_summaries), "archetypes,",
    ncol(all_summaries), "columns\n")

# -----------------------------------------------------------------------------
# Ensure motif_width is filled (re-derive from PWM file if needed)
# -----------------------------------------------------------------------------
if (any(is.na(all_summaries$motif_width)) && file.exists(opt$pwm_meme)) {
    widths <- get_archetype_pwm_widths(opt$pwm_meme)
    all_summaries[widths, motif_width := i.motif_width, on = "motif_id"]
}

# -----------------------------------------------------------------------------
# Join TF metadata + max TPM
# -----------------------------------------------------------------------------
md <- load_tf_metadata(opt$metadata)
expr <- load_gene_expression(opt$expression)
if (nrow(md) > 0 && nrow(expr) > 0) {
    tpm_per_arch <- compute_max_tpm_per_archetype(md, expr)
    all_summaries <- tpm_per_arch[all_summaries, on = "motif_id"]
    cat("Joined TF metadata: archetypes with mapped TFs =",
        sum(!is.na(all_summaries$max_tpm_across_mapped)), "/",
        nrow(all_summaries), "\n")
} else {
    all_summaries[, n_mapped_tfs := NA_integer_]
    all_summaries[, n_with_tpm := NA_integer_]
    all_summaries[, max_tpm_across_mapped := NA_real_]
    all_summaries[, max_tf_name := NA_character_]
}

# -----------------------------------------------------------------------------
# Reorder columns
# -----------------------------------------------------------------------------
canonical_first <- c("motif_id", "motif_width",
                      "n_mapped_tfs", "n_with_tpm",
                      "max_tpm_across_mapped", "max_tf_name",
                      "missing_input", "empty_input")
canonical_present <- intersect(canonical_first, names(all_summaries))
others <- setdiff(names(all_summaries), canonical_present)
setcolorder(all_summaries, c(canonical_present, others))

# Sort rows by motif_id
setorder(all_summaries, motif_id)

# -----------------------------------------------------------------------------
# Write outputs
# -----------------------------------------------------------------------------
out_main <- file.path(opt$outdir, "archetype_summary.tsv")
fwrite(all_summaries, out_main, sep = "\t", na = "NA")
cat("Wrote:", out_main, "\n")

# Compact qualifying-counts table
q_cols <- grep("^n_qualifying_", names(all_summaries), value = TRUE)
keep <- c("motif_id", "motif_width", "max_tpm_across_mapped", "max_tf_name",
           q_cols)
keep <- intersect(keep, names(all_summaries))
qual_dt <- all_summaries[, ..keep]
out_qual <- file.path(opt$outdir, "qualifying_per_archetype.tsv")
fwrite(qual_dt, out_qual, sep = "\t", na = "NA")
cat("Wrote:", out_qual, "\n")

cat("\nDone.\n")
