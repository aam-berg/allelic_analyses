# =============================================================================
# 02f_gene_expression.R — Gene expression annotation (SENSE + ANTISENSE)
# =============================================================================
#
# OUTPUT COLUMNS:
#   For each orientation in {sense, antisense}:
#     expression_n_genes_{orient}            # genes overlapped on this strand
#     expression_max_gene_id_{orient}        gene_id of the most-expressed gene
#     expression_max_gene_name_{orient}      its gene_name
#     expression_tpm_max_{orient}            tpm_mean of the most-expressed gene
#     expression_tpm_sd_at_max_{orient}      tpm_sd of that gene
#     expression_counts_sum_{orient}         sum of counts across overlapping genes
#
# PREREQUISITE:
#   Module 02b must have run first to populate
#   intragenic_gene_ids_{sense,antisense}.
#
# WHY MAX:
#   The user's stated goal is to filter motifs that aren't on expressed
#   genes. Taking the max of tpm_mean across overlapping genes is the
#   right summary: if any gene the motif sits in is highly expressed,
#   the motif is in transcribed sequence regardless of which gene
#   "wins" in the overlap. Counts sum tells you the raw read support
#   pooled across overlapping genes.
#
# JOIN MECHANICS:
#   For each motif row with intragenic_gene_ids_<orient> = "ID1|ID2|...",
#   we split on "|", look up each gene_id in the expression data.table
#   (which is keyed by gene_id), and aggregate.
#
# IF EXPRESSION DATA NOT LOADED:
#   All columns are written as NA. The pipeline still completes; downstream
#   filtering will simply skip the gene-expression criterion.
# =============================================================================

annotate_gene_expression <- function(motif_dt, resources) {

    expr_dt <- resources$expr_dt   # may be NULL
    have_expr <- !is.null(expr_dt) && nrow(expr_dt) > 0

    if (!have_expr) {
        cat("    [WARN] No RNA-seq expression table loaded. Adding NA columns.\n")
    }

    for (orient in c("sense", "antisense")) {
        col_in <- paste0("intragenic_gene_ids_", orient)

        # Initialize output columns
        motif_dt[, paste0("expression_n_genes_",         orient) := 0L]
        motif_dt[, paste0("expression_max_gene_id_",     orient) := NA_character_]
        motif_dt[, paste0("expression_max_gene_name_",   orient) := NA_character_]
        motif_dt[, paste0("expression_tpm_max_",         orient) := NA_real_]
        motif_dt[, paste0("expression_tpm_sd_at_max_",   orient) := NA_real_]
        motif_dt[, paste0("expression_counts_sum_",      orient) := NA_integer_]

        if (!have_expr) next
        if (!(col_in %in% names(motif_dt))) {
            cat("    [WARN] Column", col_in, "missing. Did 02b run? Skipping",
                orient, "\n")
            next
        }

        # Process motifs that have at least one overlapping gene on this strand
        rows_with_genes <- which(!is.na(motif_dt[[col_in]]))
        if (length(rows_with_genes) == 0) {
            cat("    No motifs with overlapping ", orient, " genes.\n", sep = "")
            next
        }

        cat("    ", orient, ": ", length(rows_with_genes), " motifs with ", orient,
            " gene overlap; joining against expression table...\n", sep = "")

        # Process in vectorized chunks: build long table of (row_id, gene_id),
        # join against expr_dt, then aggregate per row_id.
        gene_lists <- strsplit(motif_dt[rows_with_genes][[col_in]], "|", fixed = TRUE)
        rep_rows <- rep(motif_dt$row_id[rows_with_genes], lengths(gene_lists))
        all_gids <- unlist(gene_lists)

        # Join against expression table
        long_dt <- data.table(row_id = rep_rows, gene_id = all_gids)
        joined <- expr_dt[long_dt, on = "gene_id",
                           .(row_id, gene_id, gene_name,
                             tpm_mean, tpm_sd, counts_sum)]

        # Aggregate per row_id:
        #   - n_genes = number of joined gene records (including unmatched
        #     gene_ids; NAs propagate via tpm_mean)
        #   - tpm_max = max of tpm_mean
        #   - max_gene_id, max_gene_name = the gene_id/name of the row achieving
        #     the max tpm_mean
        #   - tpm_sd_at_max = tpm_sd of the max-tpm gene
        #   - counts_sum_total = sum of counts_sum
        agg <- joined[, {
            n <- .N
            valid <- !is.na(tpm_mean)
            if (any(valid)) {
                imax <- which.max(replace(tpm_mean, !valid, -Inf))
                .(n_genes        = n,
                  max_gene_id    = gene_id[imax],
                  max_gene_name  = gene_name[imax],
                  tpm_max        = tpm_mean[imax],
                  tpm_sd_at_max  = tpm_sd[imax],
                  counts_sum     = sum(counts_sum, na.rm = TRUE))
            } else {
                .(n_genes        = n,
                  max_gene_id    = NA_character_,
                  max_gene_name  = NA_character_,
                  tpm_max        = NA_real_,
                  tpm_sd_at_max  = NA_real_,
                  counts_sum     = NA_integer_)
            }
        }, by = row_id]

        # Write to motif_dt
        motif_dt[agg, on = "row_id", paste0("expression_n_genes_",         orient) := i.n_genes]
        motif_dt[agg, on = "row_id", paste0("expression_max_gene_id_",     orient) := i.max_gene_id]
        motif_dt[agg, on = "row_id", paste0("expression_max_gene_name_",   orient) := i.max_gene_name]
        motif_dt[agg, on = "row_id", paste0("expression_tpm_max_",         orient) := i.tpm_max]
        motif_dt[agg, on = "row_id", paste0("expression_tpm_sd_at_max_",   orient) := i.tpm_sd_at_max]
        motif_dt[agg, on = "row_id", paste0("expression_counts_sum_",      orient) := as.integer(i.counts_sum)]
    }

    # Brief summary
    for (orient in c("sense", "antisense")) {
        col <- paste0("expression_tpm_max_", orient)
        v <- motif_dt[[col]]
        n_with_expr <- sum(!is.na(v))
        cat("    ", col, ": ", n_with_expr, " motifs with annotated expression",
            " (", round(100*n_with_expr/nrow(motif_dt), 1), "%)\n", sep = "")
        if (n_with_expr > 0) {
            cat("      median tpm_max=", signif(median(v, na.rm = TRUE), 3),
                ", n with tpm_max>=1: ", sum(v >= 1, na.rm = TRUE), "\n", sep = "")
        }
    }
    invisible(motif_dt)
}
