# =============================================================================
# 02d_atac_overlap.R — ATAC consensus peak overlap (UNSTRANDED)
# =============================================================================
#
# OUTPUT COLUMNS:
#   atac_consensus_overlap        boolean
#   atac_consensus_n_replicates   integer (max replicate-support count among
#                                  consensus peaks the motif overlaps;
#                                  NA if no overlap)
#
# WHY UNSTRANDED:
#   ATAC peaks are not strand-typed (they reflect chromatin accessibility,
#   which doesn't have inherent directionality). So we ignore strand entirely
#   in this overlap. Same scheme as the original pipeline.
#
# CHANGE FROM OLD VERSION:
#   Previously, this looked at multiple per-sample 4DN bigBed files
#   (atac_<SAMPLE>_overlap, atac_<SAMPLE>_signalValue, etc.) after a
#   mm10->mm39 liftOver. Now we use a single F121-9-specific consensus
#   peak set computed by 03_atac/, in mm39 already, with no liftOver
#   needed and no per-sample columns.
# =============================================================================

annotate_atac_overlap <- function(motif_dt, motif_gr, resources) {

    motif_dt[, atac_consensus_overlap        := FALSE]
    motif_dt[, atac_consensus_n_replicates   := NA_integer_]

    if (is.null(resources$atac_gr) || length(resources$atac_gr) == 0) {
        cat("    [WARN] ATAC consensus GRanges not loaded; columns left NA.\n")
        return(invisible(motif_dt))
    }

    atac_gr <- resources$atac_gr

    # Build an unstranded copy of motif_gr (preserves row_id)
    motif_unstranded <- motif_gr
    strand(motif_unstranded) <- "*"

    # Align seqlevels
    shared <- intersect(seqlevels(motif_unstranded), seqlevels(atac_gr))
    if (length(shared) == 0) {
        cat("    [WARN] No shared seqlevels between motifs and ATAC; columns left NA.\n")
        return(invisible(motif_dt))
    }
    motif_unstranded <- keepSeqlevels(motif_unstranded, shared, pruning.mode = "coarse")
    atac_gr <- keepSeqlevels(atac_gr, shared, pruning.mode = "coarse")

    olap <- findOverlaps(motif_unstranded, atac_gr, ignore.strand = TRUE)
    if (length(olap) == 0) {
        cat("    No motifs overlap ATAC consensus peaks.\n")
        return(invisible(motif_dt))
    }

    # When a motif overlaps multiple peaks, take max n_replicates
    hit_dt <- data.table(
        row_id = motif_unstranded$row_id[queryHits(olap)],
        n_reps = atac_gr$n_replicates[subjectHits(olap)]
    )
    agg <- hit_dt[, .(n_reps_max = max(n_reps)), by = row_id]

    motif_dt[agg, on = "row_id", `:=`(
        atac_consensus_overlap        = TRUE,
        atac_consensus_n_replicates   = i.n_reps_max
    )]

    cat("    atac_consensus_overlap: ", sum(motif_dt$atac_consensus_overlap),
        " (", round(100*mean(motif_dt$atac_consensus_overlap), 1), "%)\n", sep = "")

    if (any(motif_dt$atac_consensus_overlap)) {
        cat("    Replicate-support distribution among overlapping motifs:\n")
        rep_table <- table(motif_dt[atac_consensus_overlap == TRUE,
                                     atac_consensus_n_replicates])
        for (k in names(rep_table)) {
            cat("      n_replicates=", k, ": ", rep_table[k], "\n", sep = "")
        }
    }
    invisible(motif_dt)
}
