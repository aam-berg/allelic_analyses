# =============================================================================
# 02c_strand_aware_tss.R — Nearest upstream/downstream TSS (SENSE + ANTISENSE)
# =============================================================================
#
# OUTPUT COLUMNS (8 distance + 24 metadata = 32 total):
#   For each direction (upstream, downstream) × strand orientation (sense, antisense):
#     {direction}_tss_dist_{orient}              integer bp
#     {direction}_tss_gene_names_{orient}        pipe-sep
#     {direction}_tss_tx_ids_{orient}            pipe-sep
#     {direction}_tss_tsls_{orient}              pipe-sep
#
# DEFINITIONS:
#   "Upstream TSS"   = nearest TSS in the 5' direction of the motif
#   "Downstream TSS" = nearest TSS in the 3' direction of the motif
#   For "+" motif: 5' = lower coord, 3' = higher coord
#   For "-" motif: 5' = higher coord, 3' = lower coord
#
#   "Sense"     = TSS on the same strand as the motif
#   "Antisense" = TSS on the opposite strand
#
# IMPLEMENTATION:
#   We use precede() / follow() which operate in coordinate space when input
#   strands are "*". We do this twice — once with the SENSE-stranded subset
#   of TSSes, once with the ANTISENSE subset — and compute distances from
#   each motif's midpoint.
#
#   Per-motif midpoint is used so that motif width doesn't bias distance
#   measurements. (TSS proximity is a regional feature; we don't need
#   single-base precision here.)
#
# COLLAPSING WHEN MULTIPLE TXs SHARE A TSS:
#   Because we're using TSS-per-transcript GRanges, multiple transcripts
#   can map to the same chrom/position/strand. We don't pre-collapse — instead
#   findNearest finds one transcript's TSS index and reports its annotations.
#   For multi-tx TSSes we'd report just one transcript_id; this is a
#   simplification but the dist is correct, and a follow-up join against
#   tx_meta could enumerate all tx at that position if needed.
# =============================================================================

# Subset TSS GRanges to a specific strand, dropping strand info (set to "*")
# so that precede/follow operate in coordinate space.
tss_by_strand_unstranded <- function(tss_gr, want_strand) {
    sel <- tss_gr[strand(tss_gr) == want_strand]
    if (length(sel) > 0) strand(sel) <- "*"
    sel
}

# Compute nearest TSS in one (motif_strand, tss_strand, direction) combination.
# Returns a data.table with row_id, dist, gene_name, transcript_id, tsl.
nearest_tss_one_combo <- function(motif_dt, motif_gr,
                                   motif_strand, tss_strand, direction,
                                   tss_gr) {
    # Filter motif_gr to the given motif_strand
    sel_motif <- motif_gr[strand(motif_gr) == motif_strand]
    if (length(sel_motif) == 0) return(NULL)

    # Use motif midpoint (single-base GRanges)
    midpoints_pos <- (start(sel_motif) + end(sel_motif)) %/% 2L
    midpoints_gr <- GRanges(
        seqnames = seqnames(sel_motif),
        ranges = IRanges(start = midpoints_pos, width = 1),
        strand = "*",
        row_id = sel_motif$row_id
    )

    # Filter TSSes to the requested strand, drop strand info
    tss_sub <- tss_by_strand_unstranded(tss_gr, tss_strand)
    if (length(tss_sub) == 0) return(NULL)

    # Align seqlevels
    shared <- intersect(seqlevels(midpoints_gr), seqlevels(tss_sub))
    if (length(shared) == 0) return(NULL)
    midpoints_gr <- keepSeqlevels(midpoints_gr, shared, pruning.mode = "coarse")
    tss_sub <- keepSeqlevels(tss_sub, shared, pruning.mode = "coarse")

    # follow() = nearest subject at LOWER coord than query
    # precede() = nearest subject at HIGHER coord than query
    if ((motif_strand == "+" && direction == "upstream") ||
        (motif_strand == "-" && direction == "downstream")) {
        idx <- follow(midpoints_gr, tss_sub)
    } else {
        idx <- precede(midpoints_gr, tss_sub)
    }

    valid <- !is.na(idx)
    if (sum(valid) == 0) return(NULL)

    rids <- midpoints_gr$row_id[valid]
    midpts <- start(midpoints_gr)[valid]
    tss_hits <- tss_sub[idx[valid]]

    data.table(
        row_id = rids,
        dist = abs(midpts - start(tss_hits)),
        gene_name = tss_hits$gene_name,
        transcript_id = tss_hits$transcript_id,
        tsl = tss_hits$transcript_support_level
    )
}

# Top-level entry point
annotate_strand_aware_tss <- function(motif_dt, motif_gr, resources) {

    tss_gr <- resources$tss_gr

    # Initialize all 32 columns with NAs
    for (orient in c("sense", "antisense")) {
        for (direction in c("upstream", "downstream")) {
            base <- paste0(direction, "_tss")
            motif_dt[, paste0(base, "_dist_", orient)         := NA_integer_]
            motif_dt[, paste0(base, "_gene_names_", orient)   := NA_character_]
            motif_dt[, paste0(base, "_tx_ids_", orient)       := NA_character_]
            motif_dt[, paste0(base, "_tsls_", orient)         := NA_character_]
        }
    }

    # 8 combinations: (motif_strand: +/-) × (orientation: sense/antisense) × (direction)
    # For each motif strand X:
    #   sense -> tss_strand = X
    #   antisense -> tss_strand = flip(X)
    for (motif_strand in c("+", "-")) {
        for (orient in c("sense", "antisense")) {
            tss_strand <- if (orient == "sense") motif_strand else if (motif_strand == "+") "-" else "+"

            for (direction in c("upstream", "downstream")) {
                hit_dt <- nearest_tss_one_combo(motif_dt, motif_gr,
                                                 motif_strand, tss_strand, direction,
                                                 tss_gr)
                if (is.null(hit_dt)) next

                base <- paste0(direction, "_tss")
                col_dist <- paste0(base, "_dist_",       orient)
                col_gnm  <- paste0(base, "_gene_names_", orient)
                col_tx   <- paste0(base, "_tx_ids_",     orient)
                col_tsl  <- paste0(base, "_tsls_",       orient)

                # Map row_id to motif_dt row positions
                idx <- match(hit_dt$row_id, motif_dt$row_id)
                set(motif_dt, i = idx, j = col_dist, value = as.integer(hit_dt$dist))
                set(motif_dt, i = idx, j = col_gnm,  value = hit_dt$gene_name)
                set(motif_dt, i = idx, j = col_tx,   value = hit_dt$transcript_id)
                set(motif_dt, i = idx, j = col_tsl,  value = hit_dt$tsl)
            }
        }
    }

    # Brief summary
    for (orient in c("sense", "antisense")) {
        for (direction in c("upstream", "downstream")) {
            col <- paste0(direction, "_tss_dist_", orient)
            v <- motif_dt[[col]]
            if (any(!is.na(v))) {
                cat(sprintf("    %s: median=%d bp, max=%d bp (n=%d)\n",
                            col,
                            as.integer(median(v, na.rm = TRUE)),
                            as.integer(max(v, na.rm = TRUE)),
                            sum(!is.na(v))))
            }
        }
    }
    invisible(motif_dt)
}
