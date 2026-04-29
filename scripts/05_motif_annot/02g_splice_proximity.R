# =============================================================================
# 02g_splice_proximity.R — Splice donor + acceptor proximity (STRAND-AWARE)
# =============================================================================
#
# OUTPUT: 8 query cells × 5 sub-columns = 40 columns total.
#
# QUERY CELLS:
#   site_type ∈ {donor, acceptor}
#   orientation ∈ {sense, antisense}      (relative to motif strand)
#   direction ∈ {upstream, downstream}    (5' or 3' of motif midpoint)
#
# COLUMNS (per cell):
#   {site_type}_{orientation}_{direction}_dist           int  (always positive bp)
#   {site_type}_{orientation}_{direction}_source         "GTF" | "SJ_only" | "GTF_and_SJ"
#   {site_type}_{orientation}_{direction}_support        int  (max n_unique reads from SJ; NA if GTF-only)
#   {site_type}_{orientation}_{direction}_motif_type     str  (e.g. "GT/AG"; "" for GTF-only)
#   {site_type}_{orientation}_{direction}_gene_names     str  (pipe-sep; "" for SJ-only with no GTF coverage)
#
# EXAMPLE COLUMN NAMES:
#   donor_sense_upstream_dist
#   acceptor_antisense_downstream_source
#   donor_sense_downstream_support
#   ...
#
# WHY 40 COLUMNS:
#   The user explicitly wants to track strand AND directionality of each
#   junction. To make all cells filterable independently, we use one column
#   per cell × per attribute. A wide schema is better than nested types here
#   because TSV is the file format and downstream code can drop columns easily.
#
# THE LOOKUP MECHANICS:
#   `resources$donors_gr` and `resources$acceptors_gr` are GRanges with
#   strand preserved AND each carrying mcols(): source, support_n_unique,
#   motif_type, gene_names. This module subsets by strand to get sense vs
#   antisense, then uses precede()/follow() in coordinate space to get the
#   nearest site in each direction from each motif midpoint.
#
# DONOR vs ACCEPTOR DEFINITION (recap):
#   Donor    = first base of intron in transcription direction
#   Acceptor = last base of intron in transcription direction
#   For "+" intron: donor = intron_start, acceptor = intron_end
#   For "-" intron: donor = intron_end, acceptor = intron_start
# =============================================================================

# Subset a junction-site GRanges by strand and drop strand info
# (so precede/follow work in coordinate space)
sites_by_strand_unstranded <- function(sites_gr, want_strand) {
    sub <- sites_gr[strand(sites_gr) == want_strand]
    if (length(sub) > 0) strand(sub) <- "*"
    sub
}

# Compute nearest site for one (motif_strand, site_strand, direction) combo
# Returns data.table: row_id, dist, source, support, motif_type, gene_names
nearest_site_one_combo <- function(motif_dt, motif_gr,
                                    motif_strand, site_strand, direction,
                                    sites_gr) {
    sel_motif <- motif_gr[strand(motif_gr) == motif_strand]
    if (length(sel_motif) == 0) return(NULL)

    # Use motif midpoint
    midpts_pos <- (start(sel_motif) + end(sel_motif)) %/% 2L
    midpts_gr <- GRanges(
        seqnames = seqnames(sel_motif),
        ranges = IRanges(start = midpts_pos, width = 1),
        strand = "*",
        row_id = sel_motif$row_id
    )

    sites_sub <- sites_by_strand_unstranded(sites_gr, site_strand)
    if (length(sites_sub) == 0) return(NULL)

    shared <- intersect(seqlevels(midpts_gr), seqlevels(sites_sub))
    if (length(shared) == 0) return(NULL)
    midpts_gr <- keepSeqlevels(midpts_gr, shared, pruning.mode = "coarse")
    sites_sub <- keepSeqlevels(sites_sub, shared, pruning.mode = "coarse")

    # follow() = nearest at LOWER coord, precede() = HIGHER coord
    if ((motif_strand == "+" && direction == "upstream") ||
        (motif_strand == "-" && direction == "downstream")) {
        idx <- follow(midpts_gr, sites_sub)
    } else {
        idx <- precede(midpts_gr, sites_sub)
    }

    valid <- !is.na(idx)
    if (sum(valid) == 0) return(NULL)

    rids <- midpts_gr$row_id[valid]
    midpts <- start(midpts_gr)[valid]
    site_hits <- sites_sub[idx[valid]]

    data.table(
        row_id     = rids,
        dist       = abs(midpts - start(site_hits)),
        source     = site_hits$source,
        support    = site_hits$support_n_unique,
        motif_type = site_hits$motif_type,
        gene_names = site_hits$gene_names
    )
}

# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------
annotate_splice_proximity <- function(motif_dt, motif_gr, resources) {

    sites_grs <- list(
        donor    = resources$donors_gr,
        acceptor = resources$acceptors_gr
    )

    # Initialize all 40 columns
    for (st in c("donor", "acceptor")) {
        for (orient in c("sense", "antisense")) {
            for (direction in c("upstream", "downstream")) {
                base <- paste0(st, "_", orient, "_", direction)
                motif_dt[, paste0(base, "_dist")       := NA_integer_]
                motif_dt[, paste0(base, "_source")     := NA_character_]
                motif_dt[, paste0(base, "_support")    := NA_integer_]
                motif_dt[, paste0(base, "_motif_type") := NA_character_]
                motif_dt[, paste0(base, "_gene_names") := NA_character_]
            }
        }
    }

    if (is.null(sites_grs$donor) || length(sites_grs$donor) == 0) {
        cat("    [WARN] No donor sites loaded; columns left NA.\n")
    }
    if (is.null(sites_grs$acceptor) || length(sites_grs$acceptor) == 0) {
        cat("    [WARN] No acceptor sites loaded; columns left NA.\n")
    }

    # Iterate over the 16 (motif_strand × site_type × orient × direction)
    # combinations. For each, look up the nearest site and write into the
    # appropriate per-cell columns (where motif rows of that strand land).
    for (motif_strand in c("+", "-")) {
        for (st in c("donor", "acceptor")) {
            sites_gr <- sites_grs[[st]]
            if (is.null(sites_gr) || length(sites_gr) == 0) next

            for (orient in c("sense", "antisense")) {
                site_strand <- if (orient == "sense")
                                    motif_strand
                                else
                                    if (motif_strand == "+") "-" else "+"

                for (direction in c("upstream", "downstream")) {
                    hit_dt <- nearest_site_one_combo(
                        motif_dt, motif_gr,
                        motif_strand, site_strand, direction,
                        sites_gr
                    )
                    if (is.null(hit_dt)) next

                    base <- paste0(st, "_", orient, "_", direction)
                    cd <- paste0(base, "_dist")
                    cs <- paste0(base, "_source")
                    cu <- paste0(base, "_support")
                    cm <- paste0(base, "_motif_type")
                    cg <- paste0(base, "_gene_names")

                    idx <- match(hit_dt$row_id, motif_dt$row_id)
                    set(motif_dt, i = idx, j = cd, value = as.integer(hit_dt$dist))
                    set(motif_dt, i = idx, j = cs, value = hit_dt$source)
                    set(motif_dt, i = idx, j = cu, value = as.integer(hit_dt$support))
                    set(motif_dt, i = idx, j = cm, value = hit_dt$motif_type)
                    set(motif_dt, i = idx, j = cg, value = hit_dt$gene_names)
                }
            }
        }
    }

    # Brief summary: median dist + source breakdown for each of the 8 cells
    for (st in c("donor", "acceptor")) {
        for (orient in c("sense", "antisense")) {
            for (direction in c("upstream", "downstream")) {
                base <- paste0(st, "_", orient, "_", direction)
                cd <- paste0(base, "_dist")
                cs <- paste0(base, "_source")
                v <- motif_dt[[cd]]
                if (any(!is.na(v))) {
                    src_table <- table(motif_dt[[cs]])
                    src_str <- paste(names(src_table), src_table, sep = "=", collapse = ", ")
                    cat(sprintf("    %s: median=%d bp, n=%d, sources={%s}\n",
                                base,
                                as.integer(median(v, na.rm = TRUE)),
                                sum(!is.na(v)),
                                src_str))
                }
            }
        }
    }
    invisible(motif_dt)
}
