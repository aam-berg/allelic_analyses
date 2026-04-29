# =============================================================================
# 02b_strand_aware_genic.R — Intragenic, promoter, genic regions (SENSE + ANTISENSE)
# =============================================================================
#
# CRITICAL DESIGN POINT:
#   Every annotation here is computed TWICE — once for "sense" (motif strand
#   matches feature strand) and once for "antisense" (motif strand opposite
#   to feature strand). This makes orientation-flipping a column relabeling
#   downstream, not a re-computation.
#
# OUTPUT COLUMNS (added to motif_dt):
#   intragenic_sense, intragenic_antisense                 (boolean)
#   intragenic_gene_ids_sense, intragenic_gene_ids_antisense          (pipe-sep)
#   intragenic_gene_names_sense, intragenic_gene_names_antisense      (pipe-sep)
#   intragenic_gene_types_sense, intragenic_gene_types_antisense      (pipe-sep)
#   is_promoter_sense, is_promoter_antisense               (boolean)
#   promoter_transcript_ids_sense, promoter_transcript_ids_antisense  (pipe-sep)
#   promoter_gene_names_sense, promoter_gene_names_antisense          (pipe-sep)
#   promoter_tsls_sense, promoter_tsls_antisense                       (pipe-sep)
#   genic_regions_sense, genic_regions_antisense                       (e.g. "5UTR|exon")
#
# STRAND-AWARE OVERLAP MECHANICS:
#   In Bioconductor, findOverlaps(query, subject, ignore.strand=FALSE) only
#   counts hits where strands AGREE (or one is "*"). That gives us the SENSE
#   relationship for free. For ANTISENSE, we flip the strand of one set:
#     query_flipped <- query
#     strand(query_flipped) <- ifelse(strand(query) == "+", "-", "+")
#     antisense_hits <- findOverlaps(query_flipped, subject, ignore.strand=FALSE)
#
#   We use this pattern throughout the module via the helpers
#   `strand_aware_overlaps()` and `flip_strand()`.
#
# PROMOTER DEFINITION:
#   For each transcript with TSS at position T on strand S, the promoter
#   region is [T - PROMOTER_UPSTREAM, T + PROMOTER_DOWNSTREAM] in the
#   transcription direction. Since we use Bioconductor's promoters() with
#   strand-aware GRanges, this is computed correctly.
#   PROMOTER_UPSTREAM/PROMOTER_DOWNSTREAM are passed in via `resources$opt`.
# =============================================================================

# ---- Helpers ----------------------------------------------------------------

flip_strand <- function(gr) {
    new_strand <- as.character(strand(gr))
    new_strand[new_strand == "+"] <- "PENDING_TO_MINUS"
    new_strand[new_strand == "-"] <- "+"
    new_strand[new_strand == "PENDING_TO_MINUS"] <- "-"
    strand(gr) <- new_strand
    gr
}

# Find overlaps between motifs and features under strict sense or antisense
# rules. Returns a Hits object using the ORIGINAL motif_gr indices.
strand_aware_overlaps <- function(motif_gr, features_gr, mode = c("sense", "antisense")) {
    mode <- match.arg(mode)
    if (length(motif_gr) == 0 || length(features_gr) == 0) {
        return(Hits(from = integer(0), to = integer(0),
                    nLnode = length(motif_gr), nRnode = length(features_gr),
                    sort.by.query = TRUE))
    }
    if (mode == "sense") {
        findOverlaps(motif_gr, features_gr, ignore.strand = FALSE)
    } else {
        findOverlaps(flip_strand(motif_gr), features_gr, ignore.strand = FALSE)
    }
}

# Ensure shared seqlevels between two GRanges before findOverlaps
align_seqlevels <- function(a, b) {
    shared <- intersect(seqlevels(a), seqlevels(b))
    if (length(shared) == 0) {
        warning("No shared seqlevels.")
        return(list(a = a[seqnames(a) %in% character(0)],
                    b = b[seqnames(b) %in% character(0)]))
    }
    a <- keepSeqlevels(a, shared, pruning.mode = "coarse")
    b <- keepSeqlevels(b, shared, pruning.mode = "coarse")
    list(a = a, b = b)
}

# =============================================================================
# Annotation 1: intragenic_*
# =============================================================================
annotate_intragenic <- function(motif_dt, motif_gr, transcripts_gr) {
    al <- align_seqlevels(motif_gr, transcripts_gr)
    motif_a <- al$a; tx_a <- al$b

    motif_dt[, `:=`(intragenic_sense           = FALSE,
                    intragenic_antisense       = FALSE,
                    intragenic_gene_ids_sense       = NA_character_,
                    intragenic_gene_ids_antisense   = NA_character_,
                    intragenic_gene_names_sense     = NA_character_,
                    intragenic_gene_names_antisense = NA_character_,
                    intragenic_gene_types_sense     = NA_character_,
                    intragenic_gene_types_antisense = NA_character_)]

    for (mode in c("sense", "antisense")) {
        hits <- strand_aware_overlaps(motif_a, tx_a, mode = mode)
        if (length(hits) == 0) next
        # Map the filtered GR's hits back to motif_dt$row_id
        rids <- motif_a$row_id[queryHits(hits)]
        sub_tx <- tx_a[subjectHits(hits)]
        # Pull gene_id; we'll map to gene_name and gene_type via tx_meta
        sub_gid <- sub_tx$gene_id

        # Build hit table
        hit_dt <- data.table(row_id = rids,
                             gene_id = sub_gid,
                             tx_id = sub_tx$tx_id)
        # Bring in gene_name and gene_type from tx_meta (via gene_id)
        # NB: tx_meta has one row per transcript; we just need a unique
        # gene_id -> (gene_name, gene_type) lookup, so collapse tx_meta first
        gid2info <- unique(tx_meta_global[, .(gene_id, gene_name, gene_type)])
        hit_dt <- merge(hit_dt, gid2info, by = "gene_id", all.x = TRUE)

        agg <- hit_dt[, .(
            gene_ids   = paste(unique(gene_id), collapse = "|"),
            gene_names = paste(unique(na.omit(gene_name)), collapse = "|"),
            gene_types = paste(unique(na.omit(gene_type)), collapse = "|")
        ), by = row_id]
        agg[gene_names == "", gene_names := NA_character_]
        agg[gene_types == "", gene_types := NA_character_]

        if (mode == "sense") {
            motif_dt[agg, on = "row_id", `:=`(
                intragenic_sense              = TRUE,
                intragenic_gene_ids_sense     = i.gene_ids,
                intragenic_gene_names_sense   = i.gene_names,
                intragenic_gene_types_sense   = i.gene_types
            )]
        } else {
            motif_dt[agg, on = "row_id", `:=`(
                intragenic_antisense              = TRUE,
                intragenic_gene_ids_antisense     = i.gene_ids,
                intragenic_gene_names_antisense   = i.gene_names,
                intragenic_gene_types_antisense   = i.gene_types
            )]
        }
    }

    cat("    intragenic_sense:     ", sum(motif_dt$intragenic_sense),
        " (", round(100*mean(motif_dt$intragenic_sense), 1), "%)\n", sep = "")
    cat("    intragenic_antisense: ", sum(motif_dt$intragenic_antisense),
        " (", round(100*mean(motif_dt$intragenic_antisense), 1), "%)\n", sep = "")
    invisible(motif_dt)
}

# =============================================================================
# Annotation 2: promoter overlap (TSS [-PROMOTER_UPSTREAM, +PROMOTER_DOWNSTREAM])
# =============================================================================
annotate_promoter <- function(motif_dt, motif_gr, tss_gr,
                               upstream, downstream, genome_fa) {
    # Build promoter GRanges with strand preserved
    prom_gr <- promoters(tss_gr, upstream = upstream, downstream = downstream)

    # Trim to chrom bounds from FASTA
    si <- seqinfo(genome_fa)
    si_std <- si[intersect(seqlevels(prom_gr), seqlevels(si))]
    seqinfo(prom_gr, new2old = match(seqlevels(si_std), seqlevels(prom_gr))) <- si_std
    prom_gr <- trim(prom_gr)

    al <- align_seqlevels(motif_gr, prom_gr)
    motif_a <- al$a; prom_a <- al$b

    motif_dt[, `:=`(is_promoter_sense              = FALSE,
                    is_promoter_antisense          = FALSE,
                    promoter_transcript_ids_sense      = NA_character_,
                    promoter_transcript_ids_antisense  = NA_character_,
                    promoter_gene_names_sense          = NA_character_,
                    promoter_gene_names_antisense      = NA_character_,
                    promoter_tsls_sense                = NA_character_,
                    promoter_tsls_antisense            = NA_character_)]

    for (mode in c("sense", "antisense")) {
        hits <- strand_aware_overlaps(motif_a, prom_a, mode = mode)
        if (length(hits) == 0) next
        rids <- motif_a$row_id[queryHits(hits)]
        prom_hit <- prom_a[subjectHits(hits)]

        hit_dt <- data.table(
            row_id = rids,
            transcript_id = prom_hit$transcript_id,
            gene_name = prom_hit$gene_name,
            tsl       = prom_hit$transcript_support_level
        )
        agg <- hit_dt[, .(
            transcript_ids = paste(unique(transcript_id), collapse = "|"),
            gene_names     = paste(unique(na.omit(gene_name)), collapse = "|"),
            tsls           = paste(unique(na.omit(tsl)), collapse = "|")
        ), by = row_id]
        agg[gene_names == "", gene_names := NA_character_]
        agg[tsls == "",       tsls := NA_character_]

        if (mode == "sense") {
            motif_dt[agg, on = "row_id", `:=`(
                is_promoter_sense                 = TRUE,
                promoter_transcript_ids_sense     = i.transcript_ids,
                promoter_gene_names_sense         = i.gene_names,
                promoter_tsls_sense               = i.tsls
            )]
        } else {
            motif_dt[agg, on = "row_id", `:=`(
                is_promoter_antisense              = TRUE,
                promoter_transcript_ids_antisense  = i.transcript_ids,
                promoter_gene_names_antisense      = i.gene_names,
                promoter_tsls_antisense            = i.tsls
            )]
        }
    }

    cat("    is_promoter_sense:     ", sum(motif_dt$is_promoter_sense),
        " (", round(100*mean(motif_dt$is_promoter_sense), 1), "%)\n", sep = "")
    cat("    is_promoter_antisense: ", sum(motif_dt$is_promoter_antisense),
        " (", round(100*mean(motif_dt$is_promoter_antisense), 1), "%)\n", sep = "")
    invisible(motif_dt)
}

# =============================================================================
# Annotation 3: genic_regions (5UTR/3UTR/CDS/exon/intron) — sense + antisense
# =============================================================================
annotate_genic_regions <- function(motif_dt, motif_gr, region_grs) {
    motif_dt[, `:=`(genic_regions_sense     = NA_character_,
                    genic_regions_antisense = NA_character_)]

    for (mode in c("sense", "antisense")) {
        # Per region, get unique row_ids that overlap
        per_region <- list()
        for (rtype in names(region_grs)) {
            gr <- region_grs[[rtype]]
            if (length(gr) == 0) next
            al <- align_seqlevels(motif_gr, gr)
            motif_a <- al$a; gr_a <- al$b
            hits <- strand_aware_overlaps(motif_a, gr_a, mode = mode)
            if (length(hits) == 0) next
            per_region[[rtype]] <- unique(motif_a$row_id[queryHits(hits)])
        }
        if (length(per_region) == 0) next

        # Build long-format then collapse
        long_dt <- rbindlist(lapply(names(per_region), function(rt) {
            data.table(row_id = per_region[[rt]], region = rt)
        }))
        agg <- long_dt[, .(
            regions = paste(sort(unique(region)), collapse = "|")
        ), by = row_id]

        col <- if (mode == "sense") "genic_regions_sense" else "genic_regions_antisense"
        motif_dt[agg, on = "row_id", (col) := i.regions]
    }

    cat("    genic_regions_sense:     ",
        sum(!is.na(motif_dt$genic_regions_sense)),
        " annotated\n", sep = "")
    cat("    genic_regions_antisense: ",
        sum(!is.na(motif_dt$genic_regions_antisense)),
        " annotated\n", sep = "")
    invisible(motif_dt)
}

# =============================================================================
# Top-level orchestrator
# =============================================================================
annotate_strand_aware_genic <- function(motif_dt, motif_gr, resources) {
    # Make tx_meta available to annotate_intragenic via global
    # (avoids passing it through every call)
    assign("tx_meta_global", resources$tx_meta, envir = .GlobalEnv)

    cat("  -- intragenic --\n")
    annotate_intragenic(motif_dt, motif_gr, resources$transcripts_gr)

    cat("  -- promoter --\n")
    annotate_promoter(motif_dt, motif_gr, resources$tss_gr,
                       upstream   = resources$opt$promoter_upstream,
                       downstream = resources$opt$promoter_downstream,
                       genome_fa  = resources$genome_fa)

    cat("  -- genic regions --\n")
    region_grs <- list(
        "5UTR"   = resources$utr5_gr,
        "3UTR"   = resources$utr3_gr,
        "CDS"    = resources$cds_gr,
        "exon"   = resources$exon_gr,
        "intron" = resources$intron_gr
    )
    annotate_genic_regions(motif_dt, motif_gr, region_grs)

    invisible(motif_dt)
}
