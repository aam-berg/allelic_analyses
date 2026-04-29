# =============================================================================
# 02e_snp_overlap.R — F1 het SNP overlap (UNSTRANDED): direct + flanking bins
# =============================================================================
#
# OUTPUT COLUMNS (per hybrid):
#   snp_<hybrid>_overlap              boolean (motif body contains a het SNP)
#   snp_<hybrid>_details              pipe-sep "chr:pos:REF>ALT"
#   snp_<hybrid>_count                integer
#
#   For each flank distance D in FLANK_DISTANCES:
#     snp_<hybrid>_flank<D>bp_overlap  boolean
#     snp_<hybrid>_flank<D>bp_details  pipe-sep
#     snp_<hybrid>_flank<D>bp_count    integer
#
# FLANK BIN GEOMETRY (unchanged from old code, well-tested):
#   For motif at [s, e] (1-based closed) and bin (d_inner, d_outer):
#     Left region:  [s - d_outer, s - d_inner - 1]
#     Right region: [e + d_inner + 1, e + d_outer]
#   Bins are concentric, NON-OVERLAPPING annular regions:
#     flank_10bp:  positions 1-10 bp from motif edge
#     flank_25bp:  positions 11-25 bp from motif edge
#     flank_50bp:  positions 26-50 bp from motif edge
#     ... and so on
#
# WHY UNSTRANDED:
#   A SNP is a position; the relevant question for our analysis is whether
#   the genomic position falls within the motif body or its flanks. SNPs
#   don't have intrinsic strand. The downstream analysis in 06_allele_pairs/
#   handles the strand-aware question of "what does the motif sequence look
#   like with each allele".
# =============================================================================

# -----------------------------------------------------------------------------
# Helper: build flanking GRanges for an annular bin (d_inner, d_outer)
# around motif regions. Returns GRanges with $row_id metadata.
# -----------------------------------------------------------------------------
make_flank_gr <- function(base_gr, d_inner, d_outer) {
    s   <- start(base_gr)
    e   <- end(base_gr)
    chr <- seqnames(base_gr)
    rid <- base_gr$row_id

    left_s  <- pmax(1L, s - as.integer(d_outer))
    left_e  <- s - as.integer(d_inner) - 1L
    right_s <- e + as.integer(d_inner) + 1L
    right_e <- e + as.integer(d_outer)

    lv <- (left_s  <= left_e)
    rv <- (right_s <= right_e)

    c(
        GRanges(seqnames = chr[lv],
                ranges = IRanges(start = left_s[lv], end = left_e[lv]),
                strand = "*", row_id = rid[lv]),
        GRanges(seqnames = chr[rv],
                ranges = IRanges(start = right_s[rv], end = right_e[rv]),
                strand = "*", row_id = rid[rv])
    )
}

# -----------------------------------------------------------------------------
# Helper: SNP overlap with aggregation by row_id, writing into motif_dt
# -----------------------------------------------------------------------------
annotate_snp_region <- function(motif_dt, query_gr, snp_gr, col_prefix) {
    col_ov <- paste0(col_prefix, "_overlap")
    col_de <- paste0(col_prefix, "_details")
    col_ct <- paste0(col_prefix, "_count")

    motif_dt[, (col_ov) := FALSE]
    motif_dt[, (col_de) := NA_character_]
    motif_dt[, (col_ct) := 0L]

    if (length(query_gr) == 0 || length(snp_gr) == 0) return(motif_dt)

    olap <- findOverlaps(query_gr, snp_gr, ignore.strand = TRUE)
    if (length(olap) == 0) return(motif_dt)

    snp_dt <- data.table(
        row_id = query_gr$row_id[queryHits(olap)],
        detail = paste0(
            as.character(seqnames(snp_gr[subjectHits(olap)])), ":",
            start(snp_gr[subjectHits(olap)]), ":",
            snp_gr$ref_allele[subjectHits(olap)], ">",
            snp_gr$alt_allele[subjectHits(olap)]
        )
    )
    snp_agg <- snp_dt[, .(
        details = paste(unique(detail), collapse = "|"),
        count   = uniqueN(detail)
    ), by = row_id]

    motif_dt[snp_agg, on = "row_id", (col_ov) := TRUE]
    motif_dt[snp_agg, on = "row_id", (col_de) := i.details]
    motif_dt[snp_agg, on = "row_id", (col_ct) := i.count]

    motif_dt
}

# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------
annotate_snp_overlap <- function(motif_dt, motif_gr, resources, hybrids,
                                  flank_dists, vcf_dir) {

    # Bin geometry: bin_i covers (FLANK[i-1], FLANK[i]] from motif edge
    flank_dists <- sort(unique(as.integer(flank_dists)))
    flank_bins <- data.table(
        d_inner = c(0L, head(flank_dists, -1)),
        d_outer = flank_dists
    )
    max_flank <- max(flank_dists)

    cat("  Flank bins:\n")
    for (i in seq_len(nrow(flank_bins))) {
        cat(sprintf("    flank_%dbp: positions %d-%d bp from motif edge\n",
                    flank_bins$d_outer[i],
                    flank_bins$d_inner[i] + 1, flank_bins$d_outer[i]))
    }

    # Make unstranded motif GR copy (preserves row_id)
    motif_unstranded <- motif_gr
    strand(motif_unstranded) <- "*"

    for (hybrid in hybrids) {
        cat("  -- hybrid:", hybrid, "--\n")

        vcf_ucsc <- file.path(vcf_dir, paste0(hybrid, "_het_snps.ucsc.vcf.gz"))
        vcf_reg  <- file.path(vcf_dir, paste0(hybrid, "_het_snps.vcf.gz"))

        base_prefix <- paste0("snp_", hybrid)

        # Initialize ALL columns for this hybrid (direct + flanks)
        motif_dt[, paste0(base_prefix, "_overlap") := FALSE]
        motif_dt[, paste0(base_prefix, "_details") := NA_character_]
        motif_dt[, paste0(base_prefix, "_count")   := 0L]
        for (d in flank_dists) {
            fp <- paste0(base_prefix, "_flank", d, "bp")
            motif_dt[, paste0(fp, "_overlap") := FALSE]
            motif_dt[, paste0(fp, "_details") := NA_character_]
            motif_dt[, paste0(fp, "_count")   := 0L]
        }

        if (file.exists(vcf_ucsc)) {
            vcf_path <- vcf_ucsc
            cat("    VCF:", basename(vcf_path), "\n")
        } else if (file.exists(vcf_reg)) {
            vcf_path <- vcf_reg
            cat("    [WARN] Using non-UCSC VCF:", basename(vcf_path), "\n")
            cat("           May use Ensembl-style chroms; few overlaps expected.\n")
        } else {
            cat("    [WARN] No VCF found. Skipping hybrid.\n")
            next
        }

        tbi_path <- paste0(vcf_path, ".tbi")
        if (!file.exists(tbi_path)) {
            cat("    [WARN] No tabix index. Run: tabix -p vcf", vcf_path, "\n")
            next
        }

        # Chrom compatibility check
        vcf_header <- scanVcfHeader(vcf_path)
        vcf_chroms <- seqlevels(vcf_header)
        motif_chroms <- unique(as.character(seqnames(motif_gr)))
        common <- intersect(motif_chroms, vcf_chroms)
        if (length(common) == 0) {
            cat("    [WARN] No common chromosomes. Skipping.\n")
            next
        }

        # Filter motif GR to common chroms
        query_gr <- motif_unstranded[seqnames(motif_unstranded) %in% common]
        query_gr <- keepSeqlevels(query_gr, common, pruning.mode = "coarse")

        # Expand by max flank for VCF reading (single I/O pass)
        query_expanded <- query_gr
        start(query_expanded) <- pmax(1L, start(query_expanded) - max_flank)
        end(query_expanded)   <- end(query_expanded) + max_flank
        query_reduced <- reduce(query_expanded, min.gapwidth = 0L)

        cat("    Reading VCF (", length(query_reduced),
            " query regions, expanded ±", max_flank, "bp)...\n", sep = "")

        tryCatch({
            svp <- ScanVcfParam(which = query_reduced, info = NA, geno = NA)
            vcf <- readVcf(vcf_path, genome = "mm39", param = svp)
            cat("    Read", nrow(vcf), "SNP records.\n")

            if (nrow(vcf) > 0) {
                snp_gr <- rowRanges(vcf)
                mcols(snp_gr)$ref_allele <- as.character(ref(vcf))
                mcols(snp_gr)$alt_allele <- vapply(alt(vcf), function(x)
                    paste(as.character(x), collapse = ","), character(1))

                # Direct overlap
                annotate_snp_region(motif_dt, query_gr, snp_gr, base_prefix)
                cat("      direct: ",
                    sum(motif_dt[[paste0(base_prefix, "_overlap")]]),
                    " (",
                    round(100*mean(motif_dt[[paste0(base_prefix, "_overlap")]]), 1),
                    "%)\n", sep = "")

                # Flanking bins
                for (i in seq_len(nrow(flank_bins))) {
                    d_inner <- flank_bins$d_inner[i]
                    d_outer <- flank_bins$d_outer[i]
                    flank_label <- paste0("flank", d_outer, "bp")
                    flank_prefix <- paste0(base_prefix, "_", flank_label)

                    flank_gr <- make_flank_gr(query_gr, d_inner, d_outer)
                    if (length(flank_gr) > 0) {
                        flank_chroms <- intersect(seqlevels(flank_gr), common)
                        if (length(flank_chroms) > 0) {
                            flank_gr <- keepSeqlevels(flank_gr, flank_chroms,
                                                       pruning.mode = "coarse")
                        }
                    }
                    annotate_snp_region(motif_dt, flank_gr, snp_gr, flank_prefix)
                    n_hit <- sum(motif_dt[[paste0(flank_prefix, "_overlap")]])
                    cat(sprintf("      %s: %d (%.1f%%)\n",
                                flank_label, n_hit,
                                100 * n_hit / nrow(motif_dt)))
                }
            }
        }, error = function(e) {
            cat("    [ERROR]:", conditionMessage(e), "\n")
        })
    }

    invisible(motif_dt)
}
