#!/usr/bin/env Rscript
# =============================================================================
# 01_build_pair_index.R — Build the master (motif_hit, het_SNP) pair index
# =============================================================================
#
# For ONE motif archetype:
#   1. Read the motif BED.
#   2. Read het-SNPs from the F121-9 VCF that fall within
#      [motif_start - MAX_PAIR_FLANK_BP, motif_end + MAX_PAIR_FLANK_BP].
#   3. Build the cross-product of (motif_hit, SNP) for SNPs in range.
#   4. Compute pair geometry — directionally aware, in transcription-direction
#      terms (NOT just absolute genomic distance).
#
# OUTPUT COLUMNS (per pair):
#   pair_id
#   motif_hit_id
#   motif_id  motif_chrom  motif_start  motif_end  motif_strand  motif_score
#   snp_chrom  snp_pos  snp_ref  snp_alt
#   snp_offset_from_motif_center_bp     SIGNED. Negative = upstream (5'),
#                                       positive = downstream (3'), in
#                                       transcription direction of the motif.
#   snp_offset_from_motif_center_abs_bp Always positive (convenience).
#   snp_distance_from_motif_edge_bp     0 if SNP is in motif body, else
#                                       distance to nearest motif edge.
#   snp_in_motif_body                   Boolean.
#   snp_upstream_of_motif               Boolean (in transcription direction).
#                                       NA if snp_in_motif_body=TRUE.
#   snp_relative_position               "body" | "upstream_flank" |
#                                       "downstream_flank"
#   snp_position_in_motif               1-based offset within motif body
#                                       (in transcription direction; NA if
#                                       SNP not in body).
#
# WHY DIRECTIONAL OFFSETS:
#   For "+" motif at [s, e] with center c, a SNP at p has offset (p - c).
#     p < c -> upstream (5') in transcription direction
#     p > c -> downstream (3') in transcription direction
#   For "-" motif: signs flip. We compute it cleanly here so downstream
#   code never has to worry about strand bookkeeping.
#
# USAGE:
#   Rscript 01_build_pair_index.R \
#     --motif_bed /path/to/AC0001.bed \
#     --vcf /path/to/F121-9_het_snps.ucsc.vcf.gz \
#     --max_flank_bp 200 \
#     --motif_id AC0001 \
#     --outdir /path/to/pair_index/
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(GenomicRanges)
    library(VariantAnnotation)
    library(data.table)
})

option_list <- list(
    make_option("--motif_bed",   type = "character"),
    make_option("--vcf",         type = "character"),
    make_option("--max_flank_bp", type = "integer", default = 200),
    make_option("--motif_id",    type = "character"),
    make_option("--outdir",      type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

cat("============================================================\n")
cat("01_build_pair_index.R — motif:", opt$motif_id, "\n")
cat("============================================================\n\n")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
STANDARD_CHROMS <- paste0("chr", c(1:19, "X", "Y"))

# ---- Read motif BED ----
motif_dt <- fread(opt$motif_bed, header = FALSE,
                  col.names = c("chrom", "start", "end", "motif_name",
                                "score", "strand"))
motif_dt <- motif_dt[chrom %in% STANDARD_CHROMS & strand %in% c("+", "-")]
if (nrow(motif_dt) == 0) {
    cat("[WARN] No motif hits on standard chromosomes. Writing empty output.\n")
    fwrite(data.table(), file.path(opt$outdir,
        paste0(opt$motif_id, "_pair_index.tsv.gz")), sep = "\t",
        compress = "gzip")
    quit(save = "no", status = 0)
}
cat("Motif hits (standard chroms):", nrow(motif_dt), "\n")

# Convert BED 0-based half-open to 1-based closed
motif_dt[, motif_start_1b := start + 1L]
motif_dt[, motif_end_1b := end]
motif_dt[, motif_center_1b := (motif_start_1b + motif_end_1b) %/% 2L]
motif_dt[, motif_hit_id := paste(opt$motif_id, chrom, motif_start_1b,
                                  strand, sep = "__")]

# Build motif GRanges expanded by flank for VCF I/O
motif_gr <- GRanges(
    seqnames = motif_dt$chrom,
    ranges = IRanges(start = motif_dt$motif_start_1b,
                     end = motif_dt$motif_end_1b),
    strand = motif_dt$strand
)
mcols(motif_gr)$motif_idx <- seq_len(nrow(motif_dt))

# Expand for VCF read
expanded <- motif_gr
start(expanded) <- pmax(1L, start(expanded) - opt$max_flank_bp)
end(expanded) <- end(expanded) + opt$max_flank_bp
expanded_reduced <- reduce(expanded, ignore.strand = TRUE)

# ---- Read VCF in this region ----
vcf_chroms <- seqlevels(scanVcfHeader(opt$vcf))
common <- intersect(seqlevels(expanded_reduced), vcf_chroms)
if (length(common) == 0) {
    stop("No common chroms between motifs and VCF. ",
         "Likely chromosome naming mismatch (UCSC vs Ensembl).")
}
expanded_reduced <- keepSeqlevels(expanded_reduced, common,
                                   pruning.mode = "coarse")
cat("Reading VCF for", length(expanded_reduced), "expanded query regions...\n")

svp <- ScanVcfParam(which = expanded_reduced, info = NA, geno = NA)
vcf <- readVcf(opt$vcf, genome = "mm39", param = svp)
cat("  SNPs in range:", nrow(vcf), "\n")

if (nrow(vcf) == 0) {
    cat("[WARN] No SNPs in range. Writing empty pair index.\n")
    fwrite(data.table(), file.path(opt$outdir,
        paste0(opt$motif_id, "_pair_index.tsv.gz")), sep = "\t",
        compress = "gzip")
    quit(save = "no", status = 0)
}

snp_gr <- rowRanges(vcf)
mcols(snp_gr)$ref_allele <- as.character(ref(vcf))
mcols(snp_gr)$alt_allele <- vapply(alt(vcf), function(x)
    paste(as.character(x), collapse = ","), character(1))

# Restrict to biallelic SNPs (single nucleotide REF and ALT)
keep <- nchar(mcols(snp_gr)$ref_allele) == 1L &
        nchar(mcols(snp_gr)$alt_allele) == 1L &
        !grepl(",", mcols(snp_gr)$alt_allele)
n_dropped <- sum(!keep)
snp_gr <- snp_gr[keep]
cat("  Biallelic SNV-only:", length(snp_gr),
    "(dropped", n_dropped, "non-SNV/multi-allelic)\n")

if (length(snp_gr) == 0) {
    fwrite(data.table(), file.path(opt$outdir,
        paste0(opt$motif_id, "_pair_index.tsv.gz")), sep = "\t",
        compress = "gzip")
    quit(save = "no", status = 0)
}

# ---- Build cross-product via overlap of motif (expanded by flank) with SNPs ----
overlap_query <- motif_gr
start(overlap_query) <- pmax(1L, start(overlap_query) - opt$max_flank_bp)
end(overlap_query) <- end(overlap_query) + opt$max_flank_bp

shared <- intersect(seqlevels(overlap_query), seqlevels(snp_gr))
overlap_query <- keepSeqlevels(overlap_query, shared, pruning.mode = "coarse")
snp_gr <- keepSeqlevels(snp_gr, shared, pruning.mode = "coarse")

hits <- findOverlaps(overlap_query, snp_gr, ignore.strand = TRUE)
if (length(hits) == 0) {
    cat("[WARN] No (motif, SNP) overlaps. Writing empty pair index.\n")
    fwrite(data.table(), file.path(opt$outdir,
        paste0(opt$motif_id, "_pair_index.tsv.gz")), sep = "\t",
        compress = "gzip")
    quit(save = "no", status = 0)
}
cat("\n(motif, SNP) pairs to process:", length(hits), "\n")

motif_idx <- queryHits(hits)
snp_idx <- subjectHits(hits)
m <- motif_dt[motif_idx]
snp_pos <- start(snp_gr)[snp_idx]
snp_ref <- mcols(snp_gr)$ref_allele[snp_idx]
snp_alt <- mcols(snp_gr)$alt_allele[snp_idx]
snp_chr <- as.character(seqnames(snp_gr))[snp_idx]

# ---- Compute directional geometry ----
# For "+" motif: direction-of-transcription matches genomic increasing
#   signed_offset = snp_pos - motif_center_1b
# For "-" motif: direction-of-transcription is genomic decreasing
#   signed_offset = motif_center_1b - snp_pos
genomic_offset <- snp_pos - m$motif_center_1b
signed_offset <- ifelse(m$strand == "+", genomic_offset, -genomic_offset)

# Body / flank classification
in_body <- (snp_pos >= m$motif_start_1b) & (snp_pos <= m$motif_end_1b)
distance_from_edge <- pmax(
    0L,
    pmax(m$motif_start_1b - snp_pos, snp_pos - m$motif_end_1b)
)

# Position within motif (in transcription direction; 1-based)
motif_width <- m$motif_end_1b - m$motif_start_1b + 1L
position_in_motif_genomic <- snp_pos - m$motif_start_1b + 1L  # 1-based, "+" oriented
position_in_motif_tx <- ifelse(
    m$strand == "+",
    position_in_motif_genomic,
    motif_width - position_in_motif_genomic + 1L
)
position_in_motif_tx[!in_body] <- NA_integer_

# Upstream / downstream in transcription direction
snp_upstream <- ifelse(in_body, NA, signed_offset < 0L)
relative_position <- ifelse(
    in_body, "body",
    ifelse(snp_upstream, "upstream_flank", "downstream_flank")
)

pair_dt <- data.table(
    pair_id = paste(m$motif_hit_id, snp_chr, snp_pos, sep = "__"),
    motif_hit_id = m$motif_hit_id,
    motif_id = opt$motif_id,
    motif_chrom = m$chrom,
    motif_start = m$motif_start_1b,
    motif_end = m$motif_end_1b,
    motif_strand = m$strand,
    motif_score = m$score,
    motif_width = motif_width,

    snp_chrom = snp_chr,
    snp_pos = snp_pos,
    snp_ref = snp_ref,
    snp_alt = snp_alt,

    snp_offset_from_motif_center_bp = signed_offset,
    snp_offset_from_motif_center_abs_bp = abs(signed_offset),
    snp_distance_from_motif_edge_bp = distance_from_edge,
    snp_in_motif_body = in_body,
    snp_upstream_of_motif = snp_upstream,
    snp_relative_position = relative_position,
    snp_position_in_motif = position_in_motif_tx
)

cat("\nPair geometry summary:\n")
cat("  Total pairs:                ", nrow(pair_dt), "\n")
cat("  Body pairs:                 ", sum(pair_dt$snp_relative_position == "body"), "\n")
cat("  Upstream-flank pairs:       ",
    sum(pair_dt$snp_relative_position == "upstream_flank"), "\n")
cat("  Downstream-flank pairs:     ",
    sum(pair_dt$snp_relative_position == "downstream_flank"), "\n")
cat("  Median |offset|:            ",
    as.integer(median(pair_dt$snp_offset_from_motif_center_abs_bp)), "bp\n")

out_path <- file.path(opt$outdir, paste0(opt$motif_id, "_pair_index.tsv.gz"))
fwrite(pair_dt, out_path, sep = "\t", na = "NA", compress = "gzip")
cat("\nWrote:", out_path, "\n")
cat("Dimensions:", nrow(pair_dt), "rows x", ncol(pair_dt), "cols\n")
