#!/usr/bin/env Rscript
# =============================================================================
# 02_annotate_motifs.R  (v3: strand-aware TSS + flanking SNP bins)
#
# Annotate a single motif archetype BED file with genomic context.
#
# For each motif hit, this script adds:
#   1.  motif_sequence          — Extracted DNA sequence from mm39
#       palindromicity_score    — How palindromic the motif is (0–1)
#   2.  is_intragenic           — Does the motif overlap any transcript?
#   3.  is_promoter             — Is the motif within [-1000,+500] of a TSS?
#       promoter_transcript_ids — Which transcript(s) define that promoter
#       promoter_tsls           — Transcript support level(s)
#       promoter_gene_names     — Gene name(s) for promoter transcripts
#   4.  intragenic_gene_types   — Gene type(s) if intragenic
#   5.  genic_regions           — 5UTR, 3UTR, CDS, exon, intron
#   6.  upstream_tss_dist       — Distance to nearest UPSTREAM TSS (strand-aware)
#       upstream_tss_gene_names, upstream_tss_tx_ids, upstream_tss_tsls
#   7.  downstream_tss_dist     — Distance to nearest DOWNSTREAM TSS (strand-aware)
#       downstream_tss_gene_names, downstream_tss_tx_ids, downstream_tss_tsls
#   8.  atac_<sample>_overlap   — ATAC peak overlap (peaks are now in mm39)
#       atac_<sample>_peak_coords, atac_<sample>_signalValue, etc.
#   9.  snp_<hybrid>_overlap    — F1 het SNP overlap (direct)
#       snp_<hybrid>_details, snp_<hybrid>_count
#  10.  snp_<hybrid>_flank<D>bp_overlap  — SNPs in non-overlapping flanking
#       snp_<hybrid>_flank<D>bp_details    bins around the motif
#       snp_<hybrid>_flank<D>bp_count      (D = 10, 25, 50, 100, 200, 400, 800)
#
# FLANKING SNP BIN LOGIC:
#   Bins are defined as concentric, non-overlapping annular regions around the
#   motif. Each bin excludes the motif itself and all smaller flanks:
#
#     flank_10bp:  positions 1–10 bp from motif edge   (but not overlapping motif)
#     flank_25bp:  positions 11–25 bp from motif edge  (excludes 0–10)
#     flank_50bp:  positions 26–50 bp from motif edge  (excludes 0–25)
#     flank_100bp: positions 51–100 bp from motif edge (excludes 0–50)
#     flank_200bp: positions 101–200 bp                (excludes 0–100)
#     flank_400bp: positions 201–400 bp                (excludes 0–200)
#     flank_800bp: positions 401–800 bp                (excludes 0–400)
#
#   For a motif at [s, e] (1-based closed), the bin (d_inner, d_outer) covers:
#     Left region:  [s - d_outer, s - d_inner - 1]
#     Right region: [e + d_inner + 1, e + d_outer]
#
# KEY CHANGE FROM V1 — STRAND-AWARE TSS DEFINITIONS:
#   "Upstream" and "downstream" are defined relative to the MOTIF'S STRAND,
#   not raw genomic coordinates. This captures the biological question of
#   how TF binding motifs interact with RNAPII transcription direction:
#
#   For a "+" strand motif:
#     5' (upstream)   = lower coordinates (genomically left)
#     3' (downstream) = higher coordinates (genomically right)
#
#   For a "-" strand motif:
#     5' (upstream)   = higher coordinates (genomically right)
#     3' (downstream) = lower coordinates (genomically left)
#
#   "Upstream TSS" = nearest TSS in the 5' direction of the motif
#   "Downstream TSS" = nearest TSS in the 3' direction of the motif
#
#   Distances are always reported as positive values (absolute distance).
#
# OTHER ASSUMPTIONS:
#   - All coordinates are mm39/GRCm39 with UCSC-style chromosome names
#   - Motif BED files have 6 columns: chrom, start, end, motif_id, score, strand
#   - BED is 0-based half-open; converted to 1-based closed for Bioconductor
#   - Genome FASTA is indexed (.fai)
#   - UCSC-style VCFs used for SNP overlap
#   - ATAC peaks already lifted to mm39 during preprocessing
#   - Palindromicity = fraction of positions matching reverse complement
#   - Multi-valued annotations are pipe-separated
#
# USAGE:
#   Rscript 02_annotate_motifs.R \
#     --motif_bed /path/to/AC0002.bed \
#     --genome_fa /path/to/mm39.fa \
#     --preprocessed_dir /path/to/preprocessed/ \
#     --vcf_dir /path/to/vcf/ \
#     --hybrids "F121-9,BL6xCAST" \
#     --outdir /path/to/annot_motifs/
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(rtracklayer)
    library(Rsamtools)
    library(Biostrings)
    library(VariantAnnotation)
    library(data.table)
    library(S4Vectors)
    library(IRanges)
})

# =============================================================================
# 0. Parse arguments
# =============================================================================

option_list <- list(
    make_option("--motif_bed", type = "character",
                help = "Path to motif BED file (one archetype, e.g. AC0002.bed)"),
    make_option("--genome_fa", type = "character",
                default = "/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/genome/mm39.fa",
                help = "Path to mm39 genome FASTA (must be indexed) [default: %default]"),
    make_option("--preprocessed_dir", type = "character",
                default = "/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/annot_motifs/preprocessed",
                help = "Directory with preprocessed resources [default: %default]"),
    make_option("--vcf_dir", type = "character",
                default = "/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/vcf",
                help = "Directory with VCF files [default: %default]"),
    make_option("--hybrids", type = "character",
                default = "F121-9,BL6xCAST",
                help = "Comma-separated hybrid names [default: %default]"),
    make_option("--promoter_upstream", type = "integer", default = 1000,
                help = "Promoter upstream distance from TSS (bp) [default: %default]"),
    make_option("--promoter_downstream", type = "integer", default = 500,
                help = "Promoter downstream distance from TSS (bp) [default: %default]"),
    make_option("--flank_distances", type = "character",
                default = "10,25,50,100,200,400,800",
                help = "Comma-separated flanking distances in bp [default: %default]"),
    make_option("--outdir", type = "character",
                default = "/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/annot_motifs",
                help = "Output directory [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$motif_bed)) stop("ERROR: --motif_bed is required.")

motif_id <- tools::file_path_sans_ext(basename(opt$motif_bed))

# Parse flank distances into sorted integer vector
FLANK_DISTS <- sort(as.integer(trimws(strsplit(opt$flank_distances, ",")[[1]])))
stopifnot("Flank distances must be positive" = all(FLANK_DISTS > 0))
MAX_FLANK <- max(FLANK_DISTS)

cat("================================================================\n")
cat("02_annotate_motifs.R (v3) — Annotating motif:", motif_id, "\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

cat("Parameters:\n")
cat("  Motif BED:        ", opt$motif_bed, "\n")
cat("  Genome FASTA:     ", opt$genome_fa, "\n")
cat("  Preprocessed dir: ", opt$preprocessed_dir, "\n")
cat("  VCF dir:          ", opt$vcf_dir, "\n")
cat("  Hybrids:          ", opt$hybrids, "\n")
cat("  Promoter window:  [-", opt$promoter_upstream, ", +",
    opt$promoter_downstream, "]\n", sep = "")
cat("  Flank distances:  ", paste(FLANK_DISTS, collapse = ", "), " bp\n")
cat("  Max flank:        ", MAX_FLANK, " bp\n")
cat("  Output dir:       ", opt$outdir, "\n\n")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
STANDARD_CHROMS <- paste0("chr", c(1:19, "X", "Y"))

# =============================================================================
# 1. Load shared resources
# =============================================================================

cat("================================================================\n")
cat("Loading shared resources\n")
cat("================================================================\n\n")

cat("[1/6] TxDb...\n")
txdb <- loadDb(file.path(opt$preprocessed_dir, "gencode_vM38_txdb.sqlite"))
cat("  ", length(transcripts(txdb)), "transcripts.\n")

cat("[2/6] Transcript metadata...\n")
tx_meta <- readRDS(file.path(opt$preprocessed_dir, "tx_metadata.rds"))
cat("  ", nrow(tx_meta), "records.\n")

cat("[3/6] TSS GRanges (per-transcript)...\n")
tss_gr <- readRDS(file.path(opt$preprocessed_dir, "tss_granges.rds"))
cat("  ", length(tss_gr), "entries.\n")

cat("[4/6] Unique-position TSS GRanges...\n")
tss_unique_gr <- readRDS(file.path(opt$preprocessed_dir, "tss_unique_granges.rds"))
cat("  ", length(tss_unique_gr), "unique positions.\n")

cat("[5/6] ATAC peaks (mm39, lifted from mm10)...\n")
atac_rds_files <- list.files(opt$preprocessed_dir, pattern = "^atac_peaks_.*\\.rds$",
                              full.names = TRUE)
atac_peaks_list <- list()
for (af in atac_rds_files) {
    sn <- sub("^atac_peaks_", "", sub("\\.rds$", "", basename(af)))
    atac_peaks_list[[sn]] <- readRDS(af)
    cat("  ", sn, ":", length(atac_peaks_list[[sn]]), "peaks (mm39).\n")
}

cat("[6/6] Genome FASTA...\n")
stopifnot("Genome FASTA not found" = file.exists(opt$genome_fa))
stopifnot("FASTA index not found" = file.exists(paste0(opt$genome_fa, ".fai")))
genome_fa <- FaFile(opt$genome_fa)
cat("  Opened:", opt$genome_fa, "\n\n")

# =============================================================================
# 2. Read motif BED file
# =============================================================================

cat("================================================================\n")
cat("Reading motif BED file\n")
cat("================================================================\n\n")

stopifnot("Motif BED not found" = file.exists(opt$motif_bed))
motif_dt <- fread(opt$motif_bed, header = FALSE,
                  col.names = c("chrom", "start", "end", "motif_id", "score", "strand"),
                  colClasses = c("character", "integer", "integer", "character",
                                 "numeric", "character"))

cat("  Read", nrow(motif_dt), "hits for:", motif_id, "\n")
motif_dt <- motif_dt[chrom %in% STANDARD_CHROMS]
cat("  After standard chrom filter:", nrow(motif_dt), "\n")
cat("  Strand distribution:\n")
print(table(motif_dt$strand))

if (nrow(motif_dt) == 0) {
    cat("WARNING: No hits on standard chromosomes. Writing empty output.\n")
    fwrite(motif_dt, file.path(opt$outdir, paste0(motif_id, "_annotated.tsv.gz")),
           sep = "\t", compress = "gzip")
    quit(save = "no", status = 0)
}

motif_dt[, row_id := .I]

# BED 0-based half-open -> 1-based closed for Bioconductor
motif_gr <- GRanges(
    seqnames = motif_dt$chrom,
    ranges   = IRanges(start = motif_dt$start + 1L, end = motif_dt$end),
    strand   = motif_dt$strand
)
motif_gr$row_id <- motif_dt$row_id

# Unstranded version for overlap queries
motif_gr_unstranded <- motif_gr
strand(motif_gr_unstranded) <- "*"

cat("  GRanges: n=", length(motif_gr), " width_range=",
    paste(range(width(motif_gr)), collapse = "-"), "\n\n")

# =============================================================================
# 3. ANNOTATION 1: Motif sequences & palindromicity
# =============================================================================

cat("================================================================\n")
cat("ANNOTATION 1: Motif sequences & palindromicity\n")
cat("================================================================\n\n")

# getSeq respects strand: "-" strand gets reverse-complemented
motif_seqs <- getSeq(genome_fa, motif_gr)
motif_dt[, motif_sequence := as.character(motif_seqs)]

# Palindromicity: fraction of positions matching reverse complement
compute_palindromicity <- function(seq_strings) {
    dna <- DNAStringSet(seq_strings)
    rc  <- reverseComplement(dna)
    vapply(seq_along(dna), function(i) {
        s <- strsplit(as.character(dna[[i]]), "")[[1]]
        r <- strsplit(as.character(rc[[i]]), "")[[1]]
        sum(s == r) / length(s)
    }, numeric(1))
}

motif_dt[, palindromicity_score := compute_palindromicity(motif_sequence)]
cat("  Palindromicity: mean=", round(mean(motif_dt$palindromicity_score), 3),
    " perfect=", sum(motif_dt$palindromicity_score == 1), "/", nrow(motif_dt), "\n\n")

# =============================================================================
# 4. ANNOTATION 2: Intergenic vs. intragenic
# =============================================================================

cat("================================================================\n")
cat("ANNOTATION 2: Intergenic vs. intragenic\n")
cat("================================================================\n\n")

all_tx <- transcripts(txdb, columns = c("tx_name"))
all_tx <- all_tx[seqnames(all_tx) %in% STANDARD_CHROMS]
all_tx <- keepSeqlevels(all_tx, STANDARD_CHROMS, pruning.mode = "coarse")

olap_tx <- findOverlaps(motif_gr_unstranded, all_tx, ignore.strand = TRUE)
motif_dt[, is_intragenic := FALSE]
motif_dt[unique(queryHits(olap_tx)), is_intragenic := TRUE]
cat("  Intragenic:", sum(motif_dt$is_intragenic),
    "(", round(100*mean(motif_dt$is_intragenic), 1), "%)\n\n")

# =============================================================================
# 5. ANNOTATION 3: Promoter overlap [-1000, +500] of TSS
# =============================================================================

cat("================================================================\n")
cat("ANNOTATION 3: Promoter overlap\n")
cat("================================================================\n\n")

tss_gr_std <- tss_gr[seqnames(tss_gr) %in% STANDARD_CHROMS]
tss_gr_std <- keepSeqlevels(tss_gr_std, STANDARD_CHROMS, pruning.mode = "coarse")

promoter_gr <- promoters(tss_gr_std,
                          upstream = opt$promoter_upstream,
                          downstream = opt$promoter_downstream)

si <- seqinfo(genome_fa)
si_std <- si[STANDARD_CHROMS]
seqinfo(promoter_gr) <- si_std
promoter_gr <- trim(promoter_gr)

olap_prom <- findOverlaps(motif_gr_unstranded, promoter_gr, ignore.strand = TRUE)

motif_dt[, `:=`(is_promoter = FALSE,
                promoter_transcript_ids = NA_character_,
                promoter_tsls = NA_character_,
                promoter_gene_names = NA_character_)]

if (length(olap_prom) > 0) {
    prom_dt <- data.table(
        row_id = motif_gr_unstranded$row_id[queryHits(olap_prom)],
        transcript_id = promoter_gr$transcript_id[subjectHits(olap_prom)],
        tsl = promoter_gr$transcript_support_level[subjectHits(olap_prom)],
        gene_name = promoter_gr$gene_name[subjectHits(olap_prom)]
    )
    prom_agg <- prom_dt[, .(
        promoter_transcript_ids = paste(unique(transcript_id), collapse = "|"),
        promoter_tsls = paste(unique(tsl), collapse = "|"),
        promoter_gene_names = paste(unique(gene_name), collapse = "|")
    ), by = row_id]

    motif_dt[prom_agg, on = "row_id", `:=`(
        is_promoter = TRUE,
        promoter_transcript_ids = i.promoter_transcript_ids,
        promoter_tsls = i.promoter_tsls,
        promoter_gene_names = i.promoter_gene_names
    )]
}
cat("  Promoter:", sum(motif_dt$is_promoter),
    "(", round(100*mean(motif_dt$is_promoter), 1), "%)\n\n")

# =============================================================================
# 6. ANNOTATION 4: Intragenic gene type
# =============================================================================

cat("================================================================\n")
cat("ANNOTATION 4: Intragenic gene type\n")
cat("================================================================\n\n")

motif_dt[, intragenic_gene_types := NA_character_]

if (length(olap_tx) > 0) {
    tx_names <- all_tx$tx_name[subjectHits(olap_tx)]
    gene_type_dt <- data.table(row_id = queryHits(olap_tx), tx_name = tx_names)
    gene_type_dt <- merge(gene_type_dt, tx_meta[, .(transcript_id, gene_type)],
                          by.x = "tx_name", by.y = "transcript_id", all.x = TRUE)
    gene_type_agg <- gene_type_dt[, .(
        intragenic_gene_types = paste(unique(na.omit(gene_type)), collapse = "|")
    ), by = row_id]
    gene_type_agg[intragenic_gene_types == "", intragenic_gene_types := NA_character_]
    motif_dt[gene_type_agg, on = "row_id", intragenic_gene_types := i.intragenic_gene_types]
}
cat("  Done.\n\n")

# =============================================================================
# 7. ANNOTATION 5: Genic region type
# =============================================================================

cat("================================================================\n")
cat("ANNOTATION 5: Genic region type (5UTR/3UTR/CDS/exon/intron)\n")
cat("================================================================\n\n")

motif_dt[, genic_regions := NA_character_]

# Extract each region type
get_regions <- function(txdb, type, standard_chroms) {
    gr <- switch(type,
        "5UTR"   = unlist(fiveUTRsByTranscript(txdb)),
        "3UTR"   = unlist(threeUTRsByTranscript(txdb)),
        "CDS"    = cds(txdb),
        "exon"   = exons(txdb),
        "intron" = unlist(intronsByTranscript(txdb))
    )
    gr <- gr[seqnames(gr) %in% standard_chroms]
    if (length(gr) > 0) gr <- keepSeqlevels(gr, standard_chroms, pruning.mode = "coarse")
    return(gr)
}

region_hits <- list()
for (rtype in c("5UTR", "3UTR", "CDS", "exon", "intron")) {
    cat("  ", rtype, "...")
    gr <- get_regions(txdb, rtype, STANDARD_CHROMS)
    if (length(gr) > 0) {
        olap <- findOverlaps(motif_gr_unstranded, gr, ignore.strand = TRUE)
        region_hits[[rtype]] <- unique(queryHits(olap))
        cat(" ", length(region_hits[[rtype]]), "motifs overlap.\n")
    } else {
        region_hits[[rtype]] <- integer(0)
        cat(" no regions.\n")
    }
}

genic_region_dt <- rbindlist(lapply(names(region_hits), function(rtype) {
    if (length(region_hits[[rtype]]) == 0) return(NULL)
    data.table(row_id = region_hits[[rtype]], region = rtype)
}))

if (nrow(genic_region_dt) > 0) {
    genic_agg <- genic_region_dt[, .(
        genic_regions = paste(sort(unique(region)), collapse = "|")
    ), by = row_id]
    motif_dt[genic_agg, on = "row_id", genic_regions := i.genic_regions]
}
cat("  Motifs with genic region:", sum(!is.na(motif_dt$genic_regions)), "\n\n")

# =============================================================================
# 8. ANNOTATION 6 & 7: Nearest UPSTREAM and DOWNSTREAM TSS
#    *** STRAND-AWARE relative to motif orientation ***
# =============================================================================

cat("================================================================\n")
cat("ANNOTATION 6 & 7: Strand-aware nearest upstream/downstream TSS\n")
cat("================================================================\n\n")

cat("  DEFINITION:\n")
cat("    'Upstream TSS'   = nearest TSS in the 5' direction of the motif\n")
cat("    'Downstream TSS' = nearest TSS in the 3' direction of the motif\n")
cat("    For + motif: 5'=left (lower coord), 3'=right (higher coord)\n")
cat("    For - motif: 5'=right (higher coord), 3'=left (lower coord)\n\n")

# Use motif midpoints for distance computation
midpoints <- (start(motif_gr) + end(motif_gr)) %/% 2L
motif_midpt_gr <- GRanges(
    seqnames = seqnames(motif_gr),
    ranges   = IRanges(start = midpoints, width = 1),
    strand   = "*",  # "*" so precede/follow work in pure coordinate space
    row_id   = motif_gr$row_id
)

# Ensure compatible seqlevels
shared_chroms <- intersect(seqlevels(motif_midpt_gr), seqlevels(tss_unique_gr))
motif_midpt_filt <- keepSeqlevels(motif_midpt_gr, shared_chroms, pruning.mode = "coarse")
tss_filt <- keepSeqlevels(tss_unique_gr, shared_chroms, pruning.mode = "coarse")

cat("  Shared chromosomes:", length(shared_chroms), "\n")
cat("  TSS positions:", length(tss_filt), "\n")

# follow(x, subject): nearest subject at LOWER coordinate than x
# precede(x, subject): nearest subject at HIGHER coordinate than x
cat("  Computing follow (nearest TSS at lower coord)...\n")
idx_lower <- follow(motif_midpt_filt, tss_filt)

cat("  Computing precede (nearest TSS at higher coord)...\n")
idx_higher <- precede(motif_midpt_filt, tss_filt)

# Map back from filtered indices to row_ids
filt_row_ids <- motif_midpt_filt$row_id
filt_midpoints <- start(motif_midpt_filt)

# Initialize columns
motif_dt[, `:=`(
    upstream_tss_dist       = NA_integer_,
    upstream_tss_gene_names = NA_character_,
    upstream_tss_tx_ids     = NA_character_,
    upstream_tss_tsls       = NA_character_,
    downstream_tss_dist       = NA_integer_,
    downstream_tss_gene_names = NA_character_,
    downstream_tss_tx_ids     = NA_character_,
    downstream_tss_tsls       = NA_character_
)]

# Helper to fill TSS info using row_id for safe indexing
fill_tss_info <- function(dt, filt_indices, filt_row_ids, filt_midpoints,
                          tss_gr, tss_idx, motif_strand, prefix) {
    # filt_indices: which elements of the filtered GRanges match this strand
    # tss_idx: result of follow/precede for ALL filtered elements
    valid <- filt_indices[!is.na(tss_idx[filt_indices])]
    if (length(valid) == 0) return(dt)

    rids <- filt_row_ids[valid]
    tss_hits <- tss_gr[tss_idx[valid]]
    dists <- abs(filt_midpoints[valid] - start(tss_hits))

    # Use match() to map row_ids back to motif_dt row positions
    dt_rows <- match(rids, dt$row_id)

    col_dist  <- paste0(prefix, "_tss_dist")
    col_genes <- paste0(prefix, "_tss_gene_names")
    col_tx    <- paste0(prefix, "_tss_tx_ids")
    col_tsls  <- paste0(prefix, "_tss_tsls")

    set(dt, i = dt_rows, j = col_dist,  value = as.integer(dists))
    set(dt, i = dt_rows, j = col_genes, value = tss_hits$gene_names)
    set(dt, i = dt_rows, j = col_tx,    value = tss_hits$transcript_ids)
    set(dt, i = dt_rows, j = col_tsls,  value = tss_hits$tsls)

    return(dt)
}

# Indices within motif_midpt_filt for each strand
plus_filt_idx  <- which(filt_row_ids %in% motif_dt[strand == "+", row_id])
minus_filt_idx <- which(filt_row_ids %in% motif_dt[strand == "-", row_id])

cat("  + strand motifs (in filtered set):", length(plus_filt_idx), "\n")
cat("  - strand motifs (in filtered set):", length(minus_filt_idx), "\n")

# For "+" strand motifs:  upstream = lower coord (follow), downstream = higher coord (precede)
motif_dt <- fill_tss_info(motif_dt, plus_filt_idx, filt_row_ids, filt_midpoints,
                          tss_filt, idx_lower, "+", "upstream")
motif_dt <- fill_tss_info(motif_dt, plus_filt_idx, filt_row_ids, filt_midpoints,
                          tss_filt, idx_higher, "+", "downstream")

# For "-" strand motifs:  upstream = higher coord (precede), downstream = lower coord (follow)
motif_dt <- fill_tss_info(motif_dt, minus_filt_idx, filt_row_ids, filt_midpoints,
                          tss_filt, idx_higher, "-", "upstream")
motif_dt <- fill_tss_info(motif_dt, minus_filt_idx, filt_row_ids, filt_midpoints,
                          tss_filt, idx_lower, "-", "downstream")

cat("\n  Upstream TSS distance (strand-aware):\n")
cat("    Min:    ", min(motif_dt$upstream_tss_dist, na.rm = TRUE), "bp\n")
cat("    Median: ", median(motif_dt$upstream_tss_dist, na.rm = TRUE), "bp\n")
cat("    Max:    ", max(motif_dt$upstream_tss_dist, na.rm = TRUE), "bp\n")
cat("  Downstream TSS distance (strand-aware):\n")
cat("    Min:    ", min(motif_dt$downstream_tss_dist, na.rm = TRUE), "bp\n")
cat("    Median: ", median(motif_dt$downstream_tss_dist, na.rm = TRUE), "bp\n")
cat("    Max:    ", max(motif_dt$downstream_tss_dist, na.rm = TRUE), "bp\n\n")

# =============================================================================
# 9. ANNOTATION 8: ATAC peak overlap (now in mm39 after liftOver)
# =============================================================================

cat("================================================================\n")
cat("ANNOTATION 8: ATAC peak overlap (mm39)\n")
cat("================================================================\n\n")

for (sample_name in names(atac_peaks_list)) {
    cat("  Sample:", sample_name, "\n")
    peaks <- atac_peaks_list[[sample_name]]

    shared <- intersect(seqlevels(motif_gr_unstranded), seqlevels(peaks))
    peaks_f <- keepSeqlevels(peaks, shared, pruning.mode = "coarse")
    motif_f <- keepSeqlevels(motif_gr_unstranded, shared, pruning.mode = "coarse")

    olap <- findOverlaps(motif_f, peaks_f, ignore.strand = TRUE)

    col_ov <- paste0("atac_", sample_name, "_overlap")
    col_co <- paste0("atac_", sample_name, "_peak_coords")
    col_si <- paste0("atac_", sample_name, "_signalValue")
    col_pv <- paste0("atac_", sample_name, "_pValue")
    col_qv <- paste0("atac_", sample_name, "_qValue")

    motif_dt[, (col_ov) := FALSE]
    motif_dt[, (col_co) := NA_character_]
    motif_dt[, (col_si) := NA_character_]
    motif_dt[, (col_pv) := NA_character_]
    motif_dt[, (col_qv) := NA_character_]

    if (length(olap) > 0) {
        # Use row_id from the filtered GRanges for safe index mapping
        peak_dt <- data.table(
            row_id = motif_f$row_id[queryHits(olap)],
            coords = paste0(as.character(seqnames(peaks_f[subjectHits(olap)])), ":",
                            start(peaks_f[subjectHits(olap)]), "-",
                            end(peaks_f[subjectHits(olap)])),
            signal = peaks_f$signalValue[subjectHits(olap)],
            pval   = peaks_f$pValue[subjectHits(olap)],
            qval   = peaks_f$qValue[subjectHits(olap)]
        )
        peak_agg <- peak_dt[, .(
            coords = paste(unique(coords), collapse = "|"),
            signal = paste(unique(round(signal, 2)), collapse = "|"),
            pval   = paste(unique(round(pval, 2)), collapse = "|"),
            qval   = paste(unique(round(qval, 2)), collapse = "|")
        ), by = row_id]

        motif_dt[peak_agg, on = "row_id", (col_ov) := TRUE]
        motif_dt[peak_agg, on = "row_id", (col_co) := i.coords]
        motif_dt[peak_agg, on = "row_id", (col_si) := i.signal]
        motif_dt[peak_agg, on = "row_id", (col_pv) := i.pval]
        motif_dt[peak_agg, on = "row_id", (col_qv) := i.qval]
    }
    cat("    Overlap:", sum(motif_dt[[col_ov]]),
        "(", round(100*mean(motif_dt[[col_ov]]), 1), "%)\n")
}
cat("\n")

# =============================================================================
# 10. ANNOTATION 9 & 10: F1 het SNP overlap + flanking bins
# =============================================================================

cat("================================================================\n")
cat("ANNOTATION 9 & 10: F1 het SNP overlap + flanking bins\n")
cat("================================================================\n\n")

cat("  Flank bins (non-overlapping):\n")
flank_bins <- data.table(
    d_inner = c(0L, head(FLANK_DISTS, -1)),
    d_outer = FLANK_DISTS
)
for (i in seq_len(nrow(flank_bins))) {
    cat("    flank_", flank_bins$d_outer[i], "bp: positions ",
        flank_bins$d_inner[i] + 1, "–", flank_bins$d_outer[i],
        " bp from motif edge\n", sep = "")
}
cat("\n")

# ---------------------------------------------------------------------------
# Helper: create flanking GRanges for an annular bin (d_inner, d_outer)
# around a set of base regions. Returns a GRanges with $row_id metadata.
#
# For motif [s, e] (1-based closed):
#   Left region:  [s - d_outer, s - d_inner - 1]
#   Right region: [e + d_inner + 1, e + d_outer]
# ---------------------------------------------------------------------------
make_flank_gr <- function(base_gr, d_inner, d_outer) {
    s   <- start(base_gr)
    e   <- end(base_gr)
    chr <- seqnames(base_gr)
    rid <- base_gr$row_id

    # Left flank
    left_s <- pmax(1L, s - as.integer(d_outer))
    left_e <- s - as.integer(d_inner) - 1L

    # Right flank
    right_s <- e + as.integer(d_inner) + 1L
    right_e <- e + as.integer(d_outer)

    # Keep only valid ranges (start <= end)
    lv <- (left_s <= left_e)
    rv <- (right_s <= right_e)

    flank_gr <- c(
        GRanges(seqnames = chr[lv],
                ranges   = IRanges(start = left_s[lv], end = left_e[lv]),
                strand   = "*",
                row_id   = rid[lv]),
        GRanges(seqnames = chr[rv],
                ranges   = IRanges(start = right_s[rv], end = right_e[rv]),
                strand   = "*",
                row_id   = rid[rv])
    )
    return(flank_gr)
}

# ---------------------------------------------------------------------------
# Helper: find SNP overlaps in a query GRanges, aggregate by row_id, and
# write results into motif_dt columns.
# ---------------------------------------------------------------------------
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

    return(motif_dt)
}

# ---------------------------------------------------------------------------
# Main loop over hybrids
# ---------------------------------------------------------------------------

hybrids <- trimws(strsplit(opt$hybrids, ",")[[1]])

for (hybrid in hybrids) {
    cat("  ── Hybrid:", hybrid, "──\n")

    vcf_ucsc <- file.path(opt$vcf_dir, paste0(hybrid, "_het_snps.ucsc.vcf.gz"))
    vcf_reg  <- file.path(opt$vcf_dir, paste0(hybrid, "_het_snps.vcf.gz"))

    # Initialize ALL columns for this hybrid (direct + flanking)
    base_prefix <- paste0("snp_", hybrid)
    motif_dt[, paste0(base_prefix, "_overlap") := FALSE]
    motif_dt[, paste0(base_prefix, "_details") := NA_character_]
    motif_dt[, paste0(base_prefix, "_count")   := 0L]
    for (d in FLANK_DISTS) {
        fp <- paste0(base_prefix, "_flank", d, "bp")
        motif_dt[, paste0(fp, "_overlap") := FALSE]
        motif_dt[, paste0(fp, "_details") := NA_character_]
        motif_dt[, paste0(fp, "_count")   := 0L]
    }

    # Find VCF file
    if (file.exists(vcf_ucsc)) {
        vcf_path <- vcf_ucsc
        cat("    Using UCSC VCF:", basename(vcf_path), "\n")
    } else if (file.exists(vcf_reg)) {
        vcf_path <- vcf_reg
        cat("    Using regular VCF:", basename(vcf_path), "\n")
        cat("    WARNING: May use Ensembl-style chroms -> no overlaps.\n")
    } else {
        cat("    WARNING: No VCF found. Skipping.\n")
        next
    }

    tbi_path <- paste0(vcf_path, ".tbi")
    if (!file.exists(tbi_path)) {
        cat("    WARNING: No tabix index. Run: tabix -p vcf", vcf_path, "\n")
        next
    }

    # Check chromosome naming compatibility
    vcf_header <- scanVcfHeader(vcf_path)
    vcf_chroms <- seqlevels(vcf_header)
    motif_chroms <- unique(as.character(seqnames(motif_gr)))
    common_chroms <- intersect(motif_chroms, vcf_chroms)
    cat("    Common chromosomes:", length(common_chroms), "\n")

    if (length(common_chroms) == 0) {
        cat("    WARNING: No common chroms! Skipping.\n")
        next
    }

    # Filter motif GRanges to common chromosomes
    query_gr <- motif_gr_unstranded[seqnames(motif_gr_unstranded) %in% common_chroms]
    query_gr <- keepSeqlevels(query_gr, common_chroms, pruning.mode = "coarse")

    # Expand query regions by MAX_FLANK for VCF reading (single I/O pass)
    query_expanded <- query_gr
    start(query_expanded) <- pmax(1L, start(query_expanded) - MAX_FLANK)
    end(query_expanded)   <- end(query_expanded) + MAX_FLANK
    query_reduced <- reduce(query_expanded, min.gapwidth = 0L)

    cat("    Reading VCF (", length(query_reduced), " query regions, ",
        "expanded ±", MAX_FLANK, "bp)...\n", sep = "")
    t0 <- Sys.time()

    tryCatch({
        svp <- ScanVcfParam(which = query_reduced, info = NA, geno = NA)
        vcf <- readVcf(vcf_path, genome = "mm39", param = svp)
        t1 <- Sys.time()
        cat("    Read", nrow(vcf), "SNP records in",
            round(difftime(t1, t0, units = "secs"), 1), "sec.\n")

        if (nrow(vcf) > 0) {
            snp_gr <- rowRanges(vcf)
            mcols(snp_gr)$ref_allele <- as.character(ref(vcf))
            mcols(snp_gr)$alt_allele <- vapply(alt(vcf), function(x)
                paste(as.character(x), collapse = ","), character(1))

            # ---- Direct motif overlap ----
            cat("    Direct overlap...\n")
            motif_dt <- annotate_snp_region(
                motif_dt, query_gr, snp_gr, base_prefix
            )
            cat("      Motifs w/ SNP:", sum(motif_dt[[paste0(base_prefix, "_overlap")]]),
                "(", round(100*mean(motif_dt[[paste0(base_prefix, "_overlap")]]), 1), "%)\n")

            # ---- Flanking bins ----
            for (i in seq_len(nrow(flank_bins))) {
                d_inner <- flank_bins$d_inner[i]
                d_outer <- flank_bins$d_outer[i]
                flank_label <- paste0("flank", d_outer, "bp")
                flank_prefix <- paste0(base_prefix, "_", flank_label)

                cat("    ", flank_label, " (", d_inner + 1, "–", d_outer,
                    "bp from edge)...\n", sep = "")

                flank_gr <- make_flank_gr(query_gr, d_inner, d_outer)

                # Filter flank regions to common chroms (should already be, but safety)
                if (length(flank_gr) > 0) {
                    flank_chroms <- intersect(seqlevels(flank_gr), common_chroms)
                    if (length(flank_chroms) > 0) {
                        flank_gr <- keepSeqlevels(flank_gr, flank_chroms,
                                                  pruning.mode = "coarse")
                    }
                }

                motif_dt <- annotate_snp_region(
                    motif_dt, flank_gr, snp_gr, flank_prefix
                )
                n_hit <- sum(motif_dt[[paste0(flank_prefix, "_overlap")]])
                cat("      Motifs w/ nearby SNP: ", n_hit,
                    " (", round(100 * n_hit / nrow(motif_dt), 1), "%)\n", sep = "")
            }
        }
    }, error = function(e) {
        cat("    ERROR:", conditionMessage(e), "\n")
    })
    cat("\n")
}

# =============================================================================
# 11. Write output
# =============================================================================

cat("================================================================\n")
cat("Writing output\n")
cat("================================================================\n\n")

motif_dt[, row_id := NULL]
out_path <- file.path(opt$outdir, paste0(motif_id, "_annotated.tsv.gz"))
cat("  File:", out_path, "\n")
cat("  Dimensions:", nrow(motif_dt), "x", ncol(motif_dt), "\n")
cat("  Columns:", paste(names(motif_dt), collapse = ", "), "\n\n")
fwrite(motif_dt, out_path, sep = "\t", na = "NA", compress = "gzip")

# =============================================================================
# 12. Summary
# =============================================================================

cat("================================================================\n")
cat("Summary for", motif_id, "\n")
cat("================================================================\n\n")

cat("  Total hits:               ", nrow(motif_dt), "\n")
cat("  Intragenic:               ", sum(motif_dt$is_intragenic),
    " (", round(100*mean(motif_dt$is_intragenic), 1), "%)\n", sep = "")
cat("  In promoter:              ", sum(motif_dt$is_promoter),
    " (", round(100*mean(motif_dt$is_promoter), 1), "%)\n", sep = "")
cat("  Mean palindromicity:      ", round(mean(motif_dt$palindromicity_score), 3), "\n")
cat("  Median upstream TSS dist: ",
    median(motif_dt$upstream_tss_dist, na.rm = TRUE), " bp\n", sep = "")
cat("  Median downstream TSS dist: ",
    median(motif_dt$downstream_tss_dist, na.rm = TRUE), " bp\n", sep = "")

for (sn in names(atac_peaks_list)) {
    col <- paste0("atac_", sn, "_overlap")
    cat("  ATAC (", sn, "): ", sum(motif_dt[[col]]),
        " (", round(100*mean(motif_dt[[col]]), 1), "%)\n", sep = "")
}
for (h in hybrids) {
    col_base <- paste0("snp_", h, "_overlap")
    if (col_base %in% names(motif_dt)) {
        cat("  SNP direct (", h, "): ", sum(motif_dt[[col_base]]),
            " (", round(100*mean(motif_dt[[col_base]]), 1), "%)\n", sep = "")
    }
    for (d in FLANK_DISTS) {
        col_flank <- paste0("snp_", h, "_flank", d, "bp_overlap")
        if (col_flank %in% names(motif_dt)) {
            cat("  SNP flank ", d, "bp (", h, "): ",
                sum(motif_dt[[col_flank]]),
                " (", round(100*mean(motif_dt[[col_flank]]), 1), "%)\n", sep = "")
        }
    }
}

cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")