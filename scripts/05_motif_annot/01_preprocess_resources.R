#!/usr/bin/env Rscript
# =============================================================================
# 01_preprocess_resources.R
#
# Pre-process shared genomic resources for the 05_motif_annot/ pipeline.
# This is run ONCE before the per-motif annotation jobs and produces:
#
#   1. TxDb SQLite database (from GENCODE GTF, mm39)
#   2. Transcript-level metadata (RDS data.table)
#   3. Gene-level metadata (RDS data.table; includes strand)
#   4. TSS GRanges per transcript (RDS, strand preserved — IMPORTANT for the
#      sense/antisense logic in 02c_strand_aware_tss.R)
#   5. Genic feature GRanges with strand: genes, transcripts, 5UTR, 3UTR,
#      CDS, exons, introns (RDS each)
#   6. ATAC consensus peaks (RDS; loaded from 03_atac/ output, mm39, unstranded)
#   7. RNA-seq gene expression table (RDS; loaded from 04c_rnaseq/, indexed
#      by gene_id for fast joins in 02f_gene_expression.R)
#   8. **Splice donor + acceptor GRanges** with `source` metadata — combined
#      from GTF annotation AND STAR SJ.out.tab files across all RNA-seq reps,
#      filtered by the support threshold in config.sh.
#
# THE SPLICE JUNCTION DATA MODEL:
#   We track splice SITES (donor positions and acceptor positions), not full
#   junctions. A donor "site" is a (chrom, position, strand) tuple; same for
#   acceptors. This is appropriate for nearest-site lookups around motifs.
#   Each unique site has an associated:
#     - source: "GTF", "SJ_only", or "GTF_and_SJ"
#     - support_n_unique: max n_unique_reads across reps (NA if GTF-only)
#     - motif_type: GT/AG | GC/AG | AT/AC | non_canonical (pipe-sep if
#                   multiple intron motifs share this site; "" for GTF-only)
#     - gene_names: pipe-sep gene names from GTF transcripts using this
#                   site (empty for SJ_only sites that don't fall in any
#                   GTF transcript)
#
#   Donor positions:
#     "+" intron:  donor = first base of intron (intron_start), strand "+"
#     "-" intron:  donor = last base of intron (intron_end), strand "-"
#   Acceptor positions:
#     "+" intron:  acceptor = last base of intron (intron_end), strand "+"
#     "-" intron:  acceptor = first base of intron (intron_start), strand "-"
#
#   These match STAR's SJ.out.tab conventions: intron_start = leftmost
#   coordinate (1-based), intron_end = rightmost coordinate (1-based).
#
# STAR SJ.out.tab MOTIF CODES:
#   0 = non_canonical
#   1 = GT/AG  (canonical, + strand)
#   2 = CT/AC  (canonical, - strand) — same intron type as 1
#   3 = GC/AG  (minor, + strand)
#   4 = CT/GC  (minor, - strand)
#   5 = AT/AC  (rare, + strand)
#   6 = GT/AT  (rare, - strand)
#
#   We collapse pairs (1,2)->"GT/AG", (3,4)->"GC/AG", (5,6)->"AT/AC",
#   0->"non_canonical".
#
# USAGE:
#   Rscript 01_preprocess_resources.R \
#     --gtf /path/to/gencode.vM38.gtf.gz \
#     --atac_consensus_bed /path/to/03_atac/consensus_peaks.bed \
#     --rnaseq_expression_tsv /path/to/04c_rnaseq/.../gene_expression_summary.tsv \
#     --rnaseq_sj_glob "/path/to/04c_rnaseq/bam/star/*_SJ.out.tab" \
#     --sj_min_unique_reads 10 \
#     --outdir /path/to/preprocessed/
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(rtracklayer)
    library(data.table)
    library(S4Vectors)
    library(IRanges)
})

# =============================================================================
# 0. Parse arguments
# =============================================================================
option_list <- list(
    make_option("--gtf", type = "character",
                help = "Path to GENCODE GTF (mm39)."),
    make_option("--atac_consensus_bed", type = "character",
                help = "Path to consensus ATAC peaks BED from 03_atac/."),
    make_option("--rnaseq_expression_tsv", type = "character",
                help = "Path to gene_expression_summary.tsv from 04c_rnaseq/."),
    make_option("--rnaseq_sj_glob", type = "character",
                help = "Glob pattern for STAR SJ.out.tab files (across reps)."),
    make_option("--sj_min_unique_reads", type = "integer", default = 10,
                help = "Min unique reads for SJ.out.tab support [default: %default]"),
    make_option("--outdir", type = "character",
                help = "Output directory for preprocessed resources."),
    make_option("--standard_chroms_only", type = "logical", default = TRUE,
                help = "Keep only chr1-19,X,Y [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

cat("==============================================================\n")
cat("01_preprocess_resources.R — 05_motif_annot/\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================\n\n")

cat("Parameters:\n")
for (n in names(opt)) cat(sprintf("  %-30s %s\n", paste0(n, ":"), opt[[n]]))
cat("\n")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

STANDARD_CHROMS <- paste0("chr", c(1:19, "X", "Y"))

# =============================================================================
# STEP 1: Build TxDb from GTF (or load if cached)
# =============================================================================
txdb_path <- file.path(opt$outdir, "gencode_txdb.sqlite")
cat("==============================================================\n")
cat("STEP 1: TxDb from GTF\n")
cat("==============================================================\n")

if (file.exists(txdb_path)) {
    cat("  TxDb cached:", txdb_path, "\n")
    txdb <- loadDb(txdb_path)
} else {
    stopifnot("GTF file not found" = file.exists(opt$gtf))
    cat("  Building TxDb (several minutes)...\n")
    t0 <- Sys.time()
    txdb <- makeTxDbFromGFF(
        file = opt$gtf, format = "gtf",
        dataSource = "GENCODE (mm39)",
        organism = "Mus musculus"
    )
    cat("  Built in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
    saveDb(txdb, file = txdb_path)
}

cat("  Genes:        ", length(genes(txdb)), "\n")
cat("  Transcripts:  ", length(transcripts(txdb)), "\n")
cat("  Exons:        ", length(exons(txdb)), "\n\n")

# =============================================================================
# STEP 2: Transcript and gene metadata + TSS GRanges (strand preserved)
# =============================================================================
cat("==============================================================\n")
cat("STEP 2: Transcript / gene metadata + TSS GRanges\n")
cat("==============================================================\n")

cat("  Importing GTF (for gene_name, gene_type, TSL)...\n")
t0 <- Sys.time()
gtf_gr <- rtracklayer::import(opt$gtf)
cat("  Imported", length(gtf_gr), "entries in",
    round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

# --- Transcript metadata ---
tx_entries <- gtf_gr[gtf_gr$type == "transcript"]
cat("  Transcript entries:", length(tx_entries), "\n")
tx_meta <- data.table(
    transcript_id            = tx_entries$transcript_id,
    gene_id                  = tx_entries$gene_id,
    gene_name                = tx_entries$gene_name,
    gene_type                = tx_entries$gene_type,
    transcript_type          = tx_entries$transcript_type,
    transcript_name          = tx_entries$transcript_name,
    transcript_support_level = tx_entries$transcript_support_level,
    seqnames                 = as.character(seqnames(tx_entries)),
    strand                   = as.character(strand(tx_entries)),
    tx_start                 = start(tx_entries),
    tx_end                   = end(tx_entries)
)
tx_meta[, tss := ifelse(strand == "+", tx_start, tx_end)]

# --- Gene metadata ---
gene_entries <- gtf_gr[gtf_gr$type == "gene"]
gene_meta <- data.table(
    gene_id    = gene_entries$gene_id,
    gene_name  = gene_entries$gene_name,
    gene_type  = gene_entries$gene_type,
    seqnames   = as.character(seqnames(gene_entries)),
    strand     = as.character(strand(gene_entries)),
    gene_start = start(gene_entries),
    gene_end   = end(gene_entries)
)
cat("  Gene entries:", nrow(gene_meta), "\n")

saveRDS(tx_meta,   file.path(opt$outdir, "tx_metadata.rds"))
saveRDS(gene_meta, file.path(opt$outdir, "gene_metadata.rds"))

# --- TSS GRanges (strand PRESERVED — important for sense/antisense queries) ---
if (opt$standard_chroms_only) {
    tx_meta_f <- tx_meta[seqnames %in% STANDARD_CHROMS & !is.na(tss) & strand %in% c("+", "-")]
    cat("  Transcripts on standard chroms with strand:", nrow(tx_meta_f), "\n")
} else {
    tx_meta_f <- tx_meta[!is.na(tss) & strand %in% c("+", "-")]
}

tss_gr <- GRanges(
    seqnames = tx_meta_f$seqnames,
    ranges   = IRanges(start = tx_meta_f$tss, width = 1),
    strand   = tx_meta_f$strand,
    transcript_id            = tx_meta_f$transcript_id,
    gene_id                  = tx_meta_f$gene_id,
    gene_name                = tx_meta_f$gene_name,
    gene_type                = tx_meta_f$gene_type,
    transcript_support_level = tx_meta_f$transcript_support_level
)
saveRDS(tss_gr, file.path(opt$outdir, "tss_per_transcript_gr.rds"))
cat("  Saved", length(tss_gr), "TSS entries (strand-preserved).\n\n")

# =============================================================================
# STEP 3: Genic feature GRanges (strand-preserved): genes, transcripts,
#         5UTR, 3UTR, CDS, exons, introns
# =============================================================================
cat("==============================================================\n")
cat("STEP 3: Genic feature GRanges\n")
cat("==============================================================\n")

filter_to_std <- function(gr) {
    if (length(gr) == 0) return(gr)
    if (opt$standard_chroms_only) {
        gr <- gr[seqnames(gr) %in% STANDARD_CHROMS]
        if (length(gr) > 0)
            gr <- keepSeqlevels(gr, STANDARD_CHROMS, pruning.mode = "coarse")
    }
    return(gr)
}

# Genes
cat("  Building genes GRanges...\n")
genes_gr <- genes(txdb)
genes_gr$gene_id <- names(genes_gr)
genes_gr <- filter_to_std(genes_gr)
saveRDS(genes_gr, file.path(opt$outdir, "genes_gr.rds"))
cat("    genes:", length(genes_gr), "(strand:",
    paste(table(strand(genes_gr)), collapse=" / "), ")\n")

# Transcripts
cat("  Building transcripts GRanges...\n")
transcripts_gr <- transcripts(txdb, columns = c("tx_name", "gene_id"))
transcripts_gr$tx_id <- transcripts_gr$tx_name
# gene_id from transcripts() can be a CharacterList; flatten
gid <- as.character(transcripts_gr$gene_id)
gid[gid == "character(0)"] <- ""
transcripts_gr$gene_id <- gid
transcripts_gr <- filter_to_std(transcripts_gr)
saveRDS(transcripts_gr, file.path(opt$outdir, "transcripts_gr.rds"))
cat("    transcripts:", length(transcripts_gr), "\n")

# 5UTR, 3UTR, CDS, exons, introns
for (rtype in c("5UTR", "3UTR", "CDS", "exon", "intron")) {
    cat("  Building", rtype, "GRanges...\n")
    gr <- switch(rtype,
        "5UTR"   = unlist(fiveUTRsByTranscript(txdb)),
        "3UTR"   = unlist(threeUTRsByTranscript(txdb)),
        "CDS"    = cds(txdb),
        "exon"   = exons(txdb),
        "intron" = unlist(intronsByTranscript(txdb))
    )
    gr <- filter_to_std(gr)
    saveRDS(gr, file.path(opt$outdir, paste0("genic_", rtype, "_gr.rds")))
    cat("    ", rtype, ":", length(gr), "\n", sep = "")
}
cat("\n")

# =============================================================================
# STEP 4: ATAC consensus peaks (mm39, F121-9-specific, from 03_atac/)
# =============================================================================
cat("==============================================================\n")
cat("STEP 4: ATAC consensus peaks\n")
cat("==============================================================\n")

if (!file.exists(opt$atac_consensus_bed)) {
    cat("  WARNING: ATAC consensus BED not found:", opt$atac_consensus_bed, "\n")
    cat("           Was 03_atac/ run? Skipping. Annotation step will produce\n")
    cat("           atac_consensus_overlap=NA columns.\n\n")
} else {
    cat("  Reading:", opt$atac_consensus_bed, "\n")
    # 03_atac consensus_peaks.bed columns:
    #   chrom, start, end, contributing_samples (comma-sep), n_replicates
    atac_dt <- fread(opt$atac_consensus_bed, header = FALSE,
                     col.names = c("chrom", "start", "end",
                                   "contributing_samples", "n_replicates"),
                     colClasses = c("character", "integer", "integer",
                                    "character", "integer"))

    atac_dt <- atac_dt[chrom %in% STANDARD_CHROMS]
    cat("  After standard chrom filter:", nrow(atac_dt), "peaks\n")

    atac_gr <- GRanges(
        seqnames = atac_dt$chrom,
        ranges = IRanges(start = atac_dt$start + 1L, end = atac_dt$end),  # BED 0-based -> 1-based
        strand = "*",
        contributing_samples = atac_dt$contributing_samples,
        n_replicates = atac_dt$n_replicates
    )
    saveRDS(atac_gr, file.path(opt$outdir, "atac_consensus_gr.rds"))
    cat("  Saved", length(atac_gr), "consensus peaks.\n")
    cat("  Replicate-support distribution:\n")
    print(table(atac_gr$n_replicates))
    cat("\n")
}

# =============================================================================
# STEP 5: RNA-seq gene expression table (from 04c_rnaseq/)
# =============================================================================
cat("==============================================================\n")
cat("STEP 5: RNA-seq gene expression table\n")
cat("==============================================================\n")

if (!file.exists(opt$rnaseq_expression_tsv)) {
    cat("  WARNING: Expression TSV not found:", opt$rnaseq_expression_tsv, "\n")
    cat("           Was 04c_rnaseq/ run? Skipping. Annotation step will\n")
    cat("           produce expression_*=NA columns.\n\n")
} else {
    cat("  Reading:", opt$rnaseq_expression_tsv, "\n")
    expr_dt <- fread(opt$rnaseq_expression_tsv, header = TRUE)
    setkey(expr_dt, gene_id)
    saveRDS(expr_dt, file.path(opt$outdir, "rnaseq_expression.rds"))
    cat("  Saved", nrow(expr_dt), "genes.\n")
    cat("  Expressed (TPM>=1):", sum(expr_dt$tpm_mean >= 1, na.rm = TRUE), "\n\n")
}

# =============================================================================
# STEP 6: Splice donor and acceptor GRanges
#         GTF + STAR SJ.out.tab (across all reps), filtered by support
# =============================================================================
cat("==============================================================\n")
cat("STEP 6: Splice donor + acceptor GRanges (GTF + STAR SJ.out.tab)\n")
cat("==============================================================\n")

# --- 6a. Extract GTF-derived donors and acceptors ---
cat("  6a. Extracting GTF-annotated splice sites from txdb introns...\n")
# Get introns as a GRangesList by transcript, then unlist with names tracking
introns_by_tx <- intronsByTranscript(txdb, use.names = TRUE)
introns_flat <- unlist(introns_by_tx)
introns_flat$transcript_id <- names(introns_flat)
introns_flat <- filter_to_std(introns_flat)

# Map transcript_id -> gene_id, gene_name (from tx_meta)
tx2gene <- tx_meta[, .(transcript_id, gene_id, gene_name)]
setkey(tx2gene, transcript_id)

intron_dt <- data.table(
    chrom = as.character(seqnames(introns_flat)),
    intron_start = start(introns_flat),
    intron_end = end(introns_flat),
    strand = as.character(strand(introns_flat)),
    transcript_id = introns_flat$transcript_id
)
intron_dt <- merge(intron_dt, tx2gene, by = "transcript_id", all.x = TRUE)
cat("    GTF introns (after merge with tx metadata):", nrow(intron_dt), "\n")

# Donor / acceptor positions (depend on strand):
#   "+" intron: donor = intron_start, acceptor = intron_end
#   "-" intron: donor = intron_end,   acceptor = intron_start
intron_dt[, donor_pos    := ifelse(strand == "+", intron_start, intron_end)]
intron_dt[, acceptor_pos := ifelse(strand == "+", intron_end,   intron_start)]

# Collapse by (chrom, position, strand) for donors and acceptors
gtf_donors_dt <- intron_dt[
    strand %in% c("+", "-"),
    .(gtf_supported     = TRUE,
      gene_names_gtf    = paste(unique(na.omit(gene_name[gene_name != ""])), collapse = "|"),
      transcript_ids_gtf = paste(unique(transcript_id), collapse = "|")),
    by = .(chrom = chrom, position = donor_pos, strand)
]
gtf_acceptors_dt <- intron_dt[
    strand %in% c("+", "-"),
    .(gtf_supported     = TRUE,
      gene_names_gtf    = paste(unique(na.omit(gene_name[gene_name != ""])), collapse = "|"),
      transcript_ids_gtf = paste(unique(transcript_id), collapse = "|")),
    by = .(chrom = chrom, position = acceptor_pos, strand)
]
cat("    Unique GTF donor positions:    ", nrow(gtf_donors_dt), "\n")
cat("    Unique GTF acceptor positions: ", nrow(gtf_acceptors_dt), "\n")

# --- 6b. Read STAR SJ.out.tab files across reps ---
cat("\n  6b. Reading STAR SJ.out.tab files...\n")

sj_files <- Sys.glob(opt$rnaseq_sj_glob)
if (length(sj_files) == 0) {
    cat("  WARNING: No SJ.out.tab files matching:", opt$rnaseq_sj_glob, "\n")
    cat("           Splice junction sources will be GTF-only. Continuing.\n")
    sj_collapsed <- data.table()
} else {
    cat("    Found", length(sj_files), "SJ files:\n")
    for (f in sj_files) cat("      ", f, "\n")

    sj_per_rep <- lapply(sj_files, function(p) {
        # SJ.out.tab columns: chrom, intron_start, intron_end, strand_code,
        # motif_code, annotated, n_unique, n_multi, max_overhang
        fread(p, header = FALSE,
              col.names = c("chrom", "intron_start", "intron_end",
                            "strand_code", "motif_code", "annotated",
                            "n_unique", "n_multi", "max_overhang"))
    })
    sj_all <- rbindlist(sj_per_rep, idcol = "rep_idx")

    # Standard chroms + valid strand only
    sj_all <- sj_all[chrom %in% STANDARD_CHROMS & strand_code %in% c(1, 2)]

    # Collapse across reps: take MAX of n_unique (junction supported if any
    # rep has it strongly), MAX of annotated flag, and consolidate motif_code
    sj_collapsed <- sj_all[, .(
        n_unique_max  = max(n_unique),
        n_multi_max   = max(n_multi),
        annotated_any = max(annotated),
        motif_code    = first(motif_code)   # motif is intrinsic to the intron
    ), by = .(chrom, intron_start, intron_end, strand_code)]

    # Filter by support
    cat("    Total junctions across reps:", nrow(sj_collapsed), "\n")
    sj_collapsed <- sj_collapsed[n_unique_max >= opt$sj_min_unique_reads]
    cat("    After n_unique >=", opt$sj_min_unique_reads, ":",
        nrow(sj_collapsed), "junctions\n")

    # Decode strand and motif
    sj_collapsed[, strand := ifelse(strand_code == 1, "+", "-")]
    sj_collapsed[, motif_type := fcase(
        motif_code == 0L, "non_canonical",
        motif_code %in% c(1L, 2L), "GT/AG",
        motif_code %in% c(3L, 4L), "GC/AG",
        motif_code %in% c(5L, 6L), "AT/AC",
        default = "unknown"
    )]

    # Donor / acceptor positions
    sj_collapsed[, donor_pos    := ifelse(strand == "+", intron_start, intron_end)]
    sj_collapsed[, acceptor_pos := ifelse(strand == "+", intron_end,   intron_start)]

    cat("    SJ motif type counts:\n")
    print(sj_collapsed[, .N, by = motif_type])
}

# --- 6c. Aggregate SJ donor and acceptor sites ---
if (nrow(sj_collapsed) > 0) {
    cat("\n  6c. Aggregating SJ donor / acceptor sites across junctions...\n")
    sj_donors_dt <- sj_collapsed[, .(
        sj_supported          = TRUE,
        sj_n_unique_max       = max(n_unique_max),
        sj_motif_types        = paste(sort(unique(motif_type)), collapse = "|"),
        sj_annotated_any      = max(annotated_any)
    ), by = .(chrom = chrom, position = donor_pos, strand)]

    sj_acceptors_dt <- sj_collapsed[, .(
        sj_supported          = TRUE,
        sj_n_unique_max       = max(n_unique_max),
        sj_motif_types        = paste(sort(unique(motif_type)), collapse = "|"),
        sj_annotated_any      = max(annotated_any)
    ), by = .(chrom = chrom, position = acceptor_pos, strand)]

    cat("    SJ donor positions (unique):    ", nrow(sj_donors_dt), "\n")
    cat("    SJ acceptor positions (unique): ", nrow(sj_acceptors_dt), "\n")
} else {
    sj_donors_dt    <- data.table(chrom = character(), position = integer(),
                                  strand = character(), sj_supported = logical(),
                                  sj_n_unique_max = integer(),
                                  sj_motif_types = character(),
                                  sj_annotated_any = integer())
    sj_acceptors_dt <- data.table(chrom = character(), position = integer(),
                                  strand = character(), sj_supported = logical(),
                                  sj_n_unique_max = integer(),
                                  sj_motif_types = character(),
                                  sj_annotated_any = integer())
}

# --- 6d. Outer-merge GTF and SJ site tables; build final GRanges ---
cat("\n  6d. Outer-merging GTF + SJ sites; building GRanges...\n")

build_site_gr <- function(gtf_dt, sj_dt, label) {
    # Outer merge
    merged <- merge(gtf_dt, sj_dt, by = c("chrom", "position", "strand"), all = TRUE)

    # Fill NAs
    merged[is.na(gtf_supported),     gtf_supported     := FALSE]
    merged[is.na(sj_supported),      sj_supported      := FALSE]
    merged[is.na(gene_names_gtf),    gene_names_gtf    := ""]
    merged[is.na(transcript_ids_gtf), transcript_ids_gtf := ""]
    merged[is.na(sj_motif_types),    sj_motif_types    := ""]

    # Source label
    merged[, source := fcase(
        gtf_supported &  sj_supported, "GTF_and_SJ",
        gtf_supported & !sj_supported, "GTF",
        !gtf_supported &  sj_supported, "SJ_only",
        default = NA_character_
    )]
    merged <- merged[!is.na(source)]

    cat("    [", label, "]\n", sep = "")
    cat("      Total sites:        ", nrow(merged), "\n")
    cat("      GTF only:           ", sum(merged$source == "GTF"), "\n")
    cat("      SJ only:            ", sum(merged$source == "SJ_only"), "\n")
    cat("      GTF and SJ:         ", sum(merged$source == "GTF_and_SJ"), "\n")

    # Build GRanges (strand-preserved!)
    gr <- GRanges(
        seqnames = merged$chrom,
        ranges = IRanges(start = merged$position, width = 1),
        strand = merged$strand,
        source = merged$source,
        support_n_unique = ifelse(merged$sj_supported, merged$sj_n_unique_max, NA_integer_),
        motif_type = merged$sj_motif_types,
        gene_names = merged$gene_names_gtf,
        transcript_ids = merged$transcript_ids_gtf
    )
    if (opt$standard_chroms_only && length(gr) > 0) {
        gr <- keepSeqlevels(gr, intersect(seqlevels(gr), STANDARD_CHROMS),
                            pruning.mode = "coarse")
    }
    return(gr)
}

donors_gr    <- build_site_gr(gtf_donors_dt,    sj_donors_dt,    "DONORS")
acceptors_gr <- build_site_gr(gtf_acceptors_dt, sj_acceptors_dt, "ACCEPTORS")

saveRDS(donors_gr,    file.path(opt$outdir, "splice_donors_gr.rds"))
saveRDS(acceptors_gr, file.path(opt$outdir, "splice_acceptors_gr.rds"))
cat("\n  Saved splice_donors_gr.rds and splice_acceptors_gr.rds\n\n")

# =============================================================================
# Final summary
# =============================================================================
cat("==============================================================\n")
cat("Preprocessing complete!\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================\n\n")

cat("Output files in:", opt$outdir, "\n")
for (f in sort(list.files(opt$outdir, pattern = "\\.(rds|sqlite)$"))) {
    sz <- file.info(file.path(opt$outdir, f))$size
    cat(sprintf("  %-45s %8.1f MB\n", f, sz / 1e6))
}
cat("\n")
