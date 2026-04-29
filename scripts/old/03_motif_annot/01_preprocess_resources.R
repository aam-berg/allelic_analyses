#!/usr/bin/env Rscript
# =============================================================================
# 01_preprocess_resources.R  (v2: adds mm10→mm39 liftOver for ATAC peaks)
#
# Pre-process shared genomic resources for the motif annotation pipeline.
#
# This script is run ONCE before the per-motif annotation jobs. It creates:
#   1. TxDb SQLite database from GENCODE GTF
#   2. TSS GRanges object with transcript metadata (RDS)
#   3. Transcript-level metadata table (RDS)
#   4. Gene-level metadata table (RDS)
#   5. ATAC peak GRanges objects (RDS, one per bigBed file) — LIFTED to mm39
#
# KEY CHANGE FROM V1:
#   The ATAC bigBed files from 4DN are in mm10 (GRCm38). We liftOver them to
#   mm39 (GRCm39) using the UCSC mm10ToMm39.over.chain file, since all other
#   data (motifs, GENCODE, genome FASTA, VCFs) are in mm39.
#
# ASSUMPTIONS:
#   - GENCODE M38 GTF is on GRCm39/mm39
#   - ATAC bigBed files from 4DN are on mm10 (GRCm38)
#   - Motif BED files use UCSC-style chromosome names (chr1, chr2, ...)
#   - The liftOver chain file mm10ToMm39.over.chain.gz will be downloaded
#     from UCSC if not already present in outdir
#   - Standard chromosomes only: chr1-chr19, chrX, chrY (chrM excluded)
#
# USAGE:
#   Rscript 01_preprocess_resources.R \
#     --gtf /path/to/gencode.vM38.annotation.gtf.gz \
#     --atac_dir /path/to/atac/bigbed/files/ \
#     --atac_files "4DNFIAEQI3RP.bb,4DNFIZNPOOZN.bb" \
#     --chain_file /path/to/mm10ToMm39.over.chain.gz \
#     --outdir /path/to/output/preprocessed/
#
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
# 0. Parse command-line arguments
# =============================================================================

option_list <- list(
    make_option("--gtf", type = "character",
                default = "/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/annotation/gencode.vM38.chr_patch_hapl_scaff.annotation.gtf.gz",
                help = "Path to GENCODE GTF file [default: %default]"),
    make_option("--atac_dir", type = "character",
                default = "/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/atac/m_musculus",
                help = "Directory containing ATAC bigBed files (mm10) [default: %default]"),
    make_option("--atac_files", type = "character",
                default = "4DNFIAEQI3RP.bb,4DNFIZNPOOZN.bb",
                help = "Comma-separated list of ATAC bigBed filenames [default: %default]"),
    make_option("--chain_file", type = "character",
                default = "",
                help = "Path to mm10ToMm39.over.chain.gz (downloaded automatically if missing)"),
    make_option("--outdir", type = "character",
                default = "/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/annot_motifs/preprocessed",
                help = "Output directory for preprocessed resources [default: %default]"),
    make_option("--standard_chroms_only", type = "logical", default = TRUE,
                help = "Keep only standard chromosomes (chr1-19,X,Y) [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("==============================================================\n")
cat("01_preprocess_resources.R  (v2: with mm10->mm39 liftOver)\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================\n\n")

cat("Parameters:\n")
cat("  GTF:                 ", opt$gtf, "\n")
cat("  ATAC dir:            ", opt$atac_dir, "\n")
cat("  ATAC files:          ", opt$atac_files, "\n")
cat("  Chain file:          ", ifelse(opt$chain_file == "", "(auto-download)", opt$chain_file), "\n")
cat("  Output dir:          ", opt$outdir, "\n")
cat("  Standard chroms only:", opt$standard_chroms_only, "\n\n")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

STANDARD_CHROMS <- paste0("chr", c(1:19, "X", "Y"))
MM10_CHR1_LEN <- 195471971L
MM39_CHR1_LEN <- 195154279L

# =============================================================================
# STEP 0: Obtain liftOver chain file (mm10 -> mm39)
# =============================================================================

cat("==============================================================\n")
cat("STEP 0: Obtaining mm10->mm39 liftOver chain file\n")
cat("==============================================================\n")

chain_path <- opt$chain_file
if (chain_path == "" || !file.exists(chain_path)) {
    chain_path <- file.path(opt$outdir, "mm10ToMm39.over.chain.gz")
}

if (file.exists(chain_path)) {
    cat("  Chain file already exists:", chain_path, "\n")
} else {
    chain_url <- "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz"
    cat("  Downloading chain file from UCSC...\n")
    cat("  URL:", chain_url, "\n")
    tryCatch({
        download.file(chain_url, chain_path, mode = "wb", quiet = FALSE)
        cat("  Downloaded successfully.\n")
    }, error = function(e) {
        stop("Could not download chain file. On a compute node without internet?\n",
             "  Download manually on a login node:\n",
             "    wget -P ", opt$outdir, " ", chain_url, "\n",
             "  Then re-run this script.\n",
             "  Error: ", conditionMessage(e))
    })
}

# Decompress if needed — some versions of rtracklayer::import.chain()
# cannot read .gz files directly
if (grepl("\\.gz$", chain_path)) {
    chain_path_unzipped <- sub("\\.gz$", "", chain_path)
    if (!file.exists(chain_path_unzipped)) {
        cat("  Decompressing chain file (import.chain needs uncompressed)...\n")
        system2("gunzip", args = c("-k", chain_path))
        # -k keeps the original .gz file
        stopifnot("Decompression failed" = file.exists(chain_path_unzipped))
        cat("  Decompressed to:", chain_path_unzipped, "\n")
    } else {
        cat("  Uncompressed chain file already exists:", chain_path_unzipped, "\n")
    }
    chain_path <- chain_path_unzipped
}

cat("  Loading chain file:", chain_path, "\n")
chain <- rtracklayer::import.chain(chain_path)
cat("  Chain file loaded:", length(chain), "chains.\n\n")

# =============================================================================
# STEP 1: Create TxDb from GENCODE GTF
# =============================================================================

txdb_path <- file.path(opt$outdir, "gencode_vM38_txdb.sqlite")

cat("==============================================================\n")
cat("STEP 1: Creating TxDb from GENCODE M38 GTF (mm39/GRCm39)\n")
cat("==============================================================\n")

if (file.exists(txdb_path)) {
    cat("  TxDb already exists:", txdb_path, "\n")
    cat("  Loading existing TxDb...\n")
    txdb <- loadDb(txdb_path)
} else {
    stopifnot("GTF file not found" = file.exists(opt$gtf))
    cat("  Creating TxDb (this may take several minutes)...\n")
    t0 <- Sys.time()
    txdb <- makeTxDbFromGFF(
        file = opt$gtf,
        format = "gtf",
        dataSource = "GENCODE vM38 (GRCm39/mm39)",
        organism = "Mus musculus"
    )
    t1 <- Sys.time()
    cat("  TxDb created in", round(difftime(t1, t0, units = "mins"), 1), "minutes.\n")
    saveDb(txdb, file = txdb_path)
    cat("  TxDb saved to:", txdb_path, "\n")
}

cat("  TxDb summary: genes=", length(genes(txdb)),
    " transcripts=", length(transcripts(txdb)),
    " exons=", length(exons(txdb)), "\n\n")

# =============================================================================
# STEP 2: Extract transcript & gene metadata from GTF
# =============================================================================

cat("==============================================================\n")
cat("STEP 2: Extracting transcript & gene metadata from GTF\n")
cat("==============================================================\n")

cat("  Importing GTF...\n")
t0 <- Sys.time()
gtf_gr <- rtracklayer::import(opt$gtf)
t1 <- Sys.time()
cat("  GTF imported in", round(difftime(t1, t0, units = "mins"), 1),
    "minutes. Total entries:", length(gtf_gr), "\n")

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

# TSS: 5'-most position. + strand: start; - strand: end.
tx_meta[, tss := ifelse(strand == "+", tx_start, tx_end)]

cat("  Gene types (top 5):\n")
print(head(tx_meta[, .N, by = gene_type][order(-N)], 5))

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

saveRDS(tx_meta, file.path(opt$outdir, "tx_metadata.rds"))
saveRDS(gene_meta, file.path(opt$outdir, "gene_metadata.rds"))
cat("  Saved tx_metadata.rds and gene_metadata.rds\n\n")

# =============================================================================
# STEP 3: Create TSS GRanges (per-transcript, for promoter annotation)
# =============================================================================

cat("==============================================================\n")
cat("STEP 3: Creating per-transcript TSS GRanges\n")
cat("==============================================================\n")

if (opt$standard_chroms_only) {
    tx_meta_f <- tx_meta[seqnames %in% STANDARD_CHROMS & !is.na(tss)]
    cat("  Filtered to standard chroms:", nrow(tx_meta_f), "transcripts.\n")
} else {
    tx_meta_f <- tx_meta[!is.na(tss)]
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

saveRDS(tss_gr, file.path(opt$outdir, "tss_granges.rds"))
cat("  Saved", length(tss_gr), "TSS entries to tss_granges.rds\n\n")

# =============================================================================
# STEP 4: Create unique-position TSS GRanges (for nearest-TSS)
# =============================================================================

cat("==============================================================\n")
cat("STEP 4: Creating unique-position TSS GRanges\n")
cat("==============================================================\n")

tss_dt <- data.table(
    seqnames      = as.character(seqnames(tss_gr)),
    tss           = start(tss_gr),
    strand        = as.character(strand(tss_gr)),
    transcript_id = tss_gr$transcript_id,
    gene_id       = tss_gr$gene_id,
    gene_name     = tss_gr$gene_name,
    gene_type     = tss_gr$gene_type,
    tsl           = tss_gr$transcript_support_level
)

# Collapse by unique (seqnames, tss) position
tss_unique_dt <- tss_dt[, .(
    strand_info    = paste(unique(strand), collapse = "|"),
    transcript_ids = paste(unique(transcript_id), collapse = "|"),
    gene_ids       = paste(unique(gene_id), collapse = "|"),
    gene_names     = paste(unique(gene_name), collapse = "|"),
    gene_types     = paste(unique(gene_type), collapse = "|"),
    tsls           = paste(unique(tsl), collapse = "|"),
    n_transcripts  = uniqueN(transcript_id)
), by = .(seqnames, tss)]

cat("  Unique TSS positions:", nrow(tss_unique_dt), "\n")

# Use "*" strand so precede/follow work in coordinate space.
# The strand-aware motif-relative logic is in 02_annotate_motifs.R.
tss_unique_gr <- GRanges(
    seqnames = tss_unique_dt$seqnames,
    ranges   = IRanges(start = tss_unique_dt$tss, width = 1),
    strand   = "*",
    strand_info    = tss_unique_dt$strand_info,
    transcript_ids = tss_unique_dt$transcript_ids,
    gene_ids       = tss_unique_dt$gene_ids,
    gene_names     = tss_unique_dt$gene_names,
    gene_types     = tss_unique_dt$gene_types,
    tsls           = tss_unique_dt$tsls,
    n_transcripts  = tss_unique_dt$n_transcripts
)

saveRDS(tss_unique_gr, file.path(opt$outdir, "tss_unique_granges.rds"))
cat("  Saved to tss_unique_granges.rds\n\n")

# =============================================================================
# STEP 5: Import ATAC peaks + liftOver mm10 -> mm39
# =============================================================================

cat("==============================================================\n")
cat("STEP 5: Importing ATAC peaks + liftOver mm10->mm39\n")
cat("==============================================================\n")

atac_files <- trimws(strsplit(opt$atac_files, ",")[[1]])
cat("  ATAC files to process:", length(atac_files), "\n")

for (af in atac_files) {
    full_path <- file.path(opt$atac_dir, af)
    sample_name <- tools::file_path_sans_ext(af)
    out_path <- file.path(opt$outdir, paste0("atac_peaks_", sample_name, ".rds"))

    cat("\n  --- Processing:", af, "---\n")

    if (!file.exists(full_path)) {
        cat("    WARNING: File not found! Skipping.\n")
        next
    }
    if (file.exists(out_path)) {
        cat("    Already preprocessed. Skipping (delete RDS to redo).\n")
        next
    }

    # Import
    cat("    Importing bigBed...\n")
    t0 <- Sys.time()
    peaks_mm10 <- rtracklayer::import(full_path)
    t1 <- Sys.time()
    cat("    Imported", length(peaks_mm10), "peaks in",
        round(difftime(t1, t0, units = "secs"), 1), "sec.\n")

    # Verify assembly
    si <- seqinfo(peaks_mm10)
    if ("chr1" %in% seqnames(si)) {
        chr1_len <- seqlengths(si)["chr1"]
        cat("    chr1 length:", chr1_len, "\n")
        if (chr1_len == MM10_CHR1_LEN) {
            cat("    CONFIRMED mm10. LiftOver will be applied.\n")
        } else if (chr1_len == MM39_CHR1_LEN) {
            cat("    Already mm39! Saving directly (no liftOver).\n")
            if (opt$standard_chroms_only) {
                peaks_mm10 <- peaks_mm10[seqnames(peaks_mm10) %in% STANDARD_CHROMS]
                peaks_mm10 <- keepSeqlevels(peaks_mm10, STANDARD_CHROMS,
                                             pruning.mode = "coarse")
            }
            saveRDS(peaks_mm10, out_path)
            next
        } else {
            cat("    WARNING: chr1 length doesn't match mm10 or mm39!\n")
            cat("    Proceeding with liftOver, but please verify assembly.\n")
        }
    }

    # Filter to standard chroms before liftOver
    if (opt$standard_chroms_only) {
        n_before <- length(peaks_mm10)
        peaks_mm10 <- peaks_mm10[seqnames(peaks_mm10) %in% STANDARD_CHROMS]
        cat("    Standard chroms:", length(peaks_mm10), "of", n_before, "\n")
    }

    # LiftOver
    cat("    Running liftOver mm10->mm39...\n")
    t0 <- Sys.time()
    lifted_list <- liftOver(peaks_mm10, chain)

    n_input  <- length(peaks_mm10)
    n_mapped <- sum(lengths(lifted_list) >= 1)
    n_split  <- sum(lengths(lifted_list) > 1)
    n_lost   <- sum(lengths(lifted_list) == 0)

    cat("    LiftOver: input=", n_input,
        " mapped=", n_mapped, " (", round(100*n_mapped/n_input, 1), "%)",
        " split=", n_split, " lost=", n_lost, "\n", sep = "")

    # For each original peak, take the first (usually only) lifted fragment.
    # ASSUMPTION: ATAC peaks are typically narrow (<1kb), so splits are rare
    # and taking the first fragment is reasonable. We log split counts above
    # so the user can assess if this is a problem.
    #
    # We use vectorized indexing rather than lapply to avoid the O(n) overhead
    # of constructing ~200k individual GRanges objects.
    keep <- lengths(lifted_list) >= 1
    unlisted <- unlist(lifted_list[keep])
    grp_lens <- lengths(lifted_list[keep])
    first_idx <- cumsum(grp_lens) - grp_lens + 1L
    peaks_mm39 <- unlisted[first_idx]
    cat("    Extracted first fragment per peak:", length(peaks_mm39), "peaks.\n")

    # Carry over metadata from original peaks
    # liftOver preserves mcols, but when we subset with first_hits we need
    # to verify this. Let's check:
    cat("    Metadata columns in lifted peaks:",
        paste(names(mcols(peaks_mm39)), collapse = ", "), "\n")

    t1 <- Sys.time()
    cat("    LiftOver done in", round(difftime(t1, t0, units = "secs"), 1), "sec.\n")

    # Filter to standard chroms in mm39
    if (opt$standard_chroms_only) {
        peaks_mm39 <- peaks_mm39[seqnames(peaks_mm39) %in% STANDARD_CHROMS]
        peaks_mm39 <- keepSeqlevels(peaks_mm39, STANDARD_CHROMS,
                                     pruning.mode = "coarse")
        cat("    After mm39 chrom filter:", length(peaks_mm39), "peaks.\n")
    }

    # Width sanity check
    w <- width(peaks_mm39)
    cat("    Width: min=", min(w), " median=", median(w),
        " mean=", round(mean(w)), " max=", max(w), "\n")

    if ("signalValue" %in% names(mcols(peaks_mm39))) {
        cat("    signalValue range:", round(range(peaks_mm39$signalValue), 2), "\n")
    }

    saveRDS(peaks_mm39, out_path)
    cat("    Saved to:", out_path, "\n")
}

# =============================================================================
# Final summary
# =============================================================================

cat("\n==============================================================\n")
cat("Preprocessing complete!\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================\n\n")

cat("Output files in:", opt$outdir, "\n")
for (f in list.files(opt$outdir, pattern = "\\.(rds|sqlite)$")) {
    sz <- file.info(file.path(opt$outdir, f))$size
    cat("  ", f, " (", round(sz / 1e6, 1), " MB)\n", sep = "")
}
cat("\nATAC peaks have been lifted from mm10 to mm39.\n")
cat("All downstream analyses use mm39 coordinates.\n\n")