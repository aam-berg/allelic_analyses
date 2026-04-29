#!/usr/bin/env Rscript
# =============================================================================
# 02_annotate_motifs.R — Per-motif annotation orchestrator
# =============================================================================
#
# Loads shared resources from preprocessed/, reads ONE motif BED file, calls
# each annotation module in order, and writes a single annotated TSV.
#
# Each module is sourced from a separate file (02a_*.R through 02g_*.R) and
# defines an `annotate_*` function that mutates motif_dt in place via
# data.table's `:=`.
#
# OUTPUT SCHEMA (~136 columns):
#   See README.md for full schema. Briefly:
#     6 BED + 2 (basic) + 12 (genic, strand-aware) + 16 (TSS, strand-aware × directional)
#     + 2 (ATAC) + 48 (SNPs × 2 hybrids × {direct + 7 flanks}) + 10 (expression, strand-aware)
#     + 40 (splice proximity) = ~136 columns
#
# USAGE:
#   Rscript 02_annotate_motifs.R \
#     --motif_bed /path/to/AC0002.bed \
#     --genome_fa /path/to/mm39.fa \
#     --preprocessed_dir /path/to/preprocessed/ \
#     --vcf_dir /path/to/vcf/ \
#     --hybrids "F121-9,BL6xCAST" \
#     --flank_distances "10,25,50,100,200,400" \
#     --promoter_upstream 1000 \
#     --promoter_downstream 500 \
#     --outdir /path/to/output/
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
# CLI
# =============================================================================
option_list <- list(
    make_option("--motif_bed", type = "character",
                help = "Path to motif BED file (one archetype)."),
    make_option("--genome_fa", type = "character",
                help = "Path to mm39 FASTA (.fai must exist)."),
    make_option("--preprocessed_dir", type = "character",
                help = "Output dir from 01_preprocess_resources.R"),
    make_option("--vcf_dir", type = "character",
                help = "Directory containing F1 het SNP VCFs (UCSC-style)."),
    make_option("--hybrids", type = "character", default = "F121-9,BL6xCAST",
                help = "Comma-sep hybrid names [default: %default]"),
    make_option("--flank_distances", type = "character",
                default = "10,25,50,100,200,400",
                help = "Comma-sep flank distances in bp [default: %default]"),
    make_option("--promoter_upstream",   type = "integer", default = 1000,
                help = "Promoter upstream window [default: %default]"),
    make_option("--promoter_downstream", type = "integer", default = 500,
                help = "Promoter downstream window [default: %default]"),
    make_option("--outdir", type = "character",
                help = "Output directory.")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$motif_bed)) stop("--motif_bed required.")

motif_id <- tools::file_path_sans_ext(basename(opt$motif_bed))

cat("================================================================\n")
cat("02_annotate_motifs.R — Motif:", motif_id, "\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

STANDARD_CHROMS <- paste0("chr", c(1:19, "X", "Y"))

# Load all annotation modules (sourced relative to this script's directory)
script_dir <- dirname(sub("--file=", "", commandArgs(trailingOnly = FALSE)[
    grepl("--file=", commandArgs(trailingOnly = FALSE))
])[1])
if (length(script_dir) == 0 || !nzchar(script_dir)) {
    # Running interactively or fallback: assume cwd
    script_dir <- "."
}
source(file.path(script_dir, "02a_motif_basic.R"))
source(file.path(script_dir, "02b_strand_aware_genic.R"))
source(file.path(script_dir, "02c_strand_aware_tss.R"))
source(file.path(script_dir, "02d_atac_overlap.R"))
source(file.path(script_dir, "02e_snp_overlap.R"))
source(file.path(script_dir, "02f_gene_expression.R"))
source(file.path(script_dir, "02g_splice_proximity.R"))

# =============================================================================
# Load shared resources from preprocessed/
# =============================================================================
cat("Loading shared resources from:", opt$preprocessed_dir, "\n\n")

load_rds_or_null <- function(path, label) {
    if (file.exists(path)) {
        cat("  [OK]    ", label, "\n", sep = "")
        readRDS(path)
    } else {
        cat("  [WARN]  ", label, " missing (", path, ")\n", sep = "")
        NULL
    }
}

resources <- list(
    opt              = opt,
    genome_fa        = NULL,
    tx_meta          = NULL,
    gene_meta        = NULL,
    transcripts_gr   = NULL,
    tss_gr           = NULL,
    utr5_gr          = NULL,
    utr3_gr          = NULL,
    cds_gr           = NULL,
    exon_gr          = NULL,
    intron_gr        = NULL,
    atac_gr          = NULL,
    expr_dt          = NULL,
    donors_gr        = NULL,
    acceptors_gr     = NULL
)

# Genome FASTA
stopifnot("Genome FASTA missing" = file.exists(opt$genome_fa))
stopifnot("FASTA index missing"  = file.exists(paste0(opt$genome_fa, ".fai")))
resources$genome_fa <- FaFile(opt$genome_fa)
cat("  [OK]    genome FASTA: ", opt$genome_fa, "\n", sep = "")

# Resource files
resources$tx_meta        <- load_rds_or_null(file.path(opt$preprocessed_dir, "tx_metadata.rds"),         "tx_metadata.rds")
resources$gene_meta      <- load_rds_or_null(file.path(opt$preprocessed_dir, "gene_metadata.rds"),       "gene_metadata.rds")
resources$transcripts_gr <- load_rds_or_null(file.path(opt$preprocessed_dir, "transcripts_gr.rds"),      "transcripts_gr.rds")
resources$tss_gr         <- load_rds_or_null(file.path(opt$preprocessed_dir, "tss_per_transcript_gr.rds"),"tss_per_transcript_gr.rds")
resources$utr5_gr        <- load_rds_or_null(file.path(opt$preprocessed_dir, "genic_5UTR_gr.rds"),       "genic_5UTR_gr.rds")
resources$utr3_gr        <- load_rds_or_null(file.path(opt$preprocessed_dir, "genic_3UTR_gr.rds"),       "genic_3UTR_gr.rds")
resources$cds_gr         <- load_rds_or_null(file.path(opt$preprocessed_dir, "genic_CDS_gr.rds"),        "genic_CDS_gr.rds")
resources$exon_gr        <- load_rds_or_null(file.path(opt$preprocessed_dir, "genic_exon_gr.rds"),       "genic_exon_gr.rds")
resources$intron_gr      <- load_rds_or_null(file.path(opt$preprocessed_dir, "genic_intron_gr.rds"),     "genic_intron_gr.rds")
resources$atac_gr        <- load_rds_or_null(file.path(opt$preprocessed_dir, "atac_consensus_gr.rds"),   "atac_consensus_gr.rds")
resources$expr_dt        <- load_rds_or_null(file.path(opt$preprocessed_dir, "rnaseq_expression.rds"),   "rnaseq_expression.rds")
resources$donors_gr      <- load_rds_or_null(file.path(opt$preprocessed_dir, "splice_donors_gr.rds"),    "splice_donors_gr.rds")
resources$acceptors_gr   <- load_rds_or_null(file.path(opt$preprocessed_dir, "splice_acceptors_gr.rds"), "splice_acceptors_gr.rds")

cat("\n")

# =============================================================================
# Read motif BED
# =============================================================================
cat("================================================================\n")
cat("Reading motif BED\n")
cat("================================================================\n\n")

stopifnot("Motif BED not found" = file.exists(opt$motif_bed))
motif_dt <- fread(opt$motif_bed, header = FALSE,
                  col.names = c("chrom", "start", "end", "motif_id", "score", "strand"),
                  colClasses = c("character", "integer", "integer", "character",
                                 "numeric", "character"))

cat("  Read", nrow(motif_dt), "hits.\n")
motif_dt <- motif_dt[chrom %in% STANDARD_CHROMS]
cat("  After standard chrom filter:", nrow(motif_dt), "\n")
cat("  Strand distribution:\n")
print(table(motif_dt$strand))

if (nrow(motif_dt) == 0) {
    cat("WARNING: No hits on standard chromosomes. Writing empty output.\n")
    # Preserve schema even on empty: motif_hit_id must always be present
    # so that downstream rbinds/joins from 06_allele_pairs/ are well-typed.
    motif_dt[, motif_hit_id := character(0)]
    setcolorder(motif_dt, c("motif_hit_id",
                             setdiff(names(motif_dt), "motif_hit_id")))
    fwrite(motif_dt, file.path(opt$outdir, paste0(motif_id, "_annotated.tsv.gz")),
           sep = "\t", compress = "gzip")
    quit(save = "no", status = 0)
}

motif_dt[, row_id := .I]   # stable identifier for joins; dropped before write

# Build motif GRanges (BED 0-based half-open -> 1-based closed)
motif_gr <- GRanges(
    seqnames = motif_dt$chrom,
    ranges   = IRanges(start = motif_dt$start + 1L, end = motif_dt$end),
    strand   = motif_dt$strand
)
motif_gr$row_id <- motif_dt$row_id

cat("  GRanges built. n=", length(motif_gr),
    " width range=", paste(range(width(motif_gr)), collapse = "-"), "\n\n", sep = "")

# =============================================================================
# Run annotation modules
# =============================================================================
flush_log <- function() if (interactive()) flush.console()
flank_distances <- as.integer(trimws(strsplit(opt$flank_distances, ",")[[1]]))
hybrids <- trimws(strsplit(opt$hybrids, ",")[[1]])

cat("================================================================\n")
cat("[02a] Motif basic (sequence + palindromicity)\n")
cat("================================================================\n")
annotate_basic(motif_dt, motif_gr, resources)
flush_log(); cat("\n")

cat("================================================================\n")
cat("[02b] Strand-aware genic annotations\n")
cat("================================================================\n")
annotate_strand_aware_genic(motif_dt, motif_gr, resources)
flush_log(); cat("\n")

cat("================================================================\n")
cat("[02c] Strand-aware nearest TSS\n")
cat("================================================================\n")
annotate_strand_aware_tss(motif_dt, motif_gr, resources)
flush_log(); cat("\n")

cat("================================================================\n")
cat("[02d] ATAC consensus peak overlap\n")
cat("================================================================\n")
annotate_atac_overlap(motif_dt, motif_gr, resources)
flush_log(); cat("\n")

cat("================================================================\n")
cat("[02e] F1 het SNP overlap (direct + flanking)\n")
cat("================================================================\n")
annotate_snp_overlap(motif_dt, motif_gr, resources, hybrids,
                      flank_dists = flank_distances, vcf_dir = opt$vcf_dir)
flush_log(); cat("\n")

cat("================================================================\n")
cat("[02f] Gene expression annotation\n")
cat("================================================================\n")
annotate_gene_expression(motif_dt, resources)
flush_log(); cat("\n")

cat("================================================================\n")
cat("[02g] Splice junction proximity\n")
cat("================================================================\n")
annotate_splice_proximity(motif_dt, motif_gr, resources)
flush_log(); cat("\n")

# =============================================================================
# Write output
# =============================================================================
cat("================================================================\n")
cat("Writing output\n")
cat("================================================================\n\n")

motif_dt[, row_id := NULL]

# Add motif_hit_id for foreign-key joining with 06_allele_pairs/.
# Colon-stripping handles the case where MEME-style names (e.g.
# "AC0001:DLX/LHX:Homeodomain") propagate into BED filenames; the
# MOODS-scanned consensus_pwms.meme has unique AC000X prefixes, so the
# stripped form is collision-free. If the BED filenames are already
# sanitized (no colons), this sub() is a no-op.
motif_id_short <- sub(":.*$", "", motif_id)
motif_dt[, motif_hit_id := paste(motif_id_short, chrom,
                                  start + 1L, strand, sep = "__")]
setcolorder(motif_dt, c("motif_hit_id",
                         setdiff(names(motif_dt), "motif_hit_id")))

out_path <- file.path(opt$outdir, paste0(motif_id, "_annotated.tsv.gz"))
fwrite(motif_dt, out_path, sep = "\t", na = "NA", compress = "gzip")
cat("  Wrote:", out_path, "\n")
cat("  Dimensions:", nrow(motif_dt), "rows x", ncol(motif_dt), "cols\n\n")

# =============================================================================
# Summary
# =============================================================================
cat("================================================================\n")
cat("Summary for", motif_id, "\n")
cat("================================================================\n\n")

cat(sprintf("  Total hits:                       %d\n", nrow(motif_dt)))
cat(sprintf("  Mean palindromicity:              %.3f\n",
            mean(motif_dt$palindromicity_score)))
cat("\n")
cat(sprintf("  intragenic_sense:                 %d (%.1f%%)\n",
            sum(motif_dt$intragenic_sense),
            100*mean(motif_dt$intragenic_sense)))
cat(sprintf("  intragenic_antisense:             %d (%.1f%%)\n",
            sum(motif_dt$intragenic_antisense),
            100*mean(motif_dt$intragenic_antisense)))
cat(sprintf("  is_promoter_sense:                %d (%.1f%%)\n",
            sum(motif_dt$is_promoter_sense),
            100*mean(motif_dt$is_promoter_sense)))
cat(sprintf("  is_promoter_antisense:            %d (%.1f%%)\n",
            sum(motif_dt$is_promoter_antisense),
            100*mean(motif_dt$is_promoter_antisense)))
cat(sprintf("  atac_consensus_overlap:           %d (%.1f%%)\n",
            sum(motif_dt$atac_consensus_overlap),
            100*mean(motif_dt$atac_consensus_overlap)))

for (h in hybrids) {
    col <- paste0("snp_", h, "_overlap")
    if (col %in% names(motif_dt)) {
        cat(sprintf("  snp_%s_overlap (direct):    %d (%.1f%%)\n",
                    h, sum(motif_dt[[col]]),
                    100*mean(motif_dt[[col]])))
    }
}

# Expression (sense): expressed gene fraction
for (orient in c("sense", "antisense")) {
    col <- paste0("expression_tpm_max_", orient)
    if (col %in% names(motif_dt)) {
        v <- motif_dt[[col]]
        cat(sprintf("  expression_tpm_max_%s expressed (>=1): %d (%.1f%% of motifs in genes)\n",
                    orient, sum(v >= 1, na.rm = TRUE),
                    if (sum(!is.na(v))>0) 100*sum(v >= 1, na.rm=TRUE)/sum(!is.na(v)) else 0))
    }
}

cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")