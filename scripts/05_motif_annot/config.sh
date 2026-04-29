#!/bin/bash
# =============================================================================
# 05_motif_annot/config.sh
#
# Shared configuration for the 05_motif_annot/ pipeline. Every script in this
# directory sources this file at the top.
# =============================================================================

# ---- Inputs ----
GENOME_FA="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/genome/mm39.fa"
MOTIF_DIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/moods_scan/merged/per_motif"
VCF_DIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/vcf"
GTF_FILE="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/annotation/gencode.vM38.chr_patch_hapl_scaff.annotation.gtf.gz"

# ---- ATAC peaks (now from 03_atac/ pipeline output, F121-9-specific, mm39) ----
# This replaces the 4DN bigBed import + mm10->mm39 liftOver from the old
# 03_motif_annot/. The peaks here are the consensus across all 4 ATAC
# replicates (≥2-of-4 per peak; threshold defined in 03_atac/config.sh).
ATAC_CONSENSUS_BED="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/atac_proc/peaks/consensus/consensus_peaks.bed"

# ---- RNA-seq inputs (from 04c_rnaseq/) ----
# Gene-level expression summary table, one row per gene_id, columns for
# per-rep counts, per-rep TPM, tpm_mean, tpm_sd, etc.
RNASEQ_EXPRESSION_TSV="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/rnaseq/expression/gene_expression_summary.tsv"

# STAR's SJ.out.tab files (one per RNA-seq replicate). The script will glob
# this pattern to find them. Format: "{prefix}_{SAMPLE}_{suffix}".
RNASEQ_SJ_GLOB="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/rnaseq/bam/star/*_SJ.out.tab"

# Minimum unique-read support for a SJ.out.tab junction to be included.
# A junction is kept if any rep has at least this many unique reads
# spanning it. (Multi-mapping reads are ignored for support counting.)
SJ_MIN_UNIQUE_READS=10

# ---- Output root ----
OUTDIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/annot_motifs_v2"
PREPROCESSED_DIR="${OUTDIR}/preprocessed"

# ---- Pipeline parameters ----
HYBRIDS="F121-9,BL6xCAST"
FLANK_DISTANCES="10,25,50,100,200,400"

# Promoter window — KEEP CONSISTENT WITH 02_motif_scanning/config.sh.
# Motifs whose midpoint falls in [TSS-PROMOTER_UPSTREAM, TSS+PROMOTER_DOWNSTREAM]
# (in transcription direction) are considered promoter-overlapping.
PROMOTER_UPSTREAM=1000
PROMOTER_DOWNSTREAM=500

# ---- SLURM defaults: per-motif annotation array ----
DEFAULT_PARTITION="medium"
DEFAULT_TIME="24:00:00"
DEFAULT_MEM="24"
DEFAULT_MAX_CONCURRENT=1000

# ---- SLURM defaults: retry of failed motifs ----
RETRY_PARTITION="medium"
RETRY_TIME="72:00:00"
RETRY_MEM="128G"
RETRY_MAX_CONCURRENT=1000

# ---- SLURM defaults: preprocess step ----
# Preprocessing now also handles RNA-seq SJ.out.tab merging + GTF junction
# extraction, so it needs slightly more time than the old version.
PREPROCESS_PARTITION="short"
PREPROCESS_TIME="02:00:00"
PREPROCESS_MEM="64G"

# ---- Environment ----
CONDA_ENV="test_m"
