#!/bin/bash
# =============================================================================
# 03_motif_annot/config.sh
#
# Shared configuration for the 03_motif_annot/ pipeline. Every script in this
# directory sources this file at the top.
#
# Edit paths below for your environment. Keep PROMOTER_UPSTREAM and
# PROMOTER_DOWNSTREAM consistent with 02_motif_scanning/config.sh.
# =============================================================================

# ---- Inputs ----
GENOME_FA="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/genome/mm39.fa"
MOTIF_DIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/moods_scan/merged/per_motif"
VCF_DIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/vcf"
GTF_FILE="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/annotation/gencode.vM38.chr_patch_hapl_scaff.annotation.gtf.gz"
CHAIN_FILE="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/genome/mm10ToMm39.over.chain"

# ---- ATAC-seq inputs (swap when changing datasets) ----
ATAC_DIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/atac/m_musculus"
ATAC_FILES="4DNFIAEQI3RP.bb,4DNFIZNPOOZN.bb"   # comma-separated bigBed filenames

# ---- Output root ----
OUTDIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/annot_motifs"
PREPROCESSED_DIR="${OUTDIR}/preprocessed"

# ---- Pipeline parameters ----
HYBRIDS="F121-9,BL6xCAST"
FLANK_DISTANCES="10,25,50,100,200,400,800"

# Promoter window — KEEP CONSISTENT WITH 02_motif_scanning/config.sh
PROMOTER_UPSTREAM=1000
PROMOTER_DOWNSTREAM=500

# ---- SLURM defaults: normal annotation array ----
DEFAULT_PARTITION="short"
DEFAULT_TIME="01:00:00"
DEFAULT_MEM="20G"
DEFAULT_MAX_CONCURRENT=1000

# ---- SLURM defaults: retry of failed motifs (bumped) ----
RETRY_PARTITION="medium"
RETRY_TIME="72:00:00"
RETRY_MEM="128G"
RETRY_MAX_CONCURRENT=1000

# ---- SLURM defaults: preprocess step ----
PREPROCESS_PARTITION="short"
PREPROCESS_TIME="00:30:00"
PREPROCESS_MEM="64G"

# ---- Environment ----
CONDA_ENV="test_m"
