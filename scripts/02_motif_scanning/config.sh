#!/bin/bash
# =============================================================================
# 02_motif_scanning/config.sh
#
# Shared configuration for the 02_motif_scanning/ pipeline. Every script in
# this directory sources this file at the top.
#
# Edit paths below for your environment. Keep PROMOTER_UPSTREAM and
# PROMOTER_DOWNSTREAM consistent with 03_motif_annot/config.sh.
# =============================================================================

# ---- Inputs ----
GENOME="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/genome/mm39.fa"
CHROM_SIZES="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/genome/mm39.chrom.sizes"
MEME_FILE="/home/alb1273/pausing_phase_project/resources/consensus_pwms.meme"
METADATA_FILE="/home/alb1273/pausing_phase_project/resources/metadata.tsv"
GTF_FILE="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/annotation/gencode.vM38.chr_patch_hapl_scaff.annotation.gtf.gz"

# ---- Output root ----
BASE_OUTDIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/moods_scan"

# ---- Pipeline parameters ----
PVALUE="1e-5"

# Promoter window — KEEP CONSISTENT WITH 03_motif_annot/config.sh
PROMOTER_UPSTREAM=1000
PROMOTER_DOWNSTREAM=500

# ---- Derived paths (rarely edited) ----
PREPDIR="${BASE_OUTDIR}/prep"
SANITYDIR="${BASE_OUTDIR}/sanity_check"
CHROMDIR="${BASE_OUTDIR}/per_chrom"
MERGEDIR="${BASE_OUTDIR}/merged"
QCDIR="${BASE_OUTDIR}/qc_report"
LOGDIR="${BASE_OUTDIR}/logs"
PALINDROME_TABLE="${SANITYDIR}/orientation_analysis.tsv"

# ---- SLURM defaults ----
SCAN_PARTITION="short"
SCAN_TIME="12:00:00"
SCAN_MEM="64G"

POST_PARTITION="short"
POST_TIME="02:00:00"
POST_MEM="64G"

QC_PARTITION="short"
QC_TIME="12:00:00"
QC_MEM="64G"

# ---- Environment ----
CONDA_ENV="motif_scanning"
