#!/bin/bash
#SBATCH --job-name=postprocess
#SBATCH --partition=short
#SBATCH --time=10:00:00
#SBATCH --mem=250G
#SBATCH --output=logs/postprocess_%j.out
#SBATCH --error=logs/postprocess_%j.err


# ---- Configuration ----
GENOME="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/genome/mm39.fa"
MEME_FILE="/home/alb1273/pausing_phase_project/resources/consensus_pwms.meme"
METADATA_FILE="/home/alb1273/pausing_phase_project/resources/metadata.tsv"
BASE_OUTDIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/moods_scan"
CHROM_SIZE="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/genome/mm39.chrom.sizes"

PVALUE="1e-5"

# Derived paths
PREPDIR="${BASE_OUTDIR}/prep"
CHROMDIR="${BASE_OUTDIR}/per_chrom"
MERGEDIR="${BASE_OUTDIR}/merged"
LOGDIR="${BASE_OUTDIR}/logs"


python 02_postprocess.py \
    --input-dir ${CHROMDIR} \
    --output-dir ${MERGEDIR} \
    --background ${PREPDIR}/background.txt \
    --meme ${MEME_FILE} \
    --pvalue ${PVALUE} \
    --chrom-sizes ${CHROM_SIZE}