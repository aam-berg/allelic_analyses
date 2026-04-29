#!/bin/bash
# =============================================================================
# 06_allele_pairs/config.sh
#
# Central configuration. Sourced by every script in 06_allele_pairs/.
# =============================================================================

# ---- Project / paths ----
PROJECT_NAME="pausing_phase_project"
SCRATCH_ROOT="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates"

# Inputs from upstream pipelines
RESOURCES_DIR="${SCRATCH_ROOT}/resources"
GENOME_FA="${RESOURCES_DIR}/genome/mm39.fa"
VCF_DIR="${RESOURCES_DIR}/vcf"

# 02_motif_scanning/ outputs (motif BEDs + PWMs)
MOTIF_BED_DIR="${SCRATCH_ROOT}/moods_scan/merged/per_motif"
PWM_MEME_FILE="/home/alb1273/pausing_phase_project/resources/consensus_pwms.meme"

# 05_motif_annot/ outputs (annotated motif tables — for joining at analysis time only)
MOTIF_ANNOT_DIR="${SCRATCH_ROOT}/annot_motifs_v2"
MOTIF_PREPROCESSED_DIR="${MOTIF_ANNOT_DIR}/preprocessed"

# 03_atac/ outputs
ATAC_ALLELE_BAM_DIR="${SCRATCH_ROOT}/atac_proc/allele_specific"
ATAC_CONSENSUS_BED="${SCRATCH_ROOT}/atac_proc/peaks/consensus/consensus_peaks.bed"

# 04_proseq/ outputs
PROSEQ_ALLELE_BAM_DIR="${SCRATCH_ROOT}/proseq_proc/allele_specific"
# Single-base bigwigs from PRO-seq (after canonical strand swap), per allele,
# per replicate, per strand. Used for the per-base profile extraction.
PROSEQ_ALLELE_BIGWIG_DIR="${SCRATCH_ROOT}/proseq_proc/bigwig/allele_specific"

# 04c_rnaseq/ outputs
RNASEQ_ALLELE_BAM_DIR="${SCRATCH_ROOT}/rnaseq/allele_specific"
RNASEQ_STAR_DIR="${SCRATCH_ROOT}/rnaseq/bam/star"

# Splice site catalog (built by 05_motif_annot/01_preprocess_resources.R)
SPLICE_DONORS_RDS="${MOTIF_PREPROCESSED_DIR}/splice_donors_gr.rds"
SPLICE_ACCEPTORS_RDS="${MOTIF_PREPROCESSED_DIR}/splice_acceptors_gr.rds"
GENES_GR_RDS="${MOTIF_PREPROCESSED_DIR}/genes_gr.rds"

# Shared lib
LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/lib"

# ---- Output root ----
OUTDIR="${SCRATCH_ROOT}/allele_pairs_v1"

# Pair index (master geometry)
PAIR_INDEX_DIR="${OUTDIR}/pair_index"          # per-motif: {motif_id}_pair_index.tsv.gz
# Per-pair annotation outputs
PWM_DIR="${OUTDIR}/pwm_scores"                  # per-motif: {motif_id}_pwm_scores.tsv.gz
ATAC_COUNTS_DIR="${OUTDIR}/atac_counts"
PROSEQ_COUNTS_DIR="${OUTDIR}/proseq_counts"
PROSEQ_PROFILE_DIR="${OUTDIR}/proseq_profile"   # ±500 bp per-base profiles
RNASEQ_JUNC_DIR="${OUTDIR}/rnaseq_junctions"
PAUSE_DIR="${OUTDIR}/pause_indices"
PSI_DIR="${OUTDIR}/psi"

# Final assembled tables
PAIR_TABLE_DIR="${OUTDIR}/pair_tables"          # per-motif master per-pair tables
LONG_TABLES_DIR="${OUTDIR}/long_tables"         # companion long-format tables

# ---- Pair geometry ----
HYBRID="F121-9"
# A pair is created if a het-SNP for HYBRID falls within the motif body OR
# within ±MAX_PAIR_FLANK_BP of either edge. Beyond this distance, we don't
# expect a SNP to affect motif function (no PWM impact, no direct ATAC
# accessibility relationship at the motif).
MAX_PAIR_FLANK_BP=200

# ---- PWM scoring ----
# Pseudocount added to PWM entries to avoid log(0). Standard MEME default.
PWM_PSEUDOCOUNT=0.01

# ---- PRO-seq pause-window definition ----
# Window is centered on the motif midpoint and has these half-widths.
# Matches 04_proseq/'s motif_window setting.
PAUSE_WINDOW_HALF_BP=10
# "Gene body" denominator for the local pause index: from the motif edge
# extending GENE_BODY_FLANK_BP outward in 3' direction, on the SAME strand.
# Excludes the pause window itself.
GENE_BODY_FLANK_BP=2000

# ---- PRO-seq aggregate profile window (for ARIMA-ready aggregations) ----
# Per-base counts are extracted in [motif_center ± PROFILE_HALF_BP].
PROFILE_HALF_BP=500

# ---- RNA-seq junction PSI ----
# When computing PSI for a splice site, count reads in
# [site - SJ_OVERLAP_BP, site + SJ_OVERLAP_BP] for the "use" denominator
# and reads spanning across the site for the "skip" denominator.
# Read counts come from STAR allele-split BAMs.
SJ_OVERLAP_BP=50

# Min reads in (ref + alt) for PSI to be reported (else NA)
PSI_MIN_TOTAL_READS=4

# ---- Per-rep replicate naming ----
# These are the sample labels used by upstream pipelines. Must match the
# allele_specific/{SAMPLE}_{ref,alt}.bam naming.
ATAC_SAMPLES=("F121-9_atac_rep1" "F121-9_atac_rep2" "F121-9_atac_rep3" "F121-9_atac_rep4")
PROSEQ_SAMPLES=("F121-9_proseq_rep1" "F121-9_proseq_rep2" "F121-9_proseq_rep3" "F121-9_proseq_rep4")
RNASEQ_SAMPLES=("WT_RNAseq_rep1" "WT_RNAseq_rep2")

# ---- SLURM defaults ----
DEFAULT_PARTITION="short"
DEFAULT_TIME="04:00:00"
DEFAULT_MEM="32G"
DEFAULT_THREADS=4
DEFAULT_MAX_CONCURRENT=200

ASSEMBLY_PARTITION="short"
ASSEMBLY_TIME="02:00:00"
ASSEMBLY_MEM="64G"

# ---- Environment ----
CONDA_ENV_R="test_m"           # for R scripts (Bioconductor stack)
CONDA_ENV_TOOLS="proseq2_env"  # for samtools, bedtools, deeptools

# ---- Helpers ----
step_header() {
    echo ""
    echo "============================================================"
    echo "$1"
    echo "Time: $(date)"
    echo "============================================================"
}

create_dirs() {
    mkdir -p "${PAIR_INDEX_DIR}" "${PWM_DIR}"
    mkdir -p "${ATAC_COUNTS_DIR}" "${PROSEQ_COUNTS_DIR}" "${PROSEQ_PROFILE_DIR}"
    mkdir -p "${RNASEQ_JUNC_DIR}" "${PAUSE_DIR}" "${PSI_DIR}"
    mkdir -p "${PAIR_TABLE_DIR}" "${LONG_TABLES_DIR}"
    mkdir -p "${OUTDIR}/logs" "${OUTDIR}/qc"
}
