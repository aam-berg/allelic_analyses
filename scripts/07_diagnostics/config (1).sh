#!/bin/bash
# =============================================================================
# 07_diagnostics/config.sh — Configuration for diagnostics pipeline
# =============================================================================

# ---- Paths ----
PROJECT_NAME="pausing_phase_project"
SCRATCH_ROOT="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates"

# Inputs from upstream pipelines
RESOURCES_DIR="${SCRATCH_ROOT}/resources"
GENOME_FA="${RESOURCES_DIR}/genome/mm39.fa"
VCF_FILE="${RESOURCES_DIR}/vcf/F121-9_het_snps.ucsc.vcf.gz"

# 02_motif_scanning/ outputs (raw motif BEDs — for total counts before any filter)
MOTIF_BED_DIR="${SCRATCH_ROOT}/moods_scan/merged/per_motif"

# 05_motif_annot/ per-motif annotated tables
MOTIF_ANNOT_DIR="${SCRATCH_ROOT}/annot_motifs_v2"

# 06_allele_pairs/ outputs (per-pair tables and long-form companions)
PAIR_TABLE_DIR="${SCRATCH_ROOT}/allele_pairs_v1/pair_tables"
PAIR_INDEX_DIR="${SCRATCH_ROOT}/allele_pairs_v1/pair_index"
ATAC_COUNTS_DIR="${SCRATCH_ROOT}/allele_pairs_v1/atac_counts"
PROSEQ_COUNTS_DIR="${SCRATCH_ROOT}/allele_pairs_v1/proseq_counts"

# Allele-specific BAMs from 03_atac/, 04_proseq/, 04c_rnaseq/
ATAC_ALLELE_BAM_DIR="${SCRATCH_ROOT}/atac_proc/allele_specific"
PROSEQ_ALLELE_BAM_DIR="${SCRATCH_ROOT}/proseq_proc/allele_specific"
RNASEQ_ALLELE_BAM_DIR="${SCRATCH_ROOT}/rnaseq/allele_specific"

# Pre-WASP BAMs (for WASP cascade diagnostics — read totals at each filter step)
ATAC_BAM_DIR="${SCRATCH_ROOT}/atac_proc/bam/wasp"          # post-alignment, pre-allele-split
PROSEQ_BAM_DIR="${SCRATCH_ROOT}/proseq_proc/bam/wasp"
RNASEQ_BAM_DIR="${SCRATCH_ROOT}/rnaseq/bam/star"

# Reference/metadata
PWM_MEME_FILE="/home/alb1273/pausing_phase_project/resources/consensus_pwms.meme"
TF_METADATA_FILE="/home/alb1273/pausing_phase_project/resources/metadata.tsv"
GENE_EXPRESSION_FILE="${SCRATCH_ROOT}/rnaseq/expression/gene_expression_summary.tsv"

# Shared lib
LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/lib"

# ---- Output root ----
OUTDIR="${SCRATCH_ROOT}/diagnostics_v1"
SUMMARY_DIR="${OUTDIR}/per_archetype_summaries"   # per-motif intermediate TSVs
TABLES_DIR="${OUTDIR}/tables"                     # aggregated tables
PLOTS_DIR="${OUTDIR}/plots"                       # all plot outputs (.pdf + .png)
WASP_DIR="${OUTDIR}/wasp"                         # WASP cascade + per-SNP coverage data
LOGS_DIR="${OUTDIR}/logs"

# ---- Filter cascade defaults ----
# Standard chromosomes — chrX kept; chrM and chrY excluded.
STANDARD_CHROMS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chrX"

# Gene type whitelist (matches downstream analysis convention).
GENE_TYPES_INCLUDE="protein_coding,lncRNA"

# Expression filter: minimum tpm_mean across replicates for the host gene.
# Host gene is taken from 05_motif_annot's expression_tpm_max_sense column
# (already mean-of-replicates; we filter on this).
EXPRESSION_TPM_MIN=1.0

# ---- Archetype-level filter ----
# Default: only consider archetypes with motif width > MIN_MOTIF_WIDTH.
# Set to 0 to include all archetypes.
MIN_MOTIF_WIDTH=8

# ---- "Qualifying TFBS" definition ----
# A TFBS is "qualifying" if it passes the full cascade AND has a het-SNP
# within QUALIFYING_DISTANCE bp of its body. Default = 0 (overlap only).
# Other valid values: 10, 25, 50, 100, 200, 400, 800 (must match the
# flank bins in 05_motif_annot/config.sh's FLANK_DISTANCES).
QUALIFYING_DISTANCE=0

# Hybrid name (must match the prefix of snp_<hybrid>_overlap columns in 05).
HYBRID="F121-9"

# ---- Power thresholds (for "powered archetype" filter) ----
# Per-pair coverage threshold per assay. A pair is "powered" if its total
# pause-window reads (PRO-seq) or motif-body reads (ATAC) meet the threshold
# in BOTH alleles.
POWER_ATAC_MIN_PER_ALLELE=10
POWER_PROSEQ_PAUSE_MIN_PER_ALLELE=5
POWER_RNASEQ_AT_JUNCTION_MIN_PER_ALLELE=4

# ---- ATAC-coupling distance diagnostic ----
# Distance bins (bp from motif edge) for the ATAC asymmetry vs SNP-distance
# diagnostic. These should mirror 05_motif_annot's flank bins.
ATAC_DIST_BIN_EDGES="0,1,25,50,100,200,400,800"

# Number of permutations for the "background" (random SNP shuffling)
# distribution in the ATAC distance diagnostic.
ATAC_DIAGNOSTIC_N_PERMUTATIONS=20

# ---- Top-archetypes panel ----
# Archetypes with at least this many qualifying TFBSs go into the top-N panel.
# (Toggleable; the script also produces N=10 / 25 / 50 versions regardless.)
TOP_ARCHETYPE_MIN_QUALIFYING=1000

# ---- Top-N selections ----
TOP_N_LIST="10,25,50"

# ---- Sample names (must match BAM file prefixes upstream) ----
ATAC_SAMPLES=("F121-9_atac_rep1" "F121-9_atac_rep2" "F121-9_atac_rep3" "F121-9_atac_rep4")
PROSEQ_SAMPLES=("F121-9_proseq_rep1" "F121-9_proseq_rep2" "F121-9_proseq_rep3" "F121-9_proseq_rep4")
RNASEQ_SAMPLES=("WT_RNAseq_rep1" "WT_RNAseq_rep2")

# ---- Strict mode ----
# Default = lenient: pipeline works with whatever subset of motifs is currently
# annotated by 05/06. Set STRICT=1 to require all motifs to be present.
STRICT="${STRICT:-0}"

# ---- SLURM defaults ----
DEFAULT_PARTITION="short"
DEFAULT_TIME="01:00:00"
DEFAULT_MEM="16G"
DEFAULT_THREADS=2
DEFAULT_MAX_CONCURRENT=200

HEAVY_PARTITION="short"
HEAVY_TIME="06:00:00"
HEAVY_MEM="32G"

# ---- Environment ----
CONDA_ENV_R="test_m"
CONDA_ENV_TOOLS="proseq2_env"

# ---- Helpers ----
step_header() {
    echo ""
    echo "============================================================"
    echo "$1"
    echo "Time: $(date)"
    echo "============================================================"
}

create_dirs() {
    mkdir -p "${SUMMARY_DIR}" "${TABLES_DIR}" "${PLOTS_DIR}" \
             "${WASP_DIR}" "${LOGS_DIR}"
}
