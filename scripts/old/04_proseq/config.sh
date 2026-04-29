#!/usr/bin/env bash
# =============================================================================
# config.sh — Central configuration for PRO-seq processing pipeline
# =============================================================================
# Source this file in every pipeline script:  source "config.sh"
# =============================================================================

# --- Conda environment ---
CONDA_ENV_NAME="proseq2_env"

# --- Thread count (adjust per SBATCH allocation) ---
THREADS=8

# --- Top-level directories ---
SCRIPT_DIR="/home/alb1273/pausing_phase_project/scripts/04_proseq"
SCRATCH_DIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/proseq"
RESOURCE_DIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources"

# --- Reference genome paths (EXISTING — from your resources directory) ---
GENOME_DIR="${RESOURCE_DIR}/genome"
MM39_FA="${GENOME_DIR}/mm39.fa"
MM39_BT2_IDX="${GENOME_DIR}/mm39"          # prefix for mm39.1.bt2, mm39.2.bt2, etc.
MM39_CHROM_SIZES="${GENOME_DIR}/mm39.chrom.sizes"
DM6_FA="${GENOME_DIR}/dm6.fa"
DM6_BT2_IDX="${GENOME_DIR}/dm6"            # prefix for dm6.1.bt2, dm6.2.bt2, etc.

# --- rRNA index (built by 00_setup_references.sh) ---
RRNA_DIR="${SCRATCH_DIR}/references/rRNA"
RRNA_BT2_IDX="${RRNA_DIR}/mouse_rRNA"

# --- VCF directory (for future WASP analysis) ---
VCF_DIR="${RESOURCE_DIR}/vcf"

# --- Working data directories ---
FASTQ_DIR="${SCRATCH_DIR}/fastq"
FASTQC_RAW_DIR="${SCRATCH_DIR}/fastqc_raw"
TRIMMED_DIR="${SCRATCH_DIR}/trimmed"
FASTQC_TRIM_DIR="${SCRATCH_DIR}/fastqc_trimmed"
BAM_DIR="${SCRATCH_DIR}/bam"
BAM_DM6_DIR="${BAM_DIR}/dm6"
BAM_RRNA_DIR="${BAM_DIR}/rRNA"
BAM_MM39_DIR="${BAM_DIR}/mm39"
BIGWIG_DIR="${SCRATCH_DIR}/bigwig"
BIGWIG_IND_DIR="${BIGWIG_DIR}/individual"
BIGWIG_MERGED_DIR="${BIGWIG_DIR}/merged"
QC_DIR="${SCRATCH_DIR}/qc"
LOG_DIR="${SCRATCH_DIR}/logs"

# --- Sample information ---
# Format: SAMPLE_NAME=SRR_ACCESSION
# These are F121-9 WT mESC PRO-seq replicates from GSE190548
declare -A SAMPLES
SAMPLES=(
    ["WT_PROseq_rep1"]="SRR17640044"
    ["WT_PROseq_rep2"]="SRR17640045"
    ["WT_PROseq_rep3"]="SRR17640046"
    ["WT_PROseq_rep4"]="SRR17640047"
)

# Ordered sample names — INDEX INTO THIS ARRAY for task arrays
SAMPLE_ORDER=("WT_PROseq_rep1" "WT_PROseq_rep2" "WT_PROseq_rep3" "WT_PROseq_rep4")

# --- PRO-seq adapter sequence ---
# Standard Illumina TruSeq Small RNA 3' adapter, used in Mahat et al. 2016.
ADAPTER_3PRIME="TGGAATTCTCGGGTGCCAAGG"

# --- Alignment parameters ---
MAPQ_THRESHOLD=10

# --- Trimming parameters ---
MIN_READ_LENGTH=20
QUALITY_CUTOFF=20

# --- Optional: 5' end deduplication ---
# PRO-seq note: normally we do NOT deduplicate because identical 5' positions
# can reflect real Pol II density. Enable this ONLY in special cases
DEDUP_5PRIME="${DEDUP_5PRIME:-true}"   #"${DEDUP_5PRIME:-false}"
DEDUP_SUFFIX="_5prime_dedup"

# =============================================================================
# Helper: create all directories
# =============================================================================
create_dirs() {
    mkdir -p "${FASTQ_DIR}" "${FASTQC_RAW_DIR}" "${TRIMMED_DIR}" "${FASTQC_TRIM_DIR}"
    mkdir -p "${BAM_DM6_DIR}" "${BAM_RRNA_DIR}" "${BAM_MM39_DIR}"
    mkdir -p "${BIGWIG_IND_DIR}" "${BIGWIG_MERGED_DIR}"
    mkdir -p "${QC_DIR}" "${LOG_DIR}" "${RRNA_DIR}"
}

# =============================================================================
# Helper: print a step header
# =============================================================================
step_header() {
    echo ""
    echo "============================================================"
    echo "$1"
    echo "Time: $(date)"
    echo "============================================================"
}

# =============================================================================
# Helper: resolve which sample(s) to process
#
# Supports three modes:
#   1. SLURM array job:  SLURM_ARRAY_TASK_ID is set → process that index
#   2. CLI argument:     $1 is a number (0-3) or sample name → process that one
#   3. No argument:      process ALL samples sequentially
#
# Sets:  RUN_SAMPLES=( list of sample names to process )
#
# Usage in scripts:
#   resolve_samples "$@"
#   for SAMPLE in "${RUN_SAMPLES[@]}"; do ... done
# =============================================================================
resolve_samples() {
    local arg="${1:-}"

    # Priority 1: SLURM_ARRAY_TASK_ID
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        local idx="${SLURM_ARRAY_TASK_ID}"
        if (( idx < 0 || idx >= ${#SAMPLE_ORDER[@]} )); then
            echo "[ERROR] SLURM_ARRAY_TASK_ID=${idx} is out of range (0-$((${#SAMPLE_ORDER[@]}-1)))"
            exit 1
        fi
        RUN_SAMPLES=("${SAMPLE_ORDER[$idx]}")
        echo "[INFO] SLURM array mode: task ${idx} → ${RUN_SAMPLES[0]}"
        return
    fi

    # Priority 2: command-line argument
    if [[ -n "${arg}" ]]; then
        # If it's a number, treat as index
        if [[ "${arg}" =~ ^[0-9]+$ ]]; then
            local idx="${arg}"
            if (( idx < 0 || idx >= ${#SAMPLE_ORDER[@]} )); then
                echo "[ERROR] Sample index ${idx} is out of range (0-$((${#SAMPLE_ORDER[@]}-1)))"
                exit 1
            fi
            RUN_SAMPLES=("${SAMPLE_ORDER[$idx]}")
            echo "[INFO] CLI mode: index ${idx} → ${RUN_SAMPLES[0]}"
        else
            # Treat as sample name
            local found=0
            for s in "${SAMPLE_ORDER[@]}"; do
                if [[ "$s" == "${arg}" ]]; then found=1; break; fi
            done
            if (( found == 0 )); then
                echo "[ERROR] Unknown sample name: ${arg}"
                echo "[ERROR] Valid names: ${SAMPLE_ORDER[*]}"
                exit 1
            fi
            RUN_SAMPLES=("${arg}")
            echo "[INFO] CLI mode: sample → ${RUN_SAMPLES[0]}"
        fi
        return
    fi

    # Priority 3: no argument — run all samples
    RUN_SAMPLES=("${SAMPLE_ORDER[@]}")
    echo "[INFO] Sequential mode: processing all ${#RUN_SAMPLES[@]} samples"
}