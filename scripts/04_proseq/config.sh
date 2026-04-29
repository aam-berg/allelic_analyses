#!/usr/bin/env bash
# =============================================================================
# 04_proseq/config.sh — Central configuration for the F121-9 PRO-seq pipeline
# =============================================================================
# Source this in every script in 04_proseq/:  source "config.sh"
#
# All paths, parameters, and helpers live here. There must be NO hardcoded
# paths inside individual step scripts.
#
# This file is the consolidated successor to the old config.sh +
# config_wasp.sh split — having two configs was fragile.
# =============================================================================

# --- Conda environment ---
# Same env as 03_atac/. Both pipelines need: bowtie2, samtools, bedtools,
# bedGraphToBigWig, bigWigMerge, fasterq-dump, fastqc, cutadapt, pigz, pysam.
# Plus entrez-direct (efetch) for the rRNA reference download in step 00.
CONDA_ENV_NAME="proseq2_env"

# --- Threads (adjust per SBATCH allocation) ---
THREADS=8

# --- Top-level directories ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Output root for PRO-seq. Kept separate from 03_atac's output root.
SCRATCH_DIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/proseq"

# Resources (REUSED from 01_resources/ and shared with 03_atac/)
RESOURCE_DIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources"

# --- Reference genomes (REUSE indices from 03_atac or earlier 04_proseq runs) ---
GENOME_DIR="${RESOURCE_DIR}/genome"
MM39_FA="${GENOME_DIR}/mm39.fa"
MM39_BT2_IDX="${GENOME_DIR}/mm39"          # prefix for mm39.{1,2,3,4,rev.1,rev.2}.bt2
MM39_CHROM_SIZES="${GENOME_DIR}/mm39.chrom.sizes"
DM6_FA="${GENOME_DIR}/dm6.fa"
DM6_BT2_IDX="${GENOME_DIR}/dm6"            # prefix for dm6 bowtie2 index

# --- rRNA index (built fresh in 00_setup_references.sh; PRO-seq specific) ---
RRNA_DIR="${SCRATCH_DIR}/references/rRNA"
RRNA_FA="${RRNA_DIR}/mouse_rRNA.fa"
RRNA_BT2_IDX="${RRNA_DIR}/mouse_rRNA"

# --- VCF (for WASP) ---
VCF_DIR="${RESOURCE_DIR}/vcf"
ALLELE_VCF="${VCF_DIR}/F121-9_het_snps.ucsc.vcf.gz"

# --- Working data directories ---
FASTQ_DIR="${SCRATCH_DIR}/fastq"
FASTQC_RAW_DIR="${SCRATCH_DIR}/fastqc_raw"
TRIMMED_DIR="${SCRATCH_DIR}/trimmed"
FASTQC_TRIM_DIR="${SCRATCH_DIR}/fastqc_trimmed"

BAM_DIR="${SCRATCH_DIR}/bam"
BAM_DM6_DIR="${BAM_DIR}/dm6"               # only used if DM6_SPIKEIN_FILTER=true
BAM_RRNA_DIR="${BAM_DIR}/rRNA"
BAM_MM39_DIR="${BAM_DIR}/mm39"

BIGWIG_DIR="${SCRATCH_DIR}/bigwig"
BIGWIG_IND_DIR="${BIGWIG_DIR}/individual"
BIGWIG_MERGED_DIR="${BIGWIG_DIR}/merged"
BIGWIG_ALLELE_DIR="${BIGWIG_DIR}/allele_specific"
BIGWIG_ALLELE_MERGED_DIR="${BIGWIG_DIR}/allele_merged"

# WASP and allele-specific outputs (consolidated under SCRATCH_DIR/wasp)
WASP_DIR="${SCRATCH_DIR}/wasp"
WASP_REPO_DIR="${SCRATCH_DIR}/tools/WASP"
WASP_FIND="${WASP_REPO_DIR}/mapping/find_intersecting_snps.py"
WASP_FILTER="${WASP_REPO_DIR}/mapping/filter_remapped_reads.py"
WASP_SNP_DIR="${WASP_DIR}/snp_files"
WASP_INTERMEDIATE_DIR="${WASP_DIR}/intermediate"
WASP_FILTERED_DIR="${WASP_DIR}/filtered_bam"

ALLELE_DIR="${SCRATCH_DIR}/allele_specific"   # split BAMs go here

QC_DIR="${SCRATCH_DIR}/qc"
LOG_DIR="${SCRATCH_DIR}/logs"

# --- Shared lib/ utilities ---
# scripts/04_proseq/.. -> scripts/  -> lib/
LIB_DIR="$(dirname "${SCRIPT_DIR}")/lib"
SPLIT_BY_ALLELE_PY="${LIB_DIR}/split_by_allele.py"
MAKE_WASP_SNPS_SH="${LIB_DIR}/make_wasp_snps.sh"

# --- Samples (loaded from samples.tsv, the canonical source) ---
SAMPLES_TSV="${SCRIPT_DIR}/samples.tsv"

declare -A SAMPLES
SAMPLE_ORDER=()
while IFS=$'\t' read -r sample_name srr condition; do
    [[ "${sample_name}" == "sample_name" || -z "${sample_name}" ]] && continue
    [[ "${sample_name}" =~ ^# ]] && continue
    SAMPLES["${sample_name}"]="${srr}"
    SAMPLE_ORDER+=("${sample_name}")
done < "${SAMPLES_TSV}"

# --- PRO-seq adapter sequence ---
# Standard Illumina TruSeq Small RNA 3' adapter (Mahat et al. 2016).
ADAPTER_3PRIME="TGGAATTCTCGGGTGCCAAGG"

# --- Trimming parameters ---
MIN_READ_LENGTH=20
QUALITY_CUTOFF=20

# --- Alignment parameters ---
# MAPQ 10: PRO-seq community standard (vs ATAC's 30). Lower because SE reads
# carry less positional information than PE, so a stricter cutoff would
# remove too much signal.
MAPQ_THRESHOLD=10

# --- Spike-in filtering ---
# The original protocol of this dataset reportedly used Drosophila S2 cells
# as a spike-in normalization control. If your data don't have a spike-in,
# set this to false to skip the dm6 alignment step (saves ~1-2 hr per sample).
#
# To verify whether spike-in is present, run the pipeline once on a single
# sample and check the dm6 alignment rate in step 03's QC output:
#   - >0.5% aligning to dm6   -> spike-in is present, leave this true
#   - <0.5% aligning to dm6   -> no spike-in, set this to false
DM6_SPIKEIN_FILTER=true

# =============================================================================
# Helpers
# =============================================================================

create_dirs() {
    mkdir -p "${FASTQ_DIR}" "${FASTQC_RAW_DIR}" "${TRIMMED_DIR}" "${FASTQC_TRIM_DIR}"
    mkdir -p "${BAM_RRNA_DIR}" "${BAM_MM39_DIR}"
    [[ "${DM6_SPIKEIN_FILTER}" == "true" ]] && mkdir -p "${BAM_DM6_DIR}"
    mkdir -p "${BIGWIG_IND_DIR}" "${BIGWIG_MERGED_DIR}"
    mkdir -p "${BIGWIG_ALLELE_DIR}" "${BIGWIG_ALLELE_MERGED_DIR}"
    mkdir -p "${WASP_SNP_DIR}" "${WASP_INTERMEDIATE_DIR}" "${WASP_FILTERED_DIR}"
    mkdir -p "${ALLELE_DIR}"
    mkdir -p "${RRNA_DIR}"
    mkdir -p "${QC_DIR}" "${LOG_DIR}"
    mkdir -p "${SCRATCH_DIR}/tmp"
}

step_header() {
    echo ""
    echo "============================================================"
    echo "$1"
    echo "Time: $(date)"
    echo "============================================================"
}

# resolve_samples — three priority levels (mirrors 03_atac/config.sh):
#   1. SLURM_ARRAY_TASK_ID is set  -> run that one sample
#   2. CLI arg is given (index or name) -> run that one sample
#   3. No arg -> run all samples sequentially
resolve_samples() {
    local arg="${1:-}"
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        local idx="${SLURM_ARRAY_TASK_ID}"
        if (( idx < 0 || idx >= ${#SAMPLE_ORDER[@]} )); then
            echo "[ERROR] SLURM_ARRAY_TASK_ID=${idx} out of range (0-$((${#SAMPLE_ORDER[@]}-1)))" >&2
            exit 1
        fi
        RUN_SAMPLES=("${SAMPLE_ORDER[$idx]}")
        echo "[INFO] SLURM array mode: task ${idx} -> ${RUN_SAMPLES[0]}"
        return
    fi
    if [[ -n "${arg}" ]]; then
        if [[ "${arg}" =~ ^[0-9]+$ ]]; then
            local idx="${arg}"
            if (( idx < 0 || idx >= ${#SAMPLE_ORDER[@]} )); then
                echo "[ERROR] Sample index ${idx} out of range (0-$((${#SAMPLE_ORDER[@]}-1)))" >&2
                exit 1
            fi
            RUN_SAMPLES=("${SAMPLE_ORDER[$idx]}")
            echo "[INFO] CLI mode: index ${idx} -> ${RUN_SAMPLES[0]}"
        else
            local found=0
            for s in "${SAMPLE_ORDER[@]}"; do [[ "$s" == "${arg}" ]] && found=1; done
            if (( found == 0 )); then
                echo "[ERROR] Unknown sample: ${arg}" >&2
                echo "[ERROR] Valid samples: ${SAMPLE_ORDER[*]}" >&2
                exit 1
            fi
            RUN_SAMPLES=("${arg}")
            echo "[INFO] CLI mode: sample -> ${RUN_SAMPLES[0]}"
        fi
        return
    fi
    RUN_SAMPLES=("${SAMPLE_ORDER[@]}")
    echo "[INFO] Sequential mode: processing all ${#RUN_SAMPLES[@]} samples"
}
