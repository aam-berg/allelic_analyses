#!/usr/bin/env bash
# =============================================================================
# 03_atac/config.sh — Central configuration for the F121-9 ATAC-seq pipeline
# =============================================================================
# Source this in every script in 03_atac/:  source "config.sh"
#
# All paths are read from this file. There must be NO hardcoded paths inside
# the per-step scripts.
# =============================================================================

# --- Conda environment ---
# Same env as 04_proseq/ — both pipelines need bowtie2, samtools, bedtools,
# bedGraphToBigWig, fasterq-dump, fastqc, cutadapt, pigz. Plus MACS2 for ATAC.
# If MACS3 is not installed:
#   pip install macs3
CONDA_ENV_NAME="proseq2_env"

# --- Threads (adjust per SBATCH allocation) ---
THREADS=8

# --- Top-level directories ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Output root for ATAC. Kept separate from 04_proseq's output root to avoid
# any cross-contamination of intermediate files.
SCRATCH_DIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/atac_proc"

# Resources (REUSED from 01_resources/ and 04_proseq/ where possible)
RESOURCE_DIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources"

# --- Reference genome (REUSE 04_proseq's mm39 bowtie2 index if present) ---
GENOME_DIR="${RESOURCE_DIR}/genome"
MM39_FA="${GENOME_DIR}/mm39.fa"
MM39_BT2_IDX="${GENOME_DIR}/mm39"           # prefix for mm39.{1,2,3,4,rev.1,rev.2}.bt2
MM39_CHROM_SIZES="${GENOME_DIR}/mm39.chrom.sizes"

# --- VCF (for WASP) ---
VCF_DIR="${RESOURCE_DIR}/vcf"
ALLELE_VCF="${VCF_DIR}/F121-9_het_snps.ucsc.vcf.gz"

# --- Working data directories ---
FASTQ_DIR="${SCRATCH_DIR}/fastq"
FASTQC_RAW_DIR="${SCRATCH_DIR}/fastqc_raw"
TRIMMED_DIR="${SCRATCH_DIR}/trimmed"
FASTQC_TRIM_DIR="${SCRATCH_DIR}/fastqc_trimmed"
BAM_DIR="${SCRATCH_DIR}/bam"
BAM_RAW_DIR="${BAM_DIR}/raw"               # post-bowtie2, before filtering
BAM_FINAL_DIR="${BAM_DIR}/final"           # MAPQ + chrM + dedup filtered

PEAKS_DIR="${SCRATCH_DIR}/peaks"
PEAKS_PER_REP_DIR="${PEAKS_DIR}/per_replicate"
PEAKS_CONSENSUS_DIR="${PEAKS_DIR}/consensus"

WASP_DIR="${SCRATCH_DIR}/wasp"
WASP_REPO_DIR="${SCRATCH_DIR}/tools/WASP"
WASP_FIND="${WASP_REPO_DIR}/mapping/find_intersecting_snps.py"
WASP_FILTER="${WASP_REPO_DIR}/mapping/filter_remapped_reads.py"
WASP_SNP_DIR="${WASP_DIR}/snp_files"
WASP_INTERMEDIATE_DIR="${WASP_DIR}/intermediate"
WASP_FILTERED_DIR="${WASP_DIR}/filtered_bam"

ALLELE_DIR="${SCRATCH_DIR}/allele_specific"

BIGWIG_DIR="${SCRATCH_DIR}/bigwig"
BIGWIG_IND_DIR="${BIGWIG_DIR}/individual"
BIGWIG_MERGED_DIR="${BIGWIG_DIR}/merged"
BIGWIG_ALLELE_DIR="${BIGWIG_DIR}/allele_specific"

QC_DIR="${SCRATCH_DIR}/qc"
LOG_DIR="${SCRATCH_DIR}/logs"

# --- Shared lib/ ---
# scripts/03_atac/.. -> scripts/  -> lib/
LIB_DIR="$(dirname "${SCRIPT_DIR}")/lib"
SPLIT_BY_ALLELE_PY="${LIB_DIR}/split_by_allele.py"
MAKE_WASP_SNPS_SH="${LIB_DIR}/make_wasp_snps.sh"

# --- Samples (loaded from samples.tsv, the canonical source) ---
SAMPLES_TSV="${SCRIPT_DIR}/samples.tsv"

# Build SAMPLES (associative: name -> SRR) and SAMPLE_ORDER (positional list)
# from samples.tsv. Skips header row, blank lines, and lines starting with #.
declare -A SAMPLES
SAMPLE_ORDER=()
while IFS=$'\t' read -r sample_name srr condition; do
    [[ "${sample_name}" == "sample_name" || -z "${sample_name}" ]] && continue
    [[ "${sample_name}" =~ ^# ]] && continue
    SAMPLES["${sample_name}"]="${srr}"
    SAMPLE_ORDER+=("${sample_name}")
done < "${SAMPLES_TSV}"

# --- ATAC trimming parameters ---
# Nextera adapters are the standard for ATAC-seq libraries. The 19-mer
# below is the conserved core that cutadapt anchors on (longer adapter
# variants share this core).
ADAPTER_R1="CTGTCTCTTATACACATCT"
ADAPTER_R2="CTGTCTCTTATACACATCT"
MIN_READ_LENGTH=20
QUALITY_CUTOFF=20

# --- Alignment parameters ---
# bowtie2 maximum insert size for PE ATAC. ATAC fragments range from ~50 bp
# (sub-nucleosomal) to ~1000 bp (di- and tri-nucleosomal). 1000 captures the
# full distribution without admitting unrealistically large inserts.
MAX_INSERT=1000

# MAPQ filter: 30 is the ATAC-seq community standard (vs 10 for PRO-seq).
# Higher because PE provides more positional info, so multi-mappers are
# more confidently identified.
MAPQ_THRESHOLD=30

# --- Peak calling parameters (MACS2) ---
MACS2_QVALUE=0.01
MACS2_GSIZE="mm"   # MACS2 shorthand for the mouse effective genome size

# --- Consensus peak parameters ---
# A peak is in the consensus set if it appears in this many replicates.
# 2 of 4 is a common threshold; tightens to 3 or 4 if you want stricter
# reproducibility.
MIN_REPS_FOR_CONSENSUS=2

# =============================================================================
# Helpers (mirror the conventions in 04_proseq/config.sh)
# =============================================================================

create_dirs() {
    mkdir -p "${FASTQ_DIR}" "${FASTQC_RAW_DIR}" "${TRIMMED_DIR}" "${FASTQC_TRIM_DIR}"
    mkdir -p "${BAM_RAW_DIR}" "${BAM_FINAL_DIR}"
    mkdir -p "${PEAKS_PER_REP_DIR}" "${PEAKS_CONSENSUS_DIR}"
    mkdir -p "${WASP_SNP_DIR}" "${WASP_INTERMEDIATE_DIR}" "${WASP_FILTERED_DIR}"
    mkdir -p "${ALLELE_DIR}"
    mkdir -p "${BIGWIG_IND_DIR}" "${BIGWIG_MERGED_DIR}" "${BIGWIG_ALLELE_DIR}"
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

# resolve_samples — three priority levels:
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
