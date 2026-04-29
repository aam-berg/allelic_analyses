#!/usr/bin/env bash
# =============================================================================
# 04c_rnaseq/config.sh — Central config for F121-9 RNA-seq pipeline
# =============================================================================
# Source this in every script in 04c_rnaseq/:  source "config.sh"
# =============================================================================

# --- Conda environment ---
# Same env as 03_atac/ and 04_proseq/ for shared tools, plus STAR + subread
# (featureCounts) which may need to be installed:
#   conda install -n proseq2_env -c bioconda star subread
# samtools, bedtools, bedGraphToBigWig, bigWigMerge, fasterq-dump, fastqc,
# cutadapt, pigz, pysam are also required.
CONDA_ENV_NAME="proseq2_env"

# --- Threads (per SBATCH allocation) ---
THREADS=8

# --- Top-level directories ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRATCH_DIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/rnaseq"
RESOURCE_DIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources"

# --- Reference genome (REUSE mm39 fasta from earlier pipelines) ---
GENOME_DIR="${RESOURCE_DIR}/genome"
MM39_FA="${GENOME_DIR}/mm39.fa"
MM39_CHROM_SIZES="${GENOME_DIR}/mm39.chrom.sizes"

# STAR index lives outside the bowtie2 index dir to avoid mixing
STAR_INDEX_DIR="${GENOME_DIR}/mm39_star"

# --- GTF annotation (vM38, mm39) ---
ANNOTATION_DIR="${RESOURCE_DIR}/annotation"
GTF_FILE="${ANNOTATION_DIR}/gencode.vM38.chr_patch_hapl_scaff.annotation.gtf.gz"

# --- VCF (for STAR's built-in WASP) ---
# The source VCF has TWO PARENTAL samples (129S1_SvImJ + CAST_EiJ), each
# scored as homozygous (1/1 or 0/0) at every site. STAR's --varVCFfile uses
# the first sample column and looks for 0/1 genotypes; with that layout it
# finds zero het records and aborts ("could not find any SNPs in VCF file").
#
# Step 00 synthesizes a STAR-compatible single-sample VCF
# (ALLELE_VCF_STAR) from the parental VCF: at every site where the parents
# are opposite homozygotes the F1 is heterozygous, so we emit one record
# per such site with a single F121_9 sample and GT=0/1. This is also where
# we drop the giant CSQ INFO blocks, shrinking the file from ~20 GB to a
# few hundred MB.
#
# ALLELE_VCF (parental, with full INFO/CSQ) remains the input for
# lib/make_wasp_snps.sh and lib/split_by_allele.py — those tools key off
# (chrom, pos, ref, alt) and don't care how genotypes are encoded.
VCF_DIR="${RESOURCE_DIR}/vcf"
ALLELE_VCF="${VCF_DIR}/F121-9_het_snps.ucsc.vcf.gz"
ALLELE_VCF_UNZIPPED="${VCF_DIR}/F121-9_het_snps.ucsc.vcf"
ALLELE_VCF_STAR="${VCF_DIR}/F121-9_het_snps.star.vcf"

# --- WASP SNP files (for lib/split_by_allele.py in step 06) ---
WASP_SNP_DIR="${SCRATCH_DIR}/wasp/snp_files"

# --- Working directories ---
FASTQ_DIR="${SCRATCH_DIR}/fastq"
FASTQC_RAW_DIR="${SCRATCH_DIR}/fastqc_raw"
TRIMMED_DIR="${SCRATCH_DIR}/trimmed"
FASTQC_TRIM_DIR="${SCRATCH_DIR}/fastqc_trimmed"

BAM_DIR="${SCRATCH_DIR}/bam"
BAM_STAR_DIR="${BAM_DIR}/star"               # raw STAR output (with vW WASP tags)
BAM_FINAL_DIR="${BAM_DIR}/final"             # WASP-passed reads only

QUANT_DIR="${SCRATCH_DIR}/quant"
QUANT_GENE_DIR="${QUANT_DIR}/gene_level"
QUANT_TX_DIR="${QUANT_DIR}/transcript_level"

BIGWIG_DIR="${SCRATCH_DIR}/bigwig"
BIGWIG_IND_DIR="${BIGWIG_DIR}/individual"
BIGWIG_MERGED_DIR="${BIGWIG_DIR}/merged"
BIGWIG_ALLELE_DIR="${BIGWIG_DIR}/allele_specific"
BIGWIG_ALLELE_MERGED_DIR="${BIGWIG_DIR}/allele_merged"

ALLELE_DIR="${SCRATCH_DIR}/allele_specific"

EXPRESSION_DIR="${SCRATCH_DIR}/expression"   # final expression tables for 05_motif_annot

QC_DIR="${SCRATCH_DIR}/qc"
LOG_DIR="${SCRATCH_DIR}/logs"

# --- Shared lib/ utilities ---
LIB_DIR="$(dirname "${SCRIPT_DIR}")/lib"
SPLIT_BY_ALLELE_PY="${LIB_DIR}/split_by_allele.py"
MAKE_WASP_SNPS_SH="${LIB_DIR}/make_wasp_snps.sh"

# --- Samples (loaded from samples.tsv) ---
SAMPLES_TSV="${SCRIPT_DIR}/samples.tsv"

# Note: SAMPLES values here are COMMA-SEPARATED SRR lists (multi-SRR per sample).
# Step 01 splits the list and concatenates downloaded FASTQs.
declare -A SAMPLES
SAMPLE_ORDER=()
while IFS=$'\t' read -r sample_name srrs condition; do
    [[ "${sample_name}" == "sample_name" || -z "${sample_name}" ]] && continue
    [[ "${sample_name}" =~ ^# ]] && continue
    SAMPLES["${sample_name}"]="${srrs}"
    SAMPLE_ORDER+=("${sample_name}")
done < "${SAMPLES_TSV}"

# --- TruSeq adapters (PE) ---
# Standard Illumina TruSeq adapters for stranded libraries.
ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

MIN_READ_LENGTH=20
QUALITY_CUTOFF=20

# --- Read length expected after trimming ---
# The original GSE200699 protocol trimmed reads to 100 nt. We'll stay close
# to that — the exact number isn't critical, but STAR's --sjdbOverhang should
# be (typical read length - 1).
READ_LENGTH=100
SJDB_OVERHANG=99

# --- STAR alignment parameters ---
# These match the parameters from the original GSE200699 paper where it
# makes sense (mismatch limit, junction filtering). We add 2-pass mode for
# better novel junction discovery and WASP for allele-specific work.
STAR_OUTFILTER_MISMATCH=4
STAR_OUTFILTER_MULTIMAP=10           # allow up to 10 multi-maps for sensitivity
STAR_OUTFILTER_TYPE="BySJout"
STAR_ALIGN_SJ_OVERHANG=8
STAR_INTRON_MOTIFS="RemoveNoncanonicalUnannotated"

# --- Library strandedness ---
# TruSeq Stranded Total RNA is REVERSE-STRANDED:
#   R1 maps to the antisense strand of the original RNA.
# Implications:
#   featureCounts -s 2 (reversely stranded)
#   STAR --quantMode GeneCounts produces a 4-col table; we read column 4
#   For PE strand-specific bigwigs:
#     RNA + strand fragment <=> R1 reverse + R2 forward
#     RNA - strand fragment <=> R1 forward + R2 reverse
LIBRARY_STRANDEDNESS="reverse"
FEATURECOUNTS_STRAND=2  # 0=unstranded, 1=forward, 2=reverse

# --- MAPQ threshold ---
# STAR's default MAPQ for unique mappers is 255. featureCounts/samtools
# treat 255 specially. For filtering, anything >=10 is essentially "unique
# in STAR's scheme" (multi-mappers get MAPQ 0/1/3 depending on count).
MAPQ_THRESHOLD=10

# =============================================================================
# Helpers
# =============================================================================

create_dirs() {
    mkdir -p "${FASTQ_DIR}" "${FASTQC_RAW_DIR}" "${TRIMMED_DIR}" "${FASTQC_TRIM_DIR}"
    mkdir -p "${BAM_STAR_DIR}" "${BAM_FINAL_DIR}"
    mkdir -p "${QUANT_GENE_DIR}" "${QUANT_TX_DIR}"
    mkdir -p "${BIGWIG_IND_DIR}" "${BIGWIG_MERGED_DIR}"
    mkdir -p "${BIGWIG_ALLELE_DIR}" "${BIGWIG_ALLELE_MERGED_DIR}"
    mkdir -p "${ALLELE_DIR}" "${EXPRESSION_DIR}"
    mkdir -p "${WASP_SNP_DIR}"
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

# Get the comma-separated SRR list as a bash array
get_srrs_for_sample() {
    local sample="$1"
    local srrs_str="${SAMPLES[${sample}]}"
    IFS=',' read -ra SRR_ARRAY <<< "${srrs_str}"
}

# resolve_samples — same pattern as 03_atac/ and 04_proseq/
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
                echo "[ERROR] Sample index ${idx} out of range" >&2
                exit 1
            fi
            RUN_SAMPLES=("${SAMPLE_ORDER[$idx]}")
            echo "[INFO] CLI mode: index ${idx} -> ${RUN_SAMPLES[0]}"
        else
            local found=0
            for s in "${SAMPLE_ORDER[@]}"; do [[ "$s" == "${arg}" ]] && found=1; done
            if (( found == 0 )); then
                echo "[ERROR] Unknown sample: ${arg}" >&2
                echo "[ERROR] Valid: ${SAMPLE_ORDER[*]}" >&2
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
