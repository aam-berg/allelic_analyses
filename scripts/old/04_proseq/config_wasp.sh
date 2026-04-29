#!/usr/bin/env bash
# =============================================================================
# config_wasp.sh — WASP-specific configuration (sources main config.sh)
# =============================================================================
# This extends config.sh with WASP-related paths.
# Source this in WASP pipeline scripts instead of config.sh directly.
# =============================================================================

source "config.sh"

# --- WASP installation ---
WASP_DIR="${SCRATCH_DIR}/tools/WASP"
WASP_FIND="${WASP_DIR}/mapping/find_intersecting_snps.py"
WASP_FILTER="${WASP_DIR}/mapping/filter_remapped_reads.py"

# --- WASP SNP directory (per-chromosome SNP files) ---
WASP_SNP_DIR="${SCRATCH_DIR}/wasp/snp_files"

# --- WASP intermediate and output directories ---
WASP_INTERMEDIATE_DIR="${SCRATCH_DIR}/wasp/intermediate"
WASP_FILTERED_DIR="${SCRATCH_DIR}/wasp/filtered_bam"
WASP_ALLELE_DIR="${SCRATCH_DIR}/wasp/allele_specific"
WASP_ALLELE_BW_DIR="${SCRATCH_DIR}/wasp/allele_bigwig"
WASP_QC_DIR="${SCRATCH_DIR}/wasp/qc"

# --- Which VCF to use for allele assignment ---
# UCSC-style (chr-prefixed) to match our BAMs
ALLELE_VCF="${VCF_DIR}/F121-9_het_snps.ucsc.vcf.gz"

# --- Helper: create WASP directories ---
create_wasp_dirs() {
    mkdir -p "${WASP_SNP_DIR}" "${WASP_INTERMEDIATE_DIR}" "${WASP_FILTERED_DIR}"
    mkdir -p "${WASP_ALLELE_DIR}" "${WASP_ALLELE_BW_DIR}" "${WASP_QC_DIR}"
}