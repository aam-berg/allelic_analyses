#!/bin/bash
#SBATCH --job-name=download_resources
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/download_resources_%j.out
#SBATCH --error=logs/download_resources_%j.err

# =============================================================================
# 01_download_resources.sbatch — Download genome (mm39) + raw MGP VCF
#
# Trimmed: this script ONLY downloads. Per-hybrid het SNP extraction is in
# 02_extract_f1_het_snps.sbatch; chrom-naming conversion is in
# 03_convert_vcf_to_ucsc.sh.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

mkdir -p logs "${GENOME_DIR}" "${VCF_DIR}"

set +u
source activate "${CONDA_ENV}"
set -u

echo "============================================"
echo "Download mm39 resources"
echo "  Genome dir:  ${GENOME_DIR}"
echo "  VCF dir:     ${VCF_DIR}"
echo "  MGP root:    ${MGP_FTP}"
echo "  Threads:     ${THREADS}"
echo "============================================"

# ---- 1. Genome FASTA (UCSC naming) ----
if [[ -f "${GENOME_FA_UCSC}" && -f "${GENOME_FA_UCSC}.fai" ]]; then
    echo "Genome already downloaded: ${GENOME_FA_UCSC}"
else
    cd "${GENOME_DIR}"
    if [[ ! -f mm39.fa.gz ]]; then
        wget -q https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
    fi
    gunzip -k mm39.fa.gz
    samtools faidx mm39.fa
    cd - >/dev/null
fi

# ---- 2. Chrom sizes (UCSC) ----
if [[ ! -f "${CHROM_SIZES_UCSC}" ]]; then
    cut -f1,2 "${GENOME_FA_UCSC}.fai" > "${CHROM_SIZES_UCSC}"
    echo "Wrote: ${CHROM_SIZES_UCSC}"
fi

# ---- 3. Raw MGP VCF + tabix index ----
if [[ -f "${MGP_VCF}" && -f "${MGP_VCF}.tbi" ]]; then
    echo "MGP VCF already downloaded: ${MGP_VCF}"
else
    cd "${VCF_DIR}"
    if [[ ! -f mgp_REL2021_snps.vcf.gz ]]; then
        wget -q "${MGP_FTP}/mgp_REL2021_snps.vcf.gz"
    fi
    if [[ ! -f mgp_REL2021_snps.vcf.gz.tbi ]]; then
        if [[ -f "${MGP_FTP}/mgp_REL2021_snps.vcf.gz.tbi" ]]; then
            wget -q "${MGP_FTP}/mgp_REL2021_snps.vcf.gz.tbi"
        else
            tabix -p vcf mgp_REL2021_snps.vcf.gz
        fi
    fi
    cd - >/dev/null
fi

echo ""
echo "Download complete. Summary:"
ls -la "${GENOME_FA_UCSC}" "${GENOME_FA_UCSC}.fai" "${CHROM_SIZES_UCSC}" "${MGP_VCF}" "${MGP_VCF}.tbi"
echo ""
echo "Next: bash 02_extract_f1_het_snps.sbatch"
