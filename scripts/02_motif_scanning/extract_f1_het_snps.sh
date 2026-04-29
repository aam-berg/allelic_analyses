#!/bin/bash
#SBATCH --job-name=extract_f1_het_snps
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/extract_f1_het_snps_%j.out
#SBATCH --error=logs/extract_f1_het_snps_%j.err

# =============================================================================
# Extract F1 Heterozygous SNPs from Mouse Genomes Project VCF
# =============================================================================
# For inbred strains, each strain is homozygous (0/0 or 1/1).
# F1 hybrids are heterozygous at sites where one parent is 0/0 and the other is 1/1.
#
# This script extracts:
#   - F121-9: 129S1_SvImJ × CAST_EiJ heterozygous sites
#   - BL6 × CAST: C57BL_6NJ × CAST_EiJ heterozygous sites
# =============================================================================

set -euo pipefail

# Configuration - adjust this path as needed
RESOURCE_DIR="${1:-/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources}"
VCF_DIR="${RESOURCE_DIR}/vcf"

# Ensure log directory exists
mkdir -p logs

# Load bcftools (adjust for your cluster)
source activate motif_scanning

cd "${VCF_DIR}"

echo "=============================================="
echo "Extracting F1 heterozygous SNPs: $(date)"
echo "Working directory: ${VCF_DIR}"
echo "=============================================="

# Verify input file exists
if [[ ! -f "mgp_REL2021_snps.vcf.gz" ]]; then
    echo "ERROR: mgp_REL2021_snps.vcf.gz not found in ${VCF_DIR}"
    echo "Please run the download script first."
    exit 1
fi

# Check available samples
echo ""
echo "Available strains in VCF:"
bcftools query -l mgp_REL2021_snps.vcf.gz | grep -E "129S1|CAST|C57BL" || true
echo ""

# =============================================================================
# F121-9: 129S1_SvImJ × CAST_EiJ
# =============================================================================
echo "[1/2] Extracting F121-9 heterozygous SNPs (129S1_SvImJ × CAST_EiJ)..."

# Remove any broken previous attempt
rm -f F121-9_het_snps.vcf.gz F121-9_het_snps.vcf.gz.tbi

# Extract SNPs where 129S1 and CAST have DIFFERENT homozygous genotypes
# This gives us sites that will be heterozygous in the F1 hybrid
bcftools view \
    -s 129S1_SvImJ,CAST_EiJ \
    -v snps \
    -m2 -M2 \
    mgp_REL2021_snps.vcf.gz \
| bcftools view \
    -i '(GT[0]="0/0" && GT[1]="1/1") || (GT[0]="1/1" && GT[1]="0/0")' \
    -Oz -o F121-9_het_snps.vcf.gz

# Index the VCF
tabix -p vcf F121-9_het_snps.vcf.gz

# Count and report
F121_COUNT=$(bcftools view -H F121-9_het_snps.vcf.gz | wc -l)
echo "  Found ${F121_COUNT} heterozygous SNPs for F121-9"

# =============================================================================
# BL6 × CAST: C57BL_6NJ × CAST_EiJ
# =============================================================================
echo ""
echo "[2/2] Extracting BL6×CAST heterozygous SNPs (C57BL_6NJ × CAST_EiJ)..."

# Remove any broken previous attempt
rm -f BL6xCAST_het_snps.vcf.gz BL6xCAST_het_snps.vcf.gz.tbi

bcftools view \
    -s C57BL_6NJ,CAST_EiJ \
    -v snps \
    -m2 -M2 \
    mgp_REL2021_snps.vcf.gz \
| bcftools view \
    -i '(GT[0]="0/0" && GT[1]="1/1") || (GT[0]="1/1" && GT[1]="0/0")' \
    -Oz -o BL6xCAST_het_snps.vcf.gz

tabix -p vcf BL6xCAST_het_snps.vcf.gz

BL6_COUNT=$(bcftools view -H BL6xCAST_het_snps.vcf.gz | wc -l)
echo "  Found ${BL6_COUNT} heterozygous SNPs for BL6×CAST"

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "=============================================="
echo "Extraction complete: $(date)"
echo "=============================================="
echo ""
echo "Output files (Ensembl chromosome naming: 1, 2, 3...):"
echo "  ${VCF_DIR}/F121-9_het_snps.vcf.gz   (${F121_COUNT} SNPs)"
echo "  ${VCF_DIR}/BL6xCAST_het_snps.vcf.gz (${BL6_COUNT} SNPs)"
echo ""
echo "Quick validation - first 5 records from F121-9:"
bcftools view -H F121-9_het_snps.vcf.gz | head -5 | cut -f1-5,10-11
echo ""
echo "To proceed with UCSC-style VCFs (chr1, chr2...), run section 4 of the main script."