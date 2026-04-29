#!/bin/bash
#SBATCH --job-name=extract_f1_het_snps
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/extract_f1_het_snps_%j.out
#SBATCH --error=logs/extract_f1_het_snps_%j.err

# =============================================================================
# 02_extract_f1_het_snps.sbatch — Per-hybrid heterozygous-SNP VCFs
#
# For each F1 hybrid (defined in config.sh), extract from the MGP VCF only
# SNPs where the two parent strains have different homozygous genotypes
# (i.e. the F1 will be heterozygous). Output is in Ensembl chromosome
# naming; convert with 03_convert_vcf_to_ucsc.sh.
#
# Output naming: ${VCF_DIR}/<HYBRID>_het_snps.vcf.gz
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

mkdir -p logs

set +u
source activate "${CONDA_ENV}"
set -u

echo "============================================"
echo "F1 het-SNP extraction"
echo "  MGP VCF:   ${MGP_VCF}"
echo "  Output:    ${VCF_DIR}"
echo "============================================"

# extract_one HYBRID_LABEL "PARENT1,PARENT2"
extract_one() {
    local label="$1"
    local parents="$2"
    local p1="${parents%,*}"
    local p2="${parents#*,}"
    local out="${VCF_DIR}/${label}_het_snps.vcf.gz"

    echo ""
    echo "[${label}]  parents: ${p1}  ×  ${p2}"

    if [[ -f "${out}" && -f "${out}.tbi" ]]; then
        echo "  ${out} exists; skipping."
        return 0
    fi

    # Restrict to SNPs where the two parents differ in homozygous genotype
    # (assumes ploidy 2; MGP gives genotypes as 0/0, 1/1, etc.)
    bcftools view \
        --threads "${THREADS}" \
        -s "${p1},${p2}" \
        -v snps \
        -m2 -M2 \
        --min-ac 1 \
        "${MGP_VCF}" \
    | bcftools view \
        -e "GT[0]=\"./.\" || GT[1]=\"./.\" || GT[0]=GT[1] || GT[0]=\"0/1\" || GT[1]=\"0/1\"" \
        -Oz -o "${out}"

    tabix -p vcf "${out}"
    n=$(bcftools view -H "${out}" | wc -l)
    echo "  Wrote ${out}  (${n} SNPs)"
}

extract_one "F121-9"   "${HYBRID_F121_PARENTS}"
extract_one "BL6xCAST" "${HYBRID_BL6CAST_PARENTS}"

echo ""
echo "Done. Convert to UCSC naming with: bash 03_convert_vcf_to_ucsc.sh"
