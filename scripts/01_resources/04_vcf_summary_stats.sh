#!/bin/bash
# =============================================================================
# 04_vcf_summary_stats.sh — Summary statistics for per-hybrid het-SNP VCFs
#
# Login-node OK. Prints chromosome distribution and total SNP count.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

set +u
source activate "${CONDA_ENV}"
set -u

summary_one() {
    local label="$1"
    local vcf="${VCF_DIR}/${label}_het_snps.vcf.gz"
    local vcf_ucsc="${VCF_DIR}/${label}_het_snps.ucsc.vcf.gz"

    echo ""
    echo "================ ${label} ================"
    if [[ -f "${vcf}" ]]; then
        echo "  Ensembl-naming VCF: ${vcf}"
        local n=$(bcftools view -H "${vcf}" | wc -l)
        echo "    Total SNPs: ${n}"
        echo "    Per-chrom (top 25):"
        bcftools view -H "${vcf}" \
            | cut -f1 \
            | sort \
            | uniq -c \
            | sort -k2,2V \
            | head -25 \
            | awk '{printf "      %-6s %s\n", $2, $1}'
    else
        echo "  (Ensembl VCF not found: ${vcf})"
    fi

    if [[ -f "${vcf_ucsc}" ]]; then
        local n=$(bcftools view -H "${vcf_ucsc}" | wc -l)
        echo "  UCSC-naming VCF:    ${vcf_ucsc}  (${n} SNPs)"
    else
        echo "  UCSC-naming VCF:    NOT FOUND — run 03_convert_vcf_to_ucsc.sh"
    fi
}

summary_one "F121-9"
summary_one "BL6xCAST"
