#!/bin/bash
# =============================================================================
# 03_convert_vcf_to_ucsc.sh — Convert per-hybrid VCFs from Ensembl → UCSC
#
# Produces ${VCF_DIR}/<HYBRID>_het_snps.ucsc.vcf.gz alongside the original.
# 03_motif_annot/02_annotate_motifs.R prefers .ucsc.vcf.gz when present.
# Run on a login node — fast.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

set +u
source activate "${CONDA_ENV}"
set -u

# Ensembl-to-UCSC mapping for the standard mouse autosomes + sex chroms.
# Mitochondrion: MT (Ensembl) → chrM (UCSC).
MAPFILE="${VCF_DIR}/.ensembl_to_ucsc.txt"
{
    for i in $(seq 1 19); do printf "%d\tchr%d\n" "${i}" "${i}"; done
    printf "X\tchrX\nY\tchrY\nMT\tchrM\n"
} > "${MAPFILE}"

convert_one() {
    local hybrid="$1"
    local in="${VCF_DIR}/${hybrid}_het_snps.vcf.gz"
    local out="${VCF_DIR}/${hybrid}_het_snps.ucsc.vcf.gz"

    if [[ ! -f "${in}" ]]; then
        echo "[${hybrid}] input ${in} not found; skipping."
        return 0
    fi
    if [[ -f "${out}" && -f "${out}.tbi" ]]; then
        echo "[${hybrid}] ${out} exists; skipping."
        return 0
    fi

    echo "[${hybrid}] ${in}  ->  ${out}"
    bcftools annotate \
        --threads "${THREADS}" \
        --rename-chrs "${MAPFILE}" \
        "${in}" \
        -Oz -o "${out}"
    tabix -p vcf "${out}"
    n=$(bcftools view -H "${out}" | wc -l)
    echo "  ${n} SNPs"
}

convert_one "F121-9"
convert_one "BL6xCAST"

echo ""
echo "UCSC-naming VCFs ready:"
ls -la "${VCF_DIR}"/*_het_snps.ucsc.vcf.gz 2>/dev/null || true
