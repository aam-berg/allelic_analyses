#!/usr/bin/env bash
# =============================================================================
# make_wasp_snps.sh — Convert a het-SNP VCF into WASP per-chromosome SNP files
# =============================================================================
#
# WASP's find_intersecting_snps.py expects a directory with one gzipped tab
# file per chromosome:
#
#     {snp_dir}/chr1.snps.txt.gz
#     {snp_dir}/chr2.snps.txt.gz
#     ...
#
# Each file has three tab-separated columns:
#
#     position(1-based)  ref_allele  alt_allele
#
# This script does that conversion in one pass:
#   1. Reads a (compressed) VCF
#   2. Filters to biallelic SNPs (single-base REF, single-base ALT)
#   3. Filters to standard chromosomes (chr1-19, chrX, chrY by default;
#      configurable via --chrom-pattern)
#   4. Splits per chromosome and writes the WASP format
#
# This script is shared between 03_atac/ and 04_proseq/ — both pipelines need
# the same SNP files for WASP.
#
# Usage:
#   bash make_wasp_snps.sh \
#       --vcf /path/to/F121-9_het_snps.ucsc.vcf.gz \
#       --output_dir /path/to/wasp_snp_files
#
# Optional:
#   --chrom-pattern PAT    awk regex for keeping chromosomes (default chr[0-9XY]+)
#   --force                regenerate even if output already exists
# =============================================================================

set -euo pipefail

VCF=""
OUTDIR=""
CHROM_PATTERN='^chr[0-9XY]+$'
FORCE=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)            VCF="$2"; shift 2 ;;
        --output_dir)     OUTDIR="$2"; shift 2 ;;
        --chrom-pattern)  CHROM_PATTERN="$2"; shift 2 ;;
        --force)          FORCE=true; shift ;;
        -h|--help)
            sed -n '4,40p' "$0"
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

if [[ -z "${VCF}" || -z "${OUTDIR}" ]]; then
    echo "Usage: $0 --vcf VCF --output_dir DIR [--chrom-pattern PAT] [--force]" >&2
    exit 1
fi

[[ -f "${VCF}" ]] || { echo "ERROR: VCF not found: ${VCF}" >&2; exit 1; }

mkdir -p "${OUTDIR}"

echo "============================================================"
echo "make_wasp_snps.sh"
echo "============================================================"
echo "  VCF:             ${VCF}"
echo "  Output dir:      ${OUTDIR}"
echo "  Chrom pattern:   ${CHROM_PATTERN}"
echo "  Force regenerate: ${FORCE}"
echo ""

# -----------------------------------------------------------------------------
# Skip if output exists, unless --force
# -----------------------------------------------------------------------------
if [[ "${FORCE}" == false ]] && ls "${OUTDIR}"/chr*.snps.txt.gz 1>/dev/null 2>&1; then
    EXISTING=$(ls "${OUTDIR}"/chr*.snps.txt.gz | wc -l)
    echo "[INFO] Found ${EXISTING} existing chromosome SNP files in ${OUTDIR}."
    echo "[INFO] Skipping. Pass --force to regenerate."
    exit 0
fi

# -----------------------------------------------------------------------------
# Quick sanity check: total het SNPs in VCF
# -----------------------------------------------------------------------------
echo "[INFO] Counting SNPs in VCF (sanity check)..."
TOTAL_SNPS=$(zcat -f "${VCF}" | grep -v -c '^#' || true)
echo "[INFO]   Total non-header records: ${TOTAL_SNPS}"

# -----------------------------------------------------------------------------
# Single-pass extraction: VCF -> temp tab file with chr/pos/ref/alt
# Filters:
#   - Standard chromosomes (regex)
#   - Biallelic SNPs only (single-base REF and single-base ALT)
# VCF is already 1-based, which is what WASP expects.
# -----------------------------------------------------------------------------
echo "[INFO] Extracting biallelic SNPs..."

TMPFILE="${OUTDIR}/_all_snps.tmp"
trap "rm -f ${TMPFILE}" EXIT

zcat -f "${VCF}" \
| grep -v '^#' \
| awk -v OFS='\t' -v pat="${CHROM_PATTERN}" '
    length($4) == 1 && length($5) == 1 && $1 ~ pat {
        print $1, $2, $4, $5
    }
' > "${TMPFILE}"

FILTERED=$(wc -l < "${TMPFILE}")
echo "[INFO]   After filter: ${FILTERED} biallelic SNPs on standard chroms"

# -----------------------------------------------------------------------------
# Split by chromosome, sort by position, gzip
# -----------------------------------------------------------------------------
echo "[INFO] Splitting by chromosome..."

# Get unique chromosomes that actually have SNPs
CHROMS=$(cut -f1 "${TMPFILE}" | sort -u)

for CHROM in ${CHROMS}; do
    OUTFILE="${OUTDIR}/${CHROM}.snps.txt.gz"
    awk -v OFS='\t' -v chrom="${CHROM}" '
        $1 == chrom { print $2, $3, $4 }
    ' "${TMPFILE}" \
    | sort -k1,1n \
    | gzip > "${OUTFILE}"

    N=$(zcat "${OUTFILE}" | wc -l)
    printf "  %-8s %12d SNPs\n" "${CHROM}" "${N}"
done

# -----------------------------------------------------------------------------
# Final sanity check
# -----------------------------------------------------------------------------
echo ""
echo "[INFO] Sanity checks:"
NFILES=$(ls "${OUTDIR}"/chr*.snps.txt.gz | wc -l)
echo "  Chromosome files written:    ${NFILES}"

TOTAL_OUT=$(zcat "${OUTDIR}"/chr*.snps.txt.gz | wc -l)
echo "  Total SNPs in WASP format:   ${TOTAL_OUT}"

if [[ "${TOTAL_OUT}" -ne "${FILTERED}" ]]; then
    echo "  [WARNING] WASP-format total (${TOTAL_OUT}) does not match filtered total (${FILTERED})."
    echo "  This shouldn't happen unless there were duplicate (chrom,pos) entries in the VCF."
fi

# Show a few example entries from chr1 if present
if [[ -f "${OUTDIR}/chr1.snps.txt.gz" ]]; then
    echo ""
    echo "  Example entries (chr1, first 3):"
    zcat "${OUTDIR}/chr1.snps.txt.gz" | head -3 | sed 's/^/    /' || true
fi

echo ""
echo "[DONE] WASP SNP files ready in: ${OUTDIR}"
