#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 03_count_atac_alleles.sh — Count ATAC ref/alt reads at each motif body
# =============================================================================
#
# For ONE motif's pair_index, produces:
#   {motif_id}_atac_counts.tsv.gz      pair_id × {ref_total, alt_total} +
#                                       per-rep counts
#   {motif_id}_atac_long.tsv.gz        long format: pair_id × rep × allele -> count
#
# COUNTS DEFINITION:
#   "Reads at motif" = number of ATAC fragments overlapping the motif body
#   (motif_start_1b - motif_end_1b inclusive). Strand-agnostic (ATAC is
#   unstranded).
#
# WHY USE ALLELE-SPLIT BAMS DIRECTLY:
#   03_atac/ produces {SAMPLE}_{ref,alt}.bam where each fragment has been
#   classified by its mate-pair SNP evidence. We just need to count fragments
#   overlapping each motif body. This is much faster than re-classifying
#   per-pair.
#
# NOTE ON PER-PAIR vs PER-MOTIF:
#   ATAC counts are a property of the MOTIF (chrom, start, end), not the pair.
#   A motif with K SNPs in/near it will have all K pairs share the same
#   ATAC counts. We compute per-motif and join back to all pairs at the end.
#
# USAGE:
#   bash 03_count_atac_alleles.sh \
#     --pair_index /path/to/pair_index/AC0001_pair_index.tsv.gz \
#     --motif_id AC0001 \
#     --outdir /path/to/atac_counts/
# =============================================================================

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

# ---- CLI ----
PAIR_INDEX=""
MOTIF_ID=""
OUT=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --pair_index) PAIR_INDEX="$2"; shift 2 ;;
        --motif_id)   MOTIF_ID="$2"; shift 2 ;;
        --outdir)     OUT="$2"; shift 2 ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

[[ -n "${PAIR_INDEX}" && -n "${MOTIF_ID}" && -n "${OUT}" ]] || \
    { echo "Usage: $0 --pair_index F --motif_id ID --outdir D" >&2; exit 1; }

mkdir -p "${OUT}"

step_header "ATAC allele counting — motif: ${MOTIF_ID}"

set +u
source activate "${CONDA_ENV_TOOLS}"
set -u

WIDE_OUT="${OUT}/${MOTIF_ID}_atac_counts.tsv.gz"
LONG_OUT="${OUT}/${MOTIF_ID}_atac_long.tsv.gz"
if [[ -f "${WIDE_OUT}" && -f "${LONG_OUT}" ]]; then
    echo "[INFO] Outputs exist; skipping."
    exit 0
fi

# Empty pair index?
N_LINES=$(zcat "${PAIR_INDEX}" | wc -l)
if (( N_LINES <= 1 )); then
    echo "[INFO] Empty pair index; writing empty outputs."
    echo -e "motif_hit_id\tatac_ref_total\tatac_alt_total" | gzip > "${WIDE_OUT}"
    echo -e "motif_hit_id\trep\tallele\tcount" | gzip > "${LONG_OUT}"
    exit 0
fi

# Build a per-motif BED (one row per unique motif_hit_id)
TMPDIR=$(mktemp -d)
trap "rm -rf ${TMPDIR}" EXIT

MOTIF_BED="${TMPDIR}/motifs.bed"
zcat "${PAIR_INDEX}" \
| awk 'BEGIN {FS=OFS="\t"} NR==1 {
    for (i=1; i<=NF; i++) col[$i]=i; next
}
{
    key = $col["motif_hit_id"]
    if (!(key in seen)) {
        seen[key] = 1
        # BED 0-based half-open
        print $col["motif_chrom"], $col["motif_start"]-1, $col["motif_end"], $col["motif_hit_id"]
    }
}' | sort -k1,1 -k2,2n > "${MOTIF_BED}"

N_MOTIFS=$(wc -l < "${MOTIF_BED}")
echo "Unique motifs:" "${N_MOTIFS}"

# Per-rep, per-allele bedtools intersect counts.
# Format: motif_hit_id <tab> rep <tab> allele <tab> count
LONG_TSV="${TMPDIR}/atac_long.tsv"
> "${LONG_TSV}"

for SAMPLE in "${ATAC_SAMPLES[@]}"; do
    for ALLELE in ref alt; do
        BAM="${ATAC_ALLELE_BAM_DIR}/${SAMPLE}_${ALLELE}.bam"
        if [[ ! -f "${BAM}" ]]; then
            echo "[WARN] missing ${BAM}; using zeros for this rep/allele."
            awk -v s="${SAMPLE}" -v a="${ALLELE}" \
                'BEGIN {OFS="\t"} {print $4, s, a, 0}' "${MOTIF_BED}" \
                >> "${LONG_TSV}"
            continue
        fi

        # bedtools coverage: -c reports count of reads in column 4
        # Use samtools view to filter to primary, properly paired (if PE)
        bedtools coverage \
            -a "${MOTIF_BED}" \
            -b "${BAM}" \
            -counts \
        | awk -v s="${SAMPLE}" -v a="${ALLELE}" \
            'BEGIN {OFS="\t"} {print $4, s, a, $5}' \
            >> "${LONG_TSV}"
    done
done

# Long output
( echo -e "motif_hit_id\trep\tallele\tcount"; cat "${LONG_TSV}" ) \
    | gzip > "${LONG_OUT}"

# Wide output: pivot to (motif_hit_id, atac_ref_total, atac_alt_total,
#                        per-rep ref/alt cols)
python3 - "${LONG_TSV}" "${WIDE_OUT}" "${ATAC_SAMPLES[@]}" << 'PYEOF'
import sys, csv, gzip
from collections import defaultdict

long_path = sys.argv[1]
wide_path = sys.argv[2]
samples = sys.argv[3:]

per_motif = defaultdict(lambda: defaultdict(int))
with open(long_path) as f:
    for line in f:
        m, s, a, c = line.rstrip('\n').split('\t')
        per_motif[m][f"{s}_{a}"] += int(c)
        per_motif[m][f"total_{a}"] += int(c)

cols = ["motif_hit_id", "atac_ref_total", "atac_alt_total"]
for s in samples:
    cols.extend([f"{s}_ref", f"{s}_alt"])

with gzip.open(wide_path, 'wt', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(cols)
    for m in sorted(per_motif):
        d = per_motif[m]
        row = [m, d.get("total_ref", 0), d.get("total_alt", 0)]
        for s in samples:
            row.extend([d.get(f"{s}_ref", 0), d.get(f"{s}_alt", 0)])
        w.writerow(row)
PYEOF

echo ""
echo "[OK] Wrote:"
echo "  ${WIDE_OUT}"
echo "  ${LONG_OUT}"
