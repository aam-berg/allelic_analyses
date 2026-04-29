#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04_count_proseq_alleles.sh — Count PRO-seq ref/alt reads in motif pause
# window and surrounding gene body, strand-aware
# =============================================================================
#
# For ONE motif's pair_index, produces:
#   {motif_id}_proseq_window_counts.tsv.gz    wide: per-motif pause + gene
#                                              body counts ref/alt with
#                                              per-rep breakdown
#   {motif_id}_proseq_window_long.tsv.gz      long: motif × rep × allele ×
#                                              window_type -> count
#
# WINDOWS (strand-aware, in TRANSCRIPTION direction):
#   PAUSE_WINDOW   = motif_center ± PAUSE_WINDOW_HALF_BP
#   GENE_BODY      = downstream of motif:
#                       "+" motif: [motif_end + 1, motif_end + GENE_BODY_FLANK_BP]
#                       "-" motif: [motif_start - GENE_BODY_FLANK_BP, motif_start - 1]
#
# READ STRAND CONVENTION:
#   PRO-seq aligns to the antisense strand of nascent RNA. The canonical
#   strand-swap (defined in lib/strand_aware_count.py and applied during
#   04_proseq/) means BAM read strand corresponds to the transcribing strand
#   in our convention. So:
#     "+" motif: count reads on "+" strand in the motif window
#     "-" motif: count reads on "-" strand in the motif window
#
#   We use samtools view with -F 16 / -f 16 flags to filter:
#     -F 16: read NOT reverse-complement -> "+" strand
#     -f 16: read IS reverse-complement -> "-" strand
#   PE handling (if PE): we count R1 only to avoid double-counting fragments,
#   via -f 64 (first in pair). For SE, no flag is needed.
#   We auto-detect SE vs PE by inspecting the first BAM.
#
# WHY USE BAMS NOT BIGWIGS:
#   For window counts (pause + gene body), BAMs let us compute exact integer
#   counts. BigWigs are needed for the per-base profile (next script).
#
# USAGE:
#   bash 04_count_proseq_alleles.sh \
#     --pair_index ... --motif_id AC0001 --outdir ...
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

step_header "PRO-seq window counting — motif: ${MOTIF_ID}"

set +u
source activate "${CONDA_ENV_TOOLS}"
set -u

WIDE_OUT="${OUT}/${MOTIF_ID}_proseq_window_counts.tsv.gz"
LONG_OUT="${OUT}/${MOTIF_ID}_proseq_window_long.tsv.gz"
if [[ -f "${WIDE_OUT}" && -f "${LONG_OUT}" ]]; then
    echo "[INFO] Outputs exist; skipping."
    exit 0
fi

# Empty pair index?
N_LINES=$(zcat "${PAIR_INDEX}" | wc -l)
if (( N_LINES <= 1 )); then
    echo "[INFO] Empty pair index; writing empty outputs."
    echo -e "motif_hit_id\tproseq_pause_ref\tproseq_pause_alt\tproseq_genebody_ref\tproseq_genebody_alt" | gzip > "${WIDE_OUT}"
    echo -e "motif_hit_id\trep\tallele\twindow_type\tcount" | gzip > "${LONG_OUT}"
    exit 0
fi

TMPDIR=$(mktemp -d)
trap "rm -rf ${TMPDIR}" EXIT

# Build per-motif BEDs for pause window and gene body, strand-aware.
# Output: 6-column BED with motif_hit_id in column 4 and motif_strand in
# column 6. We'll filter reads by strand using samtools after intersect.

PAUSE_BED="${TMPDIR}/pause.bed"
GENEBODY_BED="${TMPDIR}/genebody.bed"

zcat "${PAIR_INDEX}" \
| awk -v ph="${PAUSE_WINDOW_HALF_BP}" -v gf="${GENE_BODY_FLANK_BP}" \
'BEGIN {FS=OFS="\t"} NR==1 {
    for (i=1; i<=NF; i++) col[$i]=i; next
}
{
    key = $col["motif_hit_id"]
    if (key in seen) next
    seen[key] = 1
    chrom  = $col["motif_chrom"]
    s      = $col["motif_start"]
    e      = $col["motif_end"]
    strand = $col["motif_strand"]
    center = int((s + e) / 2)

    # Pause window (BED 0-based half-open)
    pstart = center - ph - 1
    pend   = center + ph
    if (pstart < 0) pstart = 0
    print chrom, pstart, pend, key, ".", strand >> "'"${PAUSE_BED}"'"

    # Gene body (BED 0-based half-open, in transcription direction)
    if (strand == "+") {
        gstart = e
        gend   = e + gf
    } else {
        gstart = s - 1 - gf
        if (gstart < 0) gstart = 0
        gend   = s - 1
    }
    if (gend > gstart) {
        print chrom, gstart, gend, key, ".", strand >> "'"${GENEBODY_BED}"'"
    }
}'

N_MOTIFS=$(wc -l < "${PAUSE_BED}")
echo "Unique motifs:" "${N_MOTIFS}"

# Detect PE-vs-SE from first BAM
DETECT_BAM="${PROSEQ_ALLELE_BAM_DIR}/${PROSEQ_SAMPLES[0]}_ref.bam"
IS_PE=0
if [[ -f "${DETECT_BAM}" ]]; then
    # Sample first 100 reads; if any has flag bit 1 set (paired), call PE
    PE_FLAG=$(samtools view -h "${DETECT_BAM}" 2>/dev/null \
              | head -200 \
              | awk '!/^@/ {if (and($2, 1)) print "1"; else print "0"}' \
              | sort -u | head -1)
    [[ "${PE_FLAG}" == "1" ]] && IS_PE=1
fi
echo "PE detection: IS_PE=${IS_PE}"

count_strand() {
    # $1 = BAM, $2 = BED, $3 = strand_filter ("+" or "-")
    # Output: motif_hit_id<tab>count for each motif on that strand.
    local BAM="$1" BED="$2" S="$3"
    # Filter BED by strand
    local SUB_BED="${TMPDIR}/_sub.bed"
    awk -v s="${S}" '$6==s' "${BED}" > "${SUB_BED}"
    if [[ ! -s "${SUB_BED}" ]]; then
        return 0
    fi

    local FLAG_INC=""
    local FLAG_EXC=""
    if [[ "${S}" == "+" ]]; then
        FLAG_EXC="-F 16"
    else
        FLAG_INC="-f 16"
    fi
    if (( IS_PE )); then
        # For PE: -f 64 = R1 only; combine with strand bit
        if [[ "${S}" == "+" ]]; then
            FLAG_INC="-f 64"      # require R1
            FLAG_EXC="-F 16"      # exclude reverse
        else
            FLAG_INC="-f 80"      # require R1 (64) + reverse (16)
            FLAG_EXC=""
        fi
    fi

    samtools view -bh ${FLAG_INC} ${FLAG_EXC} "${BAM}" \
    | bedtools coverage -a "${SUB_BED}" -b - -counts \
    | awk 'BEGIN {OFS="\t"} {print $4, $7}'
}

LONG_TSV="${TMPDIR}/proseq_long.tsv"
> "${LONG_TSV}"

for SAMPLE in "${PROSEQ_SAMPLES[@]}"; do
    for ALLELE in ref alt; do
        BAM="${PROSEQ_ALLELE_BAM_DIR}/${SAMPLE}_${ALLELE}.bam"
        if [[ ! -f "${BAM}" ]]; then
            echo "[WARN] missing ${BAM}; using zeros."
            for WIN_BED in "pause:${PAUSE_BED}" "genebody:${GENEBODY_BED}"; do
                WT="${WIN_BED%%:*}"
                B="${WIN_BED##*:}"
                awk -v s="${SAMPLE}" -v a="${ALLELE}" -v wt="${WT}" \
                    'BEGIN {OFS="\t"} {print $4, s, a, wt, 0}' "${B}" \
                    >> "${LONG_TSV}"
            done
            continue
        fi

        for WIN in pause genebody; do
            if [[ "${WIN}" == "pause" ]]; then
                BED="${PAUSE_BED}"
            else
                BED="${GENEBODY_BED}"
            fi

            # Combine + and - strand counts (each motif is on exactly one strand)
            { count_strand "${BAM}" "${BED}" "+"; \
              count_strand "${BAM}" "${BED}" "-"; } \
            | awk -v s="${SAMPLE}" -v a="${ALLELE}" -v wt="${WIN}" \
                'BEGIN {OFS="\t"} {print $1, s, a, wt, $2}' \
            >> "${LONG_TSV}"
        done
    done
done

# Long output
( echo -e "motif_hit_id\trep\tallele\twindow_type\tcount"; cat "${LONG_TSV}" ) \
    | gzip > "${LONG_OUT}"

# Wide output
python3 - "${LONG_TSV}" "${WIDE_OUT}" "${PROSEQ_SAMPLES[@]}" << 'PYEOF'
import sys, csv, gzip
from collections import defaultdict

long_path = sys.argv[1]
wide_path = sys.argv[2]
samples = sys.argv[3:]

per = defaultdict(lambda: defaultdict(int))
with open(long_path) as f:
    for line in f:
        m, s, a, wt, c = line.rstrip('\n').split('\t')
        c = int(c)
        per[m][f"{s}_{a}_{wt}"] += c
        per[m][f"total_{a}_{wt}"] += c

cols = ["motif_hit_id",
        "proseq_pause_ref", "proseq_pause_alt",
        "proseq_genebody_ref", "proseq_genebody_alt"]
for s in samples:
    for a in ("ref", "alt"):
        for wt in ("pause", "genebody"):
            cols.append(f"{s}_{a}_{wt}")

with gzip.open(wide_path, 'wt', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(cols)
    for m in sorted(per):
        d = per[m]
        row = [m,
               d.get("total_ref_pause", 0), d.get("total_alt_pause", 0),
               d.get("total_ref_genebody", 0), d.get("total_alt_genebody", 0)]
        for s in samples:
            for a in ("ref", "alt"):
                for wt in ("pause", "genebody"):
                    row.append(d.get(f"{s}_{a}_{wt}", 0))
        w.writerow(row)
PYEOF

echo ""
echo "[OK] Wrote:"
echo "  ${WIDE_OUT}"
echo "  ${LONG_OUT}"
