#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04b_extract_proseq_profile.sh — Per-base PRO-seq profile, ±500 bp around motif
# =============================================================================
#
# For ONE motif's pair_index, produces:
#   {motif_id}_proseq_profile.tsv.gz
#       motif_hit_id, rep, allele, strand_relative_to_motif, position_offset, count
#
# WINDOW:
#   For each motif, profile spans [center - PROFILE_HALF_BP, center + PROFILE_HALF_BP]
#   in transcription direction. position_offset ranges from -PROFILE_HALF_BP
#   to +PROFILE_HALF_BP, with 0 = motif center.
#
# STRAND_RELATIVE_TO_MOTIF:
#   "sense"     : reads on the same strand as the motif (= transcribing strand)
#   "antisense" : reads on opposite strand (background / convergent transcription)
#
# WHY PER-BASE:
#   The downstream ARIMA approach aggregates aligned per-base signal across
#   many motifs to build a high-quality average pause peak. We need single-bp
#   resolution; this is the substrate.
#
# WHY USE BIGWIGS NOT BAMS:
#   PRO-seq bigwigs (after canonical strand swap) give per-base coverage
#   directly. Much faster than re-counting from BAMs at ~1000 positions per
#   motif × ~600 motifs = ~600k extractions.
#
# IF BIGWIGS MISSING:
#   Falls back to BAM-based extraction with samtools depth. Slower but works.
#
# OUTPUT SIZE:
#   Per motif: N_motifs × 2 alleles × 4 reps × 2 strands × 1001 positions
#   For 1000 motif hits: ~16M rows. We compress aggressively and write per-motif.
#
# USAGE:
#   bash 04b_extract_proseq_profile.sh \
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

step_header "PRO-seq profile extraction — motif: ${MOTIF_ID}"

set +u
source activate "${CONDA_ENV_TOOLS}"
set -u

PROF_OUT="${OUT}/${MOTIF_ID}_proseq_profile.tsv.gz"
if [[ -f "${PROF_OUT}" ]]; then
    echo "[INFO] Profile exists; skipping."
    exit 0
fi

N_LINES=$(zcat "${PAIR_INDEX}" | wc -l)
if (( N_LINES <= 1 )); then
    echo -e "motif_hit_id\trep\tallele\tstrand_relative\tposition_offset\tcount" \
        | gzip > "${PROF_OUT}"
    exit 0
fi

TMPDIR=$(mktemp -d)
trap "rm -rf ${TMPDIR}" EXIT

# Per-motif "profile window" BED (one row per motif, ±PROFILE_HALF_BP from center)
PROFILE_BED="${TMPDIR}/profile.bed"
zcat "${PAIR_INDEX}" \
| awk -v ph="${PROFILE_HALF_BP}" \
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

    pstart = center - ph - 1
    pend   = center + ph
    if (pstart < 0) pstart = 0
    print chrom, pstart, pend, key, center, strand
}' > "${PROFILE_BED}"

N_MOTIFS=$(wc -l < "${PROFILE_BED}")
echo "Unique motifs:" "${N_MOTIFS}"

# Helper: for a given strand-specific bigwig, extract per-motif per-base coverage.
# Args: $1=bigwig, $2=output_tsv (motif_hit_id\tposition_offset\tcount),
#       $3=motif_strand_filter ("+" or "-")
# Algorithm:
#   bigWigToBedGraph -chrom=... -start=... -end=... too slow; instead use
#   bigWigAverageOverBed with -bedOut to get per-base. We use bigwig_avg as a
#   pragmatic alternative: bedtools genomecov is BAM-only.
# Cleanest path: use deeptools `multiBigwigSummary BED-file` per strand — but
# that gives summary stats, not per-base. So we use `pyBigWig` via Python.

PROFILE_OUT_TSV="${TMPDIR}/profile.tsv"
> "${PROFILE_OUT_TSV}"

# Try bigwig path first
USE_BIGWIGS=1
for SAMPLE in "${PROSEQ_SAMPLES[@]}"; do
    for ALLELE in ref alt; do
        for STRAND in plus minus; do
            BW="${PROSEQ_ALLELE_BIGWIG_DIR}/${SAMPLE}_${ALLELE}_${STRAND}.bw"
            if [[ ! -f "${BW}" ]]; then
                USE_BIGWIGS=0
                break 3
            fi
        done
    done
done

if (( USE_BIGWIGS )); then
    echo "[INFO] Using bigWigs for fast per-base extraction (pyBigWig)."

    python3 - "${PROFILE_BED}" "${PROFILE_OUT_TSV}" \
        "${PROSEQ_ALLELE_BIGWIG_DIR}" "${PROFILE_HALF_BP}" \
        "${PROSEQ_SAMPLES[@]}" << 'PYEOF'
import sys, os
try:
    import pyBigWig
except ImportError:
    sys.stderr.write("[ERROR] pyBigWig not available. Install:\n")
    sys.stderr.write("  conda install -n proseq2_env -c bioconda pyBigWig\n")
    sys.exit(2)

bed_path, out_path, bw_dir, profile_half = sys.argv[1:5]
samples = sys.argv[5:]
profile_half = int(profile_half)

bw_handles = {}
for s in samples:
    for a in ("ref", "alt"):
        for strand_label in ("plus", "minus"):
            p = os.path.join(bw_dir, f"{s}_{a}_{strand_label}.bw")
            bw_handles[(s, a, strand_label)] = pyBigWig.open(p)

# Mapping:
#   "+" motif: sense    = "plus"  bw, antisense = "minus" bw
#   "-" motif: sense    = "minus" bw, antisense = "plus"  bw

n_motifs = 0
with open(bed_path) as bed_in, open(out_path, "w") as out:
    for line in bed_in:
        chrom, bstart, bend, motif_hit_id, center, mstrand = \
            line.rstrip('\n').split('\t')
        bstart = int(bstart); bend = int(bend); center = int(center)
        n_motifs += 1

        # Build position offset array. BED is 0-based half-open.
        # We expect bend - bstart == 2*profile_half + 1 for full window.
        win_len = bend - bstart
        # Expected position offsets (length win_len): center - bstart -> 0
        # Offset for position i (0-indexed in window): (bstart + i + 1) - center
        # NB: position genomic is bstart+i (0-based) -> bstart+i+1 (1-based)
        # We use 1-based genomic for offset math.

        for s in samples:
            for a in ("ref", "alt"):
                for sense_label in ("sense", "antisense"):
                    if mstrand == "+":
                        bw_strand = "plus" if sense_label == "sense" else "minus"
                    else:
                        bw_strand = "minus" if sense_label == "sense" else "plus"
                    bw = bw_handles[(s, a, bw_strand)]

                    try:
                        vals = bw.values(chrom, bstart, bend)
                    except (RuntimeError, OSError) as e:
                        # Chromosome not in bigwig; skip
                        continue
                    if vals is None:
                        continue

                    # For "-" motif, transcription direction is genomic-decreasing.
                    # Reverse the values so position_offset is in transcription order.
                    if mstrand == "-":
                        vals = vals[::-1]

                    # Position offsets: for sense strand 0..win_len-1
                    # For "+" motif: genomic_pos[i] = bstart + i + 1
                    #   offset = genomic_pos[i] - center = (bstart + i + 1) - center
                    # For "-" motif (after reversal):
                    #   tx_pos[i] = center + (i - profile_half) (signed)
                    # In both cases offsets range from approx -profile_half to +profile_half.
                    # Compute by uniform formula on reversed array:
                    base = -(center - bstart - 1) if mstrand == "+" else -profile_half
                    # Simpler: compute offset directly:
                    if mstrand == "+":
                        offsets = [(bstart + i + 1) - center for i in range(win_len)]
                    else:
                        # Walk genome left-to-right; tx direction is right-to-left.
                        # After reversing vals, vals[0] corresponds to genomic position bend.
                        offsets = [center - (bend - i) for i in range(win_len)]

                    for off, v in zip(offsets, vals):
                        if v is None or v == 0 or (isinstance(v, float) and v != v):
                            continue   # skip NaN/zero
                        out.write(f"{motif_hit_id}\t{s}\t{a}\t{sense_label}\t{off}\t{int(v) if v == int(v) else v:g}\n")

for h in bw_handles.values():
    h.close()
sys.stderr.write(f"[INFO] Profiled {n_motifs} motifs.\n")
PYEOF

else
    echo "[WARN] bigWigs not all present; falling back to BAM-based extraction."
    echo "       This is significantly slower."

    # BAM fallback — use samtools depth across the motif window per motif.
    # Strand assignment: read flag 16 -> "-" mapping, no flag 16 -> "+" mapping.
    python3 - "${PROFILE_BED}" "${PROFILE_OUT_TSV}" \
        "${PROSEQ_ALLELE_BAM_DIR}" "${PROFILE_HALF_BP}" \
        "${PROSEQ_SAMPLES[@]}" << 'PYEOF'
import sys, os, subprocess
from collections import defaultdict

bed_path, out_path, bam_dir, profile_half = sys.argv[1:5]
samples = sys.argv[5:]
profile_half = int(profile_half)

# We'll use samtools view to stream alignments by region, then count by
# strand and position. For each motif we get per-base counts on each strand.

# Detect PE
detect_bam = os.path.join(bam_dir, f"{samples[0]}_ref.bam")
is_pe = False
try:
    out = subprocess.check_output(
        f"samtools view {detect_bam} 2>/dev/null | head -200 | awk 'NR<200' || true",
        shell=True, text=True)
    for line in out.splitlines():
        flag = int(line.split('\t')[1])
        if flag & 1:
            is_pe = True
            break
except subprocess.CalledProcessError:
    pass

with open(bed_path) as bed_in, open(out_path, 'w') as out:
    for line in bed_in:
        chrom, bstart, bend, motif_hit_id, center, mstrand = \
            line.rstrip('\n').split('\t')
        bstart = int(bstart); bend = int(bend); center = int(center)
        win_len = bend - bstart

        region = f"{chrom}:{bstart+1}-{bend}"

        for s in samples:
            for a in ("ref", "alt"):
                bam = os.path.join(bam_dir, f"{s}_{a}.bam")
                if not os.path.exists(bam):
                    continue

                # Counts per position (1-based) per genomic strand
                pos_plus = defaultdict(int)
                pos_minus = defaultdict(int)

                # We use 5'-end of each read as the position
                # (PRO-seq convention; the canonical strand-swap has
                # already been applied during 04_proseq/'s BAM construction)
                view_filter = "-F 4"
                if is_pe:
                    view_filter += " -f 64"   # R1 only
                cmd = f"samtools view {view_filter} {bam} {region}"
                try:
                    proc = subprocess.Popen(cmd, shell=True,
                                             stdout=subprocess.PIPE, text=True)
                except Exception:
                    continue

                for sam_line in proc.stdout:
                    fields = sam_line.split('\t')
                    flag = int(fields[1])
                    pos = int(fields[3])  # 1-based leftmost
                    if flag & 16:
                        # reverse strand: 5' end is at right end of alignment
                        # cigar parsing for length:
                        cigar = fields[5]
                        length = 0
                        n = 0
                        for ch in cigar:
                            if ch.isdigit():
                                n = n*10 + int(ch)
                            else:
                                if ch in "MDN=X":
                                    length += n
                                n = 0
                        end_pos = pos + length - 1
                        if bstart < end_pos <= bend:
                            pos_minus[end_pos] += 1
                    else:
                        if bstart < pos <= bend:
                            pos_plus[pos] += 1
                proc.wait()

                # Emit per-position counts for sense and antisense strands
                if mstrand == "+":
                    sense_d, antisense_d = pos_plus, pos_minus
                else:
                    sense_d, antisense_d = pos_minus, pos_plus

                for sense_label, d in (("sense", sense_d),
                                        ("antisense", antisense_d)):
                    for genomic_pos, cnt in d.items():
                        if mstrand == "+":
                            offset = genomic_pos - center
                        else:
                            offset = center - genomic_pos
                        out.write(f"{motif_hit_id}\t{s}\t{a}\t{sense_label}\t{offset}\t{cnt}\n")
PYEOF
fi

# Add header and gzip
( echo -e "motif_hit_id\trep\tallele\tstrand_relative\tposition_offset\tcount"
  cat "${PROFILE_OUT_TSV}" ) | gzip > "${PROF_OUT}"

NLINES=$(zcat "${PROF_OUT}" | wc -l)
echo ""
echo "[OK] Wrote ${PROF_OUT} (${NLINES} lines)"
