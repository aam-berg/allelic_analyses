#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 05_extract_rnaseq_junctions.sh — Per-allele junction usage at the 8 nearest
# splice sites for each motif (from 02g catalog)
# =============================================================================
#
# For ONE motif's pair_index, produces:
#   {motif_id}_rnaseq_junctions.tsv.gz
#       motif_hit_id, junction_cell, site_chrom, site_pos, site_strand,
#       site_type, source, motif_type, distance,
#       <per-rep>_use_ref, <per-rep>_use_alt, <per-rep>_skip_ref,
#       <per-rep>_skip_alt
#       (and totals across reps)
#
# The "8 directional cells" come from the matching motif annotation TSV at
# {MOTIF_ANNOT_DIR}/{motif_id}_annotated.tsv.gz, which has columns like:
#   donor_sense_upstream_dist, donor_sense_upstream_source, ...
# We use the dist + the matching site geometry to reconstruct the
# (chrom, pos, strand) of each cell's nearest site, look up that site's
# splice motif type from the donor/acceptor catalog (in preprocessed/), then
# count reads per allele at that site.
#
# READ CLASSIFICATION:
#   - "use" reads: spliced reads with one mate-end aligning across the site
#     (i.e. the read CIGAR has an N skip with one end at the site position)
#   - "skip" reads: continuous (non-spliced) reads spanning the site
#
# This script delegates the heavy lifting to a Python helper (in-line) that
# uses pysam to walk reads in each region.
#
# USAGE:
#   bash 05_extract_rnaseq_junctions.sh \
#     --pair_index ... --motif_id AC0001 \
#     --motif_annot_tsv .../AC0001_annotated.tsv.gz \
#     --outdir ...
# =============================================================================

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

PAIR_INDEX=""
MOTIF_ID=""
MOTIF_ANNOT_TSV=""
OUT=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --pair_index)        PAIR_INDEX="$2"; shift 2 ;;
        --motif_id)          MOTIF_ID="$2"; shift 2 ;;
        --motif_annot_tsv)   MOTIF_ANNOT_TSV="$2"; shift 2 ;;
        --outdir)            OUT="$2"; shift 2 ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

[[ -n "${PAIR_INDEX}" && -n "${MOTIF_ID}" && \
   -n "${MOTIF_ANNOT_TSV}" && -n "${OUT}" ]] || \
    { echo "Usage: $0 --pair_index F --motif_id ID --motif_annot_tsv F --outdir D" >&2; exit 1; }

mkdir -p "${OUT}"

step_header "RNA-seq junction extraction — motif: ${MOTIF_ID}"

set +u
source activate "${CONDA_ENV_TOOLS}"
set -u

JUNC_OUT="${OUT}/${MOTIF_ID}_rnaseq_junctions.tsv.gz"
if [[ -f "${JUNC_OUT}" ]]; then
    echo "[INFO] Output exists; skipping."
    exit 0
fi

N_LINES=$(zcat "${PAIR_INDEX}" | wc -l)
if (( N_LINES <= 1 )); then
    echo -e "motif_hit_id\tjunction_cell\tsite_chrom\tsite_pos\tsite_strand\tsite_type\tsource\tdistance" \
        | gzip > "${JUNC_OUT}"
    exit 0
fi

if [[ ! -f "${MOTIF_ANNOT_TSV}" ]]; then
    echo "[WARN] Motif annot TSV not found: ${MOTIF_ANNOT_TSV}"
    echo "       Cannot extract junction sites; writing empty file."
    echo -e "motif_hit_id\tjunction_cell\tsite_chrom\tsite_pos\tsite_strand\tsite_type\tsource\tdistance" \
        | gzip > "${JUNC_OUT}"
    exit 0
fi

# All the work happens in Python — we need pysam for read-level CIGAR parsing.
TMPDIR=$(mktemp -d)
trap "rm -rf ${TMPDIR}" EXIT

python3 - \
    "${PAIR_INDEX}" "${MOTIF_ANNOT_TSV}" "${RNASEQ_ALLELE_BAM_DIR}" \
    "${SJ_OVERLAP_BP}" "${PSI_MIN_TOTAL_READS}" "${JUNC_OUT}" \
    "${RNASEQ_SAMPLES[@]}" << 'PYEOF'
import sys
import os
import gzip
import csv
from collections import defaultdict

try:
    import pysam
except ImportError:
    sys.stderr.write("[ERROR] pysam required. conda install -c bioconda pysam\n")
    sys.exit(2)

(pair_index_path, annot_tsv_path, bam_dir,
 sj_overlap_bp, psi_min_total_reads, out_path, *samples) = sys.argv[1:]
sj_overlap_bp = int(sj_overlap_bp)
psi_min_total_reads = int(psi_min_total_reads)

# -----------------------------------------------------------------------------
# 1. Read pair index -> set of (motif_hit_id) and motif geometry
# -----------------------------------------------------------------------------
motif_geom = {}   # motif_hit_id -> (chrom, start, end, strand, center)
with gzip.open(pair_index_path, "rt") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        mhid = row["motif_hit_id"]
        if mhid in motif_geom:
            continue
        s = int(row["motif_start"])
        e = int(row["motif_end"])
        center = (s + e) // 2
        motif_geom[mhid] = (row["motif_chrom"], s, e,
                              row["motif_strand"], center)

if not motif_geom:
    sys.stderr.write("[INFO] No motifs in pair index.\n")
    sys.exit(0)

sys.stderr.write(f"Motifs to process: {len(motif_geom):,}\n")

# -----------------------------------------------------------------------------
# 2. Read motif annotation TSV → for each motif, the 8 junction cells
#    (donor/acceptor × sense/antisense × upstream/downstream)
# -----------------------------------------------------------------------------
JUNCTION_CELLS = []
for site_type in ("donor", "acceptor"):
    for orient in ("sense", "antisense"):
        for direction in ("upstream", "downstream"):
            JUNCTION_CELLS.append((site_type, orient, direction))

# We also need the motif annot TSV's BED-derived motif_hit_id to match.
# The annot TSV has chrom/start/end/motif_id/score/strand columns.
# motif_hit_id is built as: motif_id__chrom__start1b__strand
# In the annot TSV, start is BED 0-based -> start1b = start + 1.

motif_to_cells = {}   # motif_hit_id -> list of dicts (one per cell)
with gzip.open(annot_tsv_path, "rt") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        chrom = row["chrom"]
        start_1b = int(row["start"]) + 1
        end_1b = int(row["end"])
        motif_id = row["motif_id"]
        # Some motif_ids in BED have form "AC0001" or "AC0001:DLX/LHX:..." etc;
        # use the first colon-separated token to match pair_index conventions.
        motif_id_short = motif_id.split(":")[0]
        strand = row["strand"]
        mhid = f"{motif_id_short}__{chrom}__{start_1b}__{strand}"
        if mhid not in motif_geom:
            continue

        cells = []
        for st, orient, direction in JUNCTION_CELLS:
            base = f"{st}_{orient}_{direction}"
            dist_col = f"{base}_dist"
            src_col = f"{base}_source"
            mtype_col = f"{base}_motif_type"
            if dist_col not in row or row[dist_col] in ("NA", ""):
                continue
            try:
                dist = int(row[dist_col])
            except ValueError:
                continue

            # Reconstruct site (chrom, pos, strand)
            mchrom, ms, me, mstrand, center = motif_geom[mhid]
            # site_strand: sense -> motif_strand; antisense -> opposite
            site_strand = mstrand if orient == "sense" else \
                ("-" if mstrand == "+" else "+")
            # Direction in transcription:
            # If "+" motif and "upstream" -> site is at LOWER coord
            # If "+" motif and "downstream" -> site at HIGHER coord
            # If "-" motif: opposite
            if (mstrand == "+" and direction == "upstream") or \
               (mstrand == "-" and direction == "downstream"):
                site_pos = center - dist
            else:
                site_pos = center + dist

            cells.append({
                "junction_cell": base,
                "site_chrom": mchrom,
                "site_pos": site_pos,
                "site_strand": site_strand,
                "site_type": st,
                "source": row.get(src_col, ""),
                "motif_type": row.get(mtype_col, ""),
                "distance": dist,
            })
        motif_to_cells[mhid] = cells

n_with_cells = sum(1 for v in motif_to_cells.values() if v)
sys.stderr.write(
    f"Motifs with >=1 junction cell annotated: {n_with_cells:,}\n")

# -----------------------------------------------------------------------------
# 3. For each motif × cell, walk reads in [site - sj_overlap_bp, site + sj_overlap_bp]
#    in each (sample, allele) BAM and classify use vs skip.
# -----------------------------------------------------------------------------
def classify_read(read, site_pos):
    """
    Returns "use" if the read has an N skip with one boundary at site_pos
    (within ±1 bp tolerance for STAR junction position conventions).
    Returns "skip" if the read is continuous and spans site_pos (i.e., the
    site falls strictly inside the aligned reference span without an N).
    Returns None otherwise.
    """
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return None

    # Walk CIGAR and reference positions to find the introns (N operations).
    # site_pos is 1-based; pysam uses 0-based internally.
    ref_pos = read.reference_start  # 0-based
    aligned_segments = []
    cur_start = ref_pos
    n_introns = []   # list of (intron_start_0b, intron_end_0b) where end is exclusive
    for op, length in read.cigartuples or []:
        if op in (0, 7, 8):   # M, =, X
            ref_pos += length
        elif op == 2:         # D
            ref_pos += length
        elif op == 3:         # N (intron)
            intron_start = ref_pos
            ref_pos += length
            n_introns.append((intron_start, ref_pos))
        elif op == 1 or op == 4:  # I, S
            pass
        elif op == 5:         # H
            pass
        # else: P, etc. ignored

    # Reference span of the alignment
    aln_end = ref_pos
    aln_start = read.reference_start

    # Convert site_pos (1-based) to 0-based half-open reference comparisons:
    # intron_start (0-based) is the FIRST base of the intron;
    # in 1-based, intron starts at intron_start_0b + 1.
    # We treat the site as 1-based; intron-starts and intron-ends are computed
    # in 0-based for direct comparison.

    # Use: read has an N skip with one boundary == site_pos
    site_0b = site_pos - 1
    for (i_s, i_e) in n_introns:
        # i_s = first base of intron (0-based) -> the donor side (in genome coords)
        #   in 1-based: position before intron is i_s (last exonic base);
        #   first intronic base is i_s + 1 (1-based).
        # i_e = first base AFTER intron (0-based) -> the acceptor side.
        #   in 1-based: last intronic base is i_e (1-based equiv = i_e + 1 - 1 = i_e);
        #   first post-intron base is i_e + 1 (1-based).
        # STAR's SJ.out.tab convention: intron_start (column 2) and intron_end
        # (column 3) are both 1-based and refer to the FIRST and LAST base of
        # the intron respectively. So:
        #   donor (genomic strand "+"): intron_start (1-based) = i_s + 1
        #   acceptor (genomic strand "+"): intron_end (1-based) = i_e (since i_e is 0-based exclusive,
        #     last intronic base in 1-based is i_e)
        # So we accept the site_pos if it equals (i_s + 1) [donor-on-+strand]
        # or i_e [acceptor-on-+strand]. We don't need to know donor vs acceptor
        # here — both boundaries qualify as "use".
        if site_pos == (i_s + 1) or site_pos == i_e:
            return "use"

    # Skip: read aligns continuously through site (no N spanning site_pos)
    # i.e., site_pos is in [aln_start+1, aln_end] AND no intron contains site_pos.
    if (aln_start + 1) <= site_pos <= aln_end:
        for (i_s, i_e) in n_introns:
            # Intron in 1-based: first intron base = i_s+1, last = i_e
            # site falls in intron if (i_s+1) <= site_pos <= i_e
            if (i_s + 1) <= site_pos <= i_e:
                return None  # site is in intron → not "skip", not "use"
        return "skip"

    return None


# Open BAMs once
bam_handles = {}
for s in samples:
    for a in ("ref", "alt"):
        p = os.path.join(bam_dir, f"{s}_{a}.bam")
        if os.path.exists(p):
            try:
                bam_handles[(s, a)] = pysam.AlignmentFile(p, "rb")
            except Exception as e:
                sys.stderr.write(f"[WARN] cannot open {p}: {e}\n")

with gzip.open(out_path, "wt") as out_f:
    out_writer = csv.writer(out_f, delimiter='\t')
    header = ["motif_hit_id", "junction_cell", "site_chrom", "site_pos",
              "site_strand", "site_type", "source", "motif_type", "distance",
              "use_ref_total", "use_alt_total", "skip_ref_total", "skip_alt_total"]
    for s in samples:
        for a in ("ref", "alt"):
            for k in ("use", "skip"):
                header.append(f"{s}_{a}_{k}")
    out_writer.writerow(header)

    for mhid, cells in motif_to_cells.items():
        for cell in cells:
            counts = defaultdict(int)   # (sample, allele, kind) -> count
            for s in samples:
                for a in ("ref", "alt"):
                    bam = bam_handles.get((s, a))
                    if bam is None:
                        continue
                    region_start = max(0, cell["site_pos"] - sj_overlap_bp - 1)
                    region_end = cell["site_pos"] + sj_overlap_bp
                    try:
                        iter_reads = bam.fetch(cell["site_chrom"],
                                                region_start, region_end)
                    except (ValueError, KeyError):
                        continue
                    for read in iter_reads:
                        kind = classify_read(read, cell["site_pos"])
                        if kind is not None:
                            counts[(s, a, kind)] += 1

            row = [
                mhid, cell["junction_cell"], cell["site_chrom"],
                cell["site_pos"], cell["site_strand"], cell["site_type"],
                cell["source"], cell["motif_type"], cell["distance"],
                sum(counts[(s, "ref", "use")] for s in samples),
                sum(counts[(s, "alt", "use")] for s in samples),
                sum(counts[(s, "ref", "skip")] for s in samples),
                sum(counts[(s, "alt", "skip")] for s in samples),
            ]
            for s in samples:
                for a in ("ref", "alt"):
                    for k in ("use", "skip"):
                        row.append(counts[(s, a, k)])
            out_writer.writerow(row)

for h in bam_handles.values():
    h.close()

sys.stderr.write(f"[OK] Wrote {out_path}\n")
PYEOF

echo "[OK] Wrote ${JUNC_OUT}"
