#!/usr/bin/env python3
"""
strand_aware_count.py — Count PRO-seq Pol II positions inside given intervals,
                       respecting the PRO-seq strand swap.

This is the canonical implementation of the PRO-seq strand swap for
motif-centered counting. Whole-genome bigWig generation lives in
04_proseq/04_make_bigwigs.sh; this script handles the case of many small
intervals where bedtools genomecov + window-overlap would be awkward.

THE STRAND SWAP (canonical, see scripts/README.md section 7):

  PRO-seq reads come from cDNA, which is the reverse complement of the
  nascent RNA. Therefore:
    Read on BAM "+" strand (FLAG 16 NOT set)  <=>  RNA was on "-" strand
    Read on BAM "-" strand (FLAG 16 set)      <=>  RNA was on "+" strand

  The 5' end of the read corresponds to the 3' end of the RNA, which is
  where RNAPII was at run-on time (the active site). For pause analysis,
  THIS is the position that matters.

So, to count "Pol II positions where the RNA is on the + strand of the
genome inside interval [s, e]":
    1. Filter to reads with FLAG 16 SET (BAM minus strand)
    2. Take the 5' end of each read (reference_start, since it's already
       in reference orientation)
    3. Count those positions in [s, e].

Conversely for "RNA on the - strand of the genome":
    1. Filter to reads with FLAG 16 NOT set (BAM plus strand)
    2. Take the 5' end of each read in reference orientation, which is
       reference_end - 1 because pysam returns reference_start as the
       leftmost reference position regardless of read strand.

USAGE — CLI:

  Count over a BED-like intervals file. Each interval has a strand column
  that, by default, is interpreted as the RNA strand we want signal from:

    python strand_aware_count.py \\
        --bam ${SAMPLE}_final.bam \\
        --intervals motif_windows.bed \\
        --output counts.tsv

  The intervals file is expected as 6-column BED:
      chrom  start  end  name  score  strand
  where 'strand' is the RNA strand whose signal to count at this interval.
  '+' / '-' are honored; '*' or '.' counts both strands combined.

  Strand interpretation can be flipped with --strand-mode bam, which
  treats the strand column as a BAM strand (i.e., no swap applied).

USAGE — library:

  from strand_aware_count import count_at_intervals
  counts = count_at_intervals(bam_path, intervals_iter, mode="rna")
  # 'intervals_iter' yields (chrom, start, end, name, strand) tuples
  # 'counts' is a dict mapping interval index -> count

DEPENDENCIES:
  pysam
"""

import argparse
import sys

import pysam


# ---------------------------------------------------------------------------
# Core counting logic
# ---------------------------------------------------------------------------

def pol2_position_for_read(read):
    """
    Return the (chrom, position_0based, rna_strand) of the Pol II active site
    inferred from this read, or None if the read should be skipped.

    The Pol II site is the 3' end of the RNA, which is the 5' end of the
    cDNA-derived read. In reference coordinates:

        Read on BAM "+" strand: 5' of read is at reference_start
            => RNA is on "-" strand at position reference_start
        Read on BAM "-" strand: 5' of read is at reference_end - 1
            => RNA is on "+" strand at position reference_end - 1
    """
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return None

    if read.is_reverse:
        # FLAG 16 set: BAM "-" strand => RNA "+" strand
        # 5' of read in reference coords is the RIGHTMOST aligned position
        return (read.reference_name, read.reference_end - 1, "+")
    else:
        # FLAG 16 not set: BAM "+" strand => RNA "-" strand
        return (read.reference_name, read.reference_start, "-")


def count_at_intervals(bam_path, intervals, mode="rna"):
    """
    Count Pol II positions inside each interval, respecting the strand swap.

    Args:
        bam_path: path to indexed BAM.
        intervals: iterable of (chrom, start, end, name, strand) tuples.
                   start/end in BED-style (0-based, half-open).
                   strand may be '+', '-', '*', or '.'.
        mode: 'rna' (default) -> the strand column is the RNA strand, swap applied.
              'bam' -> the strand column is the BAM strand, no swap applied.

    Returns:
        list of integer counts in the same order as the intervals iterable.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    out = []

    for (chrom, start, end, name, strand) in intervals:
        if chrom not in bam.references:
            out.append(0)
            continue

        n = 0
        # We fetch reads overlapping a *padded* region. Using exactly [start, end]
        # would miss reads whose 5' end is inside the interval but whose body
        # extends outside it, depending on read orientation. A small pad is safe.
        # The default 500bp covers typical PRO-seq read lengths.
        try:
            iterator = bam.fetch(chrom, max(0, start - 500), end + 500)
        except ValueError:
            # Coordinates outside chromosome bounds; skip
            out.append(0)
            continue

        for read in iterator:
            pol2 = pol2_position_for_read(read)
            if pol2 is None:
                continue

            _, pos, rna_strand = pol2

            # In-window check (BED half-open: [start, end))
            if pos < start or pos >= end:
                continue

            # Strand check
            if strand in ("*", ".", "", None):
                # Combined: count both strands
                n += 1
            else:
                if mode == "rna":
                    # The interval's strand IS the RNA strand we want.
                    if rna_strand == strand:
                        n += 1
                elif mode == "bam":
                    # The interval's strand is a BAM strand. Convert what
                    # we know (rna_strand) back to BAM strand for comparison:
                    #   RNA "+" => BAM "-"
                    #   RNA "-" => BAM "+"
                    bam_strand = "-" if rna_strand == "+" else "+"
                    if bam_strand == strand:
                        n += 1
                else:
                    raise ValueError(f"Unknown mode: {mode}")

        out.append(n)

    bam.close()
    return out


# ---------------------------------------------------------------------------
# Optional helper: positions instead of just counts
# ---------------------------------------------------------------------------

def positions_at_interval(bam_path, chrom, start, end, rna_strand=None):
    """
    Return a list of Pol II position offsets (relative to start) for reads
    whose Pol II site falls in [start, end), optionally filtered to RNA
    strand. Useful when callers want to do their own summarization
    (e.g., distance from motif center, position-resolved kernel).
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    if chrom not in bam.references:
        bam.close()
        return []

    positions = []
    for read in bam.fetch(chrom, max(0, start - 500), end + 500):
        pol2 = pol2_position_for_read(read)
        if pol2 is None:
            continue
        _, pos, s = pol2
        if pos < start or pos >= end:
            continue
        if rna_strand is not None and rna_strand not in ("*", ".") and s != rna_strand:
            continue
        positions.append(pos - start)

    bam.close()
    return positions


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_intervals_bed(path):
    """
    Yield (chrom, start, end, name, strand) tuples from a BED-like file.
    Tolerates 3-, 4-, 5-, or 6-column input. Missing strand defaults to '*'.
    """
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3] if len(parts) >= 4 else f"{chrom}:{start}-{end}"
            strand = parts[5] if len(parts) >= 6 else "*"
            yield (chrom, start, end, name, strand)


def main():
    parser = argparse.ArgumentParser(
        description="Count PRO-seq Pol II positions inside BED intervals "
                    "with strand-swap awareness.")
    parser.add_argument("--bam", required=True, help="Indexed PRO-seq BAM.")
    parser.add_argument("--intervals", required=True,
                        help="6-column BED file (chrom start end name score strand). "
                             "Strand column interpreted per --strand-mode.")
    parser.add_argument("--output", default="-",
                        help="Output TSV path (default: stdout).")
    parser.add_argument("--strand-mode", choices=("rna", "bam"), default="rna",
                        help="Interpretation of strand column. "
                             "'rna' (default) = RNA strand; the swap is applied. "
                             "'bam' = literal BAM strand; no swap.")
    args = parser.parse_args()

    intervals = list(parse_intervals_bed(args.intervals))
    print(f"[INFO] Loaded {len(intervals):,} intervals from {args.intervals}",
          file=sys.stderr)
    print(f"[INFO] Strand mode: {args.strand_mode}", file=sys.stderr)

    counts = count_at_intervals(args.bam, intervals, mode=args.strand_mode)

    out = sys.stdout if args.output == "-" else open(args.output, "w")
    try:
        out.write("chrom\tstart\tend\tname\tstrand\tcount\n")
        for (chrom, start, end, name, strand), n in zip(intervals, counts):
            out.write(f"{chrom}\t{start}\t{end}\t{name}\t{strand}\t{n}\n")
    finally:
        if args.output != "-":
            out.close()

    total = sum(counts)
    nonzero = sum(1 for c in counts if c > 0)
    print(f"[INFO] Wrote counts for {len(counts):,} intervals "
          f"({nonzero:,} with non-zero count, {total:,} total Pol II positions)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
