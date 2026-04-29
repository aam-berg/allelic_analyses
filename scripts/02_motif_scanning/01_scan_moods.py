#!/usr/bin/env python3
"""
01_scan_moods.py — Scan one chromosome against all motif archetypes using MOODS.

This scans ALL motifs against a single chromosome in one pass, which is
efficient because MOODS uses a multi-pattern lookahead filtration algorithm.

Outputs BED6 format:
    chrom  start  end  motif_id  score  strand

Usage:
    python 01_scan_moods.py \
        --genome <genome.fa> \
        --chrom <chrN> \
        --meme <consensus_pwms.meme> \
        --background <background.txt> \
        --pvalue 1e-5 \
        --output <output.bed>
"""

import argparse
import sys
import os
import re
from collections import defaultdict, Counter
import MOODS.tools
import MOODS.scan


def parse_background(bg_path):
    """Parse background file: 'A 0.XXXX C 0.XXXX G 0.XXXX T 0.XXXX'"""
    with open(bg_path) as f:
        tokens = f.read().strip().split()
    bg = {}
    for i in range(0, len(tokens), 2):
        bg[tokens[i]] = float(tokens[i + 1])
    return [bg["A"], bg["C"], bg["G"], bg["T"]]


def parse_meme_file(meme_path):
    """
    Parse a MEME format file and extract all motif PWMs.

    Handles both 'w=6' and 'w= 6' (space after =) formats.

    Returns:
        list of (motif_id, matrix) where matrix is a list of 4 lists
        (one per nucleotide A, C, G, T) — the format MOODS expects.
    """
    motifs = []
    current_id = None
    current_matrix = None
    in_matrix = False
    expected_width = 0
    rows_read = 0

    with open(meme_path) as f:
        for line in f:
            line = line.strip()

            if line.startswith("MOTIF"):
                # Save previous motif if any
                if current_id is not None and current_matrix is not None:
                    transposed = _transpose_matrix(current_matrix)
                    motifs.append((current_id, transposed))

                parts = line.split()
                current_id = parts[1] if len(parts) > 1 else parts[0]
                current_id = current_id.split(":")[0]
                current_matrix = []
                in_matrix = False
                rows_read = 0

            elif line.startswith("letter-probability matrix"):
                in_matrix = True
                match = re.search(r'\bw=\s*(\d+)', line)
                if match:
                    expected_width = int(match.group(1))
                else:
                    print(f"WARNING: Could not parse width from: {line}",
                          file=sys.stderr)
                    expected_width = 0
                rows_read = 0

            elif in_matrix and line and not line.startswith("URL") and not line.startswith("MOTIF"):
                values = line.split()
                if len(values) == 4:
                    try:
                        row = [float(v) for v in values]
                        current_matrix.append(row)
                        rows_read += 1
                        if expected_width > 0 and rows_read >= expected_width:
                            in_matrix = False
                    except ValueError:
                        in_matrix = False
                else:
                    in_matrix = False

    # Last motif
    if current_id is not None and current_matrix is not None:
        transposed = _transpose_matrix(current_matrix)
        motifs.append((current_id, transposed))

    return motifs


def _transpose_matrix(matrix):
    """
    Transpose from position-by-nucleotide to nucleotide-by-position.
    Input:  [[pA, pC, pG, pT], ...] per position
    Output: [[A at pos0, A at pos1, ...], [C ...], [G ...], [T ...]]
    """
    result = [[], [], [], []]
    for row in matrix:
        for nuc in range(4):
            result[nuc].append(row[nuc])
    return result


def read_fasta_chrom(fasta_path, target_chrom):
    """Read a single chromosome sequence from a FASTA file."""
    seq_parts = []
    reading = False

    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if reading:
                    break  # Past our chromosome
                chrom = line[1:].strip().split()[0]
                if chrom == target_chrom:
                    reading = True
            elif reading:
                seq_parts.append(line.strip().upper())

    if not seq_parts:
        print(f"ERROR: Chromosome {target_chrom} not found in {fasta_path}",
              file=sys.stderr)
        sys.exit(1)

    return "".join(seq_parts)


def scan_chromosome(sequence, motifs, bg, pvalue, pseudocount=0.001):
    """
    Scan a chromosome sequence against all motifs using MOODS.

    Returns list of (motif_id, start, end, score, strand) tuples.
    """
    matrices = []
    matrix_ids = []
    widths = []

    for motif_id, pfm in motifs:
        log_odds = MOODS.tools.log_odds(pfm, bg, pseudocount)
        rev_comp = MOODS.tools.reverse_complement(log_odds)

        matrices.append(log_odds)
        matrices.append(rev_comp)
        matrix_ids.append((motif_id, "+"))
        matrix_ids.append((motif_id, "-"))

        w = len(pfm[0])
        widths.append(w)
        widths.append(w)

    # Compute per-matrix score thresholds from p-value
    thresholds = [
        MOODS.tools.threshold_from_p(m, bg, pvalue)
        for m in matrices
    ]

    print(f"  Scanning {len(motifs)} motifs ({len(matrices)} with rev comp) "
          f"at p={pvalue}...", file=sys.stderr)

    # Multi-pattern scan
    results = MOODS.scan.scan_dna(sequence, matrices, bg, thresholds)

    # Collect hits — MOODS returns match objects with .pos and .score attributes
    hits = []
    for i, matrix_result in enumerate(results):
        motif_id, strand = matrix_ids[i]
        w = widths[i]
        for match in matrix_result:
            pos = int(match.pos)
            if pos < 0:
                continue
            hits.append((motif_id, pos, pos + w, match.score, strand))

    return hits


def filter_overlaps_within_motif(hits):
    """
    Greedy non-overlapping filter per (motif, strand):
    sort by score descending, keep a hit only if it doesn't overlap
    any already-kept hit.

    Uses a sort-and-sweep approach: after selecting hits greedily by score,
    we check overlap efficiently by maintaining kept hits sorted by start
    and using binary search, giving O(n log n) per group instead of O(n^2).
    """
    from bisect import insort, bisect_left, bisect_right

    # Group by (motif_id, strand)
    grouped = defaultdict(list)
    for motif_id, start, end, score, strand in hits:
        grouped[(motif_id, strand)].append((score, start, end))

    filtered = []
    for (motif_id, strand), group_hits in grouped.items():
        # Sort by score descending (greedy: best score first)
        group_hits.sort(key=lambda x: -x[0])

        # Maintain sorted list of (start, end) of kept intervals
        # for efficient overlap checking
        kept_starts = []  # sorted list of start positions
        kept_intervals = {}  # start -> end mapping

        for score, start, end in group_hits:
            # Find intervals that could overlap: those with start < end
            # and end > start. We check nearby intervals via binary search.
            overlap = False

            # Check intervals starting before this end
            idx = bisect_left(kept_starts, start)

            # Check the interval just before idx (it could extend past our start)
            if idx > 0:
                prev_start = kept_starts[idx - 1]
                if kept_intervals[prev_start] > start:  # prev_end > our start
                    overlap = True

            # Check intervals starting at or after our start but before our end
            if not overlap:
                check_idx = idx
                while check_idx < len(kept_starts) and kept_starts[check_idx] < end:
                    overlap = True
                    break

            if not overlap:
                insort(kept_starts, start)
                kept_intervals[start] = end
                filtered.append((motif_id, start, end, score, strand))

    return filtered


def main():
    parser = argparse.ArgumentParser(
        description="Scan one chromosome against all motif archetypes using MOODS"
    )
    parser.add_argument("--genome", required=True, help="Genome FASTA file")
    parser.add_argument("--chrom", required=True, help="Chromosome to scan")
    parser.add_argument("--meme", required=True, help="MEME format motif file")
    parser.add_argument("--background", required=True, help="Background frequencies file")
    parser.add_argument("--pvalue", type=float, default=1e-5,
                        help="P-value threshold (default: 1e-5)")
    parser.add_argument("--pseudocount", type=float, default=0.001,
                        help="Pseudocount for PWM (default: 0.001)")
    parser.add_argument("--output", required=True, help="Output BED file")
    parser.add_argument("--no-overlap-filter", action="store_true",
                        help="Skip within-motif overlap filtering")
    args = parser.parse_args()

    # Load background
    print(f"Loading background from {args.background}", file=sys.stderr)
    bg = parse_background(args.background)
    print(f"  Background: A={bg[0]:.4f} C={bg[1]:.4f} G={bg[2]:.4f} T={bg[3]:.4f}",
          file=sys.stderr)

    # Load motifs
    print(f"Loading motifs from {args.meme}", file=sys.stderr)
    motifs = parse_meme_file(args.meme)
    print(f"  Loaded {len(motifs)} motif archetypes", file=sys.stderr)

    if len(motifs) == 0:
        print("ERROR: No motifs loaded. Check MEME file format.", file=sys.stderr)
        sys.exit(1)

    # Print motif summary
    widths = [len(m[1][0]) for m in motifs]
    print(f"  Motif widths: min={min(widths)}, max={max(widths)}, "
          f"median={sorted(widths)[len(widths)//2]}", file=sys.stderr)

    # Warn about very short motifs
    short = [(mid, len(m[0])) for mid, m in motifs if len(m[0]) <= 4]
    if short:
        print(f"  WARNING: {len(short)} motifs have width <= 4 bp — expect "
              f"very high hit counts for these", file=sys.stderr)
        for mid, w in short[:5]:
            print(f"    {mid}: {w} bp", file=sys.stderr)

    # Load chromosome sequence
    print(f"Loading chromosome {args.chrom}...", file=sys.stderr)
    sequence = read_fasta_chrom(args.genome, args.chrom)
    print(f"  Loaded {len(sequence):,} bp", file=sys.stderr)

    # Scan
    hits = scan_chromosome(sequence, motifs, bg, args.pvalue, args.pseudocount)
    print(f"  Raw hits: {len(hits):,}", file=sys.stderr)

    # Filter overlaps within each motif
    if not args.no_overlap_filter:
        hits = filter_overlaps_within_motif(hits)
        print(f"  After overlap filter: {len(hits):,}", file=sys.stderr)

    # Sort by position
    hits.sort(key=lambda x: (x[1], x[0]))

    # Write output
    outdir = os.path.dirname(os.path.abspath(args.output))
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    with open(args.output, "w") as f:
        for motif_id, start, end, score, strand in hits:
            f.write(f"{args.chrom}\t{start}\t{end}\t{motif_id}\t{score:.4f}\t{strand}\n")

    print(f"  Wrote {len(hits):,} hits to {args.output}", file=sys.stderr)

    # Per-motif summary
    motif_counts = Counter(h[0] for h in hits)
    print(f"\n  Per-motif hit counts (top 10):", file=sys.stderr)
    for mid, count in motif_counts.most_common(10):
        print(f"    {mid}: {count:,}", file=sys.stderr)
    print(f"  ...", file=sys.stderr)
    print(f"  Per-motif hit counts (bottom 5):", file=sys.stderr)
    for mid, count in motif_counts.most_common()[-5:]:
        print(f"    {mid}: {count:,}", file=sys.stderr)


if __name__ == "__main__":
    main()