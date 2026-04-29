#!/usr/bin/env python3
"""
02_postprocess.py — Merge per-chromosome MOODS results, compute QC stats,
optionally dedupe palindromic motifs.

PALINDROME DEDUPLICATION
-------------------------
MOODS scans both strands. For palindromic PWMs (forward log-odds matrix is
the same as the reverse-complement matrix), every '+' strand hit has a
duplicate '-' strand hit at the same position with ~equal score. Keeping
both is double-counting the same physical binding event.

Pass --palindrome-table /path/to/orientation_analysis.tsv (produced by
0x_sanity_check_motifs.py) to dedupe palindromic motifs to '+' strand only.
The dedup is only applied to the merged outputs (per_motif/, merged BED,
summary). The raw per_chrom/*.bed files are unchanged.

Usage:
    python 02_postprocess.py \\
        --input-dir /path/to/per_chrom \\
        --output-dir /path/to/merged \\
        --background /path/to/background.txt \\
        --meme /path/to/consensus_pwms.meme \\
        --chrom-sizes /path/to/mm39.chrom.sizes \\
        --palindrome-table /path/to/sanity_check/orientation_analysis.tsv \\
        --pvalue 1e-5
"""

import argparse
import os
import sys
import re
import glob
from collections import defaultdict, Counter


# ============================================================================
# Parsing
# ============================================================================

def parse_background(bg_path):
    with open(bg_path) as f:
        tokens = f.read().strip().split()
    bg = {}
    for i in range(0, len(tokens), 2):
        bg[tokens[i]] = float(tokens[i + 1])
    return bg


def get_motif_widths(meme_path):
    motifs = {}
    current_id = None
    with open(meme_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("MOTIF"):
                parts = line.split()
                current_id = parts[1].split(":")[0] if len(parts) > 1 else None
            elif line.startswith("letter-probability matrix") and current_id:
                m = re.search(r'\bw=\s*(\d+)', line)
                if m:
                    motifs[current_id] = int(m.group(1))
    return motifs


def parse_chrom_sizes(path):
    sizes = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


def compute_effective_genome_size(chrom_sizes):
    standard = set(f"chr{i}" for i in range(1, 20)) | {"chrX", "chrY"}
    total = 0
    used = []
    for chrom, size in sorted(chrom_sizes.items()):
        if chrom in standard:
            total += size
            used.append(chrom)
    return total, used


def parse_palindrome_table(path):
    """
    Parse the PWM-level palindromicity classification TSV from
    0x_sanity_check_motifs.py. Returns {motif_id: classification}.
    classification ∈ {'palindromic', 'near-palindromic', 'asymmetric'}.
    """
    classifications = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        mid_col = None
        for c in ("archetype_id", "motif_id"):
            if c in header:
                mid_col = header.index(c)
                break
        if mid_col is None:
            raise ValueError(f"Need 'archetype_id' or 'motif_id' column. Got: {header}")
        if "classification" not in header:
            raise ValueError(f"Need 'classification' column. Got: {header}")
        cls_col = header.index("classification")
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) > max(mid_col, cls_col):
                classifications[fields[mid_col]] = fields[cls_col]
    return classifications


# ============================================================================
# Palindrome dedup
# ============================================================================

def dedupe_palindromic_motifs(motif_hits, classifications, dedup_near=False):
    """
    Drop '-' hits for palindromic (and optionally near-palindromic) motifs.
    Returns (deduped_hits, drop_log).
    """
    deduped = {}
    drop_log = {}
    for mid, hits in motif_hits.items():
        cls = classifications.get(mid, "asymmetric")
        if cls == "palindromic" or (dedup_near and cls == "near-palindromic"):
            kept = [h for h in hits if h[4] == "+"]
            drop_log[mid] = len(hits) - len(kept)
            deduped[mid] = kept
        else:
            deduped[mid] = hits
            drop_log[mid] = 0
    return deduped, drop_log


# ============================================================================
# Main pipeline
# ============================================================================

def merge_and_analyze(input_dir, output_dir, bg, motif_widths, pvalue,
                     genome_size, classifications=None, dedup_near=False):
    os.makedirs(output_dir, exist_ok=True)
    per_motif_dir = os.path.join(output_dir, "per_motif")
    os.makedirs(per_motif_dir, exist_ok=True)

    bed_files = sorted(glob.glob(os.path.join(input_dir, "*.bed")))
    if not bed_files:
        print(f"ERROR: No BED files in {input_dir}", file=sys.stderr)
        sys.exit(1)
    print(f"Found {len(bed_files)} chromosome BED files", file=sys.stderr)

    # ---- Pass 1: collect ----
    motif_hits = defaultdict(list)
    total_hits = 0
    for bed_file in bed_files:
        n_file = 0
        with open(bed_file) as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) < 6:
                    continue
                chrom, start, end, motif_id, score, strand = fields[:6]
                motif_hits[motif_id].append(
                    (chrom, int(start), int(end), float(score), strand)
                )
                total_hits += 1
                n_file += 1
        chrom_name = os.path.basename(bed_file).replace(".bed", "")
        print(f"  {chrom_name}: {n_file:,} hits", file=sys.stderr)
    print(f"\nTotal raw hits: {total_hits:,}", file=sys.stderr)
    print(f"Motifs with hits: {len(motif_hits)}", file=sys.stderr)

    # ---- Optional: palindrome dedup ----
    drop_log = {mid: 0 for mid in motif_hits}
    if classifications is not None:
        n_known = sum(1 for mid in motif_hits if mid in classifications)
        print(f"\nPalindrome classification: {n_known}/{len(motif_hits)} "
              f"motifs have a class", file=sys.stderr)
        motif_hits, drop_log = dedupe_palindromic_motifs(
            motif_hits, classifications, dedup_near=dedup_near
        )
        n_dropped = sum(drop_log.values())
        n_pal = sum(1 for mid in motif_hits
                    if classifications.get(mid) == "palindromic")
        n_near = sum(1 for mid in motif_hits
                     if classifications.get(mid) == "near-palindromic")
        print(f"  Palindromic motifs: {n_pal}", file=sys.stderr)
        print(f"  Near-palindromic (dedup={'yes' if dedup_near else 'no'}): "
              f"{n_near}", file=sys.stderr)
        print(f"  Total '-' strand hits dropped: {n_dropped:,}", file=sys.stderr)

    total_after = sum(len(v) for v in motif_hits.values())
    print(f"\nTotal hits after dedup: {total_after:,}", file=sys.stderr)

    # ---- Per-motif BED files ----
    print(f"\nWriting per-motif BEDs to {per_motif_dir}/", file=sys.stderr)
    for mid, hits in motif_hits.items():
        hits.sort(key=lambda x: (x[0], x[1]))
        outpath = os.path.join(per_motif_dir, f"{mid}.bed")
        with open(outpath, "w") as f:
            for chrom, start, end, score, strand in hits:
                f.write(f"{chrom}\t{start}\t{end}\t{mid}\t{score:.4f}\t{strand}\n")

    # ---- Merged genome-wide BED (sorted, gzipped) ----
    merged_path = os.path.join(output_dir, "all_motifs_merged.bed.gz")
    tmp_unsorted = os.path.join(output_dir, "_tmp_all_motifs_merged.bed")
    with open(tmp_unsorted, "w") as f:
        for mid, hits in motif_hits.items():
            for chrom, start, end, score, strand in hits:
                f.write(f"{chrom}\t{start}\t{end}\t{mid}\t{score:.4f}\t{strand}\n")
    ret = os.system(f"sort -k1,1 -k2,2n {tmp_unsorted} | gzip > {merged_path}")
    if ret == 0:
        os.remove(tmp_unsorted)
        print(f"Wrote: {merged_path}", file=sys.stderr)
    else:
        print(f"WARNING: sort|gzip failed (exit {ret}); kept {tmp_unsorted}",
              file=sys.stderr)

    # ---- QC report ----
    report = generate_qc_report(motif_hits, motif_widths, bg, pvalue,
                                genome_size, classifications, drop_log,
                                dedup_near)
    with open(os.path.join(output_dir, "qc_report.txt"), "w") as f:
        f.write(report)
    print(f"Wrote QC report: qc_report.txt", file=sys.stderr)

    # ---- Summary TSV ----
    write_summary_tsv(motif_hits, motif_widths, bg, pvalue, genome_size,
                      os.path.join(output_dir, "motif_hit_summary.tsv"),
                      classifications, drop_log)
    print(f"Wrote summary TSV: motif_hit_summary.tsv", file=sys.stderr)


def generate_qc_report(motif_hits, motif_widths, bg, pvalue, genome_size,
                       classifications, drop_log, dedup_near):
    lines = []
    lines.append("=" * 80)
    lines.append("MOODS GENOME-WIDE SCAN — QC REPORT")
    lines.append("=" * 80 + "\n")
    lines.append(f"P-value: {pvalue}")
    lines.append(f"Effective genome size: {genome_size:,} bp")
    lines.append(f"Background: A={bg['A']:.4f} C={bg['C']:.4f} "
                 f"G={bg['G']:.4f} T={bg['T']:.4f}")
    lines.append(f"Motifs with hits: {len(motif_hits)}")
    total = sum(len(v) for v in motif_hits.values())
    lines.append(f"Total hits (after any dedup): {total:,}\n")

    if classifications is not None:
        n_dropped = sum(drop_log.values())
        n_pal = sum(1 for mid in motif_hits
                    if classifications.get(mid) == "palindromic")
        n_near = sum(1 for mid in motif_hits
                     if classifications.get(mid) == "near-palindromic")
        n_asym = sum(1 for mid in motif_hits
                     if classifications.get(mid, "asymmetric") == "asymmetric")
        lines.append("PALINDROME DEDUP")
        lines.append(f"  Strategy: dedupe palindromic"
                     f"{' + near-palindromic' if dedup_near else ''} → '+' only")
        lines.append(f"  Palindromic: {n_pal}")
        lines.append(f"  Near-palindromic (deduped: "
                     f"{'yes' if dedup_near else 'no'}): {n_near}")
        lines.append(f"  Asymmetric/unclassified: {n_asym}")
        lines.append(f"  Total '-' strand hits dropped: {n_dropped:,}\n")
    else:
        lines.append("PALINDROME DEDUP: not applied "
                     "(no --palindrome-table)\n")

    expected = pvalue * 2 * genome_size
    lines.append(f"Expected FP per motif (null): {expected:,.0f}")
    lines.append(f"  = p × 2 × L = {pvalue} × 2 × {genome_size:,}")
    lines.append("  NOTE: For palindromic motifs deduped to '+', "
                 "appropriate null is half this.\n")

    # Sanity 1: obs/exp
    lines.append("-" * 80)
    lines.append("SANITY 1: Observed vs. Expected")
    lines.append("-" * 80)
    lines.append(f"{'Motif':<15} {'Width':>5} {'Class':>16} "
                 f"{'Observed':>12} {'Expected':>12} {'Obs/Exp':>10}")
    lines.append("-" * 75)
    flagged = []
    all_ratios = []
    for mid in sorted(motif_hits.keys()):
        n_hits = len(motif_hits[mid])
        cls = (classifications.get(mid, "—") if classifications else "—")
        if classifications and (cls == "palindromic" or
                                (dedup_near and cls == "near-palindromic")):
            exp = expected / 2
        else:
            exp = expected
        ratio = n_hits / exp if exp > 0 else float('inf')
        all_ratios.append((mid, ratio))
        width = motif_widths.get(mid, "?")
        flag = ""
        if ratio > 10:
            flag = "HIGH"
            flagged.append((mid, ratio, "obs/exp > 10"))
        elif ratio < 0.1:
            flag = "LOW"
            flagged.append((mid, ratio, "obs/exp < 0.1"))
        lines.append(f"{mid:<15} {str(width):>5} {cls:>16} "
                     f"{n_hits:>12,} {exp:>12,.0f} {ratio:>10.2f} {flag}")
    lines.append("")
    if flagged:
        lines.append(f"  {len(flagged)} motifs flagged.")
    lines.append("")

    # Sanity 2: strand balance
    lines.append("-" * 80)
    lines.append("SANITY 2: Strand Balance")
    lines.append("-" * 80)
    if classifications is not None:
        lines.append("  NOTE: deduped palindromic motifs will show 100% '+' "
                     "(expected, not a real imbalance)")
    imbalanced = []
    for mid in sorted(motif_hits.keys()):
        cls = (classifications.get(mid, "asymmetric")
               if classifications else "asymmetric")
        if classifications and (cls == "palindromic" or
                                (dedup_near and cls == "near-palindromic")):
            continue
        hits = motif_hits[mid]
        n_plus = sum(1 for h in hits if h[4] == "+")
        if len(hits) > 100:
            frac = n_plus / len(hits)
            if frac < 0.45 or frac > 0.55:
                imbalanced.append((mid, frac, len(hits)))
    if imbalanced:
        lines.append(f"  {len(imbalanced)} non-deduped motifs imbalanced:")
        for mid, frac, total in imbalanced[:10]:
            lines.append(f"    {mid}: {frac:.1%} on '+' (n={total:,})")
    else:
        lines.append("  All non-deduped motifs balanced (45–55%)")
    lines.append("")

    # Sanity 3: chromosome distribution
    lines.append("-" * 80)
    lines.append("SANITY 3: Hits per chromosome")
    lines.append("-" * 80)
    chrom_counts = Counter()
    for hits in motif_hits.values():
        for h in hits:
            chrom_counts[h[0]] += 1
    for c in sorted(chrom_counts.keys(), key=lambda x: (
        int(x[3:]) if x[3:].isdigit() else (100 if x == "chrX" else 101)
    )):
        lines.append(f"    {c:<10} {chrom_counts[c]:>12,}")
    lines.append("")

    # Summary
    ratios = sorted(r for _, r in all_ratios)
    lines.append("=" * 80)
    lines.append("SUMMARY")
    lines.append("=" * 80)
    lines.append(f"  Motifs: {len(motif_hits)}")
    lines.append(f"  Total hits: {total:,}")
    if ratios:
        lines.append(f"  Median obs/exp: {ratios[len(ratios)//2]:.2f}")
        lines.append(f"  Motifs obs/exp > 10: {sum(1 for r in ratios if r > 10)}")
        lines.append(f"  Motifs obs/exp < 0.1: {sum(1 for r in ratios if r < 0.1)}")
    return "\n".join(lines)


def write_summary_tsv(motif_hits, motif_widths, bg, pvalue, genome_size,
                      outpath, classifications=None, drop_log=None):
    expected = pvalue * 2 * genome_size
    with open(outpath, "w") as f:
        f.write("motif_id\twidth\tpalindromic_class\tn_dropped_palindrome_dups\t"
                "n_hits\tn_plus\tn_minus\texpected_null\tobs_over_exp\t"
                "min_score\tmedian_score\tmax_score\n")
        for mid in sorted(motif_hits.keys()):
            hits = motif_hits[mid]
            n = len(hits)
            if n == 0:
                continue
            n_plus = sum(1 for h in hits if h[4] == "+")
            n_minus = n - n_plus
            scores = sorted([h[3] for h in hits])
            width = motif_widths.get(mid, -1)
            cls = (classifications.get(mid, "")
                   if classifications else "")
            n_dropped = drop_log.get(mid, 0) if drop_log else 0
            this_exp = expected / 2 if (cls == "palindromic" and n_dropped > 0) else expected
            ratio = n / this_exp if this_exp > 0 else 0
            f.write(f"{mid}\t{width}\t{cls}\t{n_dropped}\t"
                    f"{n}\t{n_plus}\t{n_minus}\t"
                    f"{this_exp:.0f}\t{ratio:.4f}\t"
                    f"{scores[0]:.4f}\t{scores[n//2]:.4f}\t{scores[-1]:.4f}\n")


# ============================================================================
# CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Merge per-chromosome MOODS BEDs, optional palindrome dedup, QC.",
    )
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--background", required=True)
    parser.add_argument("--meme", required=True)
    parser.add_argument("--pvalue", type=float, default=1e-5)
    parser.add_argument("--chrom-sizes", default=None,
                        help="UCSC chrom.sizes (RECOMMENDED)")
    parser.add_argument("--genome-size", type=int, default=None,
                        help="Used only if --chrom-sizes not provided")
    parser.add_argument("--palindrome-table", default=None,
                        help="orientation_analysis.tsv from 0x_sanity_check_motifs.py")
    parser.add_argument("--dedup-near-palindromic", action="store_true")
    args = parser.parse_args()

    # Genome size
    if args.chrom_sizes:
        chrom_sizes = parse_chrom_sizes(args.chrom_sizes)
        if chrom_sizes:
            genome_size, used = compute_effective_genome_size(chrom_sizes)
            print(f"Effective genome size from {args.chrom_sizes}: "
                  f"{genome_size:,} bp ({len(used)} chroms)", file=sys.stderr)
        elif args.genome_size:
            genome_size = args.genome_size
        else:
            print("ERROR: --chrom-sizes had no parseable sizes", file=sys.stderr)
            sys.exit(1)
    elif args.genome_size:
        genome_size = args.genome_size
    else:
        print("ERROR: need --chrom-sizes or --genome-size", file=sys.stderr)
        sys.exit(1)

    bg = parse_background(args.background)
    motif_widths = get_motif_widths(args.meme)
    print(f"Loaded {len(motif_widths)} motif widths", file=sys.stderr)

    classifications = None
    if args.palindrome_table:
        classifications = parse_palindrome_table(args.palindrome_table)
        print(f"Loaded {len(classifications)} palindrome classifications "
              f"from {args.palindrome_table}", file=sys.stderr)

    merge_and_analyze(
        args.input_dir, args.output_dir, bg, motif_widths,
        args.pvalue, genome_size,
        classifications=classifications,
        dedup_near=args.dedup_near_palindromic,
    )
    print("\nDone!", file=sys.stderr)


if __name__ == "__main__":
    main()
