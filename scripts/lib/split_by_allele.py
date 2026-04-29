#!/usr/bin/env python3
"""
split_by_allele.py — Split BAM reads by allele at heterozygous SNP positions.

Originally written for single-end PRO-seq (in 04_proseq/); generalized here to
support paired-end ATAC-seq as well. The two modes differ only in how reads
are aggregated for fragment-level classification:

  Single-end (PRO-seq):
    Each read is one fragment. Classify based on SNPs the read covers.

  Paired-end (ATAC, RNA-seq, etc.):
    A "fragment" is a read pair. We collect both mates by qname, classify
    based on the union of their SNP evidence, and write both mates to the
    output BAM corresponding to that classification.

For each read (or read pair) the classification is:

  ref       - read's allele(s) at overlapping het SNP(s) match REF
  alt       - read's allele(s) at overlapping het SNP(s) match ALT
  nosnp     - read does not overlap any het SNP
  ambiguous - read covers SNPs and shows a mix of ref and alt, OR shows only
              non-ref/non-alt bases (sequencing errors, triallelic, etc.)

Output files:
  {prefix}_ref.bam               reads supporting REF
  {prefix}_alt.bam               reads supporting ALT
  {prefix}_nosnp.bam             reads with no SNP overlap
  {prefix}_ambiguous.bam         reads with conflicting alleles
  {prefix}_allele_counts.tsv     per-SNP allele count summary

Per-SNP counts:
  In paired-end mode, we count fragment-level evidence at each SNP. If both
  mates of a fragment cover the same SNP and report the same allele, that's
  one count, not two. If the two mates of a fragment cover the same SNP and
  report DIFFERENT alleles (sequencing error or rare biology), the SNP is
  treated as "other" for that fragment and the fragment is classified as
  ambiguous overall.

Strand notes:
  This script is strand-agnostic — pysam returns query bases in reference
  orientation, so we directly compare to VCF REF / ALT alleles regardless
  of which BAM strand the read is on. Any downstream PRO-seq strand-swap
  for bigWig generation happens elsewhere (lib/strand_aware_count.py and
  04_proseq/04_make_bigwigs.sh).

Memory:
  Paired-end mode buffers reads by qname while waiting for their mates. For
  coordinate-sorted BAMs with normal insert sizes, the buffer stays small
  (~10s of MB). For very long inserts or chimeric pairs, the buffer grows.
  We emit a warning if the buffer ever exceeds 500k unmatched reads,
  suggesting name-sorting first (`samtools sort -n`).

Usage:
  # Single-end (PRO-seq), auto-detected if first read isn't paired:
  python split_by_allele.py \\
      --bam sample_wasp.bam \\
      --snp_dir /path/to/wasp_snp_files/ \\
      --output_prefix /path/to/output/sample

  # Paired-end (ATAC), explicitly:
  python split_by_allele.py \\
      --bam sample_wasp.bam \\
      --snp_dir /path/to/wasp_snp_files/ \\
      --output_prefix /path/to/output/sample \\
      --paired-end

  # Force single-end interpretation even if reads are flagged as paired:
  python split_by_allele.py ... --single-end
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict

import pysam


# ---------------------------------------------------------------------------
# SNP loading
# ---------------------------------------------------------------------------

def load_snps(snp_dir):
    """
    Load WASP-format SNP files into a nested dictionary.

    The WASP format is per-chromosome gzipped tab files with three columns:
        position(1-based)   ref_allele   alt_allele

    Returns:
        dict: {chrom: {pos_1based: (ref, alt)}}
    """
    snps = {}
    snp_files = [f for f in os.listdir(snp_dir) if f.endswith(".snps.txt.gz")]

    if not snp_files:
        raise FileNotFoundError(f"No .snps.txt.gz files found in {snp_dir}")

    for fname in sorted(snp_files):
        chrom = fname.replace(".snps.txt.gz", "")
        snps[chrom] = {}
        filepath = os.path.join(snp_dir, fname)

        with gzip.open(filepath, "rt") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    pos = int(parts[0])
                    ref = parts[1].upper()
                    alt = parts[2].upper()
                    snps[chrom][pos] = (ref, alt)

        print(f"  Loaded {len(snps[chrom]):,} SNPs for {chrom}", file=sys.stderr)

    total = sum(len(v) for v in snps.values())
    print(f"  Total SNPs loaded: {total:,}", file=sys.stderr)
    return snps


# ---------------------------------------------------------------------------
# Per-read allele extraction
# ---------------------------------------------------------------------------

def get_read_alleles(read, chrom_snps):
    """
    For a read, determine which allele it carries at each overlapping SNP.

    Args:
        read: pysam AlignedSegment.
        chrom_snps: dict {pos_1based: (ref_allele, alt_allele)} for this chrom.

    Returns:
        list of (pos, ref, alt, read_base, assignment) tuples, where
        assignment is one of 'ref', 'alt', 'other'.
    """
    if not chrom_snps:
        return []

    out = []
    # matches_only=True -> skip soft-clipped, insertion-relative, and deletion-
    # relative positions; only return aligned matches/mismatches.
    aligned_pairs = read.get_aligned_pairs(matches_only=True)

    for query_pos, ref_pos in aligned_pairs:
        snp_pos = ref_pos + 1  # 0-based -> 1-based
        if snp_pos in chrom_snps:
            ref_allele, alt_allele = chrom_snps[snp_pos]
            read_base = read.query_sequence[query_pos].upper()
            if read_base == ref_allele:
                assignment = "ref"
            elif read_base == alt_allele:
                assignment = "alt"
            else:
                assignment = "other"
            out.append((snp_pos, ref_allele, alt_allele, read_base, assignment))
    return out


def classify_assignments(assignments):
    """
    Given a list of allele assignments ('ref', 'alt', 'other') from one fragment
    (one read in SE, two reads' worth in PE), return the fragment's class.

    Rules:
      No assignments at all          -> 'nosnp'
      Only 'other' (errors)          -> 'ambiguous'
      All informative are 'ref'      -> 'ref'
      All informative are 'alt'      -> 'alt'
      Mix of 'ref' and 'alt'         -> 'ambiguous'

    The 'other' category includes sequencing errors and (rare) third alleles.
    They're filtered out before the ref/alt unanimity check.
    """
    if not assignments:
        return "nosnp"

    informative = [a for a in assignments if a in ("ref", "alt")]
    if not informative:
        return "ambiguous"

    if len(set(informative)) == 1:
        return informative[0]
    return "ambiguous"


def merge_fragment_snp_evidence(hits1, hits2):
    """
    For paired-end mode, combine per-read SNP evidence from two mates into
    fragment-level evidence.

    If both mates cover the same SNP:
      - Same allele on both reads -> count as one informative observation.
      - Different alleles on the two reads (rare; sequencing error or
        recombination) -> mark as 'other'.
    If only one mate covers a SNP -> use that mate's call.

    Returns a list of (pos, ref, alt, observed_base_or_token, assignment)
    tuples representing the FRAGMENT'S evidence at each unique SNP. The
    "observed_base_or_token" is the consensus base if both mates agree, or
    the literal string "CONFLICT" if they disagree.
    """
    by_pos = {}  # pos -> list of (pos, ref, alt, base, assignment)
    for h in hits1 + hits2:
        by_pos.setdefault(h[0], []).append(h)

    merged = []
    for pos, evidence_list in by_pos.items():
        if len(evidence_list) == 1:
            merged.append(evidence_list[0])
        else:
            # Both mates covered this SNP. Check agreement on assignment.
            assignments = {e[4] for e in evidence_list}
            ref, alt = evidence_list[0][1], evidence_list[0][2]
            if len(assignments) == 1:
                # Consistent: keep one entry.
                a = assignments.pop()
                base = evidence_list[0][3]
                merged.append((pos, ref, alt, base, a))
            else:
                # Inconsistent: mark as 'other'.
                merged.append((pos, ref, alt, "CONFLICT", "other"))
    return merged


# ---------------------------------------------------------------------------
# SE processing
# ---------------------------------------------------------------------------

def process_single_end(bam_in, snps, bam_outs, snp_counts, report_interval=5_000_000):
    """
    Iterate through bam_in treating each read as one fragment.

    bam_outs: dict mapping classification to open output BAM.
    snp_counts: defaultdict-of-dicts; updated in place with per-SNP counts.

    Returns dict of summary counts.
    """
    counts = {"ref": 0, "alt": 0, "nosnp": 0, "ambiguous": 0, "total": 0}
    for read in bam_in.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        counts["total"] += 1
        if counts["total"] % report_interval == 0:
            print(f"  Processed {counts['total']:,} reads...", file=sys.stderr)

        chrom_snps = snps.get(read.reference_name, {})
        hits = get_read_alleles(read, chrom_snps)
        cls = classify_assignments([h[4] for h in hits])
        for pos, _, _, _, a in hits:
            snp_counts[(read.reference_name, pos)][a] += 1
        counts[cls] += 1
        bam_outs[cls].write(read)
    return counts


# ---------------------------------------------------------------------------
# PE processing
# ---------------------------------------------------------------------------

def process_paired_end(bam_in, snps, bam_outs, snp_counts,
                       report_interval=5_000_000, buffer_warn_threshold=500_000):
    """
    Iterate through bam_in pairing R1/R2 by qname.

    Both mates of a pair are classified as a fragment, and BOTH mates are
    written to the output BAM corresponding to that classification (so
    downstream tools see properly-paired output).

    Buffers reads by qname while waiting for the mate. For coordinate-sorted
    BAMs with reasonable insert sizes the buffer is small. We warn if it
    exceeds buffer_warn_threshold.
    """
    pending = {}  # qname -> first-seen mate
    counts = {"ref": 0, "alt": 0, "nosnp": 0, "ambiguous": 0, "total": 0}
    orphan_count = 0
    warned = False

    def _classify_fragment(reads):
        """Classify a list of mates as one fragment. Update SNP counts; return class."""
        all_hits = []
        for r in reads:
            chrom_snps = snps.get(r.reference_name, {})
            all_hits.append(get_read_alleles(r, chrom_snps))

        if len(reads) == 2:
            merged = merge_fragment_snp_evidence(all_hits[0], all_hits[1])
        else:
            merged = all_hits[0]

        cls = classify_assignments([h[4] for h in merged])
        for pos, _, _, _, a in merged:
            # Use the chromosome of the first read; both mates of a properly
            # paired fragment must be on the same chromosome.
            snp_counts[(reads[0].reference_name, pos)][a] += 1
        return cls

    for read in bam_in.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        # Treat unpaired reads in a paired BAM as SE (rare but possible after filtering)
        if not read.is_paired or read.mate_is_unmapped:
            counts["total"] += 1
            cls = _classify_fragment([read])
            counts[cls] += 1
            bam_outs[cls].write(read)
            continue

        qname = read.query_name
        if qname in pending:
            mate = pending.pop(qname)
            counts["total"] += 1
            if counts["total"] % report_interval == 0:
                print(f"  Processed {counts['total']:,} fragments "
                      f"(buffer: {len(pending):,})...", file=sys.stderr)

            cls = _classify_fragment([mate, read])
            counts[cls] += 1
            bam_outs[cls].write(mate)
            bam_outs[cls].write(read)
        else:
            pending[qname] = read

            if len(pending) > buffer_warn_threshold and not warned:
                print(f"  WARNING: pending-mate buffer has grown to "
                      f"{len(pending):,} reads. If this BAM is coordinate-sorted "
                      f"with very long inserts or many split reads, consider "
                      f"name-sorting first (samtools sort -n) and re-running.",
                      file=sys.stderr)
                warned = True

    # Anything left in 'pending' had a mate that never appeared (orphan).
    # Process as SE.
    for read in pending.values():
        orphan_count += 1
        counts["total"] += 1
        cls = _classify_fragment([read])
        counts[cls] += 1
        bam_outs[cls].write(read)

    if orphan_count > 0:
        print(f"  NOTE: {orphan_count:,} reads had no mate in the BAM "
              f"(orphans, processed as SE).", file=sys.stderr)

    return counts


# ---------------------------------------------------------------------------
# Layout detection
# ---------------------------------------------------------------------------

def detect_layout(bam_path, n_check=1000):
    """
    Peek at the first n_check primary alignments to decide if the BAM is
    paired-end (any read with FLAG paired set) or single-end.

    Returns 'paired' or 'single'.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    n = 0
    has_paired = False
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        n += 1
        if read.is_paired:
            has_paired = True
            break
        if n >= n_check:
            break
    bam.close()
    return "paired" if has_paired else "single"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Split BAM reads by allele at het SNP positions. Supports "
                    "single-end (PRO-seq) and paired-end (ATAC) input.")
    parser.add_argument("--bam", required=True, help="Input BAM (typically WASP-filtered).")
    parser.add_argument("--snp_dir", required=True, help="WASP per-chromosome SNP directory.")
    parser.add_argument("--output_prefix", required=True,
                        help="Output prefix (e.g., /path/sample).")
    layout = parser.add_mutually_exclusive_group()
    layout.add_argument("--paired-end", action="store_true",
                        help="Force paired-end fragment-level classification.")
    layout.add_argument("--single-end", action="store_true",
                        help="Force single-end per-read classification.")
    args = parser.parse_args()

    print("=" * 60, file=sys.stderr)
    print("split_by_allele.py", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    print(f"  Input BAM:     {args.bam}", file=sys.stderr)
    print(f"  SNP dir:       {args.snp_dir}", file=sys.stderr)
    print(f"  Output prefix: {args.output_prefix}", file=sys.stderr)

    if args.paired_end:
        layout = "paired"
        print("  Layout:        PAIRED-END (forced)", file=sys.stderr)
    elif args.single_end:
        layout = "single"
        print("  Layout:        SINGLE-END (forced)", file=sys.stderr)
    else:
        layout = detect_layout(args.bam)
        print(f"  Layout:        {layout.upper()} (auto-detected)", file=sys.stderr)
    print("", file=sys.stderr)

    # --- Load SNPs ---
    print("Loading SNPs...", file=sys.stderr)
    snps = load_snps(args.snp_dir)
    print("", file=sys.stderr)

    # --- Open BAMs ---
    print("Opening BAM files...", file=sys.stderr)
    bam_in = pysam.AlignmentFile(args.bam, "rb")

    out_dir = os.path.dirname(args.output_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    bam_outs = {
        "ref":       pysam.AlignmentFile(f"{args.output_prefix}_ref.bam", "wb", template=bam_in),
        "alt":       pysam.AlignmentFile(f"{args.output_prefix}_alt.bam", "wb", template=bam_in),
        "nosnp":     pysam.AlignmentFile(f"{args.output_prefix}_nosnp.bam", "wb", template=bam_in),
        "ambiguous": pysam.AlignmentFile(f"{args.output_prefix}_ambiguous.bam", "wb", template=bam_in),
    }

    snp_counts = defaultdict(lambda: {"ref": 0, "alt": 0, "other": 0})

    # --- Process ---
    print(f"Processing reads (layout={layout})...", file=sys.stderr)
    if layout == "paired":
        counts = process_paired_end(bam_in, snps, bam_outs, snp_counts)
    else:
        counts = process_single_end(bam_in, snps, bam_outs, snp_counts)

    bam_in.close()
    for b in bam_outs.values():
        b.close()

    # --- Per-SNP counts ---
    counts_file = f"{args.output_prefix}_allele_counts.tsv"
    print(f"\nWriting per-SNP allele counts to {counts_file}...", file=sys.stderr)
    with open(counts_file, "w") as f:
        f.write("chrom\tpos\tref_allele\talt_allele\t"
                "ref_count\talt_count\tother_count\ttotal_count\tref_fraction\n")
        for (chrom, pos) in sorted(snp_counts.keys(), key=lambda x: (x[0], x[1])):
            c = snp_counts[(chrom, pos)]
            ref_allele, alt_allele = snps[chrom][pos]
            total = c["ref"] + c["alt"] + c["other"]
            ref_frac = c["ref"] / total if total > 0 else None
            ref_frac_str = f"{ref_frac:.4f}" if ref_frac is not None else "NA"
            f.write(f"{chrom}\t{pos}\t{ref_allele}\t{alt_allele}\t"
                    f"{c['ref']}\t{c['alt']}\t{c['other']}\t{total}\t"
                    f"{ref_frac_str}\n")

    covered_snps = sum(1 for c in snp_counts.values() if c["ref"] + c["alt"] > 0)

    # --- Summary ---
    unit = "fragments" if layout == "paired" else "reads"
    print("", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    print("Summary", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    print(f"  Total {unit} processed:  {counts['total']:,}", file=sys.stderr)
    for k in ("ref", "alt", "nosnp", "ambiguous"):
        pct = (100 * counts[k] / counts["total"]) if counts["total"] else 0.0
        print(f"  {k:9s}:  {counts[k]:>14,}  ({pct:5.2f}%)", file=sys.stderr)
    print("", file=sys.stderr)
    total_snps_loaded = sum(len(v) for v in snps.values())
    print(f"  Het SNPs with coverage:  {covered_snps:,} / {total_snps_loaded:,}",
          file=sys.stderr)

    allelic = counts["ref"] + counts["alt"]
    if allelic > 0:
        ref_pct = 100 * counts["ref"] / allelic
        print(f"  Global ref fraction:     {ref_pct:.2f}% "
              f"(expect ~50% after WASP)", file=sys.stderr)
        if abs(ref_pct - 50) > 5:
            print(f"  [WARNING] Ref fraction deviates from 50% by >5pp.",
                  file=sys.stderr)
        else:
            print(f"  [OK] Ref fraction within 5pp of 50%.", file=sys.stderr)

    print("", file=sys.stderr)
    print("Output files:", file=sys.stderr)
    for suffix in ("_ref.bam", "_alt.bam", "_nosnp.bam", "_ambiguous.bam",
                   "_allele_counts.tsv"):
        fpath = f"{args.output_prefix}{suffix}"
        if os.path.exists(fpath):
            sz_mb = os.path.getsize(fpath) / 1e6
            print(f"  {fpath} ({sz_mb:.1f} MB)", file=sys.stderr)


if __name__ == "__main__":
    main()
