#!/usr/bin/env python3
"""
split_by_allele.py — Split BAM reads by allele at heterozygous SNP positions.

For each read in a WASP-filtered BAM:
  1. Check if the read overlaps any het SNP
  2. Determine which allele the read carries at each overlapping SNP
  3. Assign the read to ref (allele1/129S1) or alt (allele2/CAST)
  4. Handle reads overlapping multiple SNPs (require consistency)

Output files:
  {prefix}_ref.bam       — Reads supporting the reference allele (129S1-like)
  {prefix}_alt.bam       — Reads supporting the alternative allele (CAST)
  {prefix}_nosnp.bam     — Reads not overlapping any het SNP
  {prefix}_ambiguous.bam — Reads with conflicting alleles at multiple SNPs
  {prefix}_allele_counts.tsv — Per-SNP allele counts

NOTE ON STRAND:
  Allele assignment is strand-agnostic with respect to the PRO-seq strand swap.
  In the BAM, pysam gives us aligned query bases in reference orientation, so
  we directly compare to VCF REF/ALT alleles (which are on the + strand).
  The strand swap only matters later when assigning reads to genomic strands
  for bigWig generation.

Usage:
  python split_by_allele.py \\
      --bam sample_wasp.bam \\
      --snp_dir /path/to/wasp_snp_files/ \\
      --output_prefix /path/to/output/sample
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict

import pysam


def load_snps(snp_dir):
    """
    Load WASP-format SNP files into a dictionary.

    Returns:
        dict: {chrom: {pos_1based: (ref, alt)}}
    """
    snps = {}
    snp_files = [f for f in os.listdir(snp_dir) if f.endswith('.snps.txt.gz')]

    if not snp_files:
        raise FileNotFoundError(f"No .snps.txt.gz files found in {snp_dir}")

    for fname in sorted(snp_files):
        chrom = fname.replace('.snps.txt.gz', '')
        snps[chrom] = {}
        filepath = os.path.join(snp_dir, fname)

        with gzip.open(filepath, 'rt') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    pos = int(parts[0])  # 1-based
                    ref = parts[1].upper()
                    alt = parts[2].upper()
                    snps[chrom][pos] = (ref, alt)

        print(f"  Loaded {len(snps[chrom]):,} SNPs for {chrom}")

    total = sum(len(v) for v in snps.values())
    print(f"  Total SNPs loaded: {total:,}")
    return snps


def get_read_alleles(read, chrom_snps):
    """
    For a given read, determine which allele it carries at each overlapping SNP.

    Args:
        read: pysam AlignedSegment
        chrom_snps: dict {pos_1based: (ref, alt)} for this chromosome

    Returns:
        list of tuples: [(pos, ref_allele, alt_allele, read_base, assignment)]
        assignment is 'ref', 'alt', or 'other'
    """
    if not chrom_snps:
        return []

    results = []

    # get_aligned_pairs returns (query_pos, ref_pos) where ref_pos is 0-based
    aligned_pairs = read.get_aligned_pairs(matches_only=True)

    for query_pos, ref_pos in aligned_pairs:
        snp_pos = ref_pos + 1  # convert to 1-based to match SNP dict

        if snp_pos in chrom_snps:
            ref_allele, alt_allele = chrom_snps[snp_pos]
            read_base = read.query_sequence[query_pos].upper()

            if read_base == ref_allele:
                assignment = 'ref'
            elif read_base == alt_allele:
                assignment = 'alt'
            else:
                assignment = 'other'  # sequencing error, triallelic, etc.

            results.append((snp_pos, ref_allele, alt_allele, read_base, assignment))

    return results


def classify_read(allele_hits):
    """
    Classify a read based on its allele assignments at overlapping SNPs.

    Rules:
      - If no SNPs overlap: 'nosnp'
      - If all SNP alleles are 'ref': 'ref'
      - If all SNP alleles are 'alt': 'alt'
      - If mix of 'ref' and 'alt': 'ambiguous' (possible recombination or error)
      - If any 'other' and rest are consistent: use the consistent allele
      - If only 'other': 'ambiguous'

    Returns:
        str: 'ref', 'alt', 'nosnp', or 'ambiguous'
    """
    if not allele_hits:
        return 'nosnp'

    assignments = [h[4] for h in allele_hits]

    # Filter out 'other' (sequencing errors)
    informative = [a for a in assignments if a in ('ref', 'alt')]

    if not informative:
        return 'ambiguous'

    unique_alleles = set(informative)

    if len(unique_alleles) == 1:
        return informative[0]
    else:
        # Mix of ref and alt at different SNPs — this read is ambiguous
        # (could be sequencing error, or very rarely, recombination within read)
        return 'ambiguous'


def main():
    parser = argparse.ArgumentParser(
        description='Split BAM reads by allele at het SNP positions')
    parser.add_argument('--bam', required=True, help='Input WASP-filtered BAM')
    parser.add_argument('--snp_dir', required=True, help='WASP SNP directory')
    parser.add_argument('--output_prefix', required=True,
                        help='Output prefix (e.g., /path/to/sample)')
    args = parser.parse_args()

    print("=" * 60)
    print("split_by_allele.py")
    print("=" * 60)
    print(f"  Input BAM:     {args.bam}")
    print(f"  SNP dir:       {args.snp_dir}")
    print(f"  Output prefix: {args.output_prefix}")
    print()

    # --- Load SNPs ---
    print("Loading SNPs...")
    snps = load_snps(args.snp_dir)
    print()

    # --- Open input/output BAMs ---
    print("Opening BAM files...")
    bam_in = pysam.AlignmentFile(args.bam, 'rb')

    out_dir = os.path.dirname(args.output_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    bam_ref = pysam.AlignmentFile(f"{args.output_prefix}_ref.bam", 'wb',
                                   template=bam_in)
    bam_alt = pysam.AlignmentFile(f"{args.output_prefix}_alt.bam", 'wb',
                                   template=bam_in)
    bam_nosnp = pysam.AlignmentFile(f"{args.output_prefix}_nosnp.bam", 'wb',
                                     template=bam_in)
    bam_ambig = pysam.AlignmentFile(f"{args.output_prefix}_ambiguous.bam", 'wb',
                                     template=bam_in)

    # --- Per-SNP allele counters ---
    # {(chrom, pos): {'ref': count, 'alt': count, 'other': count}}
    snp_counts = defaultdict(lambda: {'ref': 0, 'alt': 0, 'other': 0})

    # --- Process reads ---
    print("Processing reads...")
    counts = {'ref': 0, 'alt': 0, 'nosnp': 0, 'ambiguous': 0, 'total': 0}
    report_interval = 5_000_000

    for read in bam_in.fetch(until_eof=True):
        # Skip unmapped, secondary, supplementary
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        counts['total'] += 1
        if counts['total'] % report_interval == 0:
            print(f"  Processed {counts['total']:,} reads...")

        chrom = read.reference_name
        chrom_snps = snps.get(chrom, {})

        allele_hits = get_read_alleles(read, chrom_snps)
        classification = classify_read(allele_hits)

        # Update per-SNP counts
        for pos, ref_allele, alt_allele, read_base, assignment in allele_hits:
            snp_counts[(chrom, pos)][assignment] += 1

        # Write to appropriate output BAM
        counts[classification] += 1
        if classification == 'ref':
            bam_ref.write(read)
        elif classification == 'alt':
            bam_alt.write(read)
        elif classification == 'nosnp':
            bam_nosnp.write(read)
        else:
            bam_ambig.write(read)

    # --- Close BAMs ---
    bam_in.close()
    bam_ref.close()
    bam_alt.close()
    bam_nosnp.close()
    bam_ambig.close()

    # --- Write per-SNP allele counts ---
    counts_file = f"{args.output_prefix}_allele_counts.tsv"
    print(f"\nWriting per-SNP allele counts to {counts_file}...")

    with open(counts_file, 'w') as f:
        f.write("chrom\tpos\tref_allele\talt_allele\tref_count\talt_count\tother_count\ttotal_count\tref_fraction\n")
        for (chrom, pos) in sorted(snp_counts.keys(),
                                     key=lambda x: (x[0], x[1])):
            c = snp_counts[(chrom, pos)]
            ref_allele, alt_allele = snps[chrom][pos]
            total = c['ref'] + c['alt'] + c['other']
            ref_frac = c['ref'] / total if total > 0 else 'NA'
            if isinstance(ref_frac, float):
                ref_frac = f"{ref_frac:.4f}"
            f.write(f"{chrom}\t{pos}\t{ref_allele}\t{alt_allele}\t"
                    f"{c['ref']}\t{c['alt']}\t{c['other']}\t{total}\t{ref_frac}\n")

    covered_snps = sum(1 for c in snp_counts.values()
                       if c['ref'] + c['alt'] > 0)

    # --- Summary ---
    print()
    print("=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"  Total reads processed:  {counts['total']:,}")
    print(f"  Ref allele (129S1):     {counts['ref']:,} ({100*counts['ref']/max(1,counts['total']):.2f}%)")
    print(f"  Alt allele (CAST):      {counts['alt']:,} ({100*counts['alt']/max(1,counts['total']):.2f}%)")
    print(f"  No SNP overlap:         {counts['nosnp']:,} ({100*counts['nosnp']/max(1,counts['total']):.2f}%)")
    print(f"  Ambiguous:              {counts['ambiguous']:,} ({100*counts['ambiguous']/max(1,counts['total']):.2f}%)")
    print()
    print(f"  Het SNPs with coverage: {covered_snps:,} / {sum(len(v) for v in snps.values()):,}")

    ref_total = counts['ref']
    alt_total = counts['alt']
    allelic_total = ref_total + alt_total
    if allelic_total > 0:
        ref_pct = 100 * ref_total / allelic_total
        print(f"  Global ref fraction:    {ref_pct:.2f}% (expect ~50% after WASP)")
        if abs(ref_pct - 50) > 5:
            print(f"  [WARNING] Ref fraction deviates from 50% — check WASP filtering")
        else:
            print(f"  [OK] Ref fraction close to 50% — WASP filtering looks correct")

    print()
    print("Output files:")
    for suffix in ['_ref.bam', '_alt.bam', '_nosnp.bam', '_ambiguous.bam',
                    '_allele_counts.tsv']:
        fpath = f"{args.output_prefix}{suffix}"
        if os.path.exists(fpath):
            size_mb = os.path.getsize(fpath) / 1e6
            print(f"  {fpath} ({size_mb:.1f} MB)")


if __name__ == '__main__':
    main()
