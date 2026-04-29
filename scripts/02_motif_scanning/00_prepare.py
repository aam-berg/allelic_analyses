#!/usr/bin/env python3
"""
00_prepare.py — Compute genome background frequencies and create chromosome manifest.

Computes nucleotide frequencies from the genome FASTA (excluding N's and
non-standard chromosomes), and writes:
  1. background.txt  — "A 0.XXXX C 0.XXXX G 0.XXXX T 0.XXXX"
  2. chromosomes.txt  — list of chromosome names to scan (one per line)

Usage:
    python 00_prepare.py <genome.fa> <outdir>
"""

import sys
import os
from collections import Counter


def compute_background(fasta_path, outdir, include_pattern=None, exclude_pattern=None):
    """
    Compute nucleotide frequencies from a FASTA file.
    
    By default, includes only standard chromosomes (chr1-19, chrX, chrY)
    for the background calculation, but writes ALL chromosomes to the manifest
    for scanning.
    """
    counts = Counter()
    chromosomes_for_bg = []
    all_chromosomes = []
    current_chrom = None
    current_is_bg = False

    # Standard mouse chromosomes for background calculation
    standard_chroms = set(
        [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY"]
    )

    print(f"Reading genome: {fasta_path}", file=sys.stderr)

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_chrom = line[1:].split()[0]
                all_chromosomes.append(current_chrom)
                # Use only standard chromosomes for background freq
                current_is_bg = current_chrom in standard_chroms
                if current_is_bg:
                    chromosomes_for_bg.append(current_chrom)
                    print(f"  Background chrom: {current_chrom}", file=sys.stderr)
            elif current_is_bg:
                for base in line.upper():
                    if base in "ACGT":
                        counts[base] += 1

    total = sum(counts.values())
    freqs = {base: counts[base] / total for base in "ACGT"}

    print(f"\nNucleotide frequencies (from {total:,} bases):", file=sys.stderr)
    for base in "ACGT":
        print(f"  {base}: {freqs[base]:.6f}", file=sys.stderr)

    # Write background file
    os.makedirs(outdir, exist_ok=True)
    bg_path = os.path.join(outdir, "background.txt")
    with open(bg_path, "w") as f:
        f.write(" ".join(f"{base} {freqs[base]:.6f}" for base in "ACGT") + "\n")
    print(f"\nWrote background to: {bg_path}", file=sys.stderr)

    # Write chromosome manifest (standard chroms only for scanning)
    chrom_path = os.path.join(outdir, "chromosomes.txt")
    with open(chrom_path, "w") as f:
        for chrom in sorted(standard_chroms, key=lambda x: (
            # Sort numerically then X, Y
            int(x[3:]) if x[3:].isdigit() else (100 if x == "chrX" else 101)
        )):
            if chrom in all_chromosomes:
                f.write(chrom + "\n")
    print(f"Wrote chromosome list to: {chrom_path}", file=sys.stderr)

    # Also write a full chromosome list (including scaffolds) in case needed
    all_chrom_path = os.path.join(outdir, "all_chromosomes.txt")
    with open(all_chrom_path, "w") as f:
        for chrom in all_chromosomes:
            f.write(chrom + "\n")
    print(f"Wrote full chromosome list to: {all_chrom_path}", file=sys.stderr)

    return freqs, chromosomes_for_bg


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <genome.fa> <outdir>", file=sys.stderr)
        sys.exit(1)

    fasta_path = sys.argv[1]
    outdir = sys.argv[2]

    compute_background(fasta_path, outdir)
    