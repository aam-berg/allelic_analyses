#!/usr/bin/env python3
"""
03_validate_promoter_enrichment.py — Validate MOODS results by checking
enrichment of motif hits at promoter regions vs. random genomic background.

This is one of the most informative sanity checks: real TF motif hits should
be enriched near transcription start sites.

Requires bedtools in PATH.

Usage:
    python 03_validate_promoter_enrichment.py \
        --merged-dir /n/scratch/users/a/alb1273/pausing_phase_project_intermediates/moods_scan/merged \
        --gene-annotation /n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/annotation/gencode.vM38.chr_patch_hapl_scaff.annotation.gtf.gz \
        --chrom-sizes /n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/genome/mm39.chrom.sizes \
        --output promoter_enrichment.tsv
"""

import argparse
import os
import sys
import subprocess
import tempfile
import glob
from collections import defaultdict


def make_promoter_bed(gtf_path, chrom_sizes_path, outpath, upstream=1000, downstream=500):
    """
    Extract promoter regions (TSS ± window) from a GTF file.
    Uses only protein-coding genes on standard chromosomes.
    """
    print(f"Extracting promoters from {gtf_path}...", file=sys.stderr)

    # Read chrom sizes for boundary checking
    chrom_sizes = {}
    with open(chrom_sizes_path) as f:
        for line in f:
            chrom, size = line.strip().split("\t")
            chrom_sizes[chrom] = int(size)

    seen_genes = set()
    promoters = []

    # Handle gzipped GTF
    import gzip
    opener = gzip.open if gtf_path.endswith(".gz") else open

    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            if fields[2] != "gene":
                continue

            # Parse attributes
            attrs = {}
            for attr in fields[8].split(";"):
                attr = attr.strip()
                if " " in attr:
                    key, val = attr.split(" ", 1)
                    attrs[key] = val.strip('"')

            gene_name = attrs.get("gene_name", attrs.get("gene_id", ""))
            gene_type = attrs.get("gene_type", attrs.get("gene_biotype", ""))

            # Filter: protein-coding, standard chroms, unique genes
            chrom = fields[0]
            if chrom not in chrom_sizes:
                continue
            if gene_type != "protein_coding":
                continue
            if gene_name in seen_genes:
                continue
            seen_genes.add(gene_name)

            strand = fields[6]
            if strand == "+":
                tss = int(fields[3]) - 1  # GTF is 1-based
                start = max(0, tss - upstream)
                end = min(chrom_sizes[chrom], tss + downstream)
            else:
                tss = int(fields[4])  # GTF end is 1-based inclusive
                start = max(0, tss - downstream)
                end = min(chrom_sizes[chrom], tss + upstream)

            promoters.append((chrom, start, end, gene_name, "0", strand))

    with open(outpath, "w") as f:
        for p in sorted(promoters, key=lambda x: (x[0], x[1])):
            f.write("\t".join(str(x) for x in p) + "\n")

    print(f"  Extracted {len(promoters)} promoters", file=sys.stderr)
    return outpath


def compute_enrichment(motif_bed_dir, promoter_bed, chrom_sizes_path, genome_size):
    """
    For each motif, compute enrichment of hits at promoters vs. genome-wide rate.

    Enrichment = (hits_in_promoters / promoter_bp) / (total_hits / genome_bp)
    """
    results = {}

    # Get total promoter bp
    total_promoter_bp = 0
    with open(promoter_bed) as f:
        for line in f:
            fields = line.strip().split("\t")
            total_promoter_bp += int(fields[2]) - int(fields[1])

    print(f"Total promoter bp: {total_promoter_bp:,}", file=sys.stderr)

    motif_files = sorted(glob.glob(os.path.join(motif_bed_dir, "*.bed")))
    print(f"Checking {len(motif_files)} motif files...", file=sys.stderr)

    for motif_file in motif_files:
        motif_id = os.path.basename(motif_file).replace(".bed", "")

        # Count total hits
        total_hits = 0
        with open(motif_file) as f:
            for line in f:
                total_hits += 1

        if total_hits == 0:
            continue

        # Count hits overlapping promoters using bedtools
        try:
            result = subprocess.run(
                ["bedtools", "intersect", "-a", motif_file, "-b", promoter_bed, "-u", "-sorted"],
                capture_output=True, text=True, check=True
            )
            hits_in_promoters = len(result.stdout.strip().split("\n")) if result.stdout.strip() else 0
        except subprocess.CalledProcessError:
            # If sorted intersect fails (files might not be identically sorted), try without -sorted
            result = subprocess.run(
                ["bedtools", "intersect", "-a", motif_file, "-b", promoter_bed, "-u"],
                capture_output=True, text=True
            )
            hits_in_promoters = len(result.stdout.strip().split("\n")) if result.stdout.strip() else 0

        # Compute enrichment
        genome_rate = total_hits / genome_size
        promoter_rate = hits_in_promoters / total_promoter_bp if total_promoter_bp > 0 else 0
        enrichment = promoter_rate / genome_rate if genome_rate > 0 else 0

        results[motif_id] = {
            "total_hits": total_hits,
            "promoter_hits": hits_in_promoters,
            "promoter_fraction": hits_in_promoters / total_hits if total_hits > 0 else 0,
            "enrichment": enrichment,
        }

    return results


def main():
    parser = argparse.ArgumentParser(description="Validate motif hits by promoter enrichment")
    parser.add_argument("--merged-dir", required=True,
                        help="Directory with merged MOODS results (contains per_motif/)")
    parser.add_argument("--gene-annotation", required=True,
                        help="GENCODE GTF file (can be gzipped)")
    parser.add_argument("--chrom-sizes", required=True,
                        help="Chromosome sizes file (chrom<tab>size)")
    parser.add_argument("--genome-size", type=int, default=2652783500)
    parser.add_argument("--output", required=True, help="Output TSV")
    args = parser.parse_args()

    motif_bed_dir = os.path.join(args.merged_dir, "per_motif")
    if not os.path.isdir(motif_bed_dir):
        print(f"ERROR: per_motif directory not found: {motif_bed_dir}", file=sys.stderr)
        sys.exit(1)

    # Create promoter BED
    promoter_bed = os.path.join(args.merged_dir, "promoters_1kb.bed")
    if not os.path.exists(promoter_bed):
        make_promoter_bed(args.gene_annotation, args.chrom_sizes, promoter_bed)

    # Compute enrichment
    results = compute_enrichment(motif_bed_dir, promoter_bed, args.chrom_sizes, args.genome_size)

    # Write output
    with open(args.output, "w") as f:
        f.write("motif_id\ttotal_hits\tpromoter_hits\tpromoter_fraction\tenrichment_at_promoters\n")
        for motif_id in sorted(results.keys()):
            r = results[motif_id]
            f.write(f"{motif_id}\t{r['total_hits']}\t{r['promoter_hits']}\t"
                    f"{r['promoter_fraction']:.6f}\t{r['enrichment']:.4f}\n")

    print(f"\nWrote enrichment results to {args.output}", file=sys.stderr)

    # Quick summary
    enrichments = [r["enrichment"] for r in results.values() if r["enrichment"] > 0]
    if enrichments:
        enrichments.sort()
        print(f"Enrichment at promoters:", file=sys.stderr)
        print(f"  Median: {enrichments[len(enrichments)//2]:.2f}x", file=sys.stderr)
        print(f"  Range:  {enrichments[0]:.2f}x - {enrichments[-1]:.2f}x", file=sys.stderr)
        n_enriched = sum(1 for e in enrichments if e > 1.5)
        print(f"  Motifs with >1.5x enrichment: {n_enriched}/{len(enrichments)}", file=sys.stderr)
        print(f"\n  (Expect most motifs to show >1x enrichment at promoters)", file=sys.stderr)


if __name__ == "__main__":
    main()
