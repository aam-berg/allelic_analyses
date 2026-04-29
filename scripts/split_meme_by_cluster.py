#!/usr/bin/env python3
"""
Split a MEME-format motif file by cluster (archetype).

This script:
1. Reads the consensus_pwms.meme file
2. Reads the metadata.tsv to get cluster assignments
3. Creates one MEME file per cluster containing all motifs in that cluster
4. Creates a manifest file listing all clusters for array job submission

Usage:
    python split_meme_by_cluster.py \
        --meme /path/to/consensus_pwms.meme \
        --metadata /path/to/metadata.tsv \
        --outdir /path/to/output

Output:
    outdir/
        cluster_meme_files/
            AC0001.meme
            AC0002.meme
            ...
        cluster_manifest.txt  # List of cluster IDs for array jobs
"""

import argparse
import os
import sys
from collections import defaultdict
from pathlib import Path


def parse_meme_file(meme_path):
    """
    Parse a MEME file and return header and motif blocks.
    
    Returns:
        header: str - The MEME file header (everything before first MOTIF)
        motifs: dict - {motif_id: motif_block_string}
    """
    with open(meme_path, 'r') as f:
        content = f.read()
    
    # Split by MOTIF keyword
    parts = content.split('\nMOTIF ')
    
    # First part is the header
    header = parts[0]
    
    # Remaining parts are motifs
    motifs = {}
    for part in parts[1:]:
        lines = part.strip().split('\n')
        # First line contains the motif ID
        # Format: "AC0001:DLX/LHX:Homeodomain AC0001:DLX/LHX:Homeodomain"
        # or just "MOTIF_ID"
        first_line = lines[0]
        motif_id = first_line.split()[0]  # Take first whitespace-delimited token
        
        # Store the full motif block (with MOTIF prefix restored)
        motif_block = 'MOTIF ' + part
        motifs[motif_id] = motif_block
    
    return header, motifs


def load_metadata(metadata_path):
    """
    Load metadata.tsv and return motif_id -> cluster mapping.
    
    Returns:
        motif_to_cluster: dict - {motif_id: cluster_id}
        cluster_to_motifs: dict - {cluster_id: [motif_ids]}
    """
    motif_to_cluster = {}
    cluster_to_motifs = defaultdict(list)
    
    with open(metadata_path, 'r') as f:
        header = f.readline()  # Skip header
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 2:
                motif_id = fields[0]
                cluster = fields[1]
                motif_to_cluster[motif_id] = cluster
                cluster_to_motifs[cluster].append(motif_id)
    
    return motif_to_cluster, dict(cluster_to_motifs)


def extract_motif_id_from_meme(motif_block):
    """
    Extract the motif ID from a MEME motif block.
    
    The MEME file uses format like "AC0001:DLX/LHX:Homeodomain"
    but metadata might use different IDs. We need to match them.
    """
    lines = motif_block.strip().split('\n')
    first_line = lines[0].replace('MOTIF ', '')
    return first_line.split()[0]


def write_meme_file(header, motif_blocks, outpath):
    """Write a MEME file with the given header and motif blocks."""
    with open(outpath, 'w') as f:
        f.write(header)
        f.write('\n')
        for block in motif_blocks:
            f.write(block)
            if not block.endswith('\n'):
                f.write('\n')
            f.write('\n')


def main():
    parser = argparse.ArgumentParser(
        description='Split MEME file by cluster for parallel FIMO runs'
    )
    parser.add_argument('--meme', required=True, help='Input MEME file')
    parser.add_argument('--metadata', required=True, help='Metadata TSV with cluster assignments')
    parser.add_argument('--outdir', required=True, help='Output directory')
    args = parser.parse_args()
    
    # Create output directories
    outdir = Path(args.outdir)
    meme_dir = outdir / 'cluster_meme_files'
    meme_dir.mkdir(parents=True, exist_ok=True)
    
    # Parse inputs
    print(f"Parsing MEME file: {args.meme}")
    header, motifs = parse_meme_file(args.meme)
    print(f"  Found {len(motifs)} motifs")
    
    print(f"Loading metadata: {args.metadata}")
    motif_to_cluster, cluster_to_motifs = load_metadata(args.metadata)
    print(f"  Found {len(cluster_to_motifs)} clusters")
    
    # Match motifs to clusters
    # The MEME file motif IDs might be in format "AC0001:DLX/LHX:Homeodomain"
    # while metadata uses simpler IDs. We need to handle both cases.
    
    # First, try direct matching
    # If that fails, try matching the cluster ID prefix
    
    cluster_motif_blocks = defaultdict(list)
    unmatched_motifs = []
    
    for meme_motif_id, motif_block in motifs.items():
        matched = False
        
        # Try direct match first
        if meme_motif_id in motif_to_cluster:
            cluster = motif_to_cluster[meme_motif_id]
            cluster_motif_blocks[cluster].append(motif_block)
            matched = True
        else:
            # Try matching by cluster ID prefix (e.g., "AC0001:..." -> "AC0001")
            # The MEME motif ID format is "AC0001:DLX/LHX:Homeodomain"
            if ':' in meme_motif_id:
                cluster_prefix = meme_motif_id.split(':')[0]
                if cluster_prefix in cluster_to_motifs:
                    cluster_motif_blocks[cluster_prefix].append(motif_block)
                    matched = True
        
        if not matched:
            unmatched_motifs.append(meme_motif_id)
    
    if unmatched_motifs:
        print(f"  Warning: {len(unmatched_motifs)} motifs could not be matched to clusters")
        print(f"  First 5 unmatched: {unmatched_motifs[:5]}")
    
    # Write cluster-specific MEME files
    print(f"\nWriting cluster MEME files to {meme_dir}")
    cluster_list = sorted(cluster_motif_blocks.keys())
    
    for cluster in cluster_list:
        motif_blocks = cluster_motif_blocks[cluster]
        outpath = meme_dir / f'{cluster}.meme'
        write_meme_file(header, motif_blocks, outpath)
        print(f"  {cluster}: {len(motif_blocks)} motifs")
    
    # Write manifest file for array jobs
    manifest_path = outdir / 'cluster_manifest.txt'
    with open(manifest_path, 'w') as f:
        for cluster in cluster_list:
            f.write(f'{cluster}\n')
    
    print(f"\nManifest written to {manifest_path}")
    print(f"Total clusters: {len(cluster_list)}")
    
    # Write summary statistics
    summary_path = outdir / 'cluster_summary.txt'
    with open(summary_path, 'w') as f:
        f.write("cluster\tn_motifs\n")
        for cluster in cluster_list:
            f.write(f"{cluster}\t{len(cluster_motif_blocks[cluster])}\n")
    
    print(f"Summary written to {summary_path}")
    print("\nDone!")


if __name__ == '__main__':
    main()
