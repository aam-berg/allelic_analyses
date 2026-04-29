#!/usr/bin/env python3
"""
Genome-wide Motif Scanning using MOODS

Scans the mouse genome for motif occurrences using MOODS.

Usage:
    python scan_motifs_moods.py \
        --meme /path/to/consensus_pwms.meme \
        --genome /path/to/mm39.fa \
        --output /path/to/output_dir \
        --pvalue 1e-5 \
        --cluster AC0001
"""

import argparse
import os
import sys
from pathlib import Path
import math

from Bio import SeqIO
from tqdm import tqdm

import MOODS.scan
import MOODS.tools


def parse_meme_file(meme_path):
    """
    Parse a MEME-format file and extract PWMs.
    
    Returns:
        motifs: dict - {motif_id: {'pwm': list of lists, 'width': int}}
        background: list - Background frequencies [A, C, G, T]
    """
    with open(meme_path, 'r') as f:
        content = f.read()
    
    # Parse background frequencies (default to uniform)
    background = [0.25, 0.25, 0.25, 0.25]
    lines = content.split('\n')
    for i, line in enumerate(lines):
        if 'A ' in line and 'C ' in line and 'G ' in line and 'T ' in line:
            parts = line.split()
            try:
                bg_dict = {}
                for j in range(0, len(parts), 2):
                    bg_dict[parts[j]] = float(parts[j+1])
                if all(k in bg_dict for k in ['A', 'C', 'G', 'T']):
                    background = [bg_dict['A'], bg_dict['C'], bg_dict['G'], bg_dict['T']]
                    break
            except (ValueError, IndexError):
                pass
    
    # Split by MOTIF keyword and parse each
    parts = content.split('\nMOTIF ')
    motifs = {}
    
    for part in parts[1:]:
        lines = part.strip().split('\n')
        motif_id = lines[0].split()[0]
        
        # Find and parse the probability matrix
        pwm = []
        in_matrix = False
        
        for line in lines[1:]:
            line = line.strip()
            if line.startswith('letter-probability matrix'):
                in_matrix = True
                continue
            if in_matrix:
                if not line or line.startswith('URL') or line.startswith('MOTIF'):
                    break
                try:
                    probs = [float(x) for x in line.split()]
                    if len(probs) == 4:
                        pwm.append(probs)
                except ValueError:
                    break
        
        if pwm:
            motifs[motif_id] = {'pwm': pwm, 'width': len(pwm)}
    
    return motifs, background


def pwm_to_log_odds(pwm, background, pseudocount=0.0001):
    """
    Convert PWM to log-odds matrix in MOODS format.
    
    MOODS format: [[A_scores], [C_scores], [G_scores], [T_scores]]
    """
    matrix = [[], [], [], []]
    
    for pos in pwm:
        for base_idx in range(4):
            prob = max(pos[base_idx], pseudocount)
            bg = max(background[base_idx], pseudocount)
            matrix[base_idx].append(math.log(prob / bg))
    
    return matrix


def get_cluster_motifs(all_motifs, metadata_path, cluster_id):
    """Get motifs belonging to a specific cluster."""
    cluster_motifs = {}
    
    if metadata_path and os.path.exists(metadata_path):
        # Load cluster assignments from metadata
        motif_to_cluster = {}
        with open(metadata_path, 'r') as f:
            f.readline()  # skip header
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 2:
                    motif_to_cluster[fields[0]] = fields[1]
        
        for motif_id, data in all_motifs.items():
            if (motif_to_cluster.get(motif_id) == cluster_id or
                motif_id.startswith(cluster_id + ':') or
                motif_id.startswith(cluster_id + '_')):
                cluster_motifs[motif_id] = data
    else:
        # No metadata - match by prefix or exact ID
        for motif_id, data in all_motifs.items():
            if (motif_id == cluster_id or
                motif_id.startswith(cluster_id + ':') or
                motif_id.startswith(cluster_id + '_')):
                cluster_motifs[motif_id] = data
    
    return cluster_motifs


def scan_genome(genome_path, motifs_dict, background, pvalue, cluster_id, output_dir):
    """
    Scan genome for motif occurrences.
    
    Returns:
        Total number of hits found
    """
    if not motifs_dict:
        print(f"  Warning: No motifs to scan for cluster {cluster_id}")
        return 0
    
    # Convert motifs to MOODS format and calculate thresholds
    matrices = []
    thresholds = []
    motif_ids = []
    motif_widths = []
    
    for motif_id, data in motifs_dict.items():
        matrix = pwm_to_log_odds(data['pwm'], background)
        matrices.append(matrix)
        thresholds.append(MOODS.tools.threshold_from_p(matrix, background, pvalue))
        motif_ids.append(motif_id)
        motif_widths.append(data['width'])
    
    # Also add reverse complement matrices for scanning both strands
    rc_matrices = [MOODS.tools.reverse_complement(m) for m in matrices]
    all_matrices = matrices + rc_matrices
    all_thresholds = thresholds + thresholds  # same thresholds for RC
    
    print(f"  Scanning with {len(matrices)} motif(s) (+ reverse complements)")
    print(f"  P-value threshold: {pvalue}")
    
    # Setup output
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = output_dir / f"{cluster_id}.tsv"
    bed_path = output_dir / f"{cluster_id}.bed"
    
    total_hits = 0
    
    with open(tsv_path, 'w') as tsv_out, open(bed_path, 'w') as bed_out:
        tsv_out.write("cluster\tmotif_id\tchr\tstart\tend\tstrand\tscore\n")
        
        # Process each chromosome
        for record in tqdm(SeqIO.parse(genome_path, "fasta"), desc="  Scanning"):
            chrom = record.id
            
            # Skip non-standard chromosomes
            if '_' in chrom or 'random' in chrom.lower() or 'un' in chrom.lower():
                continue
            
            seq = str(record.seq).upper()
            seq_len = len(seq)
            
            # Scan with all matrices (forward + RC)
            results = MOODS.scan.scan_dna(seq, all_matrices, all_thresholds, background)
            
            # Process results
            n_matrices = len(matrices)
            for i, hits in enumerate(results):
                if i < n_matrices:
                    # Forward strand hit
                    motif_id = motif_ids[i]
                    width = motif_widths[i]
                    strand = '+'
                else:
                    # Reverse complement hit
                    motif_id = motif_ids[i - n_matrices]
                    width = motif_widths[i - n_matrices]
                    strand = '-'
                
                for hit in hits:
                    start = hit.pos
                    end = start + width
                    score = hit.score
                    
                    # Validate position
                    if start < 0 or end > seq_len:
                        continue
                    
                    total_hits += 1
                    
                    # Write TSV
                    tsv_out.write(f"{cluster_id}\t{motif_id}\t{chrom}\t{start}\t{end}\t{strand}\t{score:.4f}\n")
                    
                    # Write BED (score scaled to 0-1000)
                    bed_score = min(1000, max(0, int(score * 100)))
                    bed_out.write(f"{chrom}\t{start}\t{end}\t{cluster_id}:{motif_id}\t{bed_score}\t{strand}\n")
    
    # Compress large files
    if total_hits > 100000:
        print(f"  Compressing outputs ({total_hits:,} hits)...")
        os.system(f"gzip -f {tsv_path}")
        os.system(f"gzip -f {bed_path}")
    
    print(f"  Found {total_hits:,} motif hits")
    return total_hits


def main():
    parser = argparse.ArgumentParser(description='Scan genome for motifs using MOODS')
    parser.add_argument('--meme', required=True, help='MEME file with PWMs')
    parser.add_argument('--genome', required=True, help='Reference genome FASTA')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--pvalue', type=float, default=1e-5, help='P-value threshold')
    parser.add_argument('--metadata', help='Metadata TSV with cluster assignments')
    parser.add_argument('--cluster', help='Single cluster ID to scan')
    parser.add_argument('--manifest', help='Manifest file for array jobs')
    parser.add_argument('--array-index', type=int, help='Array job index (1-based)')
    
    args = parser.parse_args()
    
    # Determine cluster(s) to process
    if args.cluster:
        clusters = [args.cluster]
    elif args.manifest and args.array_index:
        with open(args.manifest) as f:
            all_clusters = [line.strip() for line in f if line.strip()]
        if args.array_index < 1 or args.array_index > len(all_clusters):
            print(f"Error: Array index {args.array_index} out of range (1-{len(all_clusters)})")
            sys.exit(1)
        clusters = [all_clusters[args.array_index - 1]]
    else:
        print("Error: Specify --cluster or --manifest with --array-index")
        sys.exit(1)
    
    # Parse MEME file
    print(f"Parsing MEME file: {args.meme}")
    all_motifs, background = parse_meme_file(args.meme)
    print(f"  Found {len(all_motifs)} motifs")
    print(f"  Background: A={background[0]:.3f} C={background[1]:.3f} G={background[2]:.3f} T={background[3]:.3f}")
    
    # Output directory structure
    pval_suffix = f"em{abs(int(math.log10(args.pvalue)))}"
    output_base = Path(args.output) / f"motif_archetypes_{pval_suffix}" / "motif_locations"
    
    # Process each cluster
    for cluster_id in clusters:
        print(f"\n{'='*60}")
        print(f"Processing: {cluster_id}")
        print(f"{'='*60}")
        
        # Get motifs for this cluster
        cluster_motifs = get_cluster_motifs(all_motifs, args.metadata, cluster_id)
        
        # Fallback: if cluster_id matches a motif ID directly
        if not cluster_motifs and cluster_id in all_motifs:
            cluster_motifs = {cluster_id: all_motifs[cluster_id]}
        
        if not cluster_motifs:
            print(f"  No motifs found for {cluster_id}, skipping")
            continue
        
        print(f"  Found {len(cluster_motifs)} motif(s)")
        
        # Scan
        scan_genome(
            genome_path=args.genome,
            motifs_dict=cluster_motifs,
            background=background,
            pvalue=args.pvalue,
            cluster_id=cluster_id,
            output_dir=output_base
        )
    
    print(f"\n{'='*60}")
    print("Done!")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()