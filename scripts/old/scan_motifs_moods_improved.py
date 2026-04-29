#!/usr/bin/env python3
"""
Genome-wide Motif Scanning using MOODS (CORRECTED VERSION)

CRITICAL FIX: The MOODS scan_dna API is:
    MOODS.scan.scan_dna(seq, matrices, bg, thresholds)
    
NOT:
    MOODS.scan.scan_dna(seq, matrices, thresholds, bg)  # WRONG!

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
from collections import Counter
import logging

from Bio import SeqIO
from tqdm import tqdm

import MOODS.scan
import MOODS.tools
import MOODS.parsers

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def parse_meme_file(meme_path):
    """Parse MEME file and extract PWMs."""
    with open(meme_path, 'r') as f:
        content = f.read()
    
    # Try to parse background from MEME file
    background = None
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
    
    # Parse motifs
    parts = content.split('\nMOTIF ')
    motifs = {}
    
    for part in parts[1:]:
        lines = part.strip().split('\n')
        motif_id = lines[0].split()[0]
        
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
    Convert PWM to MOODS log-odds format.
    
    MOODS format: [[A_scores], [C_scores], [G_scores], [T_scores]]
    where each inner list has length = motif width
    """
    matrix = [[], [], [], []]
    for pos in pwm:
        for base_idx in range(4):
            prob = max(pos[base_idx], pseudocount)
            bg = max(background[base_idx], pseudocount)
            matrix[base_idx].append(math.log(prob / bg))
    return matrix


def read_background_file(filepath):
    """Read background frequencies from a file."""
    bg_dict = {}
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            
            # Parse "A 0.25 C 0.25 G 0.25 T 0.25" format
            if len(parts) >= 8:
                for i in range(0, len(parts), 2):
                    if parts[i] in 'ACGT':
                        bg_dict[parts[i]] = float(parts[i+1])
            # Parse "A 0.25" format
            elif len(parts) == 2 and parts[0] in 'ACGT':
                bg_dict[parts[0]] = float(parts[1])
    
    if not all(b in bg_dict for b in 'ACGT'):
        raise ValueError(f"Background file must contain A, C, G, T. Found: {list(bg_dict.keys())}")
    
    return [bg_dict[b] for b in 'ACGT']


def compute_background_from_genome(genome_path):
    """Compute nucleotide frequencies from genome."""
    logger.info("Computing genome-wide background frequencies...")
    counts = Counter()
    
    for record in SeqIO.parse(genome_path, "fasta"):
        chrom = record.id
        if '_' in chrom or 'random' in chrom.lower() or 'un' in chrom.lower():
            continue
        seq = str(record.seq).upper()
        counts.update(c for c in seq if c in 'ACGT')
    
    total = sum(counts.values())
    if total == 0:
        logger.warning("No valid bases found, using uniform background")
        return [0.25, 0.25, 0.25, 0.25]
    
    bg = [counts[b]/total for b in 'ACGT']
    logger.info(f"  Background: A={bg[0]:.4f} C={bg[1]:.4f} G={bg[2]:.4f} T={bg[3]:.4f}")
    logger.info(f"  GC content: {(bg[1]+bg[2])*100:.1f}%")
    return bg


def get_cluster_motifs(all_motifs, metadata_path, cluster_id):
    """Get motifs belonging to a specific cluster."""
    cluster_motifs = {}
    
    if metadata_path and os.path.exists(metadata_path):
        motif_to_cluster = {}
        with open(metadata_path, 'r') as f:
            f.readline()
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
        for motif_id, data in all_motifs.items():
            if (motif_id == cluster_id or
                motif_id.startswith(cluster_id + ':') or
                motif_id.startswith(cluster_id + '_')):
                cluster_motifs[motif_id] = data
    
    return cluster_motifs


def scan_genome(genome_path, motifs_dict, background, pvalue, cluster_id, output_dir):
    """
    Scan genome for motif occurrences.
    
    CRITICAL: MOODS API is scan_dna(seq, matrices, bg, thresholds)
    """
    if not motifs_dict:
        logger.warning(f"No motifs to scan for cluster {cluster_id}")
        return 0
    
    # Convert motifs to MOODS format and calculate thresholds
    matrices = []
    thresholds = []
    motif_ids = []
    motif_widths = []
    
    for motif_id, data in motifs_dict.items():
        matrix = pwm_to_log_odds(data['pwm'], background)
        matrices.append(matrix)
        
        # Calculate threshold for this p-value
        threshold = MOODS.tools.threshold_from_p(matrix, background, pvalue)
        thresholds.append(threshold)
        
        motif_ids.append(motif_id)
        motif_widths.append(data['width'])
        
        logger.info(f"  Motif {motif_id}: width={data['width']}, threshold={threshold:.2f}")
    
    # Add reverse complements for both-strand scanning
    rc_matrices = [MOODS.tools.reverse_complement(m) for m in matrices]
    all_matrices = matrices + rc_matrices
    all_thresholds = thresholds + thresholds
    
    logger.info(f"  Scanning with {len(matrices)} motif(s) on both strands")
    logger.info(f"  P-value threshold: {pvalue}")
    
    # Prepare output
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = output_dir / f"{cluster_id}.tsv"
    bed_path = output_dir / f"{cluster_id}.bed"
    
    total_hits = 0
    
    with open(tsv_path, 'w') as tsv_out, open(bed_path, 'w') as bed_out:
        tsv_out.write("cluster\tmotif_id\tchr\tstart\tend\tstrand\tscore\n")
        
        for record in tqdm(SeqIO.parse(genome_path, "fasta"), desc="  Scanning"):
            chrom = record.id
            
            if '_' in chrom or 'random' in chrom.lower() or 'un' in chrom.lower():
                continue
            
            seq = str(record.seq).upper()
            seq_len = len(seq)
            
            # ============================================================
            # CRITICAL FIX: Correct argument order!
            # MOODS.scan.scan_dna(seq, matrices, bg, thresholds)
            # ============================================================
            results = MOODS.scan.scan_dna(seq, all_matrices, background, all_thresholds)
            
            # Process results
            n_matrices = len(matrices)
            for i, hits in enumerate(results):
                if i < n_matrices:
                    motif_id = motif_ids[i]
                    width = motif_widths[i]
                    strand = '+'
                else:
                    motif_id = motif_ids[i - n_matrices]
                    width = motif_widths[i - n_matrices]
                    strand = '-'
                
                for hit in hits:
                    start = hit.pos
                    end = start + width
                    score = hit.score
                    
                    if start < 0 or end > seq_len:
                        continue
                    
                    total_hits += 1
                    
                    tsv_out.write(f"{cluster_id}\t{motif_id}\t{chrom}\t{start}\t{end}\t{strand}\t{score:.4f}\n")
                    bed_score = min(1000, max(0, int(score * 100)))
                    bed_out.write(f"{chrom}\t{start}\t{end}\t{cluster_id}:{motif_id}\t{bed_score}\t{strand}\n")
    
    # Compress large files
    if total_hits > 100000:
        logger.info(f"  Compressing outputs ({total_hits:,} hits)...")
        os.system(f"gzip -f {tsv_path}")
        os.system(f"gzip -f {bed_path}")
    
    logger.info(f"  Found {total_hits:,} motif hits for {cluster_id}")
    return total_hits


def main():
    parser = argparse.ArgumentParser(
        description='Scan genome for motifs using MOODS',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script correctly uses the MOODS API:
    MOODS.scan.scan_dna(seq, matrices, bg, thresholds)

Background options:
  --background-file FILE   Pre-computed background file (recommended)
  --bg A C G T             Specify background directly

Example:
    python scan_motifs_moods.py --meme motifs.meme --genome mm39.fa \\
        --output results --pvalue 1e-4 --bg 0.2915 0.2085 0.2085 0.2915 \\
        --cluster AC0001
        """
    )
    parser.add_argument('--meme', required=True, help='MEME file with PWMs')
    parser.add_argument('--genome', required=True, help='Reference genome FASTA')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--pvalue', type=float, default=1e-4, help='P-value threshold (default: 1e-4)')
    parser.add_argument('--background-file', help='Pre-computed background file')
    parser.add_argument('--bg', nargs=4, type=float, metavar=('A', 'C', 'G', 'T'),
                        help='Background frequencies as 4 values')
    parser.add_argument('--metadata', help='Metadata TSV with cluster assignments')
    parser.add_argument('--cluster', help='Single cluster ID to scan')
    parser.add_argument('--manifest', help='Manifest file for array jobs')
    parser.add_argument('--array-index', type=int, help='Array job index (1-based)')
    
    args = parser.parse_args()
    
    # Determine clusters
    if args.cluster:
        clusters = [args.cluster]
    elif args.manifest and args.array_index:
        with open(args.manifest) as f:
            all_clusters = [line.strip() for line in f if line.strip()]
        if args.array_index < 1 or args.array_index > len(all_clusters):
            logger.error(f"Array index {args.array_index} out of range")
            sys.exit(1)
        clusters = [all_clusters[args.array_index - 1]]
    else:
        logger.error("Specify --cluster or --manifest with --array-index")
        sys.exit(1)
    
    # Parse MEME file
    logger.info(f"Parsing MEME file: {args.meme}")
    all_motifs, meme_background = parse_meme_file(args.meme)
    logger.info(f"  Found {len(all_motifs)} motifs")
    
    # Determine background
    if args.bg:
        background = args.bg
        logger.info(f"Using command-line background: {background}")
    elif args.background_file:
        background = read_background_file(args.background_file)
        logger.info(f"Loaded background from {args.background_file}: {background}")
    elif meme_background:
        background = meme_background
        logger.info(f"Using background from MEME file: {background}")
    else:
        background = compute_background_from_genome(args.genome)
    
    logger.info(f"Background: A={background[0]:.4f} C={background[1]:.4f} G={background[2]:.4f} T={background[3]:.4f}")
    logger.info(f"GC content: {(background[1]+background[2])*100:.1f}%")
    
    # Output directory
    pval_suffix = f"em{abs(int(math.log10(args.pvalue)))}"
    output_base = Path(args.output) / f"motif_archetypes_{pval_suffix}" / "motif_locations"
    
    # Process clusters
    for cluster_id in clusters:
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing: {cluster_id}")
        logger.info(f"{'='*60}")
        
        cluster_motifs = get_cluster_motifs(all_motifs, args.metadata, cluster_id)
        
        if not cluster_motifs and cluster_id in all_motifs:
            cluster_motifs = {cluster_id: all_motifs[cluster_id]}
        
        if not cluster_motifs:
            logger.warning(f"No motifs found for {cluster_id}, skipping")
            continue
        
        logger.info(f"  Found {len(cluster_motifs)} motif(s)")
        
        scan_genome(
            genome_path=args.genome,
            motifs_dict=cluster_motifs,
            background=background,
            pvalue=args.pvalue,
            cluster_id=cluster_id,
            output_dir=output_base
        )
    
    logger.info(f"\n{'='*60}")
    logger.info("Done!")
    logger.info(f"{'='*60}")


if __name__ == '__main__':
    main()