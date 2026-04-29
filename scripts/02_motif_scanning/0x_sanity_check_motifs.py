#!/usr/bin/env python3
"""
sanity_check_motifs.py — Motif width distribution + orientation/palindromicity analysis.

Run interactively (Jupyter) or as a script. Produces:
  1. Histogram of motif widths with cutoff annotations
  2. Per-archetype orientation similarity matrix (forward vs rev vs comp vs revcomp)
  3. Summary classification: palindromic / near-palindromic / asymmetric

Requires: numpy, matplotlib, pandas
Optional: MOODS (for log-odds comparison; falls back to PFM correlation if unavailable)

Usage:
    python sanity_check_motifs.py \
        --meme /path/to/consensus_pwms.meme \
        --metadata /path/to/metadata.tsv \
        --output-dir /path/to/output
"""

import argparse
import os
import sys
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# ============================================================================
# MEME parsing
# ============================================================================

def parse_meme_file(meme_path):
    """
    Parse MEME file → dict of archetype_id -> {
        'width': int,
        'pfm': list of lists, shape [4][width], order A/C/G/T (nucleotide × position)
    }
    """
    motifs = {}
    current_id = None
    current_matrix = []
    in_matrix = False
    expected_width = 0
    rows_read = 0

    with open(meme_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("MOTIF"):
                if current_id and current_matrix:
                    motifs[current_id] = {
                        "width": len(current_matrix),
                        "pfm": _transpose(current_matrix),
                    }
                parts = line.split()
                current_id = parts[1].split(":")[0] if len(parts) > 1 else None
                current_matrix = []
                in_matrix = False
                rows_read = 0
            elif line.startswith("letter-probability matrix"):
                in_matrix = True
                match = re.search(r"\bw=\s*(\d+)", line)
                expected_width = int(match.group(1)) if match else 0
                rows_read = 0
            elif in_matrix and line and not line.startswith("URL") and not line.startswith("MOTIF"):
                values = line.split()
                if len(values) == 4:
                    try:
                        current_matrix.append([float(v) for v in values])
                        rows_read += 1
                        if expected_width > 0 and rows_read >= expected_width:
                            in_matrix = False
                    except ValueError:
                        in_matrix = False
                else:
                    in_matrix = False

    if current_id and current_matrix:
        motifs[current_id] = {
            "width": len(current_matrix),
            "pfm": _transpose(current_matrix),
        }

    return motifs


def _transpose(matrix):
    """Position × nucleotide → Nucleotide × position."""
    result = [[], [], [], []]
    for row in matrix:
        for nuc in range(4):
            result[nuc].append(row[nuc])
    return result


# ============================================================================
# The four orientations of a PFM
# ============================================================================
#
# Given a PFM with rows [A, C, G, T] and columns = positions 0..w-1:
#
#   Forward (F):            pfm[b][j]           — the PWM as stored
#   Reverse (R):            pfm[b][w-1-j]       — read the motif backwards
#   Complement (C):         pfm[comp(b)][j]     — complement each nucleotide
#   Reverse complement (RC): pfm[comp(b)][w-1-j] — complement + reverse
#
# where comp = {A↔T, C↔G} i.e. index mapping {0↔3, 1↔2}
#
# For TF binding on dsDNA:
#   - F and RC describe the SAME physical binding event (one per strand)
#   - R and C describe the SAME physical binding event (one per strand)
#   - F/RC vs R/C are genuinely DIFFERENT binding orientations
#
# MOODS searches F and RC. A "palindromic" motif has F ≈ RC.
# If F ≈ R (or equivalently, F ≈ C), the motif has a different kind of
# internal symmetry (the motif reads similarly backwards).
#

COMP = {0: 3, 1: 2, 2: 1, 3: 0}  # A↔T, C↔G
NUC_NAMES = ["A", "C", "G", "T"]


def get_forward(pfm):
    """Return the PFM as-is (nucleotide × position)."""
    return pfm


def get_reverse(pfm):
    """Reverse: read positions right-to-left."""
    return [row[::-1] for row in pfm]


def get_complement(pfm):
    """Complement: swap A↔T, C↔G."""
    w = len(pfm[0])
    result = [[0.0] * w for _ in range(4)]
    for b in range(4):
        for j in range(w):
            result[COMP[b]][j] = pfm[b][j]
    return result


def get_reverse_complement(pfm):
    """Reverse complement: complement + reverse."""
    w = len(pfm[0])
    result = [[0.0] * w for _ in range(4)]
    for b in range(4):
        for j in range(w):
            result[COMP[b]][w - 1 - j] = pfm[b][j]
    return result


def pfm_to_array(pfm):
    """Convert PFM (list of lists) to numpy array, shape (4, width)."""
    return np.array(pfm, dtype=float)


def pfm_pearson(pfm_a, pfm_b):
    """
    Pearson correlation between two PFMs (same width).
    Flattens both matrices and computes correlation.
    """
    a = pfm_to_array(pfm_a).flatten()
    b = pfm_to_array(pfm_b).flatten()
    if np.std(a) == 0 or np.std(b) == 0:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def pfm_column_pearson(pfm_a, pfm_b):
    """
    Mean per-column Pearson correlation between two PFMs.
    More interpretable than flattened correlation for motif comparison.
    """
    a = pfm_to_array(pfm_a)  # (4, w)
    b = pfm_to_array(pfm_b)  # (4, w)
    w = a.shape[1]
    cors = []
    for j in range(w):
        col_a = a[:, j]
        col_b = b[:, j]
        if np.std(col_a) == 0 or np.std(col_b) == 0:
            # Uniform or degenerate column
            cors.append(1.0 if np.allclose(col_a, col_b) else 0.0)
        else:
            cors.append(float(np.corrcoef(col_a, col_b)[0, 1]))
    return float(np.mean(cors))


def classify_orientation_similarity(pfm):
    """
    Classify a motif's orientation properties.

    Returns dict with:
        - cor_F_RC: correlation between Forward and Reverse Complement
                    (high = palindromic; MOODS would double-count)
        - cor_F_R:  correlation between Forward and Reverse
                    (high = reads similarly backwards)
        - cor_F_C:  correlation between Forward and Complement
                    (note: cor_F_C == cor_F_R always, since R and C
                     are the same physical site on opposite strands)
        - classification: 'palindromic', 'near-palindromic', or 'asymmetric'
    """
    F = get_forward(pfm)
    R = get_reverse(pfm)
    C = get_complement(pfm)
    RC = get_reverse_complement(pfm)

    cor_F_RC = pfm_column_pearson(F, RC)
    cor_F_R = pfm_column_pearson(F, R)
    cor_F_C = pfm_column_pearson(F, C)
    cor_R_RC = pfm_column_pearson(R, RC)

    # Classification thresholds
    if cor_F_RC > 0.95:
        classification = "palindromic"
    elif cor_F_RC > 0.80:
        classification = "near-palindromic"
    else:
        classification = "asymmetric"

    return {
        "cor_F_RC": cor_F_RC,
        "cor_F_R": cor_F_R,
        "cor_F_C": cor_F_C,
        "cor_R_RC": cor_R_RC,
        "classification": classification,
    }


def consensus_from_pfm(pfm):
    """Get consensus sequence from PFM."""
    arr = pfm_to_array(pfm)
    w = arr.shape[1]
    nucs = "ACGT"
    return "".join(nucs[arr[:, j].argmax()] for j in range(w))


# ============================================================================
# Width distribution plot
# ============================================================================

def plot_width_distribution(motifs, metadata_df, cutoff=8, output_path=None):
    """
    Plot histogram of motif archetype widths with annotation lines.

    Shows:
      - Distribution of widths
      - Cutoff line at the specified width
      - Number of archetypes, TFs, families above/below cutoff
    """
    # Get width per archetype
    archetype_widths = {mid: info["width"] for mid, info in motifs.items()}
    widths = np.array(list(archetype_widths.values()))

    # Merge with metadata to get TF and family info per archetype
    # metadata maps motif_id -> cluster (=archetype), tf_name, family_name
    # We need archetype -> set of TFs, set of families
    arch_to_tfs = {}
    arch_to_families = {}
    if metadata_df is not None and len(metadata_df) > 0:
        for arch_id in archetype_widths:
            rows = metadata_df[metadata_df["cluster"] == arch_id]
            arch_to_tfs[arch_id] = set(rows["tf_name"].dropna().unique())
            arch_to_families[arch_id] = set(rows["family_name"].dropna().unique())

    # Compute stats above/below cutoff
    above_archetypes = [m for m, w in archetype_widths.items() if w >= cutoff]
    below_archetypes = [m for m, w in archetype_widths.items() if w < cutoff]

    # All unique TFs and families
    all_tfs = set()
    all_families = set()
    above_tfs = set()
    above_families = set()
    below_tfs = set()
    below_families = set()

    for arch_id in archetype_widths:
        tfs = arch_to_tfs.get(arch_id, set())
        fams = arch_to_families.get(arch_id, set())
        all_tfs |= tfs
        all_families |= fams
        if archetype_widths[arch_id] >= cutoff:
            above_tfs |= tfs
            above_families |= fams
        else:
            below_tfs |= tfs
            below_families |= fams

    n_total_arch = len(archetype_widths)
    n_above_arch = len(above_archetypes)
    n_below_arch = len(below_archetypes)
    n_total_tf = len(all_tfs)
    n_above_tf = len(above_tfs)
    n_below_tf = len(below_tfs)
    n_total_fam = len(all_families)
    n_above_fam = len(above_families)
    n_below_fam = len(below_families)

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=(10, 6))

    bins = np.arange(widths.min() - 0.5, widths.max() + 1.5, 1)
    ax.hist(widths, bins=bins, color="steelblue", edgecolor="white", alpha=0.85,
            zorder=2)

    # Cutoff line
    ax.axvline(cutoff - 0.5, color="red", linestyle="--", linewidth=2,
               zorder=3, label=f"Width = {cutoff} cutoff")

    # Annotation box — left side (below cutoff)
    below_text = (
        f"Width < {cutoff}\n"
        f"Archetypes: {n_below_arch} / {n_total_arch}\n"
        f"TFs: {n_below_tf} / {n_total_tf}\n"
        f"Families: {n_below_fam} / {n_total_fam}"
    )

    ax.text(0.02, 0.95, below_text, transform=ax.transAxes,
            fontsize=11, verticalalignment="top",  # CHANGED: was 9
            bbox=dict(boxstyle="round,pad=0.4", facecolor="lightsalmon",
                      alpha=0.8))

    # Annotation box — right side (above cutoff)
    above_text = (
        f"Width ≥ {cutoff}\n"
        f"Archetypes: {n_above_arch} / {n_total_arch}\n"
        f"TFs: {n_above_tf} / {n_total_tf}\n"
        f"Families: {n_above_fam} / {n_total_fam}"
    )
    ax.text(0.98, 0.95, above_text, transform=ax.transAxes,
            fontsize=11, verticalalignment="top", horizontalalignment="right",  # CHANGED: was 9
            bbox=dict(boxstyle="round,pad=0.4", facecolor="lightblue",
                      alpha=0.8))

    ax.tick_params(axis='both', which='major', labelsize=12)  # NEW LINE: tick labels
    ax.set_xlabel("Motif width (bp)", fontsize=14)  # CHANGED: was 12
    ax.set_ylabel("Number of motif archetypes", fontsize=14)  # CHANGED: was 12
    ax.set_title("Distribution of motif archetype widths\n"
                 "(Jeff Vierstra non-redundant archetypes, beta)", fontsize=16)  # CHANGED: was 13
    ax.legend(loc="upper center", fontsize=12)  # CHANGED: was 10

    # Add minor gridlines
    ax.grid(axis="y", alpha=0.3, zorder=0)
    ax.set_axisbelow(True)




    # ax.text(0.02, 0.95, below_text, transform=ax.transAxes,
    #         fontsize=9, verticalalignment="top",
    #         bbox=dict(boxstyle="round,pad=0.4", facecolor="lightsalmon",
    #                   alpha=0.8))

    # # Annotation box — right side (above cutoff)
    # above_text = (
    #     f"Width ≥ {cutoff}\n"
    #     f"Archetypes: {n_above_arch} / {n_total_arch}\n"
    #     f"TFs: {n_above_tf} / {n_total_tf}\n"
    #     f"Families: {n_above_fam} / {n_total_fam}"
    # )
    # ax.text(0.98, 0.95, above_text, transform=ax.transAxes,
    #         fontsize=9, verticalalignment="top", horizontalalignment="right",
    #         bbox=dict(boxstyle="round,pad=0.4", facecolor="lightblue",
    #                   alpha=0.8))

    # ax.set_xlabel("Motif width (bp)", fontsize=12)
    # ax.set_ylabel("Number of motif archetypes", fontsize=12)
    # ax.set_title("Distribution of motif archetype widths\n"
    #              "(Jeff Vierstra non-redundant archetypes, beta)", fontsize=13)
    # ax.legend(loc="upper center", fontsize=10)

    # # Add minor gridlines
    # ax.grid(axis="y", alpha=0.3, zorder=0)
    # ax.set_axisbelow(True)

    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"Saved width distribution plot: {output_path}", file=sys.stderr)

    return fig


# ============================================================================
# Orientation analysis plot
# ============================================================================

def plot_orientation_analysis(motifs, orientation_results, output_path=None):
    """
    Plot orientation similarity analysis:
      1. Scatter: cor(F, RC) vs cor(F, R) — the two key symmetry axes
      2. Histogram of cor(F, RC) — palindromicity distribution
      3. Examples of palindromic, near-palindromic, and asymmetric motifs
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    ids = list(orientation_results.keys())
    cor_f_rc = np.array([orientation_results[m]["cor_F_RC"] for m in ids])
    cor_f_r = np.array([orientation_results[m]["cor_F_R"] for m in ids])
    widths = np.array([motifs[m]["width"] for m in ids])
    classes = [orientation_results[m]["classification"] for m in ids]

    # Color by classification
    color_map = {
        "palindromic": "#e74c3c",
        "near-palindromic": "#f39c12",
        "asymmetric": "#3498db",
    }
    colors = [color_map[c] for c in classes]

    # ---- Panel A: Scatter of cor(F,RC) vs cor(F,R) ----
    ax = axes[0, 0]
    ax.scatter(cor_f_rc, cor_f_r, c=colors, s=15, alpha=0.6, edgecolors="none")
    ax.set_xlabel("cor(Forward, RevComp) — palindromicity", fontsize=10)
    ax.set_ylabel("cor(Forward, Reverse) — backward similarity", fontsize=10)
    ax.set_title("A. Orientation similarity landscape", fontsize=11)
    ax.axvline(0.95, color="red", linestyle=":", alpha=0.5)
    ax.axvline(0.80, color="orange", linestyle=":", alpha=0.5)
    ax.set_xlim(-0.3, 1.05)
    ax.set_ylim(-0.3, 1.05)

    # Legend
    for cls, clr in color_map.items():
        n = sum(1 for c in classes if c == cls)
        ax.scatter([], [], c=clr, s=30, label=f"{cls} (n={n})")
    ax.legend(fontsize=8, loc="lower right")

    # ---- Panel B: Histogram of cor(F, RC) ----
    ax = axes[0, 1]
    ax.hist(cor_f_rc, bins=50, color="steelblue", edgecolor="white", alpha=0.85)
    ax.axvline(0.95, color="red", linestyle="--", alpha=0.7, label="Palindromic (>0.95)")
    ax.axvline(0.80, color="orange", linestyle="--", alpha=0.7, label="Near-palindromic (>0.80)")
    ax.set_xlabel("cor(Forward, Reverse Complement)", fontsize=10)
    ax.set_ylabel("Count", fontsize=10)
    ax.set_title("B. Palindromicity distribution", fontsize=11)
    ax.legend(fontsize=8)

    # ---- Panel C: cor(F,RC) vs width ----
    ax = axes[1, 0]
    ax.scatter(widths, cor_f_rc, c=colors, s=15, alpha=0.6, edgecolors="none")
    ax.set_xlabel("Motif width (bp)", fontsize=10)
    ax.set_ylabel("cor(Forward, RevComp)", fontsize=10)
    ax.set_title("C. Palindromicity vs. motif width", fontsize=11)
    ax.axhline(0.95, color="red", linestyle=":", alpha=0.5)
    ax.axhline(0.80, color="orange", linestyle=":", alpha=0.5)

    # ---- Panel D: Example motif logos (as text consensus) ----
    ax = axes[1, 1]
    ax.axis("off")

    # Pick examples: most palindromic, most asymmetric, and one near-palindromic
    sorted_by_pal = sorted(ids, key=lambda m: orientation_results[m]["cor_F_RC"],
                           reverse=True)

    examples = []
    # Top palindromic (pick one with width > 6 for interest)
    for m in sorted_by_pal:
        if motifs[m]["width"] >= 8 and orientation_results[m]["classification"] == "palindromic":
            examples.append(("Palindromic", m))
            break
    if not examples:
        examples.append(("Palindromic", sorted_by_pal[0]))

    # Near-palindromic
    for m in sorted_by_pal:
        if orientation_results[m]["classification"] == "near-palindromic":
            examples.append(("Near-palindromic", m))
            break

    # Most asymmetric
    sorted_by_asym = sorted(ids, key=lambda m: orientation_results[m]["cor_F_RC"])
    for m in sorted_by_asym:
        if motifs[m]["width"] >= 8:
            examples.append(("Asymmetric", m))
            break
    if len(examples) < 3:
        examples.append(("Asymmetric", sorted_by_asym[0]))

    text_lines = ["D. Example motifs — four orientations\n"]
    text_lines.append(f"{'Type':<18} {'ID':<10} {'Width':>5}  "
                      f"{'cor(F,RC)':>9}  {'Orientations'}")
    text_lines.append("-" * 85)

    for label, mid in examples:
        pfm = motifs[mid]["pfm"]
        w = motifs[mid]["width"]
        cor_val = orientation_results[mid]["cor_F_RC"]

        cons_f = consensus_from_pfm(get_forward(pfm))
        cons_r = consensus_from_pfm(get_reverse(pfm))
        cons_c = consensus_from_pfm(get_complement(pfm))
        cons_rc = consensus_from_pfm(get_reverse_complement(pfm))

        text_lines.append(f"\n{label:<18} {mid:<10} {w:>5}  {cor_val:>9.3f}")
        text_lines.append(f"  Forward (+ strand match):    5'-{cons_f}-3'")
        text_lines.append(f"  Reverse:                     5'-{cons_r}-3'")
        text_lines.append(f"  Complement:                  5'-{cons_c}-3'")
        text_lines.append(f"  RevComp (− strand match):    5'-{cons_rc}-3'")
        text_lines.append(f"  ── MOODS searches Forward and RevComp only ──")

    ax.text(0.02, 0.98, "\n".join(text_lines), transform=ax.transAxes,
            fontsize=7.5, verticalalignment="top", fontfamily="monospace",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow",
                      alpha=0.9))

    fig.suptitle("Motif Orientation & Palindromicity Analysis", fontsize=14,
                 fontweight="bold", y=1.01)
    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"Saved orientation analysis plot: {output_path}", file=sys.stderr)

    return fig


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Sanity-check: motif width distribution and orientation analysis"
    )
    parser.add_argument("--meme", required=True,
                        help="MEME format motif file (consensus_pwms.meme)")
    parser.add_argument("--metadata", required=True,
                        help="metadata.tsv with columns: motif_id, cluster, "
                             "tf_name, family_name, ...")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for plots and tables")
    parser.add_argument("--width-cutoff", type=int, default=8,
                        help="Width cutoff for annotation (default: 8)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # ---- Load data ----
    print("Loading MEME file...", file=sys.stderr)
    motifs = parse_meme_file(args.meme)
    print(f"  Loaded {len(motifs)} motif archetypes", file=sys.stderr)

    print("Loading metadata...", file=sys.stderr)
    metadata_df = pd.read_csv(args.metadata, sep="\t")
    print(f"  Loaded {len(metadata_df)} rows, "
          f"{metadata_df['cluster'].nunique()} unique clusters, "
          f"{metadata_df['tf_name'].nunique()} unique TFs, "
          f"{metadata_df['family_name'].nunique()} unique families",
          file=sys.stderr)

    # Check which archetypes in MEME are in metadata and vice versa
    meme_ids = set(motifs.keys())
    meta_clusters = set(metadata_df["cluster"].unique())
    in_both = meme_ids & meta_clusters
    only_meme = meme_ids - meta_clusters
    only_meta = meta_clusters - meme_ids
    print(f"\n  Archetypes in MEME: {len(meme_ids)}", file=sys.stderr)
    print(f"  Clusters in metadata: {len(meta_clusters)}", file=sys.stderr)
    print(f"  In both: {len(in_both)}", file=sys.stderr)
    if only_meme:
        print(f"  In MEME only: {len(only_meme)} "
              f"(first 5: {sorted(only_meme)[:5]})", file=sys.stderr)
    if only_meta:
        print(f"  In metadata only: {len(only_meta)} "
              f"(first 5: {sorted(only_meta)[:5]})", file=sys.stderr)

    # ---- Width distribution ----
    print("\nPlotting width distribution...", file=sys.stderr)
    width_plot_path = os.path.join(args.output_dir, "motif_width_distribution.png")
    plot_width_distribution(motifs, metadata_df, cutoff=args.width_cutoff,
                            output_path=width_plot_path)

    # ---- Orientation / palindromicity analysis ----
    print("Analyzing motif orientations...", file=sys.stderr)
    orientation_results = {}
    for mid, info in motifs.items():
        orientation_results[mid] = classify_orientation_similarity(info["pfm"])

    # Summary
    n_pal = sum(1 for v in orientation_results.values()
                if v["classification"] == "palindromic")
    n_near = sum(1 for v in orientation_results.values()
                 if v["classification"] == "near-palindromic")
    n_asym = sum(1 for v in orientation_results.values()
                 if v["classification"] == "asymmetric")
    print(f"\n  Classification summary:", file=sys.stderr)
    print(f"    Palindromic (cor > 0.95):        {n_pal}", file=sys.stderr)
    print(f"    Near-palindromic (0.80-0.95):     {n_near}", file=sys.stderr)
    print(f"    Asymmetric (cor < 0.80):          {n_asym}", file=sys.stderr)

    orient_plot_path = os.path.join(args.output_dir,
                                     "motif_orientation_analysis.png")
    plot_orientation_analysis(motifs, orientation_results,
                              output_path=orient_plot_path)

    # ---- Write orientation summary table ----
    rows = []
    for mid in sorted(orientation_results.keys()):
        info = motifs[mid]
        ori = orientation_results[mid]
        pfm = info["pfm"]
        rows.append({
            "archetype_id": mid,
            "width": info["width"],
            "consensus_forward": consensus_from_pfm(get_forward(pfm)),
            "consensus_revcomp": consensus_from_pfm(get_reverse_complement(pfm)),
            "consensus_reverse": consensus_from_pfm(get_reverse(pfm)),
            "consensus_complement": consensus_from_pfm(get_complement(pfm)),
            "cor_F_RC": round(ori["cor_F_RC"], 4),
            "cor_F_R": round(ori["cor_F_R"], 4),
            "cor_F_C": round(ori["cor_F_C"], 4),
            "cor_R_RC": round(ori["cor_R_RC"], 4),
            "classification": ori["classification"],
        })

    orient_df = pd.DataFrame(rows)
    orient_tsv_path = os.path.join(args.output_dir, "orientation_analysis.tsv")
    orient_df.to_csv(orient_tsv_path, sep="\t", index=False)
    print(f"\n  Wrote orientation table: {orient_tsv_path}", file=sys.stderr)

    # ---- Print notable findings ----
    print("\n" + "=" * 70, file=sys.stderr)
    print("NOTABLE FINDINGS", file=sys.stderr)
    print("=" * 70, file=sys.stderr)

    # Palindromic motifs
    pal_df = orient_df[orient_df["classification"] == "palindromic"].sort_values(
        "cor_F_RC", ascending=False)
    print(f"\nPalindromic motifs (n={len(pal_df)}):", file=sys.stderr)
    print(f"  These produce DUPLICATE hits at the same position on + and − strands.",
          file=sys.stderr)
    print(f"  For downstream analysis, consider deduplicating by keeping only + hits.",
          file=sys.stderr)
    if len(pal_df) > 0:
        print(f"\n  {'ID':<10} {'Width':>5} {'cor(F,RC)':>10} "
              f"{'Forward':<20} {'RevComp':<20}", file=sys.stderr)
        for _, row in pal_df.head(10).iterrows():
            print(f"  {row['archetype_id']:<10} {row['width']:>5} "
                  f"{row['cor_F_RC']:>10.4f} "
                  f"{row['consensus_forward']:<20} "
                  f"{row['consensus_revcomp']:<20}", file=sys.stderr)

    # Motifs where Forward ≈ Reverse (reads same backwards)
    # This is rare and interesting — suggests internal palindromic structure
    high_f_r = orient_df[(orient_df["cor_F_R"] > 0.80) &
                          (orient_df["classification"] != "palindromic")]
    if len(high_f_r) > 0:
        print(f"\nMotifs with high Forward-Reverse similarity (reads same backwards, "
              f"n={len(high_f_r)}):", file=sys.stderr)
        for _, row in high_f_r.head(5).iterrows():
            print(f"  {row['archetype_id']}: cor(F,R)={row['cor_F_R']:.3f}, "
                  f"consensus={row['consensus_forward']}", file=sys.stderr)

    # Potentially duplicated motifs (same consensus)
    dup_consensus = orient_df.groupby("consensus_forward")["archetype_id"].apply(list)
    dups = dup_consensus[dup_consensus.apply(len) > 1]
    if len(dups) > 0:
        print(f"\nArchetypes sharing identical consensus (n={len(dups)} groups):",
              file=sys.stderr)
        for cons, ids in dups.items():
            if len(cons) <= 30:  # Don't print super long ones
                print(f"  {cons}: {', '.join(ids)}", file=sys.stderr)

    # Check if any motif's RC consensus == another motif's forward consensus
    fwd_to_id = {}
    for _, row in orient_df.iterrows():
        fwd_to_id.setdefault(row["consensus_forward"], []).append(row["archetype_id"])

    rc_matches = []
    for _, row in orient_df.iterrows():
        rc_cons = row["consensus_revcomp"]
        if rc_cons in fwd_to_id:
            matching_ids = [m for m in fwd_to_id[rc_cons]
                            if m != row["archetype_id"]]
            if matching_ids:
                rc_matches.append((row["archetype_id"], matching_ids, rc_cons))

    if rc_matches:
        # Deduplicate (A->B and B->A are the same pair)
        seen_pairs = set()
        unique_matches = []
        for mid, matches, cons in rc_matches:
            pair = tuple(sorted([mid, matches[0]]))
            if pair not in seen_pairs:
                seen_pairs.add(pair)
                unique_matches.append((mid, matches, cons))

        print(f"\nMotif pairs where one's forward = another's revcomp "
              f"(n={len(unique_matches)}):", file=sys.stderr)
        print(f"  These are effectively the same motif entered twice "
              f"(one per strand).", file=sys.stderr)
        for mid, matches, cons in unique_matches[:10]:
            print(f"  {mid} (fwd) ↔ {', '.join(matches)} (RC): {cons}",
                  file=sys.stderr)

    print("\nDone!", file=sys.stderr)


if __name__ == "__main__":
    main()
    