#!/usr/bin/env python3
"""
0x_apply_promoter_patch.py — one-shot in-place patch for 03_qc_report.py

Replaces the hardcoded promoter window (TSS-1000 / TSS+200) with two new
CLI arguments (--promoter-upstream, --promoter-downstream), defaulting to
1000 and 500 respectively. Idempotent: re-running on an already-patched
file is a no-op.

Run from inside 02_motif_scanning/ after copying in the new files:

    python 0x_apply_promoter_patch.py 03_qc_report.py

Then verify with:

    grep -n 'promoter-upstream\\|promoter_upstream' 03_qc_report.py

You should see the new argparse entries and the function-signature update.
"""

import sys
import os
import re


PATCHES = [
    # --- (1) parse_gtf_for_annotations signature ---
    {
        "name": "Function signature: parse_gtf_for_annotations",
        "old": "def parse_gtf_for_annotations(gtf_path, tmpdir):",
        "new": ("def parse_gtf_for_annotations("
                "gtf_path, tmpdir, "
                "promoter_upstream=1000, promoter_downstream=500):"),
    },
    # --- (2) Replacement of hardcoded 1000/200 inside that function ---
    {
        "name": "Promoter BED window (1000/200 → parameter)",
        "old": (
            '        if strand == "+":\n'
            '                prom_start = max(0, pos - 1000)\n'
            '                prom_end = pos + 200\n'
            '            else:\n'
            '                prom_start = max(0, pos - 200)\n'
            '                prom_end = pos + 1000\n'
        ),
        "new": (
            '        if strand == "+":\n'
            '                prom_start = max(0, pos - promoter_upstream)\n'
            '                prom_end = pos + promoter_downstream\n'
            '            else:\n'
            '                prom_start = max(0, pos - promoter_downstream)\n'
            '                prom_end = pos + promoter_upstream\n'
        ),
    },
    # --- (3) Call site in main(): pass the args through ---
    {
        "name": "main(): pass promoter args into parse_gtf_for_annotations",
        "old": "annotation_beds = parse_gtf_for_annotations(args.gtf, tmpdir)",
        "new": ("annotation_beds = parse_gtf_for_annotations("
                "args.gtf, tmpdir, "
                "args.promoter_upstream, args.promoter_downstream)"),
    },
]


# argparse insertion: add two new arguments. We insert them right after the
# --pvalue argument, which we locate by string match.
ARGPARSE_INSERT_ANCHOR = (
    'parser.add_argument("--pvalue", type=float, default=1e-5,\n'
    '                        help="P-value used for scanning (default: 1e-5)")'
)
ARGPARSE_NEW_ARGS = (
    '\n'
    '    parser.add_argument("--promoter-upstream", type=int, default=1000,\n'
    '                        help="Promoter bp upstream of TSS (default: 1000)")\n'
    '    parser.add_argument("--promoter-downstream", type=int, default=500,\n'
    '                        help="Promoter bp downstream of TSS (default: 500)")'
)


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    path = sys.argv[1]
    if not os.path.isfile(path):
        print(f"ERROR: file not found: {path}", file=sys.stderr)
        sys.exit(1)

    with open(path) as f:
        src = f.read()

    # --- Idempotency guard: detect already-patched files ---
    already_patched = ("--promoter-upstream" in src and
                       "promoter_upstream=1000" in src)
    if already_patched:
        print(f"{path} already appears patched; nothing to do.")
        sys.exit(0)

    n_applied = 0

    for p in PATCHES:
        if p["new"] in src:
            print(f"  Skipped (already applied): {p['name']}")
            continue
        if p["old"] not in src:
            print(f"  WARNING: anchor not found for: {p['name']}",
                  file=sys.stderr)
            print(f"           expected:\n{p['old']}", file=sys.stderr)
            continue
        src = src.replace(p["old"], p["new"], 1)
        print(f"  Applied: {p['name']}")
        n_applied += 1

    # --- Argparse insertion ---
    if ARGPARSE_INSERT_ANCHOR in src and "--promoter-upstream" not in src:
        src = src.replace(
            ARGPARSE_INSERT_ANCHOR,
            ARGPARSE_INSERT_ANCHOR + ARGPARSE_NEW_ARGS,
            1,
        )
        print("  Applied: argparse: --promoter-upstream / --promoter-downstream")
        n_applied += 1
    elif "--promoter-upstream" in src:
        print("  Skipped (already applied): argparse insertion")
    else:
        print("  WARNING: argparse anchor not found — please add the two "
              "args manually:", file=sys.stderr)
        print(ARGPARSE_NEW_ARGS, file=sys.stderr)

    if n_applied == 0:
        print("\nNo changes applied (file already up to date or anchors moved).")
        sys.exit(0)

    backup = path + ".bak"
    print(f"\nBacking up original to {backup}")
    os.rename(path, backup)
    with open(path, "w") as f:
        f.write(src)
    print(f"Wrote patched: {path}")
    print(f"\nVerify with: grep -n 'promoter-upstream\\|promoter_upstream' {path}")


if __name__ == "__main__":
    main()
