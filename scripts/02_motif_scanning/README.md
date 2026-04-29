# `02_motif_scanning/` — Genome-wide MOODS motif scan

Runs MOODS over `mm39.fa` for every PWM in `consensus_pwms.meme`, optionally
de-duplicates palindromic motifs to a single strand, and produces a QC
report (genomic context, TSS distances, score distributions).

Inputs come from `01_resources/` (genome) and the project resources dir
(MEME file, metadata, GTF).

---

## Pipeline at a glance

```
0x_run_sanity_check.sh           (PRE-FLIGHT — produces palindrome table)
            │
            ▼
00_prepare.py                    (background freqs + chromosome list)
            │
            ▼
01_scan_moods.py                 (per-chrom array; via 01_scan_moods.sbatch)
            │
            ▼
02_postprocess.py                (merge, dedup palindromes, summary)
            │
            ▼
03_qc_report.py                  (genomic context, TSS, plots)
```

`run_pipeline.sh` is the orchestrator that runs steps 00 → 03 with
SLURM dependencies.

---

## Quickstart

```bash
# 1. Pre-flight sanity check (login node, ~1 min). Produces palindrome table.
bash 0x_run_sanity_check.sh

# 2. Dry-run the pipeline (prints what it would do)
bash run_pipeline.sh

# 3. Submit
bash run_pipeline.sh --submit

# 4. Or, to skip palindrome dedup
bash run_pipeline.sh --submit --no-palindrome-dedup
```

The sanity check is optional in the strict sense — `02_postprocess.py`
falls through gracefully without `--palindrome-table` — but you almost
always want it. Without it, palindromic motifs (like the ZF dimer that
reads `ACGTACGT`) get counted twice (once on each strand).

---

## Outputs

```
${BASE_OUTDIR}/
├── prep/
│   ├── background.txt
│   └── chromosomes.txt
├── sanity_check/
│   ├── motif_width_distribution.png
│   ├── motif_orientation_analysis.png
│   └── orientation_analysis.tsv     # palindrome table (consumed by step 2)
├── per_chrom/                       # raw MOODS output, one BED per chrom
│   ├── chr1.bed
│   ├── ...
│   └── chrY.bed
├── merged/
│   ├── per_motif/                   # one BED per motif (after dedup)
│   ├── all_motifs_merged.bed.gz     # sorted, gzipped
│   ├── motif_hit_summary.tsv        # per-motif n_hits, scores, dedup count
│   └── qc_report.txt                # text QC summary
├── qc_report/                       # comprehensive HTML/PDF QC
└── logs/
```

The per-motif BEDs in `merged/per_motif/` are what
`03_motif_annot/02_annotate_motifs.R` consumes.

---

## Configuration

All paths and parameters live in `config.sh`. Every script in this
directory sources it:

```bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
```

Key variables:

| Variable | Default | Purpose |
|---|---|---|
| `GENOME` | `…/mm39.fa` | Path to indexed FASTA |
| `CHROM_SIZES` | `…/mm39.chrom.sizes` | UCSC chrom.sizes (drives effective genome size) |
| `MEME_FILE` | `…/consensus_pwms.meme` | Vierstra archetypes, MEME format |
| `METADATA_FILE` | `…/metadata.tsv` | Vierstra archetype metadata |
| `GTF_FILE` | `…/gencode.vM38…gtf.gz` | For the QC report's TSS distances and promoter overlap |
| `BASE_OUTDIR` | `…/moods_scan` | Pipeline output root |
| `PVALUE` | `1e-5` | MOODS scanning threshold |
| `PROMOTER_UPSTREAM` | `1000` | bp upstream of TSS for promoter window |
| `PROMOTER_DOWNSTREAM` | `500` | bp downstream of TSS for promoter window |
| `CONDA_ENV` | `motif_scanning` | python+MOODS+bedtools |

**Promoter window**: keep `PROMOTER_UPSTREAM` and `PROMOTER_DOWNSTREAM`
identical to the values in `03_motif_annot/config.sh`. The annotation
pipeline assumes [-1000, +500] by default.

---

## Scripts

| Script | Type | Purpose |
|---|---|---|
| `00_prepare.py` | python | Compute genome ACGT freqs and write `chromosomes.txt`. |
| `01_scan_moods.py` | python | Run MOODS for one chromosome against all motifs. |
| `01_scan_moods.sbatch` | SLURM array | Wrapper for `01_scan_moods.py` (one task per chromosome). |
| `02_postprocess.py` | python | Merge per-chrom BEDs, palindrome dedup, summary. |
| `02_postprocess.sbatch` | SLURM | Wrapper for `02_postprocess.py`. |
| `03_qc_report.py` | python | Comprehensive QC: genomic context, TSS, plots. |
| `03_qc_report.sbatch` | SLURM | Wrapper for `03_qc_report.py`. |
| `0x_sanity_check_motifs.py` | python | Pre-flight: width hist, orientation analysis, palindrome table. |
| `0x_run_sanity_check.sh` | bash | Wrapper for the sanity check (login-node OK). |
| `0x_apply_promoter_patch.py` | python | One-shot in-place patcher for `03_qc_report.py` (parameterizes promoter window). |
| `run_pipeline.sh` | bash | Master orchestrator (steps 0 → 3, with SLURM dependencies). |

The `0x_*` prefix marks utilities/pre-flight scripts that aren't part of
the main numerical sequence.

---

## Palindrome dedup — what it does

Some PWMs are palindromic: forward log-odds matrix ≈ reverse-complement
matrix (e.g. ZF dimers reading `ACGTACGT`). MOODS scans both strands and
reports both hits at every position, even though only one binding event
is happening. Without dedup, palindromic motifs show up with 2× their
true hit count, with a 50/50 strand split that's an artifact, not biology.

`0x_sanity_check_motifs.py` classifies each PWM by `cor(F, RC)`:

| `cor(F, RC)` | Class | What `02_postprocess.py` does |
|---|---|---|
| > 0.95 | `palindromic` | Drops `−` strand hits (keeps `+` only). |
| 0.80–0.95 | `near-palindromic` | Keeps both by default. Drop with `--dedup-near-palindromic`. |
| ≤ 0.80 | `asymmetric` | Keeps both — strand IS the orientation signal. |

The dedup is non-destructive: only `merged/per_motif/`, `merged/*.bed.gz`,
and the summary are affected. Raw `per_chrom/*.bed` are untouched, so
you can re-run with different settings without re-scanning.

---

## Resource estimates

| Step | Time | Memory |
|---|---|---|
| `0x_run_sanity_check.sh` | ~1 min | 4 GB (login-node OK) |
| `00_prepare.py` | ~1 min | 4 GB |
| `01_scan_moods.py` (per chrom) | 2–8 hr | 32–64 GB |
| `02_postprocess.py` | 30–60 min | 32–64 GB |
| `03_qc_report.py` | 1–4 hr | 32–64 GB |

The array job for step 1 runs 21 chromosomes (1–19, X, Y). Skip chrM.

---

## Known issues / things to watch

1. **`03_qc_report.py` was originally hardcoded to [-1000, +200]** and
   has been parameterized via `--promoter-upstream`/`--promoter-downstream`.
   `0x_apply_promoter_patch.py` does this in place; run it once after
   copying the new files in:
   ```bash
   python 0x_apply_promoter_patch.py 03_qc_report.py
   ```
   The patch is idempotent.

2. **`03_validate_promoter_enrichment.py` is removed** — `03_qc_report.py`
   does the same checks better with consistent parameters.

3. **MEME parser strips `:` in motif IDs**: motif IDs like
   `AC0001:Tbx5_etc` become `AC0001`. This is intentional (matches the
   palindrome-table key format) but worth knowing if you see surprising
   key collisions.

4. **`set -euo pipefail`** is set in every wrapper. If you see weird
   failures, check the wrapper isn't being invoked through a chain that
   strips it.

5. **No version pinning** for the `motif_scanning` conda env. Add an
   `environment.yml`.

6. **chrM is excluded** from effective-genome-size calculation (set in
   `02_postprocess.py::compute_effective_genome_size`). If you ever
   want to include it, edit that function.
