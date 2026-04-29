# `07_diagnostics/` — Overview & sanity-check plots

This pipeline produces the descriptive overview figures we need before any
hypothesis testing. It does **not** modify any upstream pipeline outputs —
it only reads the per-motif annotated tables from `05_motif_annot/`, the
per-pair tables from `06_allele_pairs/`, and the allele-specific BAMs from
the upstream assay pipelines.

The full statistical / inferential analyses (allelic regressions, ARIMA peak
detection, ΔPWM × Δpause × ΔPSI tests) live in `08_metaprofiles/`,
`09_arima/`, and `10_main_analyses/`. This pipeline answers four questions
without testing any hypotheses:

1. **Which TFBSs survive our standard filtering cascade, and how many?**
2. **Which motif archetypes have enough SNP-overlapping TFBSs to power
   downstream analyses?**
3. **Over what genomic distances can a SNP plausibly affect chromatin
   accessibility?** (Empirical answer to the SNP-coupling-distance question.)
4. **Is WASP working correctly, and how many archetypes are jointly powered
   for the full ATAC × PRO-seq × RNA-seq analysis?**

## Filter cascade (the central design choice)

In transcription-direction order:

```
0. All scanned TFBSs (raw 02_motif_scanning output)
1. Standard chromosomes (chr1–19, X; exclude M, Y)
2. Intragenic on the motif's sense strand
3. Host gene type ∈ {protein_coding, lncRNA}
4. Host gene expressed (TPM ≥ 1, mean across reps)        [toggleable]
5. Outside promoter-proximal region (-1000 to +500 from any TSS)
6. ATAC-accessible (overlaps any ATAC consensus peak)
```

All steps use `_sense` columns from `05_motif_annot/`. Step 4 is toggleable
via `EXPRESSION_TPM_MIN` in `config.sh`.

A "**qualifying TFBS**" passes the full cascade AND has at least one het-SNP
within `QUALIFYING_DISTANCE` bp (default 0 = body overlap; can include
flank bins up to 800 bp matching `05_motif_annot`'s flank schema).

A "**powered TFBS**" is a qualifying TFBS where at least one of its
paired SNPs has sufficient allele-specific reads (toggleable per assay).

## File listing

```
07_diagnostics/
├── README.md                       (this file)
├── config.sh                       paths, parameters, defaults
├── lib_diagnostics.R               shared R utilities
├── 00_setup.sh                     pre-flight input validation
├── 01_build_summary.R              per-motif: cascade + qualifying counts
├── 02_aggregate_summary.R          merge per-motif rows into master tables
├── 03_plot_filter_cascade.R        plot (i): cascade bars + heatmap
├── 04_plot_qualifying_histogram.R  plot (ii): histogram + ECDF
├── 05_plot_top_archetypes.R        plot (iii): bars + logos + max-TPM color
├── 06_plot_pwm_score_dists.R       plot (iv): PWM score dists (per-arch + global)
├── 07_plot_snp_distance_diag.R     ATAC-coupling distance diagnostic
├── 08_compute_wasp_balance.sh      WASP cascade + per-SNP allelic depth
├── 09_plot_wasp_balance.R          WASP cascade + balance + coverage plots
├── 10_plot_powered_archetypes.R    joint power filter histogram
├── run_pipeline.sh                 SLURM orchestrator
└── 06_allele_pairs_max_flank_extension.patch   one-time patch for upstream
```

## Usage

```bash
# 1. (One time only) extend 06_allele_pairs/ to max flank 800 bp.
#    This is needed for the SNP-distance ATAC diagnostic (plot 07).
cd ../06_allele_pairs
patch -p1 < ../07_diagnostics/06_allele_pairs_max_flank_extension.patch
bash run_pipeline.sh --apply        # re-run 06 with the new flank
cd ../07_diagnostics

# 2. Validate inputs (lenient by default).
bash 00_setup.sh
# Strict mode:
STRICT=1 bash 00_setup.sh

# 3. Dry-run / submit pipeline.
bash run_pipeline.sh                # dry run
bash run_pipeline.sh --apply        # actually submit

# Selective re-runs:
bash run_pipeline.sh --apply --start-from 03   # skip 01/02 if already done
bash run_pipeline.sh --apply --skip 08,09       # don't recompute WASP
```

## Plot guide

### (i) Filter cascade — `plots/cascade_*.pdf`
Bar plot with one bar per cascade stage, showing the number of TFBSs (or
unique bp) passing through each filter. Two archetype subsets:
- `_all`: all archetypes
- `_widthgt8`: only archetypes with motif width > 8 bp (default for
  downstream analyses; toggleable via `MIN_MOTIF_WIDTH` in `config.sh`)

Plus a per-archetype heatmap showing survival fraction at each stage,
sorted by final-stage retention. Useful for spotting archetypes that drop
out unexpectedly.

### (ii) Qualifying TFBSs per archetype — `plots/qualifying_hist_*.pdf`
Histogram + ECDF of qualifying-TFBS counts per archetype. Vertical lines:
- Red dashed: median
- Green-gradient dashed: cumulative thresholds (≥250, ≥500, ≥1000, ≥2000),
  with archetype counts annotated

One PDF per qualifying-distance threshold (overlap, flank10, flank25, ...,
flank800), each with `_all` and `_widthgt8` variants.

### (iii) Top archetypes panel — `plots/top_archetypes_n*_overlap.pdf`
For top-N archetypes by qualifying-TFBS count: horizontal bar (count) +
fill color = `log2(max TPM + 1)` across mapped TFs + ggseqlogo of the PWM
on the right. Generated for N ∈ `TOP_N_LIST` (default `10,25,50`).

"Max TPM" = mean of replicate TPMs from `gene_expression_summary.tsv`,
maximum across TFs mapped to the archetype (via `metadata.tsv`'s `cluster`
column joined to `tf_name` → `gene_name`).

### (iv) PWM score distributions — `plots/pwm_score_dists_*.pdf`
Two views:
- **Per-archetype faceted** (`pwm_score_dists_per_archetype_top25.pdf`):
  violins for each top archetype, comparing all TFBSs vs qualifying TFBSs.
  Tests whether SNPs are enriched/depleted at high-scoring TFBS positions.
- **Global pooled** (`pwm_score_dists_global_zscore.pdf`):
  within-archetype z-normalized score density, pooled across all
  width-passed archetypes. Cross-archetype-comparable.

### SNP-distance ATAC diagnostic — `plots/snp_distance_atac_*.pdf`
Empirical answer to "over what distance can a SNP affect chromatin
accessibility?". Per (motif, SNP) pair, computes:
- `|log2((atac_ref+1)/(atac_alt+1))|` at the motif body
- Bin by `|SNP-to-motif-edge distance|` (default bins 0, 1-25, 25-50,
  50-100, 100-200, 200-400, 400-800)

Compared against a permutation background (random SNP–motif assignments).
Where the observed bin median merges with background = the practical
effective coupling distance for ATAC asymmetry.

**Hard-coded distances are now empirical.** The previous suggestion of
"~250 bp for ATAC, ~50 bp for cooperative effects" gets replaced by what
this plot says. Adjust `ATAC_DIST_BIN_EDGES` in `config.sh` to refine
binning.

### WASP cascade — `plots/wasp_cascade.pdf`
Per assay × replicate, dodged bars showing reads at each WASP stage:
Total input → WASP-passed (vW=1) → Allele-split total. Healthy WASP
output: minimal drop from input to WASP-passed; allele-split ≈ WASP-passed.

### Allele balance — `plots/allele_balance.pdf`
Two-panel plot. Top: per-replicate ref/alt bars (post-WASP). Bottom:
reference fraction with red dashed line at 0.5. Strong departures from 0.5
indicate residual reference bias.

### Per-SNP coverage — `plots/per_snp_coverage_hist.pdf`
Histogram of total reads at each unique SNP position (summed across
replicates and alleles within each assay). Vertical thresholds at
10/25/50 reads, with annotated counts. Tells you what fraction of het-SNPs
have enough depth to call alleles confidently.

### Powered archetypes — `plots/powered_archetypes_*.pdf`
**The headline diagnostic.** For each power definition (ATAC-only,
PRO-seq-only, RNA-seq-only, ATAC + PRO-seq, all three), histogram of
"powered motif hits per archetype" in the same aesthetic as plot (ii).
Tells you exactly how many archetypes you can power for the full
allelic-pair analysis.

A pair is "powered" under each definition if:
- ATAC: `min(atac_ref_total, atac_alt_total) ≥ POWER_ATAC_MIN_PER_ALLELE` (default 10)
- PRO-seq: `pause_min_reads_per_allele ≥ POWER_PROSEQ_PAUSE_MIN_PER_ALLELE` (default 5)
- RNA-seq: at least one junction cell has `min(ref_total, alt_total) ≥
  POWER_RNASEQ_AT_JUNCTION_MIN_PER_ALLELE` (default 4)

A motif hit is "powered" if at least one of its paired SNPs is powered.

## Output directory layout

```
${SCRATCH_ROOT}/diagnostics_v1/
├── per_archetype_summaries/        intermediate (one TSV per motif)
│   ├── AC0001_summary.tsv          per-motif cascade row
│   ├── AC0001_scores.tsv.gz        per-TFBS scores at post-cascade stage
│   └── ...
├── tables/
│   ├── archetype_summary.tsv          master table (~50 cols × ~700 rows)
│   ├── qualifying_per_archetype.tsv    compact (motif × qualifying counts)
│   ├── snp_distance_atac_diagnostic.tsv  bin summaries from plot 07
│   └── powered_per_archetype.tsv         from plot 10
├── plots/
│   ├── cascade_*.pdf / .png
│   ├── qualifying_hist_*.pdf / .png
│   ├── top_archetypes_*.pdf / .png
│   ├── pwm_score_dists_*.pdf / .png
│   ├── snp_distance_atac_*.pdf / .png
│   ├── wasp_cascade.pdf / .png
│   ├── allele_balance.pdf / .png
│   ├── per_snp_coverage_hist.pdf / .png
│   └── powered_archetypes_*.pdf / .png
├── wasp/
│   ├── wasp_cascade.tsv
│   ├── allele_balance.tsv
│   ├── unique_snps.bed
│   └── per_snp_coverage.tsv.gz
└── logs/
```

All sized to fit in scratch; intermediate per-motif summaries are tiny
(~1 KB each) and per-TFBS score files are typically <1 MB each.

## Lenient vs strict mode

By default (lenient), the pipeline runs on whatever subset of motifs is
currently present in `05_motif_annot/` and `06_allele_pairs/`. Missing
motifs simply produce no rows in the summary table. Set `STRICT=1` to
require all motifs from `02_motif_scanning/` to have annotations and pair
tables.

This is critical because `06_allele_pairs/` is the slowest upstream
pipeline; we typically want to start producing 07 plots from the partial
05+06 results before everything finishes.

## Cross-pipeline dependencies

Reads:
- `02_motif_scanning/per_motif/*.bed` (counted in setup, not actually read)
- `05_motif_annot/{motif_id}_annotated.tsv.gz` (per-motif annotated tables)
- `06_allele_pairs/pair_tables/{motif_id}_pair_table.tsv.gz` (per-pair tables)
- `06_allele_pairs/pair_index/{motif_id}_pair_index.tsv.gz` (for unique SNPs)
- Allele-specific BAMs from `03_atac/`, `04_proseq/`, `04c_rnaseq/`
- `consensus_pwms.meme`, `metadata.tsv`, `gene_expression_summary.tsv`

Writes:
- `${SCRATCH_ROOT}/diagnostics_v1/` (all outputs)

## Required environment

R packages (in `test_m`): `optparse`, `data.table`, `ggplot2`, `scales`,
`patchwork`, `ggseqlogo`, `GenomicRanges`. All standard.

Tool packages (in `proseq2_env`): `samtools`, `bcftools`. Used by stage 08
only.

If `ggseqlogo` is not installed:
```r
install.packages("ggseqlogo")  # CRAN
```

## What this pipeline deliberately does NOT do

- It does not run any allele-specific signal extraction at single-bp
  resolution (that's `08_metaprofiles/`).
- It does not fit any models (that's `09_arima/`).
- It does not perform any hypothesis tests (that's `10_main_analyses/`).
- It does not modify upstream outputs (the one exception is the
  `06_allele_pairs/` patch, applied separately).

## Troubleshooting

**`[SKIP] No pair tables found`** — `06_allele_pairs/` hasn't run, or
ran to a different output directory. Check `PAIR_TABLE_DIR` in `config.sh`.

**Empty cascade plot** — all motifs are missing `_sense` columns. Ensure
`05_motif_annot/` was patched and re-run with the strand-aware refactor
(see 05's README §"Compatibility break").

**`expression_tpm_max_sense` all NA** — the gene expression summary file
is missing or has different column names. Check `GENE_EXPRESSION_FILE` in
`config.sh`.

**Per-SNP coverage step takes hours** — bcftools mpileup at ~1 M unique
SNPs across 20 BAMs can take time. Run via SLURM (the orchestrator does
this); single-node bash runs may time out.

**`ggseqlogo` errors in plot (iii)** — the package is missing or the PWM
matrix has wrong shape. Falls back to bars-only if logo construction
fails.
