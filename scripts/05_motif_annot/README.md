# `05_motif_annot/` — Per-motif genomic annotation

Annotates each motif hit (one row per (motif_archetype, genomic_position) pair) with strand-aware genomic context, ATAC accessibility, F1 het SNP overlap, RNA-seq gene expression, and splice junction proximity. Produces one `${motif_id}_annotated.tsv.gz` file per motif archetype.

This is the refactor of the old `03_motif_annot/`. The major changes are:

1. **Strand-aware everything** — every motif row gets two values for each genic feature (`_sense` for motif's strand, `_antisense` for the opposite). Orientation flipping = column relabeling, not recomputation.
2. **ATAC peaks from new `03_atac/`** — F121-9-specific consensus peaks in mm39 directly. No more 4DN bigBed import + mm10→mm39 liftOver.
3. **Gene expression** — joins against `04c_rnaseq/`'s expression table to produce 10 `expression_*` columns.
4. **Splice junction proximity (NEW)** — 40 columns covering nearest splice donor and acceptor in each (sense/antisense × upstream/downstream) cell, sourced from **GTF + STAR's `SJ.out.tab` discovered junctions** (filtered by ≥10 unique reads in any rep), tracking source/support/motif type.

## Pipeline layout

```
05_motif_annot/
├── README.md                           ← you are here
├── config.sh                           paths, params, SLURM defaults
├── 01_preprocess_resources.R           builds shared GRanges and tables
├── 01_run_preprocess_resources.sh      SLURM wrapper for the above
├── 02_annotate_motifs.R                orchestrator, sources 02a-02g
├── 02a_motif_basic.R                   sequence + palindromicity
├── 02b_strand_aware_genic.R            intragenic, promoter, genic regions
├── 02c_strand_aware_tss.R              nearest TSS distances
├── 02d_atac_overlap.R                  ATAC consensus peak overlap
├── 02e_snp_overlap.R                   F1 het SNP overlap (direct + flanking)
├── 02f_gene_expression.R               RNA-seq expression annotation
├── 02g_splice_proximity.R              splice donor + acceptor proximity
├── 02_annotate_motifs.sh               SLURM array submission
├── 0x_find_missing_motifs.sh           find motifs without output
├── 0x_resubmit_missing_motifs.sh       resubmit on a higher-resource queue
└── run_pipeline.sh                     end-to-end orchestrator
```

## Usage

```bash
# Stage 1: preprocess (single SLURM job, ~1 hour)
sbatch 01_run_preprocess_resources.sh

# Stage 2: per-motif array (after Stage 1 finishes)
bash 02_annotate_motifs.sh             # dry-run
bash 02_annotate_motifs.sh --apply     # actually submit

# Optional Stage 3: retry any failed motifs on a bigger queue
bash 0x_find_missing_motifs.sh
bash 0x_resubmit_missing_motifs.sh --apply

# OR: top-level orchestrator
bash run_pipeline.sh --apply           # submits Stage 1; then run Stage 2 manually
bash run_pipeline.sh --apply --start-from 02  # skip Stage 1
```

## Output schema (~136 columns per motif row)

The output TSV has this column structure. All `*_sense` columns describe features on the same strand as the motif row; all `*_antisense` columns describe features on the opposite strand. Unstranded annotations (ATAC, SNPs, sequence) have no suffix.

### Original BED columns (6)
| Column | Description |
|--------|-------------|
| `chrom`, `start`, `end` | BED 0-based half-open coordinates |
| `motif_id` | Motif archetype identifier (e.g. `AC0002`) |
| `score` | MOODS score from motif scanning |
| `strand` | `+` or `-` (motif orientation) |

### `02a` — Motif basic (2 columns)
| Column | Description |
|--------|-------------|
| `motif_sequence` | DNA sequence in motif strand orientation (already reverse-complemented for `-` strand) |
| `palindromicity_score` | Fraction of positions matching reverse complement (0=none, 1=fully palindromic) |

### `02b` — Strand-aware genic (12 columns)
For each `_orient` ∈ {`_sense`, `_antisense`}:

| Column | Description |
|--------|-------------|
| `intragenic{orient}` | Boolean — motif overlaps a transcript on this strand |
| `intragenic_gene_ids{orient}` | Pipe-separated GENCODE gene IDs |
| `intragenic_gene_names{orient}` | Pipe-separated gene symbols |
| `intragenic_gene_types{orient}` | Pipe-separated biotypes (e.g. `protein_coding\|lncRNA`) |
| `is_promoter{orient}` | Boolean — motif falls in [TSS-1000, TSS+500] of any transcript on this strand |
| `promoter_transcript_ids{orient}` | Pipe-separated transcript IDs |
| `promoter_gene_names{orient}` | Pipe-separated gene names |
| `promoter_tsls{orient}` | Pipe-separated transcript support levels |
| `genic_regions{orient}` | Pipe-separated regions overlapped (e.g. `5UTR\|exon` or `intron\|exon`) |

Promoter window matches `02_motif_scanning/` (currently `-1000 / +500`; configured in `config.sh`).

### `02c` — Strand-aware nearest TSS (16 columns)
For each `_direction` ∈ {`upstream_tss_`, `downstream_tss_`} × `_orient` ∈ {`_sense`, `_antisense`}:

| Column | Description |
|--------|-------------|
| `{direction}dist{orient}` | Distance from motif midpoint to nearest TSS in this direction (always positive bp) |
| `{direction}gene_names{orient}` | Pipe-separated gene names of TSS-bearing transcripts |
| `{direction}tx_ids{orient}` | Pipe-separated transcript IDs |
| `{direction}tsls{orient}` | Pipe-separated transcript support levels |

`upstream` = 5' direction in transcription sense; `downstream` = 3' direction.

### `02d` — ATAC consensus peak overlap (2 columns)
| Column | Description |
|--------|-------------|
| `atac_consensus_overlap` | Boolean — motif overlaps any ATAC consensus peak |
| `atac_consensus_n_replicates` | Max replicate-support count (1-4) among overlapping peaks; NA if no overlap |

Source: `03_atac/peaks/consensus/consensus_peaks.bed`, F121-9-specific, mm39, no liftOver. Unstranded.

### `02e` — F1 het SNP overlap (48 columns: 2 hybrids × {direct + 7 flanks} × 3)
For each hybrid in {`F121-9`, `BL6xCAST`}:

| Column pattern | Description |
|----------------|-------------|
| `snp_<hybrid>_overlap` | Direct overlap with motif body |
| `snp_<hybrid>_details` | Pipe-separated `chr:pos:REF>ALT` |
| `snp_<hybrid>_count` | Number of unique SNPs |
| `snp_<hybrid>_flank<D>bp_overlap` | Same, for annular bin D bp from motif edge |
| `snp_<hybrid>_flank<D>bp_details` | (same) |
| `snp_<hybrid>_flank<D>bp_count` | (same) |

Flank bins are concentric and **non-overlapping**:

| Bin | Distance from motif edge |
|-----|--------------------------|
| `flank10bp`  | positions 1–10 |
| `flank25bp`  | positions 11–25 |
| `flank50bp`  | positions 26–50 |
| `flank100bp` | positions 51–100 |
| `flank200bp` | positions 101–200 |
| `flank400bp` | positions 201–400 |
| `flank800bp` | positions 401–800 |

(Configurable via `FLANK_DISTANCES` in `config.sh`.) Unstranded.

### `02f` — Gene expression (10 columns)
For each `_orient` ∈ {`_sense`, `_antisense`}:

| Column | Description |
|--------|-------------|
| `expression_n_genes{orient}` | Number of genes overlapping the motif on this strand |
| `expression_max_gene_id{orient}` | Gene ID of the most-expressed overlapping gene |
| `expression_max_gene_name{orient}` | Its gene symbol |
| `expression_tpm_max{orient}` | Its `tpm_mean` (across RNA-seq replicates) |
| `expression_tpm_sd_at_max{orient}` | Its `tpm_sd` |
| `expression_counts_sum{orient}` | Sum of `counts_sum` across all overlapping genes |

Source: `04c_rnaseq/expression/gene_expression_summary.tsv`. The "max" choice for the summary statistic reflects the goal of filtering out motifs that aren't on transcribed genes — if any overlapping gene is highly expressed, the motif is in transcribed sequence.

### `02g` — Splice junction proximity (40 columns)
**This is the largest annotation block by column count.** For each `_site` ∈ {`donor_`, `acceptor_`} × `_orient` ∈ {`sense`, `antisense`} × `_direction` ∈ {`upstream`, `downstream`}, we report 5 sub-columns:

| Column pattern | Description |
|----------------|-------------|
| `{site}{orient}_{direction}_dist` | Distance from motif midpoint to nearest splice site of this type/strand in this direction (positive bp) |
| `{site}{orient}_{direction}_source` | `GTF` (annotated only), `SJ_only` (STAR-discovered, not in GTF), or `GTF_and_SJ` (both — strongest evidence) |
| `{site}{orient}_{direction}_support` | Max `n_unique_reads` across reps for SJ-supported sites; NA for `GTF`-only |
| `{site}{orient}_{direction}_motif_type` | `GT/AG`, `GC/AG`, `AT/AC`, or `non_canonical` (pipe-separated if multiple); empty for GTF-only |
| `{site}{orient}_{direction}_gene_names` | Pipe-separated gene names of GTF transcripts using this site; empty for SJ-only sites without GTF coverage |

#### Splice-site definitions
- **Donor** = first base of intron in transcription direction (the GT of GT/AG)
- **Acceptor** = last base of intron in transcription direction (the AG of GT/AG)

In genomic coordinates:
| Strand | Donor position | Acceptor position |
|--------|----------------|-------------------|
| `+` | `intron_start` (lower coord) | `intron_end` (higher coord) |
| `-` | `intron_end` (higher coord) | `intron_start` (lower coord) |

#### Source tracking — three states

| `source` value | Meaning |
|----------------|---------|
| `GTF` | This donor/acceptor position appears in at least one GTF-annotated intron, but no STAR-discovered junction with ≥10 unique reads supports it. Annotated but not detectable in our RNA-seq. |
| `SJ_only` | This position appears in STAR's `SJ.out.tab` with sufficient support, but is not in any GTF intron. **Novel/unannotated junction.** |
| `GTF_and_SJ` | Both — annotated and supported by RNA-seq reads. Highest confidence. |

#### STAR `SJ.out.tab` filter
A junction is kept if:
- `n_unique` (column 7) ≥ `SJ_MIN_UNIQUE_READS` (default 10) in any replicate
- Strand code (column 4) ∈ {1, 2} (i.e., not `0`/undefined)

Multi-mapping reads (column 8) are **not** counted toward support.

#### Example column names (8 cells × 5 sub-cols = 40 columns)
```
donor_sense_upstream_dist        donor_sense_upstream_source
donor_sense_upstream_support     donor_sense_upstream_motif_type
donor_sense_upstream_gene_names

donor_sense_downstream_dist       (and 4 more sub-cols)
donor_antisense_upstream_dist     (and 4 more sub-cols)
donor_antisense_downstream_dist   (and 4 more sub-cols)

acceptor_sense_upstream_dist      (and 4 more sub-cols)
acceptor_sense_downstream_dist    (and 4 more sub-cols)
acceptor_antisense_upstream_dist  (and 4 more sub-cols)
acceptor_antisense_downstream_dist (and 4 more sub-cols)
```

## Strand-aware overlap mechanics (the `_sense` / `_antisense` convention)

For every genic annotation, we compute two values:

- **Sense:** motif strand = feature strand
- **Antisense:** motif strand opposite to feature strand

Implementation uses `findOverlaps(query, subject, ignore.strand=FALSE)`, which only counts hits where strands agree. For antisense queries we flip the motif strand first, then call `findOverlaps` with the same `ignore.strand=FALSE` flag. This is encapsulated in `lib/strand_aware_overlaps()` (defined inline within `02b`).

If you later want to flip a motif's orientation (e.g., to ask "what if the motif is actually on the other strand"), you simply rename the columns: `_sense` ↔ `_antisense`. **No recomputation needed.** This was the central design constraint of this refactor.

## Cross-pipeline dependencies

This pipeline reads from:
- `01_resources/`: mm39 FASTA, GENCODE vM38 GTF, F121-9 + BL6xCAST het VCFs (UCSC-style)
- `03_atac/peaks/consensus/consensus_peaks.bed` — required for `atac_consensus_*` columns; if missing, those columns are NA
- `04c_rnaseq/expression/gene_expression_summary.tsv` — required for `expression_*_sense/antisense` columns; if missing, those columns are NA
- `04c_rnaseq/bam/star/*_SJ.out.tab` — required for `*_source = SJ_only / GTF_and_SJ` enrichment; if missing, splice-site sources reduce to GTF-only

`01_preprocess_resources.R` warns if any optional input is missing but continues; the per-motif annotation (`02_annotate_motifs.R`) similarly warns and writes NAs for missing categories.

## Resource expectations

On the user's HPC (typical):

| Stage | Time | Memory | Notes |
|-------|------|--------|-------|
| 01 preprocess | ~1 hour | 64 GB | TxDb build dominates; SJ merge is fast |
| 02 per-motif (each) | ~30–60 min | 16–24 GB | VCF I/O and the 40 splice queries are the heaviest |

The per-motif step is embarrassingly parallel; submit all motifs as a SLURM array with `--max-concurrent 1000` (default). Total wall time: 1–2 hours for ~600 motifs at 1000-wide concurrency, ~6–10 hours at 100-wide.

## Compatibility break with old `03_motif_annot/`

Old column names that no longer exist (you'll need to update downstream code):

| Old | New |
|-----|-----|
| `is_intragenic` | Split into `intragenic_sense` and `intragenic_antisense` |
| `intragenic_gene_types` | `intragenic_gene_types_sense`, `intragenic_gene_types_antisense` |
| `is_promoter` | `is_promoter_sense`, `is_promoter_antisense` |
| `genic_regions` | `genic_regions_sense`, `genic_regions_antisense` |
| `upstream_tss_*`, `downstream_tss_*` | All gain `_sense` / `_antisense` suffix |
| `atac_<sample>_overlap`, `atac_<sample>_signalValue`, etc. | Replaced by `atac_consensus_overlap` + `atac_consensus_n_replicates` |

NEW (didn't exist in old version):
- `intragenic_gene_ids_sense/_antisense` (10 expression cols depend on this)
- `expression_*` (10 cols)
- All 40 splice junction columns

## Troubleshooting

- **"GTF not found"** at preprocess time → check `GTF_FILE` in `config.sh`. Path should resolve from a compute node.
- **All splice columns `_source = "GTF"`** → no STAR `SJ.out.tab` files matched `RNASEQ_SJ_GLOB`, OR all junctions had < `SJ_MIN_UNIQUE_READS`. Lower the threshold in `config.sh` or check the glob.
- **All `expression_*_sense = NA`** → either `04c_rnaseq/expression/gene_expression_summary.tsv` doesn't exist, or `intragenic_gene_ids_sense` is itself NA (i.e., motif isn't intragenic on that strand). Both are expected behaviors; check whether the underlying motif is actually intragenic.
- **Per-motif task OOM at 20 GB** → bump `DEFAULT_MEM` in `config.sh` to 32 GB, or use the retry path (`0x_resubmit_missing_motifs.sh`) which submits to `medium` partition with 128 GB.
- **`No common chromosomes` for SNP overlap** → VCF is using Ensembl-style chrom names (`1`, `2`, ...) instead of UCSC (`chr1`, `chr2`, ...). The pipeline expects UCSC-style; convert with `bcftools annotate --rename-chrs` if needed.
