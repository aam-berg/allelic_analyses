# 06_allele_pairs/

The allele-pair pipeline. This is where the central biological unit of the
pausing-phase project lives: the **(motif_hit, het_SNP) pair** in F121-9
(129S1×CAST) F1 hybrid mESCs. Every pair gets per-allele PWM scores,
per-allele ATAC accessibility, per-allele PRO-seq pause signal, and per-allele
junction-usage PSI for all 8 directional splice-site cells.

Outputs from this pipeline feed `07_analysis/` (regression, ARIMA aggregation,
plots, statistical tests).

---

## 1. What is an allele pair?

An allele pair is a tuple `(motif_hit, het_SNP)` where the SNP falls in or
near the motif body. "Near" = within `MAX_PAIR_FLANK_BP` (default 200 bp) of
either motif edge. SNPs further away than this are not paired with the motif
because they don't plausibly affect motif function.

A single SNP can pair with multiple motifs (if it falls in or near several);
a single motif can pair with multiple SNPs. The pair is the analysis unit;
the motif is a grouping key.

The **central data model** is the per-motif master pair table at
`${OUTDIR}/pair_tables/{motif_id}_pair_table.tsv.gz`. One row per pair.

---

## 2. Architecture (normalized schema)

This pipeline is **normalized**: it stores only pair-level columns plus a
`motif_hit_id` foreign key. Static genomic features (gene context, distance
to TSS, splice-site distances, etc., ~136 columns) live in
`05_motif_annot/{motif_id}_annotated.tsv.gz` and are joined at analysis time.

Why: motif features are 1:1 with `motif_hit_id`, so duplicating them across
all pairs of a motif would multiply file size by ~8-50× without adding
information. Joining at analysis time is one line in R/Python.

**The contract that makes this work** is the `motif_hit_id` key, defined as

```
motif_hit_id = paste(motif_id, chrom, motif_start_1based, strand, sep="__")
```

Both `06_allele_pairs/` and `05_motif_annot/` outputs include this column.
A small patch to `05_motif_annot/02_annotate_motifs.R` adds the column there
(see "patches" below).

---

## 3. Pipeline stages

```
01 build_pair_index        per-motif: cross-product of motif hits × het-SNPs in range
02 score_pwm               per-motif: PWM ref/alt scoring for body-SNPs
03 count_atac              per-motif: ATAC ref/alt counts at motif body
04 count_proseq_window     per-motif: PRO-seq pause + gene body counts (strand-aware)
04b extract_proseq_profile per-motif: ±500 bp per-base PRO-seq profile (ARIMA substrate)
04c compute_pause_indices  per-motif: local pause-index columns
05 extract_rnaseq_junc     per-motif: junction use/skip counts at the 8 cells
05b compute_psi            per-motif: PSI per allele per cell
06 assemble_pair_table     per-motif: master per-pair table
07 pair_qc                 single:    cross-motif QC report
```

Run via SLURM array fan-out (one task per motif):

```bash
bash 00_setup.sh              # validate upstream inputs
bash run_pipeline.sh          # dry run
bash run_pipeline.sh --apply  # submit
```

---

## 4. Master pair table — schema

`${OUTDIR}/pair_tables/{motif_id}_pair_table.tsv.gz`. One row per pair.

### 4.1 Pair identity & geometry (~20 cols)

| Column                                | Type   | Notes |
|---------------------------------------|--------|-------|
| `pair_id`                             | str    | `{motif_hit_id}__{snp_chrom}__{snp_pos}` |
| `motif_hit_id`                        | str    | foreign key to `05_motif_annot/`. `{motif_id}__{chrom}__{start_1b}__{strand}`. |
| `motif_id`                            | str    | e.g. `AC0001` |
| `motif_chrom`, `motif_start`, `motif_end`, `motif_strand`, `motif_score`, `motif_width` | mixed | 1-based inclusive coords |
| `snp_chrom`, `snp_pos`, `snp_ref`, `snp_alt` | mixed | SNP coords (1-based); biallelic SNVs only |
| `snp_offset_from_motif_center_bp`     | int    | **SIGNED**. Negative = upstream (5'); positive = downstream (3'); in transcription direction. |
| `snp_offset_from_motif_center_abs_bp` | int    | always positive |
| `snp_distance_from_motif_edge_bp`     | int    | 0 if SNP is in body |
| `snp_in_motif_body`                   | bool   | |
| `snp_upstream_of_motif`               | bool   | NA if SNP is in body. In transcription direction. |
| `snp_relative_position`               | factor | `body` / `upstream_flank` / `downstream_flank` |
| `snp_position_in_motif`               | int    | 1-based position within motif, in transcription direction (NA if not in body) |

**Why directional offsets?** Pausing has direction. A SNP at +100 bp
(downstream) of a motif on the "+" strand is mechanistically very different
from a SNP at -100 bp (upstream). All offsets are computed in the
transcription direction so downstream code never has to worry about strand
bookkeeping.

### 4.2 PWM scoring (~6 cols, body-SNPs only)

| Column                  | Type    | Notes |
|-------------------------|---------|-------|
| `pwm_scoring_applicable` | bool   | TRUE only if `snp_in_motif_body && motif_width == PWM_width` |
| `motif_sequence_ref`    | str     | motif span as REF, in transcription direction |
| `motif_sequence_alt`    | str     | same with SNP swapped to ALT |
| `pwm_score_ref`         | float   | log-odds against uniform background |
| `pwm_score_alt`         | float   | same |
| `delta_pwm_score`       | float   | `pwm_score_alt - pwm_score_ref` |

**Score interpretation:** higher = better motif match. `delta_pwm_score > 0`
means ALT allele binds the TF more tightly than REF (according to the PWM).
Pseudocount = 0.01 to avoid log(0).

For flank SNPs, all five PWM columns are NA. ATAC accessibility asymmetry is
the proxy for the affinity effect there.

### 4.3 ATAC allele-specific (~10 cols)

Per-motif (broadcast to all pairs of that motif).

| Column            | Notes                                                      |
|-------------------|------------------------------------------------------------|
| `atac_ref_total`, `atac_alt_total` | summed across reps                            |
| `{rep}_ref`, `{rep}_alt`           | per-rep counts                                |

Long-format companion at `${OUTDIR}/atac_counts/{motif_id}_atac_long.tsv.gz`
(`motif_hit_id × rep × allele → count`) for proper stat tests with replicate
structure.

### 4.4 PRO-seq window counts (~20 cols)

Per-motif (broadcast). Strand-aware — counts only reads on the motif's
transcribing strand.

| Column                            | Notes |
|-----------------------------------|-------|
| `proseq_pause_ref/alt`            | reads in `motif_center ± PAUSE_WINDOW_HALF_BP` (default ±10 bp) |
| `proseq_genebody_ref/alt`         | reads in `GENE_BODY_FLANK_BP` (default 2000 bp) downstream of motif on same strand |
| `{rep}_{ref/alt}_{pause/genebody}` | per-rep counts |

Long-format companion at
`${OUTDIR}/proseq_counts/{motif_id}_proseq_window_long.tsv.gz`.

### 4.5 Pause indices (~12 cols)

| Column                                 | Notes |
|----------------------------------------|-------|
| `pause_density_ref/alt`                | reads/bp in pause window |
| `genebody_density_ref/alt`             | reads/bp in gene body window |
| `pause_index_ref/alt`                  | `pause_density / max(genebody_density, eps)` |
| `log2_pause_index_ratio_alt_over_ref`  | log2 ratio |
| `pause_total_reads`                    | both alleles, both windows of the pause |
| `genebody_total_reads`                 | both alleles, gene body |
| `pause_min_reads_per_allele`           | `min(pause_ref, pause_alt)` — power filter |
| `genebody_min_reads_per_allele`        | `min(genebody_ref, genebody_alt)` |

**Important:** these "local pause indices" are mathematically analogous to
classical promoter-proximal pause indices, but they're applied at intragenic
TFBSs — a different mechanistic phenomenon. The goal is per-site
quantification when reads happen to be deep enough; for most sites we'll
rely on the aggregate ARIMA approach (see §6) instead.

### 4.6 PSI per allele × 8 directional cells (~~~130 cols)

For each motif, 05_motif_annot annotates 8 nearest splice sites (donor /
acceptor × sense / antisense × upstream / downstream of motif). For each cell
we compute per-allele PSI from RNA-seq.

For each of the 8 cells `{cell}` in
```
donor_sense_upstream, donor_sense_downstream,
donor_antisense_upstream, donor_antisense_downstream,
acceptor_sense_upstream, acceptor_sense_downstream,
acceptor_antisense_upstream, acceptor_antisense_downstream
```
the table has ~16 columns prefixed with the cell name:

| Column suffix              | Notes |
|----------------------------|-------|
| `_site_chrom`, `_site_pos`, `_site_strand` | reconstructed from motif geometry + dist + direction |
| `_site_type`               | `donor` or `acceptor` (redundant with cell name) |
| `_source`                  | `GTF` / `SJ_only` / `GTF_and_SJ` (from STAR SJ filtering) |
| `_motif_type`              | splice motif (GT-AG, etc.) |
| `_distance`                | bp from motif center to site, in transcription direction |
| `_total_use_ref/alt`       | reads spliced through this site (one CIGAR-N boundary at site_pos) |
| `_total_skip_ref/alt`      | reads continuous across this site (no N intron at site_pos) |
| `_ref_total`, `_alt_total` | use + skip per allele (denominator) |
| `_psi_ref`, `_psi_alt`     | `use / (use + skip)` if ≥ `PSI_MIN_TOTAL_READS` reads, else NA |
| `_delta_psi`               | `psi_alt - psi_ref` |

**Interpretation:** PSI here is junction-centric — "given reads that cover
this splice site, what fraction splices through it?" — which is the right
question for asking how pausing affects co-transcriptional splice choice.
For sense-strand cells, this is the primary signal. Antisense cells serve
as background / convergent-transcription controls.

A long-form variant `(motif_hit_id × cell × ...)` is at
`${OUTDIR}/psi/{motif_id}_psi_long.tsv.gz` (~20 cols, 8 rows per motif).

### 4.7 Total master-table column count

About 175-185 columns per motif, dominated by the 130 PSI columns
(8 cells × 16 sub-columns). The static 136-column motif annotation from
`05_motif_annot/` is joined at analysis time.

---

## 5. Companion long-format tables

Beyond the master per-pair tables, the pipeline writes long-format
companions designed for replicate-aware statistical testing and aggregation.

| Path                                                              | Schema                                    |
|-------------------------------------------------------------------|-------------------------------------------|
| `${OUTDIR}/atac_counts/{motif_id}_atac_long.tsv.gz`               | `motif_hit_id, rep, allele, count`        |
| `${OUTDIR}/proseq_counts/{motif_id}_proseq_window_long.tsv.gz`    | `motif_hit_id, rep, allele, window_type, count` |
| `${OUTDIR}/proseq_profile/{motif_id}_proseq_profile.tsv.gz`       | `motif_hit_id, rep, allele, strand_relative, position_offset, count` |
| `${OUTDIR}/rnaseq_junctions/{motif_id}_rnaseq_junctions.tsv.gz`   | per-motif × cell with all per-rep columns |

The PRO-seq per-base profile (`04b`) is the substrate for ARIMA aggregation
— see next section.

---

## 6. Downstream analysis hooks

This pipeline is the **substrate**. Real biological tests live in
`07_analysis/`. But the pipeline is designed around the analyses we'll need.

### 6.1 Aggregate ARIMA pause peaks (the key downstream method)

PRO-seq drop-out is high enough that single-site pause peaks are usually
invisible at intragenic TFBSs. The standard remedy is to **aggregate many
sites**: average the per-base normalized PRO-seq signal across all motifs in
a partition, then apply ARIMA outlier detection to the average curve.

The substrate is `${OUTDIR}/proseq_profile/{motif_id}_proseq_profile.tsv.gz`
which gives, for each motif, single-bp PRO-seq counts in
`[motif_center − 500, motif_center + 500]` (transcription-direction signed
offsets), per rep, per allele, per strand-relative-to-motif (sense /
antisense).

Typical analysis recipe:

1. Define a partition of pairs (e.g. quartiles by `delta_pwm_score`,
   or quartiles by `log2(atac_ref_total/atac_alt_total)`).
2. For each partition, average per-base PRO-seq counts (with library-size
   normalization) across all motifs × pairs in that partition. Do this
   separately per allele.
3. Apply ARIMA outlier detection at the motif-center vicinity to score
   "peak strength" per partition per allele.
4. Compare ARIMA scores across partitions and across alleles. The hypothesis
   predicts a monotonic relationship between affinity (or accessibility)
   delta and ARIMA peak strength.

The aggregation logic, library-size normalization, and ARIMA fitting all
live in `07_analysis/` — they're parameter choices, not pipeline outputs.
The pipeline guarantees that for any partition you can build (e.g. via a
data.table groupby on the master pair table), you have the per-base profile
ready to aggregate.

### 6.2 Site-specific pause index analyses

For motif × SNP combinations with deep enough coverage:

```r
# delta-PWM vs delta-pause (log-ratio):
plot(pair_table$delta_pwm_score,
     pair_table$log2_pause_index_ratio_alt_over_ref,
     subset = pause_min_reads_per_allele >= 10 & pwm_scoring_applicable)
```

Filter on `pause_min_reads_per_allele` (≥10 is conservative for these
narrow windows; you can sweep this threshold).

### 6.3 ATAC asymmetry vs delta-pause

```r
pair_table[, atac_log2_ref_alt := log2((atac_ref_total + 1) / (atac_alt_total + 1))]
plot(pair_table$atac_log2_ref_alt,
     pair_table$log2_pause_index_ratio_alt_over_ref)
```

This is the affinity-as-accessibility test. ATAC measures total binding,
which (independent of PWM) reflects in vivo TF occupancy.

### 6.4 Pause asymmetry → PSI asymmetry

The headline pause-affects-splicing test, focused on the most biologically
motivated cell:

```r
plot(pair_table$log2_pause_index_ratio_alt_over_ref,
     pair_table$donor_sense_downstream_delta_psi)
```

Subset on pairs where:
- `pause_min_reads_per_allele >= 10`
- `donor_sense_downstream_psi_min_reads_per_allele >= 4`

For sense-direction donors specifically, downstream donors are the
co-transcriptional splice acceptor for the closest exon-exon boundary
downstream of the pause site. Predicted direction: more pausing → more time
for splicing decisions → could go either way (more or less inclusion
depending on splice site strength), so the test is two-sided.

### 6.5 Stratified context analyses

The static genomic-context columns from `05_motif_annot/` (gene region,
TSS distance, exon/intron, distance to nearest junction) facet every plot
above. Standard recipe:

```r
master <- fread("...pair_tables/AC0001_pair_table.tsv.gz")
annot  <- fread("...annot_motifs_v2/AC0001_annotated.tsv.gz")
master <- annot[master, on = "motif_hit_id"]
# Now master has all 175 + 136 = ~310 columns.
ggplot(master) +
  geom_point(aes(delta_pwm_score, log2_pause_index_ratio_alt_over_ref)) +
  facet_wrap(~ gene_region)   # exon vs intron vs intergenic
```

---

## 7. Resource expectations

For ~700 motif archetypes with ~1k motif hits each (~700k motif hits,
~3-5M pairs) at typical F121-9 het-SNP density:

| Stage                  | Time/motif | Memory  | Notes                            |
|------------------------|-----------|---------|----------------------------------|
| 01 build_pair_index    | 5-30 min  | 4-8 GB  | I/O: VCF region reads            |
| 02 score_pwm           | 1-3 min   | 2 GB    | Genome FA reads via Rsamtools    |
| 03 count_atac          | 5-15 min  | 4 GB    | bedtools coverage on 4 reps × 2 alleles |
| 04 count_proseq        | 10-30 min | 4 GB    | strand-aware samtools view filter |
| 04b extract_profile    | 30-90 min | 4 GB    | pyBigWig path; 2-4× slower if BAM fallback |
| 04c pause_indices      | <30 sec   | 2 GB    | trivial reduction                 |
| 05 extract_junctions   | 30-90 min | 4 GB    | pysam read-by-read                |
| 05b compute_psi        | <30 sec   | 2 GB    | trivial reduction                 |
| 06 assemble            | 1-2 min   | 4 GB    | data.table joins                  |
| 07 qc                  | 10-30 min (single)| 16 GB | reads all per-motif tables    |

Total wall-clock for the full pipeline at ~700 motifs with reasonable SLURM
parallelism (`MAX_CONCURRENT=200`): 2-6 hours, dominated by 04b and 05.

---

## 8. Cross-pipeline dependencies

This pipeline reads outputs from:

- **`02_motif_scanning/`**: per-motif BED files (motif hits)
- **`03_atac/`**: `${SCRATCH}/atac_proc/allele_specific/{SAMPLE}_{ref,alt}.bam`
- **`04_proseq/`**:
  - `${SCRATCH}/proseq_proc/allele_specific/{SAMPLE}_{ref,alt}.bam`
  - `${SCRATCH}/proseq_proc/bigwig/allele_specific/{SAMPLE}_{ref,alt}_{plus,minus}.bw`
- **`04c_rnaseq/`**: `${SCRATCH}/rnaseq/allele_specific/{SAMPLE}_{ref,alt}.bam`
- **`05_motif_annot/`**:
  - `${SCRATCH}/annot_motifs_v2/{motif_id}_annotated.tsv.gz` (for the 8
    junction cells per motif)
  - `${SCRATCH}/annot_motifs_v2/preprocessed/splice_donors_gr.rds` (optional)
- **External**: `/home/alb1273/pausing_phase_project/resources/consensus_pwms.meme`

Run `bash 00_setup.sh` to validate that all inputs exist before submitting
any heavy jobs.

---

## 9. Required environment additions

Beyond what `proseq2_env` and `test_m` already have:

```bash
# pyBigWig (for fast per-base profile extraction in 04b)
conda install -n proseq2_env -c bioconda pyBigWig

# pysam (for junction-read classification in 05; usually already present)
conda install -n proseq2_env -c bioconda pysam
```

If `pyBigWig` is missing, `04b` falls back to BAM-based extraction
(2-4× slower but functionally equivalent).

R packages (in `test_m`): `optparse`, `GenomicRanges`, `VariantAnnotation`,
`Biostrings`, `Rsamtools`, `data.table`. All standard Bioconductor.

---

## 10. Patches required for upstream pipelines

### 10.1 `05_motif_annot/02_annotate_motifs.R`: add `motif_hit_id` column

The normalized join contract requires `motif_hit_id` in both the pair
tables and the annotated-motif tables. The pair tables already have it
(produced by `01_build_pair_index.R` here). The annotated-motif tables need
a one-line backward-compatible addition.

The patch is delivered alongside this pipeline as
`05_motif_annot_motif_hit_id.patch` — apply with:

```bash
cd scripts/05_motif_annot
patch < ../06_allele_pairs/05_motif_annot_motif_hit_id.patch
```

After the patch, re-run `05_motif_annot/run_pipeline.sh --apply`. This will
overwrite `${SCRATCH}/annot_motifs_v2/*.tsv.gz` with the new column added —
the rest of the schema is unchanged.

If you don't apply the patch, `06_allele_pairs/` still works (it
reconstructs `motif_hit_id` from `chrom`/`start`/`strand` columns
on-the-fly in `05_extract_rnaseq_junctions.sh`). But analysis-time joins
against the annotated tables become more awkward.

---

## 11. Troubleshooting

**"No common chroms between motifs and VCF" error in stage 01**
  Chromosome naming mismatch (UCSC `chr1` vs Ensembl `1`). The F121-9 het
  VCF should already be UCSC-named (`F121-9_het_snps.ucsc.vcf.gz`). If not,
  rename via `bcftools annotate --rename-chrs`.

**`pyBigWig` import error in stage 04b**
  Falls back to BAM extraction automatically. Install `pyBigWig` in
  `proseq2_env` for ~3× speedup.

**Stage 04 reports zero PRO-seq counts at obvious motifs**
  Check (a) PE detection (printed at start of stage 04): if you have PE data
  but `IS_PE=0`, the strand filtering is wrong. (b) Canonical strand swap:
  the PRO-seq BAMs from `04_proseq/` should already have reads on the
  transcribing strand (not the originally-sequenced antisense strand). If
  they don't, `04_proseq/`'s `06_swap_strand` step didn't run.

**Stage 05 reports zero PSI signal**
  RNA-seq is much sparser than PRO-seq for splicing reads. Check
  `psi_min_reads_per_allele` — if most cells have <4 reads/allele, you may
  need to either lower `PSI_MIN_TOTAL_READS` or accept that single-site PSI
  isn't powered for most pairs.

**Stage 06 produces tables with NA in unexpected columns**
  Check that all upstream stages completed successfully. The orchestrator's
  dependency chain should prevent this, but in case of interrupted runs
  (e.g., `--start-from 06`), some inputs may be missing. Re-run from
  the earliest interrupted stage.

**Disk usage**
  Profile tables (04b) are the largest output; ~50 GB total at ~700 motifs
  × ~1k hits/motif × 1001 positions × 4 reps × 2 alleles × 2 strands. Plan
  scratch space accordingly.

---

## 12. First-run verification checklist

1. **Stage 00**: `bash 00_setup.sh` shows `[OK]` for all required inputs.
2. **Stage 01 dry-run on one motif**:
   `Rscript 01_build_pair_index.R --motif_bed AC0001.bed --vcf ... --motif_id AC0001 --outdir /tmp/test_pi/`
   Inspect output for sane offset distributions.
3. **Stage 02 dry-run**: same — confirm `pwm_score_ref` looks reasonable
   (>0 for body matches; PWM-width mismatch warnings expected only for
   archetype-set inconsistencies).
4. **Stage 04 strand sanity**: for a motif with known strong intragenic
   pausing, manually verify counts on transcribing strand >> antisense.
5. **Stage 04b first motif**: spot-check that `position_offset = 0` reads
   align with motif center, and that "sense" strand has more signal than
   "antisense" for true motif hits.
6. **Stage 05 first motif**: verify that pairs with high `_total_use_ref` +
   `_total_use_alt` correspond to known constitutive junctions in IGV.
7. **Stage 07 QC**: `pair_qc_global.txt` reports a sane fraction of pairs
   with measurable PWM, ATAC, PRO-seq, and PSI signal — typically:
   - 20-40% with body-SNP PWM scoring
   - 60-80% with non-zero ATAC reads
   - 30-50% with non-zero pause-window reads
   - 5-15% with PSI in at least one sense-direction cell

---

## 13. File listing

```
06_allele_pairs/
├── README.md                       (this file)
├── config.sh                       paths, sample names, parameter defaults
├── 00_setup.sh                     pre-flight input validation
├── 01_build_pair_index.R           per-motif: pair geometry
├── 02_score_pwm_per_allele.R       per-motif: PWM ref/alt for body-SNPs
├── 03_count_atac_alleles.sh        per-motif: ATAC ref/alt counts
├── 04_count_proseq_alleles.sh      per-motif: PRO-seq window counts
├── 04b_extract_proseq_profile.sh   per-motif: per-base profile (ARIMA substrate)
├── 04c_compute_pause_indices.R     per-motif: local pause indices
├── 05_extract_rnaseq_junctions.sh  per-motif: junction use/skip counts
├── 05b_compute_psi_per_allele.R    per-motif: PSI per cell
├── 06_assemble_pair_table.R        per-motif: master pair table
├── 07_pair_qc.R                    cross-motif QC report
├── run_pipeline.sh                 SLURM orchestrator
└── 05_motif_annot_motif_hit_id.patch   one-line patch for upstream
```
