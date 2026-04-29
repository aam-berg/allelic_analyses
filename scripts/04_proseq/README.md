# `04_proseq/` — F121-9 PRO-seq Processing Pipeline

Processes raw single-end PRO-seq from GSE190548 (4 replicates of F121-9
mESCs, untreated) into:

1. **Strand-specific bigWigs** — for browser visualization and motif-aggregate plots.
2. **Allele-specific BAMs and bigWigs** — used by `06_allele_pairs/` to
   quantify allele-specific RNAPII pausing in motif windows.

## Input

GEO series GSE190548, 4 F121-9 WT mESC PRO-seq replicates:

| Replicate    | SRR          |
|--------------|--------------|
| rep1         | SRR17640044  |
| rep2         | SRR17640045  |
| rep3         | SRR17640046  |
| rep4         | SRR17640047  |

The canonical sample-to-SRR mapping is in `samples.tsv`. Edit it to add or
swap samples — nothing else needs changing.

## Output

```
${SCRATCH_DIR}/
├── fastq/                          # raw + downloaded FASTQ
├── trimmed/                        # 3'-adapter trimmed FASTQ
├── bam/
│   ├── rRNA/                       # discarded rRNA reads (not kept)
│   ├── dm6/                        # only if DM6_SPIKEIN_FILTER=true
│   └── mm39/   {SAMPLE}_final.bam   # MAPQ + chrM filtered, primary BAM
├── bigwig/
│   ├── individual/   {SAMPLE}_{3,5}prime_{plus,minus}.bw   # per-rep, 4 files each
│   ├── merged/       merged_{3,5}prime_{plus,minus}.bw     # sum across all reps
│   ├── allele_specific/  {SAMPLE}_{ref,alt}_3prime_{plus,minus}.bw
│   └── allele_merged/    merged_{ref,alt}_3prime_{plus,minus}.bw
├── wasp/
│   ├── snp_files/   chr*.snps.txt.gz   (built by lib/make_wasp_snps.sh)
│   ├── filtered_bam/{SAMPLE}_wasp.bam
│   └── intermediate/...
├── allele_specific/{SAMPLE}_{ref,alt,nosnp,ambiguous}.bam
└── qc/
    ├── pipeline_qc_report.txt      # final compiled QC
    └── allele_counts_merged.tsv    # per-SNP allelic read counts across all reps
```

## Pipeline steps

| Step  | Script                          | What it does                                         | Mode             |
|-------|---------------------------------|------------------------------------------------------|------------------|
| 00    | `00_setup_references.sh`        | Reuse mm39/dm6 if present, build rRNA index         | single           |
| 01    | `01_download_fastq.sh`          | fasterq-dump (SE), FastQC raw                        | array (4)        |
| 02    | `02_trim_qc.sh`                 | cutadapt 3' adapter, FastQC trimmed                  | array (4)        |
| 03    | `03_align_filter.sh`            | rRNA → dm6* → mm39 → MAPQ + chrM filter             | array (4)        |
| 04    | `04_make_bigwigs.sh`            | Per-rep strand-specific bigWigs (3'/5' × +/−)       | array (4)        |
| 05    | `05_wasp_setup.sh`              | Clone WASP, build SNP files                          | single           |
| 06    | `06_wasp_filter.sh`             | WASP mapping bias removal (SE)                       | array (4)        |
| 07    | `07_allele_specific.sh`         | Split by allele via lib/split_by_allele.py          | array (4)        |
| 08    | `08_make_bigwigs_allele.sh`     | Per-allele 3'-end bigWigs (ref/alt × +/−)           | array (4)        |
| 09    | `09_merge_qc.sh`                | Cross-rep merging + final QC + merged allele counts  | single           |

\* Step 03's dm6 stage is skipped if `DM6_SPIKEIN_FILTER=false` in `config.sh`.

## Running the pipeline

**Dry run (default):**
```bash
cd 04_proseq
bash run_pipeline.sh
```

**Submit:**
```bash
bash run_pipeline.sh --submit
```

**Restart from a specific stage:**
```bash
bash run_pipeline.sh --submit --start-from 05
```

Every step is **idempotent**: it checks for its expected output and skips if
present. To force a re-run of a stage, delete the relevant output files first.

## Key design decisions

**No deduplication.** PRO-seq measures Pol II positions at single-nucleotide
resolution. At strong pause sites, many Pol II molecules sit at the *exact
same* nucleotide. Position-based duplicate removal would collapse real
biological signal into one read. This pipeline therefore does no dedup —
neither MarkDuplicates nor `samtools markdup -r`. (The earlier version of
this pipeline had an optional 5'-dedup flag that was removed in the refactor;
it was never appropriate for this analysis.)

**MAPQ ≥ 10.** PRO-seq community standard, distinct from ATAC's MAPQ ≥ 30.
Lower because SE reads carry less positional information than PE, and a
stricter cutoff would remove too much signal at multi-copy regions like
microsatellites and pericentromeric repeats.

**rRNA filter is `--sensitive` (bowtie2 default), not `--very-fast`.** For a
filter, we want high sensitivity so rRNA reads don't leak through to the
mm39 alignment. The earlier pipeline used `--very-fast`; the refactor
switched to defaults (essentially `--sensitive`) for safety. Runtime hit
is small (a few minutes per sample).

**Drosophila spike-in filter is optional.** Controlled by
`DM6_SPIKEIN_FILTER` in `config.sh`. The `03_align_filter.sh` script aligns
to dm6 first and discards aligned reads only if this is `true`. To verify
spike-in presence in your data, run on one sample and check the dm6
alignment rate in `qc/alignment_stats_*.tsv`. Below 0.5% means there's
effectively no spike-in present and the filter can be safely disabled.

**Strand swap is canonical.** PRO-seq reads are sequenced from cDNA, the
reverse complement of the nascent RNA. The 5' end of the read corresponds
to the 3' end of the RNA, which is the RNAPII active site. The full
reasoning lives in `scripts/README.md` section 7 (canonical reference).
The bigWig generation in step 04 and the allele-specific quantification in
`06_allele_pairs/` both follow this convention.

**chrM removal.** Mitochondrial transcription is heavy and skews
normalization. Removed in step 03 alongside the MAPQ filter.

**No 5'/3' UMI handling.** This dataset's protocol does not use UMIs
(verified from GEO metadata). If a future PRO-seq dataset does, UMI
extraction would slot in at the start of step 02.

**Shared `lib/` utilities.** `split_by_allele.py` and `make_wasp_snps.sh`
live in `scripts/lib/` and are shared with `03_atac/`. Both pipelines call
the same utility, with the SE/PE difference resolved at the
`split_by_allele.py` invocation site (`--single-end` here, `--paired-end`
in ATAC).

## Cross-pipeline dependencies

This pipeline reads from:

- `01_resources/`: `mm39.fa`, `mm39.chrom.sizes`, `F121-9_het_snps.ucsc.vcf.gz`
- `lib/`: `make_wasp_snps.sh`, `split_by_allele.py`
- `03_atac/` (optional): the mm39 bowtie2 index, if already built (reused if
  present; built fresh if not)

Outputs consumed by:

- `06_allele_pairs/04_proseq_alleles_at_motifs.R`: reads
  `allele_specific/{SAMPLE}_{ref,alt}.bam` and the per-allele bigWigs

## Resource expectations

On a typical HPC node (8 cores, 32 GB):

| Step  | Wall time per sample | Memory   | Notes                                      |
|-------|----------------------|----------|--------------------------------------------|
| 00    | 0–60 min             | 16 GB    | Most time is rRNA fetch + index build      |
| 01    | 30–60 min            | 8 GB     | I/O bound (SRA download)                   |
| 02    | 5–15 min             | 4 GB     | cutadapt is fast                           |
| 03    | 2–4 hr               | 16 GB    | Three sequential bowtie2 runs (rRNA, dm6, mm39) |
| 04    | 30 min               | 8 GB     | bedtools + bedGraphToBigWig                |
| 05    | 5 min                | 4 GB     | One-off setup                              |
| 06    | 2–4 hr               | 32 GB    | Re-alignment doubles bowtie2 work          |
| 07    | 30 min               | 8 GB     | Pure I/O + classification                  |
| 08    | 30 min               | 8 GB     | bedtools + bedGraphToBigWig per allele     |
| 09    | 30–60 min            | 16 GB    | bigWig merging + Python allele-count merge |

## Troubleshooting

- **dm6 alignment rate is ~0%.** No spike-in present. Set
  `DM6_SPIKEIN_FILTER=false` in `config.sh` and re-run from step 03.
- **rRNA fraction is suspiciously high (>30%).** Check the rRNA reference
  in `${RRNA_FA}`; ensure it includes 28S/18S/5.8S/5S and 45S pre-rRNA.
- **WASP removes a very high fraction of reads (>25%).** Suggests highly
  divergent regions. Inspect the WASP QC stats in `qc/wasp_stats_*.tsv`
  and the strand balance after WASP.
- **Strand ratio (+/−) is far from 1.0.** Could indicate library prep
  issues. Healthy PRO-seq has ratio 0.8–1.2 genome-wide.
- **Step 04 bigWig file sizes for plus and minus are very different.** Same
  cause as strand-ratio issue above — investigate.

## Future extensions

- UMI handling for protocols that use them.
- Spike-in normalization factors (currently bigWigs are raw counts; a
  spike-in-derived scaling factor could be applied at the merge stage).
- Strand-specific consensus peak calls if a peak-calling step is ever needed
  (PRO-seq peak callers like dREG could slot in alongside step 04).
