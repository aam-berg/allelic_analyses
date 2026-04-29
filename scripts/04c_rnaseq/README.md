# `04c_rnaseq/` — F121-9 RNA-seq Processing Pipeline

Processes paired-end stranded RNA-seq from GSE200699 (2 replicates of F121-9
mESCs, untreated WT) into:

1. **Gene expression table** — TPM mean/SD + raw counts per gene, used by
   `05_motif_annot/02f_gene_expression.R` for the `expression_*` columns.
2. **Strand-specific bigWigs** — for browser visualization and motif-aggregate
   plots.
3. **Allele-specific BAMs and counts** — for future allele-specific expression
   and PSI work in `06_allele_pairs/`.
4. **STAR splice-junction tables (`SJ.out.tab`)** — usable later as RNA-seq
   evidence for splice-junction proximity in `05_motif_annot/`. (For now,
   `05_motif_annot/` uses GTF-based junction coordinates only.)

## Input

GEO series GSE200699 (subset of GSE200702, F121-9 WT untreated reps):

| Replicate           | GSM         | SRRs (lane-split, concatenated)      |
|---------------------|-------------|---------------------------------------|
| WT_RNAseq_rep1      | GSM6042438  | SRR18735736, SRR18735737              |
| WT_RNAseq_rep2      | GSM6042439  | SRR18735734, SRR18735735              |

Multiple SRRs per sample = lane-split technical replicates. Step 01
concatenates them at the FASTQ level into a single `_R1.fastq.gz` /
`_R2.fastq.gz` pair per sample.

The canonical sample mapping is in `samples.tsv`. SRRs are comma-separated.

## Library protocol

**TruSeq Stranded Total RNA, paired-end, 100 nt reads.**
- Reverse-stranded: R1 maps to the antisense strand of the original RNA.
- featureCounts: `-s 2`
- STAR `--quantMode GeneCounts` 4-col table: read **column 4** (reverse).
- Strand-specific bigWigs:
  - RNA + strand fragment ⇔ R1 reverse + R2 forward (BAM)
  - RNA − strand fragment ⇔ R1 forward + R2 reverse (BAM)

## Output

```
${SCRATCH_DIR}/
├── fastq/                       # concatenated PE fastqs
├── trimmed/                     # adapter-trimmed PE fastqs
├── bam/
│   ├── star/                    # raw STAR output, with vW WASP tags
│   │   {SAMPLE}_Aligned.sortedByCoord.out.bam
│   │   {SAMPLE}_SJ.out.tab     # splice junctions
│   │   {SAMPLE}_ReadsPerGene.out.tab  # built-in stranded gene counts
│   └── final/{SAMPLE}_final.bam # WASP-passed reads only
├── quant/
│   ├── gene_level/{SAMPLE}_gene_counts.tsv
│   └── transcript_level/{SAMPLE}_tx_counts.tsv
├── bigwig/
│   ├── individual/{SAMPLE}_{plus,minus}.bw     # stranded RNA signal
│   ├── merged/merged_{plus,minus}.bw
│   ├── allele_specific/{SAMPLE}_{ref,alt}_{plus,minus}.bw
│   └── allele_merged/merged_{ref,alt}_{plus,minus}.bw
├── allele_specific/{SAMPLE}_{ref,alt,nosnp,ambiguous}.bam
├── expression/
│   └── gene_expression_summary.tsv          # ← used by 05_motif_annot/
└── qc/
    ├── pipeline_qc_report.txt
    └── allele_counts_merged.tsv
```

## Pipeline steps

| Step | Script                            | What it does                                              |
|------|-----------------------------------|-----------------------------------------------------------|
| 00   | `00_setup_references.sh`          | STAR index from mm39 + vM38 GTF (reuses if present)        |
| 01   | `01_download_fastq.sh`            | fasterq-dump + concat multi-SRRs per sample (PE)           |
| 02   | `02_trim_qc.sh`                   | cutadapt PE TruSeq, FastQC                                 |
| 03   | `03_align_star.sh`                | STAR + WASP via `vW` tag + 2-pass + GeneCounts             |
| 04   | `04_quantify.sh`                  | featureCounts gene + tx, stranded; per-sample TPM           |
| 05   | `05_make_bigwigs.sh`              | Strand-specific PE bigWigs (RNA + and − strand)             |
| 06   | `06_allele_specific.sh`           | Filter by WASP tag, split via lib/split_by_allele.py PE     |
| 07   | `07_make_bigwigs_allele.sh`       | Per-allele strand-specific bigWigs                          |
| 08   | `08_merge_qc.sh`                  | Cross-rep merging + final QC + expression summary table     |

## Running the pipeline

```bash
cd 04c_rnaseq
bash run_pipeline.sh                 # dry run (default)
bash run_pipeline.sh --submit        # actually submit
bash run_pipeline.sh --submit --start-from 03   # skip 00-02
```

Every step is idempotent — re-running skips completed work.

## Key design decisions

**STAR aligner.** RNA-seq must use a splice-aware aligner; STAR is the
field standard for general-purpose bulk RNA-seq and the alignment provides
exactly what we need for both gene quantification and the future splicing
analysis (BAMs + SJ.out.tab). bowtie2 is *wrong* for RNA-seq because
nascent reads span splice junctions.

**STAR's native WASP.** Instead of the explicit find/remap/filter cycle
used in `03_atac/` and `04_proseq/`, we use STAR's built-in WASP via
`--varVCFfile` + `--waspOutputMode SAMtag`. STAR adds a `vW:i:N` tag to
each read indicating WASP filter status:

| Tag value | Meaning |
|-----------|---------|
| `vW:i:1`  | Passed WASP test (kept for allele-specific analysis) |
| `vW:i:2-7`| Failed (different alignment when allele swapped, etc.) |
| no tag    | Read does not overlap any SNP (kept) |

For allele-specific analysis, step 06 keeps reads with `vW:i:1` or no
`vW` tag. For gene-level quantification (step 04), we use the unfiltered
BAM — since gene counting is not allele-specific, WASP filtering would
reduce statistical power without benefit.

**2-pass alignment (`--twopassMode Basic`).** STAR's first pass discovers
novel splice junctions; the second pass realigns reads using those
junctions plus the GTF-supplied junctions. Better sensitivity for novel
splicing events. The original GSE200699 paper used `--outFilterType
BySJout` which similarly favors junction-supported alignments.

**MAPQ ≥ 10 in step 04.** STAR's MAPQ scheme: 255 = unique, 3 = 2 multi,
1 = ≤9 multi, 0 = ≥10 multi. ≥10 keeps unique alignments only — strict
for allele-specific work; reasonable for stranded gene counting.

**Strand-specific bigWig logic for PE reverse-stranded library.** For
each fragment, we determine which strand of the original RNA it came from
based on the BAM flags of R1 and R2:

- RNA + strand fragments: R1 has flag 16 (reverse) AND R2 has flag NOT 16
- RNA − strand fragments: R1 has flag NOT 16 AND R2 has flag 16

In step 05/07 we extract these fragment sets with `samtools view` flag
filters, then convert to fragment-coverage bigWigs.

**Gene expression metrics.** Step 08 produces `gene_expression_summary.tsv`
with one row per gene and columns:
- `gene_id`, `gene_name`, `gene_strand`
- per-rep raw counts (one column per replicate)
- per-rep TPM (one column per replicate)
- `tpm_mean`, `tpm_sd`, `counts_sum`
- `tpm_log2_mean` (log2(tpm_mean + 1) for filtering)

Used by `05_motif_annot/02f_gene_expression.R` to add `expression_*_sense`
and `expression_*_antisense` columns to the per-motif annotation tables.

**Allele splitting in PE mode.** `lib/split_by_allele.py --paired-end` is
fragment-aware: classifies a fragment based on combined SNP evidence from
both mates. Identical to the ATAC pipeline's allele-split step.

## Cross-pipeline dependencies

This pipeline reads from:
- `01_resources/`: `mm39.fa`, `gencode.vM38.annotation.gtf.gz`,
  `F121-9_het_snps.ucsc.vcf.gz`
- `lib/`: `split_by_allele.py`, `make_wasp_snps.sh`

Outputs consumed by:
- `05_motif_annot/02f_gene_expression.R`: reads `expression/gene_expression_summary.tsv`
- `06_allele_pairs/05_splice_jx_alleles.R` (future): reads
  `allele_specific/{SAMPLE}_{ref,alt}.bam` and STAR's `SJ.out.tab` files

## Resource expectations

On a typical HPC node (8 cores, 64 GB):

| Step | Wall time | Memory | Notes                                        |
|------|-----------|--------|----------------------------------------------|
| 00   | 1-2 hr    | 64 GB  | STAR index gen for mm39 with sjdb            |
| 01   | 1-2 hr    | 16 GB  | Multi-SRR download + concat (I/O bound)      |
| 02   | 1-2 hr    | 8 GB   | cutadapt is fast                             |
| 03   | 4-8 hr    | 64 GB  | STAR with 2-pass + WASP                      |
| 04   | 30 min    | 8 GB   | featureCounts                                 |
| 05   | 1 hr      | 16 GB  | Strand-split + bigWig conversion             |
| 06   | 1 hr      | 16 GB  | WASP tag filter + allele split               |
| 07   | 1 hr      | 16 GB  | Per-allele strand-specific bigWigs           |
| 08   | 30 min    | 16 GB  | Cross-rep merging + expression summary       |

## Troubleshooting

- **STAR runs out of memory.** Increase `--mem` to 96G. mm39 index loading
  is ~30 GB; 2-pass + WASP can push peak usage higher.
- **`vW:i:6` tags everywhere.** STAR couldn't find genotype info for SNPs.
  The VCF needs a `GT` field (which our F121-9 het VCF has — `0|1` or `1|0`
  per SNP).
- **All reads tagged `vW:i:1`.** Means no SNPs are being detected by STAR.
  Verify the VCF is the chr-prefixed UCSC version, matching the BAM.
- **Strand check on bigWigs.** For housekeeping genes with known strand,
  the plus.bw and minus.bw should clearly differ. If not, the strand swap
  may be wrong — verify against IGV.
- **STAR not found.** `conda install -n proseq2_env -c bioconda star`
- **featureCounts not found.** `conda install -n proseq2_env -c bioconda subread`

## Future extensions

- Salmon/RSEM tx-level quantification with bias correction (currently
  using STAR + featureCounts)
- DESeq2 normalized expression (TPM is fine for filtering; DESeq2 is
  better for differential analysis)
- rMATS or DEXSeq splice junction usage (currently raw SJ.out.tab only)
- Allele-specific expression from `allele_specific/*.bam` (the BAMs are
  produced; analysis script lives in `06_allele_pairs/` when needed)
