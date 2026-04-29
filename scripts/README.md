# Pausing Phase Project — Pipeline Architecture

This repository processes raw sequencing data and motif annotations into a final,
analysis-ready table of **(motif hit, heterozygous SNP) pairs** annotated with
allele-specific chromatin accessibility, RNAPII pausing, and splice junction
proximity, in F121-9 (129S1 × CAST/EiJ) F1 hybrid mouse embryonic stem cells.

---

## 1. Scientific goal

We want to use F1 hybrid heterozygous SNPs as **natural perturbations** to ask
whether allelic differences in TFBS chromatin accessibility or motif affinity
*cause* allelic differences in:

(a) RNAPII pausing at intragenic transcription factor binding sites, and
(b) the usage of splice junctions a few hundred base pairs downstream.

Because the two alleles in an F1 share *cellular context, trans factors, and
distal regulatory state*, a within-cell allelic comparison controls away most
confounders that plague between-sample comparisons. The unit of analysis is
the **allele pair**: a (motif hit, het SNP) tuple where the two alleles differ
at exactly one position.

## 2. The "allele pair" concept

Every downstream statistical question is phrased over allele pairs. An allele
pair is one row in the master table produced by `06_allele_pairs/`, with the
schema (described in detail in `06_allele_pairs/README.md`):

```
motif_id    chrom    motif_start    motif_end    motif_strand
gene_strand    gene_id              # the host gene's strand and ID (for intragenic)
snp_chrom   snp_pos  snp_ref  snp_alt   # the het SNP defining the pair
score_ref   score_alt    delta_score   # affinity proxy (PWM rescore by default)
atac_ref_count   atac_alt_count        # F121-9 ATAC reads per allele in window
proseq_pause_ref proseq_pause_alt      # F121-9 PRO-seq pause-window counts per allele
nearest_3pSS_dist  nearest_3pSS_id     # splice junction proximity (annotation only)
nearest_3pSS_psi_ref  nearest_3pSS_psi_alt   # PLACEHOLDER — populated when RNA-seq arrives
... etc.
```

**Two SNP-to-motif relations matter and the pipeline supports both:**

- *Direct overlap*: the het SNP sits inside the motif. This is the strongest
  candidate for affinity-driven causation (and is where allelic motif rescoring
  changes the score).
- *Nearby SNP*: the het SNP is within a configurable flanking window
  (10/25/50/100 bp etc.) of the motif. Useful for chromatin-accessibility
  questions, where a SNP that disrupts a neighboring TF site or a nucleosome
  positioning element can change accessibility at our motif of interest. This
  flanking machinery is already implemented in `02_annotate_motifs.R` for SNP
  *annotation*; it is repurposed at the *allele-quantification* stage in
  `06_allele_pairs/03_atac_alleles_at_motifs.R`.

## 3. Directory layout

```
pausing_phase_project/scripts/
├── README.md                      # ← this file
├── lib/                           # shared utilities (see lib/README.md)
│   ├── README.md
│   ├── split_by_allele.py         # used by 03_atac/, 04_proseq/
│   ├── make_wasp_snps.sh          # used by 03_atac/, 04_proseq/
│   └── strand_aware_count.py      # used by 06_allele_pairs/
├── 01_resources/                  # genome, VCF, GTF prep                          [DONE, frozen]
├── 02_motif_scanning/             # MOODS genome-wide scan                          [DONE, frozen]
├── 03_atac/                       # F121-9 ATAC processing → per-allele BAMs       [PLANNED]
├── 04_proseq/                     # F121-9 PRO-seq processing → per-allele BAMs    [REFACTOR]
├── 05_motif_annot/                # positional annotations only (no signal counts) [REFACTOR from 03_motif_annot/]
├── 06_allele_pairs/               # the join layer — produces the master table     [PLANNED]
└── 07_analysis/                   # final stats and plots                           [PLANNED]
```

A note on numbering: the rename from `03_motif_annot/ → 05_motif_annot/` and the
insertion of `03_atac/` are deliberate. Pipelines that ingest sequencing reads
(`03_atac/`, `04_proseq/`) come before the annotation join (`05_motif_annot/`,
`06_allele_pairs/`). The old `03_motif_annot/` directory can be archived once
the refactor is complete.

## 4. Cross-stage dependencies (DAG)

```
                                  ┌────────────────────────────────────────────┐
                                  │  02_motif_scanning                         │
                                  │  per-archetype BED files, one row per hit  │
                                  └─────────────────────┬──────────────────────┘
                                                        │
                                                        ▼
  ┌──────────────────────────────┐         ┌──────────────────────────────────┐
  │ 01_resources                 │         │ 05_motif_annot                   │
  │   • mm39 genome FASTA        │────────▶│   • intragenic / promoter / UTR  │
  │   • F121-9 het VCF (UCSC)    │         │   • TSS distances (strand-aware) │
  │   • GENCODE vM38 GTF         │         │   • splice junction proximity    │
  │   • mm10→mm39 chain          │         │   • ATAC peak overlap (boolean)  │
  └─────────────┬────────────────┘         │   • SNP overlap (direct + flank) │
                │                          └──────────────────┬───────────────┘
                │                                             │
                ▼                                             │
  ┌──────────────────────────────┐                            │
  │ 03_atac    (PE reads)        │                            │
  │   trim → bowtie2 → MAPQ →    │                            │
  │   chrM filter → MACS2 peaks  │                            │
  │   → WASP → split_by_allele   │─────┐                      │
  │     ref/alt/nosnp/ambig BAMs │     │                      │
  │     allele-specific bigWigs  │     │                      │
  └──────────────────────────────┘     │                      │
                                       ▼                      ▼
  ┌──────────────────────────────┐  ┌─────────────────────────────────┐
  │ 04_proseq  (SE reads)        │  │ 06_allele_pairs                 │
  │   trim → bowtie2 → MAPQ →    │  │   • build (motif × het SNP)     │
  │   chrM filter →              │  │     pair table                  │
  │   WASP → split_by_allele     │─▶│   • score motifs under ref/alt  │
  │     ref/alt/nosnp/ambig BAMs │  │   • per-allele ATAC counts      │
  │     allele-specific bigWigs  │  │   • per-allele PRO-seq pause    │
  └──────────────────────────────┘  │   • nearest splice junction(s)  │
                                    │   • [PLACEHOLDER] AS splicing   │
                                    └────────────────┬────────────────┘
                                                     │
                                                     ▼
                                    ┌────────────────────────────────┐
                                    │ 07_analysis                    │
                                    │   • pair-level statistics      │
                                    │   • figures                    │
                                    └────────────────────────────────┘
```

## 5. Stage summary

Each numbered directory has its own `README.md` with the full description of
inputs, outputs, parameters, and how to run it. Brief overviews:

**`01_resources/`** — Downloads and prepares mm39 genome, F121-9 het SNP VCF
(both UCSC and Ensembl-style chromosome naming), and GENCODE vM38 annotation.
Produces the canonical resource files referenced everywhere downstream.

**`02_motif_scanning/`** — Runs MOODS over mm39 with the consensus PWM library.
Produces one BED per motif archetype (e.g., `AC0007.bed`) plus a sanity-check
report identifying palindromic motifs.

**`03_atac/`** — F121-9 ATAC-seq pipeline (GSE198517, 4 replicates). Mirrors
`04_proseq/` structure. Produces consensus peaks (used by `05_motif_annot/`)
and per-replicate / per-allele BAMs (used by `06_allele_pairs/`).

**`04_proseq/`** — F121-9 PRO-seq pipeline. Produces per-replicate / per-allele
BAMs and bigWigs.

**`05_motif_annot/`** — Decomposed from the previous monolithic
`02_annotate_motifs.R` into modular per-annotation-family stages
(`02a_basic.R`, `02b_promoter_genic.R`, etc.). Each stage adds columns to a
shared per-archetype table. **No signal counts here** — only positional
features.

**`06_allele_pairs/`** — The join layer. Builds the (motif × het SNP) pair
table, rescores motifs under each allele, extracts allele-specific ATAC and
PRO-seq counts in motif windows, and computes splice junction proximity.
Designed with placeholder columns for allele-specific splicing PSI to be
populated when an RNA-seq dataset is acquired.

**`07_analysis/`** — Statistics and figures over the master pair table.

## 6. Conventions

These hold across every directory; deviations are bugs.

**Genome and naming.** mm39 / GRCm39 throughout. UCSC-style chromosome names
(`chr1`, `chrX`, ...). The Ensembl-style F121-9 VCF is converted to UCSC-style
in `01_resources/03_convert_vcf_to_ucsc.sh`; downstream pipelines use only the
UCSC-style version.

**Standard chromosomes.** chr1–19, chrX, chrY. chrM is excluded everywhere.

**Path management.** All paths in `config.sh` (or `config_wasp.sh`,
`config_atac.sh`, etc.); never hardcoded inside scripts. Each pipeline directory
sources its own `config.sh` at the top.

**SLURM submission.** Every pipeline orchestrator supports `--dry-run`
(default) and `--submit` (or `--apply`) modes. Dependencies between stages
use `sbatch --dependency=afterok:JOBID`.

**Per-sample QC, parallel-safe.** Per-sample pipelines (steps run as SLURM
arrays) write per-sample QC files (e.g., `alignment_stats_${SAMPLE}.tsv`) that
are merged in a final non-array step. This prevents race conditions on shared
output files.

**Idempotency.** Every step checks for existence of its expected output and
skips if present. Re-running a pipeline is always safe.

**Conda environment.** All scripts source `conda activate ${CONDA_ENV}` from
config. The current environment is `proseq2_env` for the read-processing
pipelines and `test_m` for the R-based annotation pipelines. (TODO: unify
these once stable.)

## 7. PRO-seq strand swap convention

This is documented in extensive prose inside `04_proseq/04_make_bigwigs.sh`,
but stated once here as the canonical reference because every downstream
PRO-seq counting step depends on it.

PRO-seq reads come from cDNA, which is the reverse complement of the original
nascent RNA. Therefore:

- Read mapped to **+** strand of the genome (FLAG 16 *not* set) ⟺
  RNA was on the **−** strand (so this read contributes to *minus-strand* signal)
- Read mapped to **−** strand of the genome (FLAG 16 set) ⟺
  RNA was on the **+** strand (so this read contributes to *plus-strand* signal)

For the **3' end of the RNA** (= RNAPII active site):

- Use the **5' end of the read** (`bedtools genomecov -5`)
- Plus-strand RNA signal ← reads with FLAG 16 set
- Minus-strand RNA signal ← reads with FLAG 16 not set

For an intragenic motif inside a gene transcribed on the **+** strand, you want
PRO-seq reads where the RNA was on the **+** strand → reads with FLAG 16 set.
For an intragenic motif on a **−** strand gene, you want reads with FLAG 16
not set.

`lib/strand_aware_count.py` is the single source of truth for this logic in
the `06_allele_pairs/` stage. Do not re-derive the strand swap inside other
scripts.

## 8. Status

| Directory             | Status                              | Notes                                |
|-----------------------|-------------------------------------|--------------------------------------|
| `01_resources/`       | DONE, frozen                        | Do not modify.                       |
| `02_motif_scanning/`  | DONE, frozen (currently running)    | Do not modify.                       |
| `lib/`                | Skeleton present                    | Will be filled in as stages need.    |
| `03_atac/`            | Planned                             | Next implementation target.          |
| `04_proseq/`          | Refactor planned                    | Existing code in `04_proseq/`.       |
| `05_motif_annot/`     | Refactor planned                    | Existing code in `03_motif_annot/`.  |
| `06_allele_pairs/`    | Planned                             | Includes RNA-seq placeholders.       |
| `07_analysis/`        | Planned                             | Final stats and plots.               |

## 9. Common operations

```bash
# Submit a pipeline (after dry-run inspection):
cd 03_atac && bash run_pipeline.sh --submit

# Dry-run only (default): print what would be submitted
cd 03_atac && bash run_pipeline.sh

# Re-run just one stage of motif annotation (e.g., re-annotate ATAC overlap
# after switching to a new ATAC peak set):
cd 05_motif_annot && bash run_pipeline.sh --stages 02e --apply

# Inspect pipeline DAG status:
squeue -u $(whoami)
```

## 10. Where things are *not*

A few things are intentionally not in this repository:

- Raw sequencing data lives in `/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/` (not version controlled).
- WASP itself is cloned by `06_wasp_setup.sh` (or `03_atac/05_wasp_setup.sh`) into the scratch tree.
- The `consensus_pwms.meme` motif library lives at `/home/alb1273/pausing_phase_project/resources/consensus_pwms.meme`.
- Final analysis figures will go to `07_analysis/output/` (still planned).
