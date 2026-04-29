# `lib/` — Shared Pipeline Utilities

Utilities used by **two or more** pipeline directories. They live here, rather
than being copied into each pipeline, to prevent silent drift between copies.

## Contents

| File                         | Used by                                  | Purpose                                                                 |
|------------------------------|------------------------------------------|-------------------------------------------------------------------------|
| `split_by_allele.py`         | `03_atac/`, `04_proseq/`                 | Assign reads in a BAM to ref / alt / nosnp / ambiguous based on overlap with het SNP positions. Supports both single-end (PRO-seq) and paired-end (ATAC) data. |
| `make_wasp_snps.sh`          | `03_atac/`, `04_proseq/`                 | Convert a heterozygous SNP VCF into the per-chromosome SNP-file format that WASP requires. |
| `strand_aware_count.py`      | `06_allele_pairs/`                       | Single source of truth for the PRO-seq strand-swap logic. Counts Pol II positions (3' end of RNA = 5' end of read on the opposite strand) inside given intervals. |

## Why these are shared

**`split_by_allele.py`**: the read-classification logic — "for each read, look
at every het SNP it covers; assign the read to ref / alt / nosnp / ambiguous"
— is identical for PRO-seq and ATAC at the per-base level. The *only*
difference is that ATAC is paired-end, so the script supports both modes
(autodetected, with an explicit `--paired-end` / `--single-end` override).
Putting this in one place means a bug fix or reclassification rule change
benefits both pipelines simultaneously.

**`make_wasp_snps.sh`**: both ATAC and PRO-seq pipelines need the F121-9 het
SNPs in the per-chromosome `chrN.snps.txt.gz` format that WASP's
`find_intersecting_snps.py` consumes. The conversion is a few lines but easy
to get subtly wrong (1-based vs 0-based, biallelic filtering, chromosome
naming). Centralizing avoids two slightly different versions.

**`strand_aware_count.py`**: the PRO-seq strand swap
(read on BAM "+" → RNA was on "−") is one of the most error-prone pieces of
the project. There must be exactly one place where this swap is implemented.
`04_proseq/04_make_bigwigs.sh` does it for whole-genome bigWigs (using
`samtools view -f 16` / `-F 16` flags); `06_allele_pairs/` needs it for
narrow motif-centered windows (which is awkward to do with bedtools when
intervals are tiny). The Python implementation here is what `06` uses.

## Conventions

- Each script has a `--help` and works as a standalone CLI tool.
- Scripts prefer arguments over environment variables, but accept paths from
  the calling pipeline's `config.sh` via the caller passing them on the
  command line.
- All scripts produce informative stderr/stdout that tells you what was done
  and what passed sanity checks.
- Every script is idempotent where reasonable (skip if output exists).

## Adding to lib/

A utility belongs here if and only if **at least two** pipeline directories
need to call it with substantively the same logic. Per-pipeline helpers stay
inside their pipeline.
