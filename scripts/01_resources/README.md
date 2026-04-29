# `01_resources/` — Genome and VCF resource preparation

Downloads the mouse genome (mm39) and Mouse Genomes Project (MGP) VCF,
extracts per-hybrid heterozygous SNPs for each F1 of interest, and
converts chromosome naming for downstream UCSC-style consumers.

These scripts are run **once per project** (or whenever the source data
needs refreshing). Outputs land in `${RESOURCE_DIR}` defined in `config.sh`.

---

## Pipeline at a glance

```
01_download_resources.sbatch     (mm39 FASTA + MGP VCF)
            │
            ▼
02_extract_f1_het_snps.sbatch    (per-hybrid het VCFs, Ensembl naming)
            │
            ▼
03_convert_vcf_to_ucsc.sh        (parallel UCSC-naming VCFs)
            │
            ▼
04_vcf_summary_stats.sh          (sanity stats — optional)
```

---

## Quickstart

```bash
# Step 1 — download genome and raw MGP VCF
sbatch 01_download_resources.sbatch

# Step 2 — extract per-hybrid het SNPs
sbatch 02_extract_f1_het_snps.sbatch

# Step 3 — convert to UCSC chrom naming
bash 03_convert_vcf_to_ucsc.sh

# Step 4 — summary stats (optional, login-node OK)
bash 04_vcf_summary_stats.sh
```

---

## Configuration

All paths are in `config.sh`. Every script in this directory sources it:

```bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
```

To change the output root, the strain hybrids, or the MGP version, edit
`config.sh`. Don't edit individual scripts.

Key variables:

| Variable | Default | Purpose |
|---|---|---|
| `RESOURCE_DIR` | `/n/scratch/users/a/alb1273/.../resources` | Output root |
| `MGP_FTP` | EBI URL | MGP release base URL |
| `HYBRID_F121_PARENTS` | `129S1_SvImJ,CAST_EiJ` | Sample IDs in MGP VCF for F121-9 parents |
| `HYBRID_BL6CAST_PARENTS` | `C57BL_6NJ,CAST_EiJ` | Sample IDs for BL6×CAST parents |
| `THREADS` | 4 | Passed to bcftools |
| `CONDA_ENV` | `motif_scanning` | Provides bcftools, tabix, samtools, wget |

---

## Outputs

```
${RESOURCE_DIR}/
├── genome/
│   ├── mm39.fa
│   ├── mm39.fa.fai
│   └── mm39.chrom.sizes
└── vcf/
    ├── mgp_REL2021_snps.vcf.gz             # raw MGP, Ensembl naming
    ├── mgp_REL2021_snps.vcf.gz.tbi
    ├── F121-9_het_snps.vcf.gz              # per-hybrid, Ensembl
    ├── F121-9_het_snps.vcf.gz.tbi
    ├── F121-9_het_snps.ucsc.vcf.gz         # per-hybrid, UCSC
    ├── F121-9_het_snps.ucsc.vcf.gz.tbi
    ├── BL6xCAST_het_snps.vcf.gz
    ├── BL6xCAST_het_snps.vcf.gz.tbi
    ├── BL6xCAST_het_snps.ucsc.vcf.gz
    └── BL6xCAST_het_snps.ucsc.vcf.gz.tbi
```

`03_motif_annot/02_annotate_motifs.R` prefers `<HYBRID>_het_snps.ucsc.vcf.gz`
when present (since motif BEDs use UCSC naming) and falls back to the
Ensembl-naming version with a warning. **Always run step 3.**

---

## Adding a new hybrid

1. Add a new `HYBRID_FOO_PARENTS="<sample1>,<sample2>"` to `config.sh`.
   The sample IDs must match the column names in the MGP VCF — check
   with `bcftools view -h ${MGP_VCF} | tail -1`.
2. Add a corresponding `extract_one "FOO" "${HYBRID_FOO_PARENTS}"` line
   to `02_extract_f1_het_snps.sbatch`.
3. Add a corresponding `convert_one "FOO"` line to `03_convert_vcf_to_ucsc.sh`.
4. Add to `summary_one` in `04_vcf_summary_stats.sh` if you want stats.
5. Re-run steps 2–3.

---

## Notes

- `01_download_resources.sbatch` does **only** download — it doesn't
  extract or convert. That's intentional, so the slow internet step is
  separate from the fast data-munging steps.
- The chain file `mm10ToMm39.over.chain.gz` (used by `03_motif_annot/`)
  is not downloaded by these scripts; the R preprocess script in
  `03_motif_annot/` will fetch it on demand.
- `set -euo pipefail` is set in every script; conda activation is
  bracketed by `set +u; source activate; set -u` (conda activation
  scripts reference unset variables).
