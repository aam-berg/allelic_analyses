#!/bin/bash
# =============================================================================
# 01_resources/config.sh
#
# Shared configuration for the 01_resources/ pipeline. Every script in this
# directory sources this file at the top:
#
#   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#   source "${SCRIPT_DIR}/config.sh"
#
# Edit paths below for your environment.
# =============================================================================

# ---- Output root ----
RESOURCE_DIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources"

# ---- Derived directories ----
GENOME_DIR="${RESOURCE_DIR}/genome"
VCF_DIR="${RESOURCE_DIR}/vcf"
ANNOT_DIR="${RESOURCE_DIR}/annotation"

# ---- Genome ----
GENOME_FA_UCSC="${GENOME_DIR}/mm39.fa"               # chr1, chr2, ...
GENOME_FA_ENSEMBL="${GENOME_DIR}/mm39.ensembl.fa"    # 1, 2, ...
CHROM_SIZES_UCSC="${GENOME_DIR}/mm39.chrom.sizes"
CHROM_SIZES_ENSEMBL="${GENOME_DIR}/mm39.chrom.sizes.ensembl"

# ---- mm10 → mm39 liftOver chain (used by 03_motif_annot) ----
CHAIN_FILE="${GENOME_DIR}/mm10ToMm39.over.chain"
CHAIN_FILE_GZ="${CHAIN_FILE}.gz"

# ---- Mouse Genomes Project ----
MGP_FTP="https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels"
MGP_VCF="${VCF_DIR}/mgp_REL2021_snps.vcf.gz"

# ---- Strain hybrid definitions (parents for the F1 het-SNP extraction) ----
HYBRID_F121_PARENTS="129S1_SvImJ,CAST_EiJ"
HYBRID_BL6CAST_PARENTS="C57BL_6NJ,CAST_EiJ"

# ---- Compute settings ----
THREADS=4

# ---- Environment ----
CONDA_ENV="motif_scanning"
