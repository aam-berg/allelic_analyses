#!/usr/bin/env bash
#set -euo pipefail
# =============================================================================
# 00_setup_references.sh — Download reference genomes and build bowtie2 indices
# =============================================================================
#
# WHAT THIS DOES:
#   1. Downloads mm39 (mouse) genome FASTA from UCSC
#   2. Downloads dm6 (Drosophila) genome FASTA from UCSC — needed for spike-in filtering
#   3. Builds bowtie2 indices for both genomes
#   4. Generates mm39 chromosome sizes file (needed for bedGraphToBigWig)
#
# WHY:
#   - We need mm39 for the primary alignment (user's requested genome build)
#   - We need dm6 because the original PRO-seq experiment included Drosophila
#     spike-in cells for normalization. We align to dm6 first and DISCARD those
#     reads so they don't contaminate our mouse signal.
#   - bowtie2 is our aligner of choice for PRO-seq: it handles the short,
#     single-end reads well and is the standard in the field.
#
# RUNTIME: ~45 min for mm39 index, ~15 min for dm6 index (8 threads)
# MEMORY:  ~8 GB peak for mm39 index building
# sbatch --partition=short --time=12:00:00 --mem=64G --cpus-per-task=8 ./00_setup_references.sh
# USAGE:
#   bash 00_setup_references.sh
# =============================================================================

source "config.sh"

echo "============================================================"
echo "PRO-seq Pipeline: Reference Genome Setup"
echo "Started: $(date)"
echo "============================================================"

# --- Activate conda environment ---
echo ""
source activate ${CONDA_ENV_NAME}
echo "[INFO] Using bowtie2 version: $(bowtie2 --version | head -1)"
echo "[INFO] Using samtools version: $(samtools --version | head -1)"

# --- Create directories ---
echo ""
echo "[INFO] Creating directory structure..."
mkdir -p "${MM39_DIR}" "${DM6_DIR}"
echo "  mm39 dir: ${MM39_DIR}"
echo "  dm6 dir:  ${DM6_DIR}"

# =============================================================================
# 1. Download mm39 genome
# =============================================================================
echo ""
echo "============================================================"
echo "STEP 1: Download mm39 (mouse) genome"
echo "============================================================"

MM39_FA="${MM39_DIR}/mm39.fa"
if [[ -f "${MM39_FA}" ]]; then
    echo "[INFO] mm39 FASTA already exists: ${MM39_FA}"
    echo "[INFO] Skipping download. Delete file to re-download."
else
    echo "[INFO] Downloading mm39 from UCSC..."
    wget -q --show-progress -O "${MM39_DIR}/mm39.fa.gz" \
        "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz"
    echo "[INFO] Decompressing..."
    gunzip "${MM39_DIR}/mm39.fa.gz"
    echo "[INFO] mm39 FASTA ready: ${MM39_FA}"
fi

# Sanity check
echo "[SANITY CHECK] mm39 FASTA:"
echo "  File size: $(du -h "${MM39_FA}" | cut -f1)"
echo "  Number of chromosomes: $(grep -c '^>' "${MM39_FA}")"
echo "  First 5 chromosome names:"
grep '^>' "${MM39_FA}" | head -5 | sed 's/^/    /'

# =============================================================================
# 2. Generate mm39 chromosome sizes
# =============================================================================
echo ""
echo "============================================================"
echo "STEP 2: Generate mm39 chromosome sizes"
echo "============================================================"

if [[ -f "${MM39_CHROM_SIZES}" ]]; then
    echo "[INFO] Chrom sizes file already exists: ${MM39_CHROM_SIZES}"
else
    echo "[INFO] Generating chromosome sizes with samtools faidx..."
    samtools faidx "${MM39_FA}"
    # Extract chr name and length (columns 1 and 2 of .fai)
    cut -f1,2 "${MM39_FA}.fai" > "${MM39_CHROM_SIZES}"
    echo "[INFO] Chrom sizes saved: ${MM39_CHROM_SIZES}"
fi

echo "[SANITY CHECK] First 10 entries of chrom sizes:"
head -10 "${MM39_CHROM_SIZES}" | sed 's/^/    /'

# =============================================================================
# 3. Build mm39 bowtie2 index
# =============================================================================
echo ""
echo "============================================================"
echo "STEP 3: Build mm39 bowtie2 index"
echo "============================================================"

if [[ -f "${MM39_BT2_IDX}.1.bt2" ]]; then
    echo "[INFO] mm39 bowtie2 index already exists: ${MM39_BT2_IDX}.*"
    echo "[INFO] Skipping. Delete index files to rebuild."
else
    echo "[INFO] Building bowtie2 index for mm39 (this takes ~30-45 min)..."
    echo "[INFO] Using ${THREADS} threads"
    bowtie2-build --threads "${THREADS}" "${MM39_FA}" "${MM39_BT2_IDX}"
    echo "[INFO] mm39 bowtie2 index complete."
fi

echo "[SANITY CHECK] mm39 index files:"
ls -lh "${MM39_BT2_IDX}"*.bt2 2>/dev/null | sed 's/^/    /' || echo "    ERROR: index files not found!"

# =============================================================================
# 4. Download dm6 genome
# =============================================================================
echo ""
echo "============================================================"
echo "STEP 4: Download dm6 (Drosophila) spike-in genome"
echo "============================================================"

DM6_FA="${DM6_DIR}/dm6.fa"
if [[ -f "${DM6_FA}" ]]; then
    echo "[INFO] dm6 FASTA already exists: ${DM6_FA}"
else
    echo "[INFO] Downloading dm6 from UCSC..."
    wget -q --show-progress -O "${DM6_DIR}/dm6.fa.gz" \
        "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz"
    echo "[INFO] Decompressing..."
    gunzip "${DM6_DIR}/dm6.fa.gz"
    echo "[INFO] dm6 FASTA ready: ${DM6_FA}"
fi

echo "[SANITY CHECK] dm6 FASTA:"
echo "  File size: $(du -h "${DM6_FA}" | cut -f1)"
echo "  Number of chromosomes: $(grep -c '^>' "${DM6_FA}")"

# =============================================================================
# 5. Build dm6 bowtie2 index
# =============================================================================
echo ""
echo "============================================================"
echo "STEP 5: Build dm6 bowtie2 index"
echo "============================================================"

if [[ -f "${DM6_BT2_IDX}.1.bt2" ]]; then
    echo "[INFO] dm6 bowtie2 index already exists: ${DM6_BT2_IDX}.*"
else
    echo "[INFO] Building bowtie2 index for dm6 (this takes ~10-15 min)..."
    bowtie2-build --threads "${THREADS}" "${DM6_FA}" "${DM6_BT2_IDX}"
    echo "[INFO] dm6 bowtie2 index complete."
fi

echo "[SANITY CHECK] dm6 index files:"
ls -lh "${DM6_BT2_IDX}"*.bt2 2>/dev/null | sed 's/^/    /' || echo "    ERROR: index files not found!"

# =============================================================================
# Done
# =============================================================================
echo ""
echo "============================================================"
echo "Reference genome setup COMPLETE"
echo "Finished: $(date)"
echo "============================================================"
echo ""
echo "Summary:"
echo "  mm39 FASTA:       ${MM39_FA}"
echo "  mm39 bt2 index:   ${MM39_BT2_IDX}"
echo "  mm39 chrom sizes: ${MM39_CHROM_SIZES}"
echo "  dm6 FASTA:        ${DM6_FA}"
echo "  dm6 bt2 index:    ${DM6_BT2_IDX}"
echo ""
echo "Next step: Run 01_download_and_process.sh"