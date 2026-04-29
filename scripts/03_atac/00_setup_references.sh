#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 00_setup_references.sh — mm39 bowtie2 index for ATAC alignment
# =============================================================================
#
# WHAT THIS DOES:
#   Ensures the mm39 bowtie2 index, FASTA, and chromosome sizes file exist.
#
# DESIGN CHOICE — REUSE EXISTING INDEX:
#   If 04_proseq/00_setup_references.sh has already run, the mm39 index lives
#   at ${MM39_BT2_IDX} and is fully usable here. We check for it and skip
#   rebuilding if present (~30-45 min saved per project setup).
#
# WHAT WE DO NOT BUILD HERE:
#   - dm6 spike-in index (not used in ATAC)
#   - rRNA index (not used in ATAC)
#   These are 04_proseq-specific.
#
# SBATCH RESOURCES (only if building from scratch):
#   sbatch --partition=short --time=2:00:00 --mem=16G --cpus-per-task=8 \
#       -o logs/00_%j.out -e logs/00_%j.err 00_setup_references.sh
# =============================================================================

source "config.sh"

step_header "ATAC Pipeline Step 00: Reference Setup"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] bowtie2 version: $(bowtie2 --version 2>&1 | head -1)"
echo "[INFO] samtools version: $(samtools --version | head -1)"

create_dirs
mkdir -p "${GENOME_DIR}"

# -----------------------------------------------------------------------------
# 1. mm39 bowtie2 index
# -----------------------------------------------------------------------------
if [[ -f "${MM39_BT2_IDX}.1.bt2" ]]; then
    echo ""
    echo "[INFO] mm39 bowtie2 index already present:"
    ls -lh "${MM39_BT2_IDX}".*.bt2 2>/dev/null | sed 's/^/    /'
    echo "[INFO] Reusing it (likely from 04_proseq/00_setup_references.sh)."
    echo "[INFO] To rebuild, delete the .bt2 files and re-run."
else
    echo ""
    echo "[INFO] mm39 bowtie2 index not found. Will build from scratch."

    # FASTA
    if [[ ! -f "${MM39_FA}" ]]; then
        echo "[INFO] Downloading mm39.fa from UCSC..."
        wget -q --show-progress -O "${MM39_FA}.gz" \
            "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz"
        gunzip "${MM39_FA}.gz"
    fi

    # FASTA index
    if [[ ! -f "${MM39_FA}.fai" ]]; then
        echo "[INFO] samtools faidx..."
        samtools faidx "${MM39_FA}"
    fi

    # Bowtie2 index
    echo "[INFO] Building bowtie2 index (this takes 30-45 min)..."
    bowtie2-build --threads "${THREADS}" "${MM39_FA}" "${MM39_BT2_IDX}"
fi

# -----------------------------------------------------------------------------
# 2. Chromosome sizes (used by bedGraphToBigWig)
# -----------------------------------------------------------------------------
if [[ ! -f "${MM39_CHROM_SIZES}" ]]; then
    echo ""
    echo "[INFO] Generating chromosome sizes..."
    if [[ ! -f "${MM39_FA}.fai" ]]; then
        samtools faidx "${MM39_FA}"
    fi
    cut -f1,2 "${MM39_FA}.fai" > "${MM39_CHROM_SIZES}"
fi

# -----------------------------------------------------------------------------
# 3. Sanity check
# -----------------------------------------------------------------------------
echo ""
echo "[SANITY CHECK]"
echo "  mm39 bowtie2 index files:"
ls -lh "${MM39_BT2_IDX}".*.bt2 | sed 's/^/    /'
echo "  Chromosome sizes file: ${MM39_CHROM_SIZES}"
head -5 "${MM39_CHROM_SIZES}" | sed 's/^/    /'

# -----------------------------------------------------------------------------
# 4. Verify ATAC-specific tools are present
# -----------------------------------------------------------------------------
echo ""
echo "[INFO] Verifying ATAC pipeline tools..."
for tool in bowtie2 samtools bedtools bedGraphToBigWig fasterq-dump fastqc cutadapt pigz macs3; do
    if command -v "${tool}" &>/dev/null; then
        echo "  [OK]    ${tool}: $(command -v ${tool})"
    else
        echo "  [MISSING] ${tool}"
        if [[ "${tool}" == "macs3" ]]; then
            echo "          Install with: pip install macs3"
        fi
    fi
done

step_header "Step 00 COMPLETE"
echo "Next: bash 01_download_fastq.sh   (or sbatch as array of 4)"
