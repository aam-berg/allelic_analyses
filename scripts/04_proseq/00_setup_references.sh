#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04_proseq/00_setup_references.sh — Reference setup for PRO-seq
# =============================================================================
#
# WHAT THIS DOES:
#   1. mm39 bowtie2 index (REUSE if present, else download + build)
#   2. dm6 bowtie2 index (REUSE if present; only build if DM6_SPIKEIN_FILTER=true)
#   3. mm39 chromosome sizes file
#   4. mouse rRNA reference + bowtie2 index (PRO-seq specific; built fresh)
#
# DESIGN — REUSE FIRST:
#   The mm39 and dm6 indices are typically built once during the project's
#   first sequencing pipeline run (in either 03_atac/ or this directory). We
#   check for their presence and skip rebuilding (~30-45 min and ~10-15 min
#   saved respectively). Re-runs are cheap; only the rRNA index is unique
#   to this pipeline.
#
# WHAT'S DIFFERENT FROM THE OLDER VERSION OF THIS SCRIPT:
#   The previous 04_proseq had two competing setup scripts (one for
#   mm39+dm6, one for rRNA). They've been merged here into a single
#   idempotent script.
#
# SBATCH RESOURCES:
#   sbatch --partition=short --time=2:00:00 --mem=16G --cpus-per-task=8 \
#       -o logs/00_%j.out -e logs/00_%j.err 00_setup_references.sh
# =============================================================================

source "config.sh"

step_header "PRO-seq Pipeline Step 00: Reference Setup"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] bowtie2:  $(bowtie2 --version 2>&1 | head -1)"
echo "[INFO] samtools: $(samtools --version | head -1)"

create_dirs
mkdir -p "${GENOME_DIR}"

# -----------------------------------------------------------------------------
# 1. mm39 bowtie2 index
# -----------------------------------------------------------------------------
if [[ -f "${MM39_BT2_IDX}.1.bt2" ]]; then
    echo ""
    echo "[INFO] mm39 bowtie2 index already present:"
    ls -lh "${MM39_BT2_IDX}".*.bt2 2>/dev/null | sed 's/^/    /' | head -3
    echo "[INFO] Reusing it (likely from 03_atac/ or earlier 04_proseq/ run)."
else
    echo ""
    echo "[INFO] mm39 bowtie2 index not found. Building..."

    if [[ ! -f "${MM39_FA}" ]]; then
        echo "[INFO] Downloading mm39.fa from UCSC..."
        wget -q --show-progress -O "${MM39_FA}.gz" \
            "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz"
        gunzip "${MM39_FA}.gz"
    fi
    if [[ ! -f "${MM39_FA}.fai" ]]; then
        samtools faidx "${MM39_FA}"
    fi

    echo "[INFO] Building bowtie2 index (30-45 min)..."
    bowtie2-build --threads "${THREADS}" "${MM39_FA}" "${MM39_BT2_IDX}"
fi

# -----------------------------------------------------------------------------
# 2. mm39 chromosome sizes
# -----------------------------------------------------------------------------
if [[ ! -f "${MM39_CHROM_SIZES}" ]]; then
    echo ""
    echo "[INFO] Generating mm39 chromosome sizes..."
    [[ -f "${MM39_FA}.fai" ]] || samtools faidx "${MM39_FA}"
    cut -f1,2 "${MM39_FA}.fai" > "${MM39_CHROM_SIZES}"
fi

# -----------------------------------------------------------------------------
# 3. dm6 bowtie2 index (only if spike-in filter is enabled)
# -----------------------------------------------------------------------------
if [[ "${DM6_SPIKEIN_FILTER}" == "true" ]]; then
    if [[ -f "${DM6_BT2_IDX}.1.bt2" ]]; then
        echo ""
        echo "[INFO] dm6 bowtie2 index already present. Reusing."
    else
        echo ""
        echo "[INFO] dm6 bowtie2 index not found. Building..."
        if [[ ! -f "${DM6_FA}" ]]; then
            wget -q --show-progress -O "${DM6_FA}.gz" \
                "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz"
            gunzip "${DM6_FA}.gz"
        fi
        bowtie2-build --threads "${THREADS}" "${DM6_FA}" "${DM6_BT2_IDX}"
    fi
else
    echo ""
    echo "[INFO] DM6_SPIKEIN_FILTER=false; skipping dm6 index setup."
fi

# -----------------------------------------------------------------------------
# 4. Mouse rRNA reference and bowtie2 index
# -----------------------------------------------------------------------------
echo ""
if [[ -f "${RRNA_FA}" ]]; then
    echo "[INFO] rRNA FASTA already present: ${RRNA_FA}"
else
    echo "[INFO] Downloading mouse rRNA sequences from NCBI..."
    # Mouse rRNA accessions:
    #   NR_003279.1 — 28S ribosomal RNA (Rn28s1)
    #   NR_003278.3 — 18S ribosomal RNA (Rn18s)
    #   NR_003280.2 — 5.8S ribosomal RNA (Rn5-8s1)
    #   NR_030686.1 — 5S ribosomal RNA (Rn5s)
    #   BK000964.3  — 45S pre-ribosomal RNA (full ~13kb rDNA repeat unit, including
    #                                         ITS/ETS spacers transcribed by Pol I)
    ACCESSIONS=("NR_003279.1" "NR_003278.3" "NR_003280.2" "NR_030686.1" "BK000964.3")

    : > "${RRNA_FA}"
    for ACC in "${ACCESSIONS[@]}"; do
        echo "  Fetching ${ACC}..."
        efetch -db nucleotide -id "${ACC}" -format fasta >> "${RRNA_FA}"
    done
fi

echo "[SANITY CHECK] rRNA FASTA:"
echo "  Size:        $(du -h "${RRNA_FA}" | cut -f1)"
echo "  Sequences:   $(grep -c '^>' "${RRNA_FA}")"
echo "  Names:"
grep '^>' "${RRNA_FA}" | head -10 | sed 's/^/    /'

if [[ -f "${RRNA_BT2_IDX}.1.bt2" ]]; then
    echo ""
    echo "[INFO] rRNA bowtie2 index already present."
else
    echo ""
    echo "[INFO] Building rRNA bowtie2 index..."
    bowtie2-build "${RRNA_FA}" "${RRNA_BT2_IDX}"
fi

# -----------------------------------------------------------------------------
# 5. Verify all required references exist
# -----------------------------------------------------------------------------
echo ""
echo "[PRE-FLIGHT CHECK]"
ALL_OK=1
declare -a CHECKS=("${MM39_BT2_IDX}.1.bt2" "${MM39_CHROM_SIZES}" "${RRNA_BT2_IDX}.1.bt2")
[[ "${DM6_SPIKEIN_FILTER}" == "true" ]] && CHECKS+=("${DM6_BT2_IDX}.1.bt2")
for F in "${CHECKS[@]}"; do
    if [[ -f "${F}" ]]; then
        echo "  [OK]    ${F}"
    else
        echo "  [MISSING] ${F}"
        ALL_OK=0
    fi
done

# Verify essential tools
echo ""
echo "[INFO] Verifying tools..."
for tool in bowtie2 samtools bedtools bedGraphToBigWig bigWigMerge fasterq-dump fastqc cutadapt pigz efetch; do
    if command -v "${tool}" &>/dev/null; then
        echo "  [OK]    ${tool}"
    else
        echo "  [MISSING] ${tool}"
        ALL_OK=0
    fi
done

(( ALL_OK )) || { echo ""; echo "[ERROR] Pre-flight check failed."; exit 1; }

step_header "Step 00 COMPLETE"
echo "Next: bash 01_download_fastq.sh   (or sbatch as array of ${#SAMPLE_ORDER[@]})"
