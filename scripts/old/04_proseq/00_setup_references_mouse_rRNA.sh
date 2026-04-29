#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 00_setup_references.sh — Build rRNA bowtie2 index for PRO-seq filtering
# =============================================================================
#
# WHAT THIS DOES:
#   1. Downloads mouse rRNA sequences (28S, 18S, 5.8S, 5S, 45S pre-rRNA)
#   2. Builds a bowtie2 index for rRNA filtering
#
# WHY:
#   rDNA loci are among the most heavily transcribed regions in the genome.
#   In PRO-seq, a substantial fraction of reads (sometimes 5-20%) come from
#   active Pol I transcription of rDNA. These reads:
#     - Consume sequencing depth without contributing to Pol II analysis
#     - Some map uniquely and pass MAPQ filters (so MAPQ alone doesn't remove them)
#     - Can skew normalization
#   Filtering against an rRNA reference before genome alignment is standard in
#   dedicated PRO-seq pipelines (e.g., proseq2.0 from the Danko lab).
#
# NOTE:
#   This script assumes you already have mm39 and dm6 bowtie2 indices
#   (from your existing resources directory). It only builds the rRNA index.
#
# SBATCH RESOURCES:
#   sbatch --partition=short --time=10:00:00 --mem=16G --cpus-per-task=4 00_setup_references_mouse_rRNA.sh
# =============================================================================

source "config.sh"

step_header "PRO-seq Pipeline Step 00: Build rRNA Reference Index"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] bowtie2 version: $(bowtie2 --version 2>&1 | head -1)"

create_dirs

# =============================================================================
# 1. Create a mouse rRNA FASTA
# =============================================================================
RRNA_FA="${RRNA_DIR}/mouse_rRNA.fa"

if [[ -f "${RRNA_FA}" ]]; then
    echo "[INFO] rRNA FASTA already exists: ${RRNA_FA}"
else
    echo "[INFO] Downloading mouse rRNA sequences from NCBI..."

    # Mouse rRNA accessions:
    #   NR_003279.1 — Mus musculus 28S ribosomal RNA (Rn28s1)
    #   NR_003278.3 — Mus musculus 18S ribosomal RNA (Rn18s)
    #   NR_003280.2 — Mus musculus 5.8S ribosomal RNA (Rn5-8s1)
    #   NR_030686.1 — Mus musculus 5S ribosomal RNA (Rn5s)
    #   BK000964.3  — Mus musculus 45S pre-ribosomal RNA (complete rDNA repeat unit)
    #
    # The 45S pre-rRNA (BK000964.3) is the full ~13kb rDNA repeat unit
    # transcribed by RNA Pol I. It contains 18S, 5.8S, and 28S plus spacers.
    # Including this ensures we catch reads from internal/external transcribed
    # spacers (ITS/ETS) that are part of the nascent Pol I transcript.

    ACCESSIONS=("NR_003279.1" "NR_003278.3" "NR_003280.2" "NR_030686.1" "BK000964.3")

    > "${RRNA_FA}"  # empty the file

    for ACC in "${ACCESSIONS[@]}"; do
        echo "  Fetching ${ACC}..."
        # Use efetch from entrez-direct (already in your conda env)
        efetch -db nucleotide -id "${ACC}" -format fasta >> "${RRNA_FA}"
    done

    echo "[INFO] rRNA FASTA created: ${RRNA_FA}"
fi

# Sanity check
echo "[SANITY CHECK] rRNA FASTA:"
echo "  File size: $(du -h "${RRNA_FA}" | cut -f1)"
echo "  Number of sequences: $(grep -c '^>' "${RRNA_FA}")"
echo "  Sequence names:"
grep '^>' "${RRNA_FA}" | sed 's/^/    /'

# =============================================================================
# 2. Build bowtie2 index for rRNA
# =============================================================================
if [[ -f "${RRNA_BT2_IDX}.1.bt2" ]]; then
    echo "[INFO] rRNA bowtie2 index already exists: ${RRNA_BT2_IDX}.*"
else
    echo "[INFO] Building bowtie2 index for rRNA..."
    bowtie2-build "${RRNA_FA}" "${RRNA_BT2_IDX}"
    echo "[INFO] rRNA bowtie2 index complete."
fi

echo "[SANITY CHECK] rRNA index files:"
ls -lh "${RRNA_BT2_IDX}"*.bt2 2>/dev/null | sed 's/^/    /' || echo "    ERROR: index files not found!"

# =============================================================================
# 3. Verify all required references exist
# =============================================================================
step_header "Verifying all reference indices"

ALL_OK=1
for F in "${MM39_BT2_IDX}.1.bt2" "${DM6_BT2_IDX}.1.bt2" "${RRNA_BT2_IDX}.1.bt2" "${MM39_CHROM_SIZES}"; do
    if [[ -f "${F}" ]]; then
        echo "  [OK]    ${F}"
    else
        echo "  [MISSING] ${F}"
        ALL_OK=0
    fi
done

if (( ALL_OK )); then
    echo ""
    echo "[INFO] All references are ready."
else
    echo ""
    echo "[WARNING] Some references are missing. Check paths in config.sh."
fi

step_header "Step 00 COMPLETE"
echo "rRNA index: ${RRNA_BT2_IDX}"
echo ""
echo "Next: Run 01_download_fastq.sh"