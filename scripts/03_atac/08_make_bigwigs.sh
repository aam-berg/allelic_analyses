#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 08_make_bigwigs.sh — Per-sample fragment-coverage bigWigs
# =============================================================================
#
# WHAT THIS DOES (per sample):
#   1. Per-replicate bigWig from the final (pre-WASP) BAM
#   2. Allele-specific bigWigs from the {SAMPLE}_ref.bam and {SAMPLE}_alt.bam
#
# Cross-sample merging (combined and per-allele) is done in step 09, which
# runs as a single non-array job.
#
# WHY FRAGMENT COVERAGE (-bg -pc) AND NOT 5'-END COVERAGE:
#   For ATAC visualization and motif-window quantification, the fragment span
#   carries the relevant biological signal: an open region produces fragments
#   that span it, regardless of which end of the fragment is closer to the
#   actual Tn5 cut. This matches the ENCODE ATAC convention.
#
#   For single-base Tn5 cut-site precision (footprinting), apply the +4/-5
#   shift and use -5 to extract cut sites. Not done here; can be added as a
#   parallel bigWig type later if 07_analysis/ needs it.
#
# WHY NO STRAND SPLIT:
#   Unlike PRO-seq, ATAC reads on + and - strands at the same accessible
#   region carry the same biological signal (both ends of the same fragment).
#   We produce a single combined-strand bigWig.
#
# IDEMPOTENCY:
#   Skip generation if the output bigWig already exists.
#
# PARALLEL EXECUTION:
#   sbatch --array=0-3 --partition=short --time=2:00:00 --mem=16G \
#       --cpus-per-task=4 -o logs/08_%A_%a.out -e logs/08_%A_%a.err \
#       08_make_bigwigs.sh
# =============================================================================

source "config.sh"

step_header "ATAC Pipeline Step 08: Per-sample BigWigs"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] bedtools:         $(bedtools --version)"
echo "[INFO] bedGraphToBigWig: $(which bedGraphToBigWig)"

create_dirs
resolve_samples "${1:-}"

# Verify chrom sizes
if [[ ! -f "${MM39_CHROM_SIZES}" ]]; then
    echo "[ERROR] mm39 chromosome sizes file missing: ${MM39_CHROM_SIZES}" >&2
    echo "[ERROR] Run 00_setup_references.sh first." >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Helper: BAM -> fragment-coverage bigWig
# -----------------------------------------------------------------------------
make_atac_bigwig() {
    local BAM="$1"
    local OUTBW="$2"
    local LABEL="$3"

    if [[ -f "${OUTBW}" ]]; then
        echo "  [SKIP] ${LABEL}: bigWig already exists."
        return 0
    fi
    if [[ ! -f "${BAM}" ]]; then
        echo "  [SKIP] ${LABEL}: input BAM missing (${BAM})"
        return 0
    fi

    local TMPBG="${OUTBW%.bw}.bedGraph"

    # bedtools genomecov:
    #   -ibam  input BAM
    #   -bg    bedGraph output (vs -bga which would emit zero-coverage rows)
    #   -pc    PE coverage from inserts (R1 leftmost to R2 rightmost), not
    #          per-read coverage. The "-pc" flag is what makes this fragment-
    #          level rather than read-level.
    # We then sort and gate to chromosomes present in chrom_sizes (drops
    # any contigs that bedtools might emit but bedGraphToBigWig doesn't know).
    bedtools genomecov -ibam "${BAM}" -bg -pc \
    | sort -k1,1 -k2,2n \
    | awk 'NR==FNR {chr[$1]=1; next} ($1 in chr)' "${MM39_CHROM_SIZES}" - \
    > "${TMPBG}"

    if [[ ! -s "${TMPBG}" ]]; then
        echo "  [WARN] ${LABEL}: bedGraph is empty; skipping bigWig conversion."
        rm -f "${TMPBG}"
        return 0
    fi

    bedGraphToBigWig "${TMPBG}" "${MM39_CHROM_SIZES}" "${OUTBW}"
    rm -f "${TMPBG}"

    local SIZE_MB
    SIZE_MB=$(awk "BEGIN {printf \"%.1f\", $(stat -c%s "${OUTBW}")/1048576}")
    echo "  [OK]   ${LABEL}: ${OUTBW} (${SIZE_MB} MB)"
}

# -----------------------------------------------------------------------------
# Loop over samples
# -----------------------------------------------------------------------------
for SAMPLE in "${RUN_SAMPLES[@]}"; do
    step_header "BigWigs for ${SAMPLE}"

    # ---- Per-replicate bigWig (from pre-WASP final BAM) ----
    FINAL_BAM="${BAM_FINAL_DIR}/${SAMPLE}_final.bam"
    PER_REP_BW="${BIGWIG_IND_DIR}/${SAMPLE}.bw"

    if [[ ! -f "${FINAL_BAM}" ]]; then
        echo "[ERROR] Final BAM missing for ${SAMPLE}: ${FINAL_BAM}" >&2
        exit 1
    fi
    make_atac_bigwig "${FINAL_BAM}" "${PER_REP_BW}" "${SAMPLE} (per-rep)"

    # ---- Allele-specific bigWigs (from WASP-filtered allele-split BAMs) ----
    for ALLELE in ref alt; do
        BAM="${ALLELE_DIR}/${SAMPLE}_${ALLELE}.bam"
        OUTBW="${BIGWIG_ALLELE_DIR}/${SAMPLE}_${ALLELE}.bw"
        make_atac_bigwig "${BAM}" "${OUTBW}" "${SAMPLE} ${ALLELE}"
    done
done

step_header "Step 08 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output:"
echo "  Per-replicate bigWigs:    ${BIGWIG_IND_DIR}/"
echo "  Allele-specific bigWigs:  ${BIGWIG_ALLELE_DIR}/"
echo ""
echo "Next: bash 09_merge_qc.sh   (single non-array job; produces merged bigWigs + final QC report)"
