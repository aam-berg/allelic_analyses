#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04_proseq/04_make_bigwigs.sh — Per-rep strand-specific bigWigs
# =============================================================================
#
# WHAT THIS DOES (per sample):
#   Generates 4 single-nucleotide-resolution bigWigs:
#     {SAMPLE}_3prime_plus.bw   - 3' end of nascent RNA, plus strand
#     {SAMPLE}_3prime_minus.bw  - 3' end of nascent RNA, minus strand
#     {SAMPLE}_5prime_plus.bw   - 5' end of nascent RNA, plus strand
#     {SAMPLE}_5prime_minus.bw  - 5' end of nascent RNA, minus strand
#
#   These are RAW COUNTS (no normalization). Cross-rep merging happens in
#   step 09; allele-specific bigWigs come from step 08.
#
# THE PRO-seq STRAND SWAP (canonical, see scripts/README.md section 7):
#   PRO-seq reads = cDNA = REVERSE COMPLEMENT of RNA. So:
#     Read on BAM "+" strand (FLAG 16 NOT set) <=> RNA on "-" strand
#     Read on BAM "-" strand (FLAG 16 set)     <=> RNA on "+" strand
#
#   Position assignment:
#     5' of read = 3' of RNA = Pol II active site  -> bedtools -5
#     3' of read = 5' of RNA                       -> bedtools -3
#
#   Therefore for the 3' end of RNA (Pol II position):
#     plus_3prime.bw    : reads with FLAG 16 set     (`-f 16`), bedtools -5
#     minus_3prime.bw   : reads without FLAG 16     (`-F 16`), bedtools -5
#   And for the 5' end of RNA:
#     plus_5prime.bw    : reads with FLAG 16 set     (`-f 16`), bedtools -3
#     minus_5prime.bw   : reads without FLAG 16     (`-F 16`), bedtools -3
#
# CHANGES FROM PREVIOUS VERSION:
#   - All DEDUP_5PRIME flag handling removed (dedup itself is gone in step 03).
#   - No cross-rep merging here. Merging moved to step 09 (single non-array
#     job that runs after both step 04 and step 08 complete).
#
# PARALLEL EXECUTION:
#   sbatch --array=0-3 --partition=short --time=2:00:00 --mem=8G \
#       --cpus-per-task=4 -o logs/04_%A_%a.out -e logs/04_%A_%a.err \
#       04_make_bigwigs.sh
# =============================================================================

source "config.sh"

step_header "PRO-seq Pipeline Step 04: Per-rep BigWigs"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] bedtools:         $(bedtools --version)"
echo "[INFO] bedGraphToBigWig: $(which bedGraphToBigWig)"

create_dirs
resolve_samples "${1:-}"

[[ -f "${MM39_CHROM_SIZES}" ]] || \
    { echo "[ERROR] mm39 chrom sizes missing: ${MM39_CHROM_SIZES}" >&2; exit 1; }

# -----------------------------------------------------------------------------
# Helper: BAM -> single-nt bigWig with strand and end filtering
# -----------------------------------------------------------------------------
#   Args:
#     $1 = BAM
#     $2 = samtools strand flag ("-f 16" or "-F 16")
#     $3 = bedtools end flag ("-5" or "-3")
#     $4 = output bigWig
#     $5 = label for logging
make_bigwig() {
    local BAM="$1"
    local SAM_FLAG="$2"
    local BT_END_FLAG="$3"
    local OUTBW="$4"
    local LABEL="$5"

    if [[ -f "${OUTBW}" ]]; then
        echo "  [SKIP] ${LABEL}: bigWig already exists."
        return 0
    fi

    local TMPBG="${OUTBW%.bw}.bedGraph"

    samtools view -b ${SAM_FLAG} "${BAM}" \
    | bedtools genomecov -ibam stdin -bg ${BT_END_FLAG} \
    | sort -k1,1 -k2,2n \
    | awk 'NR==FNR {chr[$1]=1; next} ($1 in chr)' "${MM39_CHROM_SIZES}" - \
    > "${TMPBG}"

    if [[ ! -s "${TMPBG}" ]]; then
        echo "  [WARN] ${LABEL}: bedGraph empty; skipping bigWig."
        rm -f "${TMPBG}"
        return 0
    fi

    local N_POS TOTAL_COUNTS
    N_POS=$(wc -l < "${TMPBG}")
    TOTAL_COUNTS=$(awk '{s+=$4} END {print s+0}' "${TMPBG}")

    bedGraphToBigWig "${TMPBG}" "${MM39_CHROM_SIZES}" "${OUTBW}"
    rm -f "${TMPBG}"
    echo "  [OK]   ${LABEL}: ${N_POS} positions, ${TOTAL_COUNTS} total counts -> ${OUTBW}"
}

# -----------------------------------------------------------------------------
# Loop over samples
# -----------------------------------------------------------------------------
for SAMPLE in "${RUN_SAMPLES[@]}"; do
    step_header "BigWigs for ${SAMPLE}"

    BAM="${BAM_MM39_DIR}/${SAMPLE}_final.bam"
    [[ -f "${BAM}" ]] || { echo "[ERROR] Final BAM missing: ${BAM}" >&2; exit 1; }

    # 3' end of RNA (Pol II active site)
    echo "[INFO] 3' end of RNA (Pol II position):"
    make_bigwig "${BAM}" "-f 16" "-5" \
        "${BIGWIG_IND_DIR}/${SAMPLE}_3prime_plus.bw"  "3prime_plus"
    make_bigwig "${BAM}" "-F 16" "-5" \
        "${BIGWIG_IND_DIR}/${SAMPLE}_3prime_minus.bw" "3prime_minus"

    # 5' end of RNA
    echo "[INFO] 5' end of RNA:"
    make_bigwig "${BAM}" "-f 16" "-3" \
        "${BIGWIG_IND_DIR}/${SAMPLE}_5prime_plus.bw"  "5prime_plus"
    make_bigwig "${BAM}" "-F 16" "-3" \
        "${BIGWIG_IND_DIR}/${SAMPLE}_5prime_minus.bw" "5prime_minus"

    # Sanity: file-size strand balance (should be roughly comparable)
    echo ""
    echo "[SANITY CHECK] ${SAMPLE} bigWig file sizes:"
    for END in 3prime 5prime; do
        P_BW="${BIGWIG_IND_DIR}/${SAMPLE}_${END}_plus.bw"
        M_BW="${BIGWIG_IND_DIR}/${SAMPLE}_${END}_minus.bw"
        if [[ -f "${P_BW}" && -f "${M_BW}" ]]; then
            PSIZE=$(du -h "${P_BW}" | cut -f1)
            MSIZE=$(du -h "${M_BW}" | cut -f1)
            echo "  ${END}: plus=${PSIZE}, minus=${MSIZE}"
        fi
    done
done

step_header "Step 04 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output: ${BIGWIG_IND_DIR}/"
echo ""
echo "Next: bash 05_wasp_setup.sh"
echo "      (cross-rep merging happens in step 09)"
