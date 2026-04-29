#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04_proseq/08_make_bigwigs_allele.sh — Per-allele strand-specific bigWigs
# =============================================================================
#
# WHAT THIS DOES (per sample):
#   Generates 4 single-nt-resolution bigWigs from the allele-split BAMs:
#     {SAMPLE}_ref_3prime_plus.bw    REF allele, plus-strand RNA, 3' end
#     {SAMPLE}_ref_3prime_minus.bw   REF allele, minus-strand RNA, 3' end
#     {SAMPLE}_alt_3prime_plus.bw    ALT allele, plus-strand RNA, 3' end
#     {SAMPLE}_alt_3prime_minus.bw   ALT allele, minus-strand RNA, 3' end
#
#   These are 3'-end (Pol II position) only — that's the analysis-relevant
#   signal for pause biology. 5'-end allele bigWigs are not produced because
#   they're not used downstream; the per-rep 5'-end bigWigs in step 04 cover
#   genome-wide 5' visualization needs.
#
# WHY 3' END ONLY (not 5'):
#   The motif-pause analysis in 06_allele_pairs/ uses 3' end (Pol II active
#   site) signal in narrow windows around motifs. The 5' end signal exists
#   but isn't part of the pause-window analysis. Producing per-allele 5'
#   bigWigs would double the file count without supporting the planned use.
#   If 5' end allele tracks become useful later, this script can easily
#   add them.
#
# STRAND SWAP (canonical, see scripts/README.md section 7):
#   Plus-strand RNA signal:    BAM "-" reads (FLAG 16 SET, `-f 16`), bedtools -5
#   Minus-strand RNA signal:   BAM "+" reads (FLAG 16 NOT SET, `-F 16`), bedtools -5
#
# PARALLEL EXECUTION:
#   sbatch --array=0-3 --partition=short --time=2:00:00 --mem=8G \
#       --cpus-per-task=4 -o logs/08_%A_%a.out -e logs/08_%A_%a.err \
#       08_make_bigwigs_allele.sh
# =============================================================================

source "config.sh"

step_header "PRO-seq Pipeline Step 08: Per-Allele BigWigs"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] bedtools:         $(bedtools --version)"
echo "[INFO] bedGraphToBigWig: $(which bedGraphToBigWig)"

create_dirs
resolve_samples "${1:-}"

[[ -f "${MM39_CHROM_SIZES}" ]] || \
    { echo "[ERROR] mm39 chrom sizes missing: ${MM39_CHROM_SIZES}" >&2; exit 1; }

# -----------------------------------------------------------------------------
# Helper (same as step 04 but kept local for this script's self-containment)
# -----------------------------------------------------------------------------
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

    if [[ ! -f "${BAM}" ]]; then
        echo "  [SKIP] ${LABEL}: input BAM missing (${BAM})"
        return 0
    fi

    local READ_COUNT
    READ_COUNT=$(samtools view -c ${SAM_FLAG} "${BAM}")
    if (( READ_COUNT == 0 )); then
        echo "  [WARN] ${LABEL}: 0 reads with flag ${SAM_FLAG}; skipping bigWig."
        return 0
    fi

    local TMPBG="${OUTBW%.bw}.bedGraph"
    samtools view -b ${SAM_FLAG} "${BAM}" \
    | bedtools genomecov -ibam stdin -bg ${BT_END_FLAG} \
    | sort -k1,1 -k2,2n \
    | awk 'NR==FNR {chr[$1]=1; next} ($1 in chr)' "${MM39_CHROM_SIZES}" - \
    > "${TMPBG}"

    if [[ ! -s "${TMPBG}" ]]; then
        echo "  [WARN] ${LABEL}: empty bedGraph; no bigWig created."
        rm -f "${TMPBG}"
        return 0
    fi

    local TOTAL_COUNTS
    TOTAL_COUNTS=$(awk '{s+=$4} END {print s+0}' "${TMPBG}")
    bedGraphToBigWig "${TMPBG}" "${MM39_CHROM_SIZES}" "${OUTBW}"
    rm -f "${TMPBG}"
    echo "  [OK]   ${LABEL}: ${READ_COUNT} reads, ${TOTAL_COUNTS} total counts -> ${OUTBW}"
}

# -----------------------------------------------------------------------------
# Loop over samples and alleles
# -----------------------------------------------------------------------------
for SAMPLE in "${RUN_SAMPLES[@]}"; do
    step_header "Per-allele bigWigs: ${SAMPLE}"

    for ALLELE in ref alt; do
        BAM="${ALLELE_DIR}/${SAMPLE}_${ALLELE}.bam"
        echo ""
        echo "  --- ${ALLELE} allele ---"

        # Plus-strand RNA: BAM minus reads, bedtools -5
        make_bigwig "${BAM}" "-f 16" "-5" \
            "${BIGWIG_ALLELE_DIR}/${SAMPLE}_${ALLELE}_3prime_plus.bw"  \
            "${ALLELE}_3prime_plus"

        # Minus-strand RNA: BAM plus reads, bedtools -5
        make_bigwig "${BAM}" "-F 16" "-5" \
            "${BIGWIG_ALLELE_DIR}/${SAMPLE}_${ALLELE}_3prime_minus.bw" \
            "${ALLELE}_3prime_minus"
    done
done

step_header "Step 08 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output: ${BIGWIG_ALLELE_DIR}/"
echo ""
echo "Next: bash 09_merge_qc.sh"
