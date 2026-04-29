#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04c_rnaseq/07_make_bigwigs_allele.sh — Per-allele strand-specific bigWigs
# =============================================================================
#
# WHAT THIS DOES (per sample):
#   Generates 4 strand-specific bigWigs from the allele-split BAMs:
#     {SAMPLE}_ref_plus.bw     REF allele, RNA + strand
#     {SAMPLE}_ref_minus.bw    REF allele, RNA - strand
#     {SAMPLE}_alt_plus.bw     ALT allele, RNA + strand
#     {SAMPLE}_alt_minus.bw    ALT allele, RNA - strand
#
#   Same strand-flag logic as step 05 (TruSeq Stranded reverse-stranded
#   library) applied to the ref.bam and alt.bam files from step 06.
#
# PARALLEL EXECUTION:
#   sbatch --array=0-1 --partition=short --time=2:00:00 --mem=16G \
#       --cpus-per-task=4 -o logs/07_%A_%a.out -e logs/07_%A_%a.err \
#       07_make_bigwigs_allele.sh
# =============================================================================

source "config.sh"

step_header "RNA-seq Pipeline Step 07: Per-Allele Strand-Specific BigWigs"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] samtools:         $(samtools --version | head -1)"
echo "[INFO] bedtools:         $(bedtools --version)"
echo "[INFO] bedGraphToBigWig: $(which bedGraphToBigWig)"

create_dirs
resolve_samples "${1:-}"

[[ -f "${MM39_CHROM_SIZES}" ]] || \
    { echo "[ERROR] mm39 chrom sizes missing: ${MM39_CHROM_SIZES}" >&2; exit 1; }

# -----------------------------------------------------------------------------
# Helper: same as step 05 — extract RNA-strand fragments and convert to bigwig
# -----------------------------------------------------------------------------
make_strand_bigwig() {
    local BAM="$1"
    local STRAND="$2"
    local OUTBW="$3"
    local LABEL="$4"

    if [[ -f "${OUTBW}" ]]; then
        echo "  [SKIP] ${LABEL}: bigWig already exists."
        return 0
    fi
    if [[ ! -f "${BAM}" ]]; then
        echo "  [SKIP] ${LABEL}: input BAM missing."
        return 0
    fi

    local TMP="${OUTBW%.bw}"
    local R1_BAM="${TMP}_r1.bam"
    local R2_BAM="${TMP}_r2.bam"
    local FRAG_BAM="${TMP}_frag.bam"
    local TMPBG="${TMP}.bedGraph"

    if [[ "${STRAND}" == "plus" ]]; then
        samtools view -bh -f 80 "${BAM}"             > "${R1_BAM}"
        samtools view -bh -f 128 -F 16 "${BAM}"      > "${R2_BAM}"
    elif [[ "${STRAND}" == "minus" ]]; then
        samtools view -bh -f 64 -F 16 "${BAM}"       > "${R1_BAM}"
        samtools view -bh -f 144 "${BAM}"            > "${R2_BAM}"
    else
        echo "[ERROR] Unknown strand: ${STRAND}" >&2
        return 1
    fi

    samtools merge -f -@ "${THREADS}" "${FRAG_BAM}" "${R1_BAM}" "${R2_BAM}"
    samtools sort -@ "${THREADS}" -o "${FRAG_BAM}.sorted" "${FRAG_BAM}"
    mv "${FRAG_BAM}.sorted" "${FRAG_BAM}"
    samtools index "${FRAG_BAM}"

    local READ_COUNT
    READ_COUNT=$(samtools view -c "${FRAG_BAM}")

    if (( READ_COUNT == 0 )); then
        echo "  [WARN] ${LABEL}: 0 reads after strand split"
        rm -f "${R1_BAM}" "${R2_BAM}" "${FRAG_BAM}" "${FRAG_BAM}.bai"
        return 0
    fi

    bedtools genomecov -ibam "${FRAG_BAM}" -bg -pc \
    | sort -k1,1 -k2,2n \
    | awk 'NR==FNR {chr[$1]=1; next} ($1 in chr)' "${MM39_CHROM_SIZES}" - \
    > "${TMPBG}"

    if [[ ! -s "${TMPBG}" ]]; then
        echo "  [WARN] ${LABEL}: empty bedGraph"
        rm -f "${R1_BAM}" "${R2_BAM}" "${FRAG_BAM}" "${FRAG_BAM}.bai" "${TMPBG}"
        return 0
    fi

    bedGraphToBigWig "${TMPBG}" "${MM39_CHROM_SIZES}" "${OUTBW}"
    rm -f "${R1_BAM}" "${R2_BAM}" "${FRAG_BAM}" "${FRAG_BAM}.bai" "${TMPBG}"

    echo "  [OK]   ${LABEL}: ${READ_COUNT} reads -> ${OUTBW}"
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
        make_strand_bigwig "${BAM}" "plus" \
            "${BIGWIG_ALLELE_DIR}/${SAMPLE}_${ALLELE}_plus.bw"  "${ALLELE}_plus"
        make_strand_bigwig "${BAM}" "minus" \
            "${BIGWIG_ALLELE_DIR}/${SAMPLE}_${ALLELE}_minus.bw" "${ALLELE}_minus"
    done
done

step_header "Step 07 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output: ${BIGWIG_ALLELE_DIR}/"
echo ""
echo "Next: bash 08_merge_qc.sh"
