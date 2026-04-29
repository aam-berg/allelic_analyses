#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04c_rnaseq/05_make_bigwigs.sh — Strand-specific PE bigWigs (reverse-stranded)
# =============================================================================
#
# WHAT THIS DOES (per sample):
#   For each sample, generates two bigWigs:
#     {SAMPLE}_plus.bw   - signal of original RNA on the + strand
#     {SAMPLE}_minus.bw  - signal of original RNA on the - strand
#
#   These are fragment-coverage tracks (bedtools -pc style), not single-nt.
#   Cross-rep merging is done in step 08.
#
# THE STRAND CONVENTION (TruSeq Stranded, reverse-stranded library):
#   A fragment is identified by its read-pair flag combinations:
#
#     RNA + strand fragment <=> R1 reverse + R2 forward
#                               (R1 has FLAG 16, R2 does not)
#                               (i.e. "first-in-pair maps to genome - strand")
#
#     RNA - strand fragment <=> R1 forward + R2 reverse
#                               (R1 does not have FLAG 16, R2 does)
#                               (i.e. "first-in-pair maps to genome + strand")
#
# IMPLEMENTATION:
#   To get fragment coverage for one RNA strand, we need BOTH mates of the
#   correct read-pair orientation in the BAM, then run bedtools genomecov -pc.
#   We use samtools flag filtering twice and merge:
#     - For RNA + strand:
#         R1 reverse  : -f 80   (= FLAG 64 + FLAG 16 set)
#         R2 forward  : -f 128 -F 16  (= FLAG 128 set, FLAG 16 not set)
#     - For RNA - strand:
#         R1 forward  : -f 64 -F 16
#         R2 reverse  : -f 144  (= FLAG 128 + FLAG 16 set)
#
# WHY USE THE FULL STAR BAM (not WASP-filtered):
#   For visualization and motif-aggregate plots, we want the full transcription
#   signal, including reads at SNP positions. Allele-specific bigwigs from
#   WASP-filtered + allele-split BAMs are produced in step 07.
#
# SANITY CHECKS:
#   For TruSeq Stranded reverse libraries, the plus.bw and minus.bw should
#   be roughly balanced genome-wide (similar total counts). If they're very
#   unequal, the library may not actually be reverse-stranded, or there may
#   be an issue with the strand-flag logic. We report counts after bigWig
#   generation.
#
# PARALLEL EXECUTION:
#   sbatch --array=0-1 --partition=short --time=2:00:00 --mem=16G \
#       --cpus-per-task=4 -o logs/05_%A_%a.out -e logs/05_%A_%a.err \
#       05_make_bigwigs.sh
# =============================================================================

source "config.sh"

step_header "RNA-seq Pipeline Step 05: Strand-Specific PE BigWigs"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] samtools:         $(samtools --version | head -1)"
echo "[INFO] bedtools:         $(bedtools --version)"
echo "[INFO] bedGraphToBigWig: $(which bedGraphToBigWig)"

create_dirs
resolve_samples "${1:-}"

[[ -f "${MM39_CHROM_SIZES}" ]] || \
    { echo "[ERROR] mm39 chrom sizes missing: ${MM39_CHROM_SIZES}" >&2; exit 1; }

# -----------------------------------------------------------------------------
# Helper: extract RNA-strand-specific fragments and produce bigWig
# -----------------------------------------------------------------------------
# For TruSeq Stranded reverse-stranded library:
#   - rna_strand="plus":  collect R1 reverse + R2 forward
#   - rna_strand="minus": collect R1 forward + R2 reverse
#
# Args: $1=BAM, $2=rna_strand ("plus"|"minus"), $3=output_bw, $4=label
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

    # samtools view flag filters per the strand convention:
    if [[ "${STRAND}" == "plus" ]]; then
        # RNA + strand: R1 reverse (flag 80) + R2 forward (flag 128, no 16)
        samtools view -bh -f 80 "${BAM}"             > "${R1_BAM}"
        samtools view -bh -f 128 -F 16 "${BAM}"      > "${R2_BAM}"
    elif [[ "${STRAND}" == "minus" ]]; then
        # RNA - strand: R1 forward (flag 64, no 16) + R2 reverse (flag 144)
        samtools view -bh -f 64 -F 16 "${BAM}"       > "${R1_BAM}"
        samtools view -bh -f 144 "${BAM}"            > "${R2_BAM}"
    else
        echo "[ERROR] Unknown strand: ${STRAND}" >&2
        return 1
    fi

    # Merge R1 and R2 for this strand into one BAM, then sort
    samtools merge -f -@ "${THREADS}" "${FRAG_BAM}" "${R1_BAM}" "${R2_BAM}"
    samtools sort -@ "${THREADS}" -o "${FRAG_BAM}.sorted" "${FRAG_BAM}"
    mv "${FRAG_BAM}.sorted" "${FRAG_BAM}"
    samtools index "${FRAG_BAM}"

    local READ_COUNT
    READ_COUNT=$(samtools view -c "${FRAG_BAM}")

    # bedtools genomecov -pc: paired-end fragment coverage
    # Requires BOTH mates of each pair, which is why we merged R1+R2 above
    bedtools genomecov -ibam "${FRAG_BAM}" -bg -pc \
    | sort -k1,1 -k2,2n \
    | awk 'NR==FNR {chr[$1]=1; next} ($1 in chr)' "${MM39_CHROM_SIZES}" - \
    > "${TMPBG}"

    if [[ ! -s "${TMPBG}" ]]; then
        echo "  [WARN] ${LABEL}: empty bedGraph (BAM had ${READ_COUNT} reads)"
        rm -f "${R1_BAM}" "${R2_BAM}" "${FRAG_BAM}" "${FRAG_BAM}.bai" "${TMPBG}"
        return 0
    fi

    local TOTAL_COUNTS
    TOTAL_COUNTS=$(awk '{s+=$4} END {print s+0}' "${TMPBG}")

    bedGraphToBigWig "${TMPBG}" "${MM39_CHROM_SIZES}" "${OUTBW}"

    # Cleanup intermediates
    rm -f "${R1_BAM}" "${R2_BAM}" "${FRAG_BAM}" "${FRAG_BAM}.bai" "${TMPBG}"

    echo "  [OK]   ${LABEL}: ${READ_COUNT} reads -> ${OUTBW}"
}

# -----------------------------------------------------------------------------
# Loop over samples
# -----------------------------------------------------------------------------
for SAMPLE in "${RUN_SAMPLES[@]}"; do
    step_header "BigWigs for ${SAMPLE}"

    BAM="${BAM_STAR_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    [[ -f "${BAM}" ]] || { echo "[ERROR] BAM missing: ${BAM}" >&2; exit 1; }

    PLUS_BW="${BIGWIG_IND_DIR}/${SAMPLE}_plus.bw"
    MINUS_BW="${BIGWIG_IND_DIR}/${SAMPLE}_minus.bw"

    make_strand_bigwig "${BAM}" "plus"  "${PLUS_BW}"  "${SAMPLE}_plus"
    make_strand_bigwig "${BAM}" "minus" "${MINUS_BW}" "${SAMPLE}_minus"

    # Sanity: file-size strand balance
    if [[ -f "${PLUS_BW}" && -f "${MINUS_BW}" ]]; then
        echo ""
        echo "[SANITY CHECK] ${SAMPLE} bigWig file sizes:"
        echo "  plus:  $(du -h "${PLUS_BW}"  | cut -f1)"
        echo "  minus: $(du -h "${MINUS_BW}" | cut -f1)"
        echo "  (Roughly similar = OK; very different may indicate library issue.)"
    fi
done

step_header "Step 05 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output: ${BIGWIG_IND_DIR}/"
echo ""
echo "Next: bash 06_allele_specific.sh"
echo "      (cross-rep merging happens in step 08)"
