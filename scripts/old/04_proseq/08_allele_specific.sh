#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 08_allele_specific.sh — Split reads by allele + allele-specific bigWigs
# =============================================================================
#
# WHAT THIS DOES (per sample):
#   1. Runs split_by_allele.py on the WASP-filtered BAM:
#      - Assigns each read to ref (129S1), alt (CAST), nosnp, or ambiguous
#      - Generates per-SNP allele count tables
#   2. Sorts and indexes allele-specific BAMs
#   3. Generates strand-specific 3' end bigWigs for each allele
#      (same strand-swap logic as step 04, applied to allele-split BAMs)
#
# OUTPUT (per sample):
#   BAMs:
#     {sample}_ref.bam        — reads carrying the 129S1 (reference) allele
#     {sample}_alt.bam        — reads carrying the CAST (alternative) allele
#     {sample}_nosnp.bam      — reads not overlapping any het SNP
#     {sample}_ambiguous.bam  — reads with conflicting alleles
#
#   BigWigs (3' end = Pol II position, for each allele):
#     {sample}_ref_3prime_plus.bw / minus.bw
#     {sample}_alt_3prime_plus.bw / minus.bw
#
#   Per-SNP allele counts:
#     {sample}_allele_counts.tsv
#
# WHY PER-ALLELE bigWigs?
#   These let you visualize allele-specific transcription in a genome browser.
#   At loci with strong allelic differences, you'll see the signal differ
#   between ref and alt tracks. This is useful for validation and exploration.
#   The per-SNP allele counts table is what you'll use for the systematic
#   TFBS analysis downstream.
#
# PARALLEL EXECUTION:
#   sbatch --array=0-3 --partition=short --time=12:00:00 --mem=32G --cpus-per-task=8 \
#       -o logs/08_%A_%a.out -e logs/08_%A_%a.err 08_allele_specific.sh
# =============================================================================

source "config_wasp.sh"

step_header "PRO-seq Pipeline Step 08: Allele-Specific Processing"

source activate "${CONDA_ENV_NAME}"
create_wasp_dirs
resolve_samples "${1:-}"

SPLIT_SCRIPT="split_by_allele.py"
if [[ ! -f "${SPLIT_SCRIPT}" ]]; then
    echo "[ERROR] split_by_allele.py not found: ${SPLIT_SCRIPT}"
    exit 1
fi

# BigWig helper (same as step 04)
make_bigwig() {
    local BAM="$1"
    local SAM_FLAG="$2"
    local BT_END_FLAG="$3"
    local OUTBW="$4"
    local LABEL="$5"

    local TMPBG="${OUTBW%.bw}.bedGraph"

    if [[ -f "${OUTBW}" ]]; then
        echo "  [INFO] ${LABEL}: exists. Skipping."
        return 0
    fi

    local READ_COUNT
    READ_COUNT=$(samtools view -c ${SAM_FLAG} "${BAM}")

    if [[ "${READ_COUNT}" -eq 0 ]]; then
        echo "  [WARN] ${LABEL}: 0 reads with flag ${SAM_FLAG}. Skipping bigWig."
        return 0
    fi

    echo "  [INFO] ${LABEL}: ${READ_COUNT} reads → bigWig..."

    samtools view -b ${SAM_FLAG} "${BAM}" \
    | bedtools genomecov -ibam stdin -bg ${BT_END_FLAG} \
    | sort -k1,1 -k2,2n \
    | awk -v OFS='\t' 'NR==FNR {chr[$1]=1; next} ($1 in chr)' "${MM39_CHROM_SIZES}" - \
    > "${TMPBG}"

    local TOTAL_COUNTS
    TOTAL_COUNTS=$(awk '{s+=$4} END {print s+0}' "${TMPBG}")
    echo "  [INFO] ${LABEL}: ${TOTAL_COUNTS} total counts"

    if [[ -s "${TMPBG}" ]]; then
        bedGraphToBigWig "${TMPBG}" "${MM39_CHROM_SIZES}" "${OUTBW}"
    else
        echo "  [WARN] ${LABEL}: empty bedGraph. No bigWig created."
    fi
    rm -f "${TMPBG}"
}

for SAMPLE in "${RUN_SAMPLES[@]}"; do
    step_header "Allele-specific processing: ${SAMPLE}"

    WASP_BAM="${WASP_FILTERED_DIR}/${SAMPLE}_wasp.bam"
    OUT_PREFIX="${WASP_ALLELE_DIR}/${SAMPLE}"

    if [[ ! -f "${WASP_BAM}" ]]; then
        echo "[ERROR] WASP-filtered BAM not found: ${WASP_BAM}"
        exit 1
    fi

    # =========================================================================
    # STEP 1: Split reads by allele
    # =========================================================================
    echo "[STEP 1] Splitting reads by allele..."

    if [[ -f "${OUT_PREFIX}_ref.bam" && -f "${OUT_PREFIX}_allele_counts.tsv" ]]; then
        echo "[INFO] Allele-split files already exist. Skipping split."
    else
        python "${SPLIT_SCRIPT}" \
            --bam "${WASP_BAM}" \
            --snp_dir "${WASP_SNP_DIR}" \
            --output_prefix "${OUT_PREFIX}"
    fi

    # =========================================================================
    # STEP 2: Sort and index allele-specific BAMs
    # =========================================================================
    echo ""
    echo "[STEP 2] Sorting and indexing allele-specific BAMs..."

    for ALLELE in ref alt; do
        RAW="${OUT_PREFIX}_${ALLELE}.bam"
        SORTED="${OUT_PREFIX}_${ALLELE}_sorted.bam"

        if [[ -f "${SORTED}" && -f "${SORTED}.bai" ]]; then
            echo "  [INFO] ${ALLELE}: already sorted."
            continue
        fi

        if [[ ! -f "${RAW}" ]]; then
            echo "  [WARN] ${ALLELE}: BAM not found. Skipping."
            continue
        fi

        samtools sort -@ "${THREADS}" -o "${SORTED}" "${RAW}"
        samtools index "${SORTED}"

        # Replace unsorted with sorted
        mv "${SORTED}" "${RAW}"
        mv "${SORTED}.bai" "${RAW}.bai"

        COUNT=$(samtools view -c "${RAW}")
        echo "  [INFO] ${ALLELE}: ${COUNT} reads (sorted + indexed)"
    done

    # =========================================================================
    # STEP 3: Generate allele-specific bigWigs (3' end of nascent RNA)
    # =========================================================================
    echo ""
    echo "[STEP 3] Generating allele-specific bigWigs (3' end = Pol II site)..."

    BWDIR="${WASP_ALLELE_BW_DIR}"

    for ALLELE in ref alt; do
        BAM="${OUT_PREFIX}_${ALLELE}.bam"

        if [[ ! -f "${BAM}" ]]; then
            echo "  [WARN] ${ALLELE}: BAM not found. Skipping bigWigs."
            continue
        fi

        echo ""
        echo "  --- ${ALLELE} allele ---"

        # Strand swap + 3' end logic (identical to step 04):
        # Plus strand signal: BAM − strand reads (flag 16), bedtools -5
        make_bigwig "${BAM}" "-f 16" "-5" \
            "${BWDIR}/${SAMPLE}_${ALLELE}_3prime_plus.bw" \
            "${ALLELE}_3prime_plus"

        # Minus strand signal: BAM + strand reads (NOT flag 16), bedtools -5
        make_bigwig "${BAM}" "-F 16" "-5" \
            "${BWDIR}/${SAMPLE}_${ALLELE}_3prime_minus.bw" \
            "${ALLELE}_3prime_minus"
    done

    # =========================================================================
    # STEP 4: Per-sample QC summary
    # =========================================================================
    echo ""
    echo "[STEP 4] QC summary..."

    COUNTS_FILE="${OUT_PREFIX}_allele_counts.tsv"
    if [[ -f "${COUNTS_FILE}" ]]; then
        TOTAL_SNPS_COVERED=$(awk 'NR>1 && ($5+$6) > 0' "${COUNTS_FILE}" | wc -l)
        TOTAL_REF_READS=$(awk 'NR>1 {s+=$5} END {print s+0}' "${COUNTS_FILE}")
        TOTAL_ALT_READS=$(awk 'NR>1 {s+=$6} END {print s+0}' "${COUNTS_FILE}")
        TOTAL_ALLELIC=$((TOTAL_REF_READS + TOTAL_ALT_READS))

        echo "[SANITY CHECK] Allele-specific counts: ${SAMPLE}"
        echo "  Het SNPs with coverage: ${TOTAL_SNPS_COVERED}"
        echo "  Total ref reads at SNPs: ${TOTAL_REF_READS}"
        echo "  Total alt reads at SNPs: ${TOTAL_ALT_READS}"

        if [[ "${TOTAL_ALLELIC}" -gt 0 ]]; then
            REF_PCT=$(awk "BEGIN {printf \"%.2f\", ${TOTAL_REF_READS}/${TOTAL_ALLELIC}*100}")
            echo "  Global ref fraction: ${REF_PCT}%"
            echo "  (Should be close to 50% after WASP — confirms bias removal)"
        fi

        # Distribution of coverage at SNPs
        echo ""
        echo "  Coverage distribution at het SNPs:"
        awk 'NR>1 {total=$5+$6; if(total>=1) print total}' "${COUNTS_FILE}" \
        | sort -n \
        | awk '
            BEGIN {n=0}
            {vals[n++]=$1}
            END {
                if (n==0) {print "    No covered SNPs"; exit}
                printf "    Min:    %d\n", vals[0]
                printf "    Q1:     %d\n", vals[int(n*0.25)]
                printf "    Median: %d\n", vals[int(n*0.5)]
                printf "    Q3:     %d\n", vals[int(n*0.75)]
                printf "    Max:    %d\n", vals[n-1]
                printf "    SNPs with >=5 reads: %d\n", 0
            }
        '
        # Count SNPs with decent coverage
        SNPS_GE5=$(awk 'NR>1 && ($5+$6) >= 5' "${COUNTS_FILE}" | wc -l)
        echo "    SNPs with >=5 allelic reads: ${SNPS_GE5}"
    fi

    echo ""
    echo "[DONE] ${SAMPLE}"
    echo "  Allele BAMs:     ${OUT_PREFIX}_{ref,alt}.bam"
    echo "  Allele bigWigs:  ${BWDIR}/${SAMPLE}_{ref,alt}_3prime_{plus,minus}.bw"
    echo "  Allele counts:   ${COUNTS_FILE}"
done

step_header "Step 08 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output:"
echo "  Allele BAMs:    ${WASP_ALLELE_DIR}/"
echo "  Allele bigWigs: ${WASP_ALLELE_BW_DIR}/"
echo "  Allele counts:  ${WASP_ALLELE_DIR}/*_allele_counts.tsv"
echo ""
echo "Next: Run 09_wasp_qc.sh to compile cross-sample QC"