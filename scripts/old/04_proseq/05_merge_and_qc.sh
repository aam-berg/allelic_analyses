#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 05_merge_and_qc.sh — Merge replicate bigWigs + compile QC summary
# =============================================================================
#
# WHAT THIS DOES:
#   1. Collects per-sample QC files from parallel array jobs into unified tables
#   2. Merges bigWig files across all 4 replicates (sum of counts)
#   3. If DEDUP_5PRIME=true, also merges deduplicated bigWigs and includes
#      deduplication stats in the QC report
#   4. Compiles a comprehensive QC report
#
# NOTE: This script always runs on ALL samples. Do not use with --array.
#   It should be run AFTER all step 04 array tasks have completed.
#
# SBATCH RESOURCES:
#   sbatch --partition=short --time=2:00:00 --mem=16G --cpus-per-task=4 \
#       -o logs/05_%j.out -e logs/05_%j.err 05_merge_and_qc.sh
#
#   # With dedup:
#   DEDUP_5PRIME=true sbatch ... 05_merge_and_qc.sh
# =============================================================================

source "config.sh"

step_header "PRO-seq Pipeline Step 05: Merge Replicates & Final QC"

source activate "${CONDA_ENV_NAME}"
create_dirs

if [[ "${DEDUP_5PRIME}" == "true" ]]; then
    echo "[INFO] 5' end deduplication: ENABLED — will also merge dedup bigWigs"
else
    echo "[INFO] 5' end deduplication: DISABLED (default)"
fi

# =============================================================================
# 1. Consolidate per-sample QC files from array jobs
# =============================================================================
step_header "STEP 1: Consolidate QC statistics"

# --- Trimming stats ---
TRIM_COMBINED="${QC_DIR}/trimming_stats.tsv"
echo -e "sample\traw_reads\treads_with_adapter\treads_too_short\treads_after_trim\tpct_adapter\tpct_surviving" \
    > "${TRIM_COMBINED}"
for SAMPLE in "${SAMPLE_ORDER[@]}"; do
    F="${QC_DIR}/trimming_stats_${SAMPLE}.tsv"
    if [[ -f "${F}" ]]; then
        cat "${F}" >> "${TRIM_COMBINED}"
    else
        echo "[WARNING] Missing: ${F}"
    fi
done
echo "[INFO] Combined trimming stats:"
column -t "${TRIM_COMBINED}"

# --- Alignment stats ---
ALIGN_COMBINED="${QC_DIR}/alignment_stats.tsv"
# Take header from first per-sample file
FIRST_SAMPLE="${SAMPLE_ORDER[0]}"
if [[ -f "${QC_DIR}/alignment_stats_${FIRST_SAMPLE}.tsv" ]]; then
    head -1 "${QC_DIR}/alignment_stats_${FIRST_SAMPLE}.tsv" > "${ALIGN_COMBINED}"
else
    echo "[WARNING] No alignment stats found. Was step 03 completed?"
fi
for SAMPLE in "${SAMPLE_ORDER[@]}"; do
    F="${QC_DIR}/alignment_stats_${SAMPLE}.tsv"
    if [[ -f "${F}" ]]; then
        tail -1 "${F}" >> "${ALIGN_COMBINED}"
    else
        echo "[WARNING] Missing: ${F}"
    fi
done
echo ""
echo "[INFO] Combined alignment stats:"
column -t "${ALIGN_COMBINED}"

# --- Dedup stats (if applicable) ---
if [[ "${DEDUP_5PRIME}" == "true" ]]; then
    DEDUP_COMBINED="${QC_DIR}/dedup_stats.tsv"
    echo -e "sample\tbefore_dedup\tafter_dedup\tremoved\tpct_duplicates" \
        > "${DEDUP_COMBINED}"
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        F="${QC_DIR}/dedup_stats_${SAMPLE}.tsv"
        if [[ -f "${F}" ]]; then
            cat "${F}" >> "${DEDUP_COMBINED}"
        else
            echo "[WARNING] Missing dedup stats: ${F}"
        fi
    done
    echo ""
    echo "[INFO] Combined dedup stats:"
    column -t "${DEDUP_COMBINED}"
fi

# =============================================================================
# 2. Verify all individual bigWigs exist before merging
# =============================================================================
step_header "STEP 2: Verify individual bigWigs"

ALL_BW_OK=1
for SAMPLE in "${SAMPLE_ORDER[@]}"; do
    for END in 3prime 5prime; do
        for STRAND in plus minus; do
            BW="${BIGWIG_IND_DIR}/${SAMPLE}_${END}_${STRAND}.bw"
            if [[ -f "${BW}" ]]; then
                echo "  [OK]    ${BW}"
            else
                echo "  [MISSING] ${BW}"
                ALL_BW_OK=0
            fi
        done
    done
done

if (( ! ALL_BW_OK )); then
    echo ""
    echo "[ERROR] Some bigWig files are missing. Ensure all step 04 array tasks completed."
    echo "[ERROR] Cannot proceed with merge."
    exit 1
fi
echo ""
echo "[INFO] All 16 individual bigWig files present (4 samples × 4 files)."

# Also verify dedup bigWigs if applicable
if [[ "${DEDUP_5PRIME}" == "true" ]]; then
    echo ""
    echo "[INFO] Verifying deduplicated bigWigs..."
    ALL_DEDUP_BW_OK=1
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        for END in 3prime 5prime; do
            for STRAND in plus minus; do
                BW="${BIGWIG_IND_DIR}/${SAMPLE}${DEDUP_SUFFIX}_${END}_${STRAND}.bw"
                if [[ -f "${BW}" ]]; then
                    echo "  [OK]    ${BW}"
                else
                    echo "  [MISSING] ${BW}"
                    ALL_DEDUP_BW_OK=0
                fi
            done
        done
    done

    if (( ! ALL_DEDUP_BW_OK )); then
        echo ""
        echo "[ERROR] Some deduplicated bigWig files are missing."
        echo "[ERROR] Ensure all step 04 array tasks completed with DEDUP_5PRIME=true."
        exit 1
    fi
    echo ""
    echo "[INFO] All 16 deduplicated bigWig files also present."
fi

# =============================================================================
# 3. Merge replicate bigWigs (standard)
# =============================================================================
step_header "STEP 3: Merge replicate bigWigs"

for END in 3prime 5prime; do
    for STRAND in plus minus; do
        echo ""
        echo "--- Merging ${END}_${STRAND} ---"

        MERGED_BW="${BIGWIG_MERGED_DIR}/WT_PROseq_merged_${END}_${STRAND}.bw"
        MERGED_BG="${BIGWIG_MERGED_DIR}/WT_PROseq_merged_${END}_${STRAND}.bedGraph"

        if [[ -f "${MERGED_BW}" ]]; then
            echo "[INFO] Already exists: ${MERGED_BW}. Skipping."
            continue
        fi

        BW_FILES=()
        for SAMPLE in "${SAMPLE_ORDER[@]}"; do
            BW_FILES+=("${BIGWIG_IND_DIR}/${SAMPLE}_${END}_${STRAND}.bw")
        done

        echo "  Input: ${BW_FILES[*]}"

        # bigWigMerge: sums counts across all input bigWig files at each position
        bigWigMerge "${BW_FILES[@]}" "${MERGED_BG}"

        sort -k1,1 -k2,2n "${MERGED_BG}" -o "${MERGED_BG}"

        # Filter to valid chromosomes
        awk 'NR==FNR {chr[$1]=1; next} ($1 in chr)' "${MM39_CHROM_SIZES}" "${MERGED_BG}" \
            > "${MERGED_BG}.tmp"
        mv "${MERGED_BG}.tmp" "${MERGED_BG}"

        bedGraphToBigWig "${MERGED_BG}" "${MM39_CHROM_SIZES}" "${MERGED_BW}"
        rm -f "${MERGED_BG}"

        echo "  [DONE] ${MERGED_BW} ($(du -h "${MERGED_BW}" | cut -f1))"
    done
done

echo ""
echo "[INFO] Merged bigWig files:"
ls -lh "${BIGWIG_MERGED_DIR}"/WT_PROseq_merged_*.bw 2>/dev/null | grep -v "${DEDUP_SUFFIX}" | sed 's/^/    /'

# =============================================================================
# 3b. Merge deduplicated replicate bigWigs (optional)
# =============================================================================
if [[ "${DEDUP_5PRIME}" == "true" ]]; then
    step_header "STEP 3b: Merge deduplicated replicate bigWigs"

    for END in 3prime 5prime; do
        for STRAND in plus minus; do
            echo ""
            echo "--- Merging ${END}_${STRAND} (dedup) ---"

            MERGED_BW="${BIGWIG_MERGED_DIR}/WT_PROseq_merged${DEDUP_SUFFIX}_${END}_${STRAND}.bw"
            MERGED_BG="${BIGWIG_MERGED_DIR}/WT_PROseq_merged${DEDUP_SUFFIX}_${END}_${STRAND}.bedGraph"

            if [[ -f "${MERGED_BW}" ]]; then
                echo "[INFO] Already exists: ${MERGED_BW}. Skipping."
                continue
            fi

            BW_FILES=()
            for SAMPLE in "${SAMPLE_ORDER[@]}"; do
                BW_FILES+=("${BIGWIG_IND_DIR}/${SAMPLE}${DEDUP_SUFFIX}_${END}_${STRAND}.bw")
            done

            echo "  Input: ${BW_FILES[*]}"

            bigWigMerge "${BW_FILES[@]}" "${MERGED_BG}"

            sort -k1,1 -k2,2n "${MERGED_BG}" -o "${MERGED_BG}"

            awk 'NR==FNR {chr[$1]=1; next} ($1 in chr)' "${MM39_CHROM_SIZES}" "${MERGED_BG}" \
                > "${MERGED_BG}.tmp"
            mv "${MERGED_BG}.tmp" "${MERGED_BG}"

            bedGraphToBigWig "${MERGED_BG}" "${MM39_CHROM_SIZES}" "${MERGED_BW}"
            rm -f "${MERGED_BG}"

            echo "  [DONE] ${MERGED_BW} ($(du -h "${MERGED_BW}" | cut -f1))"
        done
    done

    echo ""
    echo "[INFO] Merged deduplicated bigWig files:"
    ls -lh "${BIGWIG_MERGED_DIR}"/*"${DEDUP_SUFFIX}"*.bw | sed 's/^/    /'
fi

# =============================================================================
# 4. Compile QC report
# =============================================================================
step_header "STEP 4: Compile QC Summary Report"

QC_REPORT="${QC_DIR}/pipeline_qc_report.txt"

cat > "${QC_REPORT}" << 'HEADER'
================================================================================
PRO-seq Processing Pipeline — QC Summary Report
================================================================================
HEADER

{
    echo "Generated: $(date)"
    echo "Genome build: mm39"
    echo "Spike-in genome: dm6"
    echo "Pipeline: rRNA filter → dm6 filter → mm39 align → MAPQ/chrM filter"
    echo "Deduplication: NONE (standard for PRO-seq without UMIs)"
    if [[ "${DEDUP_5PRIME}" == "true" ]]; then
        echo ""
        echo "*** OPTIONAL 5' END DEDUPLICATION: ENABLED ***"
        echo "  Method: samtools markdup -r -s"
        echo "  Logic:  Reads sharing the same 5' mapping position + strand are"
        echo "          grouped. The read with the highest base-quality sum is KEPT;"
        echo "          all others are REMOVED."
        echo "  Suffix: ${DEDUP_SUFFIX}"
        echo "  Note:   Standard (non-dedup) outputs are ALSO produced."
    fi
    echo ""

    echo "--- TRIMMING STATISTICS ---"
    column -t "${TRIM_COMBINED}"
    echo ""

    echo "--- ALIGNMENT STATISTICS ---"
    column -t "${ALIGN_COMBINED}"
    echo ""

    if [[ "${DEDUP_5PRIME}" == "true" && -f "${DEDUP_COMBINED}" ]]; then
        echo "--- 5' END DEDUPLICATION STATISTICS ---"
        column -t "${DEDUP_COMBINED}"
        echo ""
        echo "Interpretation:"
        echo "  <10% removed:  Low duplication (typical for good PRO-seq libraries)"
        echo "  10-30%:        Moderate — may indicate some PCR over-amplification"
        echo "  >30%:          High — library complexity may be low"
        echo "  REMINDER: At strong Pol II pause sites, identical 5' positions are"
        echo "  expected biological signal. High 'duplication' may be real biology."
        echo ""
    fi

    echo "--- FINAL READ COUNTS ---"
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        BAM="${BAM_MM39_DIR}/${SAMPLE}_final.bam"
        if [[ -f "${BAM}" ]]; then
            COUNT=$(samtools view -c "${BAM}")
            PLUS=$(samtools view -c -F 16 "${BAM}")
            MINUS=$(samtools view -c -f 16 "${BAM}")
            RATIO=$(awk "BEGIN {printf \"%.3f\", ${PLUS}/${MINUS:-1}}")
            echo "  ${SAMPLE}: ${COUNT} total  (+: ${PLUS}, -: ${MINUS}, ratio: ${RATIO})"
        fi
    done
    echo ""

    if [[ "${DEDUP_5PRIME}" == "true" ]]; then
        echo "--- FINAL READ COUNTS (5' DEDUPLICATED) ---"
        for SAMPLE in "${SAMPLE_ORDER[@]}"; do
            DEDUP_BAM="${BAM_MM39_DIR}/${SAMPLE}_final${DEDUP_SUFFIX}.bam"
            if [[ -f "${DEDUP_BAM}" ]]; then
                COUNT=$(samtools view -c "${DEDUP_BAM}")
                PLUS=$(samtools view -c -F 16 "${DEDUP_BAM}")
                MINUS=$(samtools view -c -f 16 "${DEDUP_BAM}")
                RATIO=$(awk "BEGIN {printf \"%.3f\", ${PLUS}/${MINUS:-1}}")
                echo "  ${SAMPLE} (dedup): ${COUNT} total  (+: ${PLUS}, -: ${MINUS}, ratio: ${RATIO})"
            fi
        done
        echo ""
    fi

    echo "--- OUTPUT FILES ---"
    echo ""
    echo "Individual bigWigs (${BIGWIG_IND_DIR}/):"
    ls -lh "${BIGWIG_IND_DIR}"/*.bw 2>/dev/null | grep -v "${DEDUP_SUFFIX}" | awk '{print "  "$NF" ("$5")"}' || echo "  (none)"
    echo ""
    if [[ "${DEDUP_5PRIME}" == "true" ]]; then
        echo "Individual bigWigs — deduplicated (${BIGWIG_IND_DIR}/):"
        ls -lh "${BIGWIG_IND_DIR}"/*"${DEDUP_SUFFIX}"*.bw 2>/dev/null | awk '{print "  "$NF" ("$5")"}' || echo "  (none)"
        echo ""
    fi
    echo "Merged bigWigs (${BIGWIG_MERGED_DIR}/):"
    ls -lh "${BIGWIG_MERGED_DIR}"/*.bw 2>/dev/null | grep -v "${DEDUP_SUFFIX}" | awk '{print "  "$NF" ("$5")"}' || echo "  (none)"
    echo ""
    if [[ "${DEDUP_5PRIME}" == "true" ]]; then
        echo "Merged bigWigs — deduplicated (${BIGWIG_MERGED_DIR}/):"
        ls -lh "${BIGWIG_MERGED_DIR}"/*"${DEDUP_SUFFIX}"*.bw 2>/dev/null | awk '{print "  "$NF" ("$5")"}' || echo "  (none)"
        echo ""
    fi
    echo "Final BAMs (${BAM_MM39_DIR}/):"
    ls -lh "${BAM_MM39_DIR}"/*_final.bam 2>/dev/null | awk '{print "  "$NF" ("$5")"}' || echo "  (none)"
    if [[ "${DEDUP_5PRIME}" == "true" ]]; then
        echo ""
        echo "Deduplicated BAMs (${BAM_MM39_DIR}/):"
        ls -lh "${BAM_MM39_DIR}"/*"${DEDUP_SUFFIX}".bam 2>/dev/null | awk '{print "  "$NF" ("$5")"}' || echo "  (none)"
    fi
    echo ""

    cat << 'GUIDE'
--- QUALITY ASSESSMENT GUIDE ---

1. ADAPTER CONTENT: 30-60% is normal for PRO-seq
2. rRNA RATE: 5-20% is typical; Pol I transcription produces many reads
3. SPIKE-IN RATE: 1-10% is typical
4. MM39 ALIGNMENT: >70% is good for PRO-seq
5. STRAND BALANCE: +/- ratio should be ~1.0 (0.8-1.2)
6. FINAL READ COUNT: >5M per replicate is adequate; >10M is good; >20M is excellent
7. REPLICATE CONCORDANCE: Compare bigWigs in IGV; patterns should be similar

For WASP allele-specific analysis:
  Use *_final.bam files as input to WASP find_intersecting_snps.py
  with F121-9_het_snps.vcf.gz from the VCF directory.
GUIDE

} >> "${QC_REPORT}"

echo ""
echo "--- QC Report ---"
cat "${QC_REPORT}"

# =============================================================================
# Done
# =============================================================================
step_header "PIPELINE COMPLETE"
echo ""
echo "FINAL OUTPUT SUMMARY"
echo ""
echo "Per-replicate bigWigs (raw counts, single-nt resolution):"
echo "  ${BIGWIG_IND_DIR}/"
for SAMPLE in "${SAMPLE_ORDER[@]}"; do
    echo "    ${SAMPLE}_3prime_{plus,minus}.bw"
    echo "    ${SAMPLE}_5prime_{plus,minus}.bw"
done
echo ""
if [[ "${DEDUP_5PRIME}" == "true" ]]; then
    echo "Per-replicate bigWigs (5' deduplicated):"
    echo "  ${BIGWIG_IND_DIR}/"
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        echo "    ${SAMPLE}${DEDUP_SUFFIX}_3prime_{plus,minus}.bw"
        echo "    ${SAMPLE}${DEDUP_SUFFIX}_5prime_{plus,minus}.bw"
    done
    echo ""
fi
echo "Merged bigWigs (sum of 4 replicates):"
echo "  ${BIGWIG_MERGED_DIR}/"
echo "    WT_PROseq_merged_{3prime,5prime}_{plus,minus}.bw"
if [[ "${DEDUP_5PRIME}" == "true" ]]; then
    echo "    WT_PROseq_merged${DEDUP_SUFFIX}_{3prime,5prime}_{plus,minus}.bw"
fi
echo ""
echo "Final BAMs (for WASP):"
echo "  ${BAM_MM39_DIR}/*_final.bam"
if [[ "${DEDUP_5PRIME}" == "true" ]]; then
    echo "  ${BAM_MM39_DIR}/*_final${DEDUP_SUFFIX}.bam"
fi
echo ""
echo "QC report: ${QC_REPORT}"