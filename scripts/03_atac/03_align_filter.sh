#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 03_align_filter.sh — bowtie2 PE align → filter (MAPQ + chrM + dedup)
# =============================================================================
#
# WHAT THIS DOES (per sample):
#   1. Aligns trimmed PE reads to mm39 with bowtie2 (--very-sensitive, -X 1000)
#   2. Filters:
#        - properly paired only (-f 2)
#        - MAPQ >= 30
#        - removes unmapped, secondary, supplementary, QC-fail, duplicates
#          (samtools -F 1804)
#        - removes chrM
#   3. Marks and removes PCR duplicates with samtools markdup
#   4. Indexes the final BAM
#
# WHY THESE FILTERS:
#   --very-sensitive:    standard ATAC; precision matters for narrow peaks
#   -X 1000:             ATAC fragments range to ~1 kb (di- and tri-nucleosomal)
#   --no-mixed --no-discordant: only keep proper pairs; ATAC analyses depend on
#                        fragment-level information
#   MAPQ >= 30:          ATAC community standard, removes multi-mappers
#   -F 1804:             ENCODE ATAC standard exclusion mask
#   chrM removal:        mt-DNA is open and Tn5-accessible, often 30-60% of
#                        an ATAC library; removing here gets us a cleaner
#                        nuclear signal
#   markdup -r:          ATAC dups are PCR artifacts (unlike PRO-seq pause
#                        positions); removed before peak calling
#
# WHY samtools markdup INSTEAD OF picard MarkDuplicates:
#   Faster, no Java, equivalent results. Requires fixmate first to add the
#   ms (mate score) and MC tags it needs.
#
# PARALLEL EXECUTION:
#   sbatch --array=0-3 --partition=short --time=12:00:00 --mem=32G \
#       --cpus-per-task=8 -o logs/03_%A_%a.out -e logs/03_%A_%a.err \
#       03_align_filter.sh
# =============================================================================

source "config.sh"

step_header "ATAC Pipeline Step 03: Align, Filter, Dedup"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] bowtie2:  $(bowtie2 --version 2>&1 | head -1)"
echo "[INFO] samtools: $(samtools --version | head -1)"

create_dirs
resolve_samples "${1:-}"

# Verify reference
if [[ ! -f "${MM39_BT2_IDX}.1.bt2" ]]; then
    echo "[ERROR] mm39 bowtie2 index missing at ${MM39_BT2_IDX}" >&2
    echo "[ERROR] Run 00_setup_references.sh first." >&2
    exit 1
fi

for SAMPLE in "${RUN_SAMPLES[@]}"; do
    step_header "Processing ${SAMPLE}"

    R1="${TRIMMED_DIR}/${SAMPLE}_R1_trimmed.fastq.gz"
    R2="${TRIMMED_DIR}/${SAMPLE}_R2_trimmed.fastq.gz"

    RAW_BAM="${BAM_RAW_DIR}/${SAMPLE}_raw.bam"
    FILT_BAM="${BAM_FINAL_DIR}/${SAMPLE}_filtered.bam"   # intermediate
    FINAL_BAM="${BAM_FINAL_DIR}/${SAMPLE}_final.bam"     # deduplicated, indexed

    if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
        echo "[ERROR] Trimmed FASTQs missing for ${SAMPLE}." >&2
        exit 1
    fi

    # ===========================================================================
    # STEP A — Align with bowtie2 PE
    # ===========================================================================
    BOWTIE2_LOG="${LOG_DIR}/${SAMPLE}_bowtie2.log"

    if [[ -f "${RAW_BAM}" && -f "${RAW_BAM}.bai" ]]; then
        echo "[STEP A] Raw BAM exists. Skipping alignment."
    else
        echo "[STEP A] Aligning ${SAMPLE} with bowtie2 PE..."
        bowtie2 \
            --very-sensitive \
            -X "${MAX_INSERT}" \
            --no-mixed \
            --no-discordant \
            --no-unal \
            -x "${MM39_BT2_IDX}" \
            -1 "${R1}" -2 "${R2}" \
            -p "${THREADS}" \
            2> "${BOWTIE2_LOG}" \
        | samtools view -bS -@ 2 - \
        | samtools sort -@ "${THREADS}" -o "${RAW_BAM}"

        samtools index "${RAW_BAM}"
        echo "[STEP A] Done. Raw BAM: ${RAW_BAM}"
    fi

    # Parse bowtie2 stats
    TOTAL_PAIRS=$(grep "reads; of these:" "${BOWTIE2_LOG}" | awk '{print $1}' | head -1)
    CONCORDANT_1=$(grep "aligned concordantly exactly 1 time" "${BOWTIE2_LOG}" | awk '{print $1}')
    CONCORDANT_M=$(grep "aligned concordantly >1 times" "${BOWTIE2_LOG}" | awk '{print $1}')
    OVERALL=$(grep "overall alignment rate" "${BOWTIE2_LOG}" | awk '{print $1}')

    echo ""
    echo "[SANITY CHECK] bowtie2 alignment:"
    echo "  Input pairs:                ${TOTAL_PAIRS}"
    echo "  Concordant unique:          ${CONCORDANT_1:-0}"
    echo "  Concordant multi:           ${CONCORDANT_M:-0}"
    echo "  Overall alignment rate:     ${OVERALL}"
    echo "[INTERPRETATION]"
    echo "  Overall rate >85%: good for ATAC"
    echo "  Overall rate <70%: investigate (contamination? wrong genome?)"

    # ===========================================================================
    # STEP B — Filter: MAPQ + properly paired + remove chrM
    # ===========================================================================
    if [[ -f "${FINAL_BAM}" && -f "${FINAL_BAM}.bai" ]]; then
        echo ""
        echo "[STEP B/C] Final BAM exists. Skipping filter + dedup."
    else
        echo ""
        echo "[STEP B] Filtering: MAPQ>=${MAPQ_THRESHOLD}, properly paired, no chrM..."

        # -F 1804 = exclude unmapped (4) + mate unmapped (8) + secondary (256)
        #         + QC fail (512) + duplicate (1024)
        # -f 2    = require properly paired
        # -q THR  = MAPQ threshold
        # awk     = filter out chrM (header passes through unchanged)
        samtools view -b -F 1804 -f 2 -q "${MAPQ_THRESHOLD}" "${RAW_BAM}" \
        | samtools view -h - \
        | awk 'BEGIN{OFS="\t"} /^@/ {print; next} $3 != "chrM" {print}' \
        | samtools view -bS - \
        | samtools sort -@ "${THREADS}" -o "${FILT_BAM}"

        samtools index "${FILT_BAM}"

        # ===========================================================================
        # STEP C — Deduplicate with samtools markdup
        # ===========================================================================
        DEDUP_LOG="${LOG_DIR}/${SAMPLE}_markdup.log"
        echo "[STEP C] Deduplicating with samtools markdup..."

        # markdup needs: name-sorted -> fixmate (-m adds ms/MC tags) ->
        # coord-sort -> markdup (-r removes dups)
        samtools sort -n -@ "${THREADS}" "${FILT_BAM}" \
        | samtools fixmate -m - - \
        | samtools sort -@ "${THREADS}" - \
        | samtools markdup -r -s - "${FINAL_BAM}" \
            2> "${DEDUP_LOG}"

        samtools index "${FINAL_BAM}"

        # Clean up the intermediate filtered BAM (we keep only the final)
        rm -f "${FILT_BAM}" "${FILT_BAM}.bai"
    fi

    # ===========================================================================
    # Sanity stats
    # ===========================================================================
    AFTER_ALIGN=$(samtools view -c "${RAW_BAM}")
    AFTER_FILT_STAGE=$(samtools view -c -F 1804 -f 2 -q "${MAPQ_THRESHOLD}" "${RAW_BAM}")
    CHRM_RAW=$(samtools view -c "${RAW_BAM}" chrM 2>/dev/null || echo "0")
    FINAL_COUNT=$(samtools view -c "${FINAL_BAM}")

    DEDUP_LOG="${LOG_DIR}/${SAMPLE}_markdup.log"
    DUPS_REMOVED=$(grep -E "DUPLICATE PAIR:|DUPLICATE TOTAL:" "${DEDUP_LOG}" 2>/dev/null \
                   | head -1 | awk '{print $NF}' || echo "?")
    DUPS_PCT="?"
    if [[ "${DUPS_REMOVED}" =~ ^[0-9]+$ ]] && (( AFTER_FILT_STAGE > 0 )); then
        DUPS_PCT=$(awk "BEGIN {printf \"%.2f\", ${DUPS_REMOVED}/${AFTER_FILT_STAGE}*100}")
    fi

    CHRM_PCT=$(awk "BEGIN {printf \"%.2f\", ${CHRM_RAW}/${AFTER_ALIGN:-1}*100}")

    echo ""
    echo "[SANITY CHECK] ${SAMPLE} filtering cascade:"
    printf "  After bowtie2:           %12d\n" "${AFTER_ALIGN}"
    printf "  chrM in raw:             %12d  (%s%%)\n" "${CHRM_RAW}" "${CHRM_PCT}"
    printf "  After MAPQ+paired+chrM:  %12d\n" "${AFTER_FILT_STAGE}"
    printf "  Duplicates removed:      %12s  (%s%%)\n" "${DUPS_REMOVED}" "${DUPS_PCT}"
    printf "  Final reads:             %12d\n" "${FINAL_COUNT}"
    echo ""
    echo "[INTERPRETATION]"
    echo "  chrM <30%: typical good prep"
    echo "  chrM >50%: poor mitochondrial depletion; library still usable"
    echo "  chrM >70%: very high; effective depth significantly reduced"
    echo "  Dup rate <30%: typical good library"
    echo "  Dup rate 30-50%: moderate; common with limited input"
    echo "  Dup rate >50%: low complexity; downstream stats may be weak"

    # Strand balance sanity check (should be ~50/50 for ATAC)
    PLUS=$(samtools view -c -F 16 "${FINAL_BAM}")
    MINUS=$(samtools view -c -f 16 "${FINAL_BAM}")
    if [[ "${MINUS}" -gt 0 ]]; then
        STRAND_RATIO=$(awk "BEGIN {printf \"%.3f\", ${PLUS}/${MINUS}}")
    else
        STRAND_RATIO="NA"
    fi
    echo "  + strand: ${PLUS}, - strand: ${MINUS}, ratio: ${STRAND_RATIO} (expect ~1.0)"

    # ---- Per-sample QC (parallel-safe) ----
    cat > "${QC_DIR}/alignment_stats_${SAMPLE}.tsv" << EOF
sample	total_pairs	overall_rate	after_align	chrM_raw	chrM_pct	after_filter	dups_removed	dups_pct	final	plus	minus	strand_ratio
${SAMPLE}	${TOTAL_PAIRS}	${OVERALL}	${AFTER_ALIGN}	${CHRM_RAW}	${CHRM_PCT}	${AFTER_FILT_STAGE}	${DUPS_REMOVED}	${DUPS_PCT}	${FINAL_COUNT}	${PLUS}	${MINUS}	${STRAND_RATIO}
EOF

done

step_header "Step 03 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output:"
echo "  Final BAMs:       ${BAM_FINAL_DIR}/*_final.bam"
echo "  Per-sample stats: ${QC_DIR}/alignment_stats_*.tsv"
echo ""
echo "Next: bash 04_call_peaks.sh per-rep   (then 'consensus' after array completes)"
