#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04_proseq/03_align_filter.sh — rRNA filter → dm6 filter → mm39 → MAPQ + chrM
# =============================================================================
#
# WHAT THIS DOES (per sample):
#   1. Align to rRNA reference, DISCARD reads that map to rRNA
#   2. Align surviving reads to dm6 (Drosophila spike-in), DISCARD
#      reads that map to dm6 — but ONLY IF DM6_SPIKEIN_FILTER=true
#   3. Align surviving reads to mm39 with bowtie2 --very-sensitive
#   4. Filter: MAPQ >= MAPQ_THRESHOLD, remove secondary/supplementary,
#      remove unmapped, remove chrM
#   5. Generate per-sample QC stats
#
# CRITICAL: NO DEDUPLICATION
#   The previous version of this script had an optional 5'-end dedup flag.
#   It has been REMOVED. PRO-seq measures Pol II positions at single-nt
#   resolution; at strong pause sites many Pol II molecules sit at the same
#   nucleotide, and dedup would collapse real biology into one read.
#
# CHANGES FROM PREVIOUS VERSION:
#   - All DEDUP_5PRIME flag handling removed.
#   - rRNA bowtie2 changed from --very-fast to default (--sensitive). For a
#     filter we want to catch ALL rRNA reads; --very-fast risks letting some
#     leak through to mm39 alignment.
#   - dm6 step is conditional on DM6_SPIKEIN_FILTER. If false, the rRNA
#     unmapped FASTQ is fed straight into mm39.
#
# WHY BOWTIE2 (NOT STAR/HISAT2) FOR PRO-seq:
#   PRO-seq measures NASCENT RNA which has not been spliced. Splice-aware
#   aligners would incorrectly split reads across nonexistent splice
#   junctions. bowtie2 (non-splice-aware) is the field standard.
#
# PARALLEL EXECUTION:
#   sbatch --array=0-3 --partition=short --time=12:00:00 --mem=32G \
#       --cpus-per-task=8 -o logs/03_%A_%a.out -e logs/03_%A_%a.err \
#       03_align_filter.sh
# =============================================================================

source "config.sh"

step_header "PRO-seq Pipeline Step 03: Align + Filter (NO DEDUP)"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] bowtie2:  $(bowtie2 --version 2>&1 | head -1)"
echo "[INFO] samtools: $(samtools --version | head -1)"
echo "[INFO] dm6 spike-in filter: ${DM6_SPIKEIN_FILTER}"

create_dirs
resolve_samples "${1:-}"

# Verify reference indices
for IDX in "${RRNA_BT2_IDX}" "${MM39_BT2_IDX}"; do
    if [[ ! -f "${IDX}.1.bt2" ]]; then
        echo "[ERROR] Missing bowtie2 index: ${IDX}.1.bt2" >&2
        echo "[ERROR] Run 00_setup_references.sh first." >&2
        exit 1
    fi
done
if [[ "${DM6_SPIKEIN_FILTER}" == "true" && ! -f "${DM6_BT2_IDX}.1.bt2" ]]; then
    echo "[ERROR] DM6_SPIKEIN_FILTER=true but dm6 index missing." >&2
    exit 1
fi

for SAMPLE in "${RUN_SAMPLES[@]}"; do
    step_header "Processing ${SAMPLE}"

    INPUT="${TRIMMED_DIR}/${SAMPLE}_trimmed.fastq.gz"
    [[ -f "${INPUT}" ]] || { echo "[ERROR] Trimmed FASTQ missing: ${INPUT}" >&2; exit 1; }

    # --- Initialize stats ---
    TRIMMED_READS=""; RRNA_ALIGNED=0; RRNA_PCT=0
    DM6_ALIGNED=0; DM6_PCT=0
    MM39_INPUT=""; MM39_TOTAL=0; MM39_PCT=0
    MM39_ALIGNED_1=0; MM39_ALIGNED_MULTI=0

    # =========================================================================
    # STEP A — Filter rRNA reads
    # =========================================================================
    echo "[STEP A] Filtering rRNA reads..."

    RRNA_UNMAPPED="${BAM_RRNA_DIR}/${SAMPLE}_rRNA_unmapped.fastq.gz"
    RRNA_LOG="${LOG_DIR}/${SAMPLE}_rRNA_bowtie2.log"

    if [[ -f "${RRNA_UNMAPPED}" ]]; then
        echo "[INFO] rRNA-filtered reads exist. Skipping."
    else
        # NOTE: Switched from --very-fast (previous version) to default presets.
        # For a filter we want sensitivity; missing rRNA reads here means they
        # leak through to mm39 and contaminate the signal.
        bowtie2 \
            --no-unal \
            --un-gz "${RRNA_UNMAPPED}" \
            -x "${RRNA_BT2_IDX}" \
            -U "${INPUT}" \
            -p "${THREADS}" \
            2> "${RRNA_LOG}" \
        | samtools view -bS -@ 2 - > /dev/null   # discard the rRNA BAM
    fi

    TRIMMED_READS=$(grep "reads; of these:" "${RRNA_LOG}" | awk '{print $1}' | head -1)
    RRNA_1=$(grep "aligned exactly 1 time" "${RRNA_LOG}" | awk '{print $1}' | head -1)
    RRNA_M=$(grep "aligned >1 times" "${RRNA_LOG}" | awk '{print $1}' | head -1)
    RRNA_ALIGNED=$(( ${RRNA_1:-0} + ${RRNA_M:-0} ))
    RRNA_PCT=$(awk "BEGIN {printf \"%.2f\", ${RRNA_ALIGNED}/${TRIMMED_READS:-1}*100}")
    POST_RRNA=$(( ${TRIMMED_READS:-0} - ${RRNA_ALIGNED} ))

    echo "[SANITY CHECK] rRNA filter:"
    echo "  Input:           ${TRIMMED_READS}"
    echo "  rRNA reads:      ${RRNA_ALIGNED} (${RRNA_PCT}%)"
    echo "  Surviving:       ${POST_RRNA}"
    echo "[INTERPRETATION] rRNA rates 5-20% typical; >30% high but not unusual"

    # =========================================================================
    # STEP B — Filter dm6 spike-in reads (conditional)
    # =========================================================================
    if [[ "${DM6_SPIKEIN_FILTER}" == "true" ]]; then
        echo ""
        echo "[STEP B] Filtering dm6 spike-in reads..."

        DM6_UNMAPPED="${BAM_DM6_DIR}/${SAMPLE}_dm6_unmapped.fastq.gz"
        DM6_LOG="${LOG_DIR}/${SAMPLE}_dm6_bowtie2.log"

        if [[ -f "${DM6_UNMAPPED}" ]]; then
            echo "[INFO] dm6-filtered reads exist. Skipping."
        else
            bowtie2 \
                --very-sensitive \
                --no-unal \
                --un-gz "${DM6_UNMAPPED}" \
                -x "${DM6_BT2_IDX}" \
                -U "${RRNA_UNMAPPED}" \
                -p "${THREADS}" \
                2> "${DM6_LOG}" \
            | samtools view -bS -@ 2 - > /dev/null   # discard the dm6 BAM too
        fi

        DM6_INPUT=$(grep "reads; of these:" "${DM6_LOG}" | awk '{print $1}' | head -1)
        DM6_1=$(grep "aligned exactly 1 time" "${DM6_LOG}" | awk '{print $1}' | head -1)
        DM6_M=$(grep "aligned >1 times" "${DM6_LOG}" | awk '{print $1}' | head -1)
        DM6_ALIGNED=$(( ${DM6_1:-0} + ${DM6_M:-0} ))
        DM6_PCT=$(awk "BEGIN {printf \"%.2f\", ${DM6_ALIGNED}/${DM6_INPUT:-1}*100}")
        POST_DM6=$(( ${DM6_INPUT:-0} - ${DM6_ALIGNED} ))

        echo "[SANITY CHECK] dm6 spike-in:"
        echo "  Input:           ${DM6_INPUT}"
        echo "  dm6 reads:       ${DM6_ALIGNED} (${DM6_PCT}%)"
        echo "  Passed to mm39:  ${POST_DM6}"
        echo "[INTERPRETATION]"
        echo "  1-10%:  Typical when spike-in is present"
        echo "  <0.5%:  No spike-in present in this dataset; consider"
        echo "          DM6_SPIKEIN_FILTER=false in config to skip this step."

        MM39_INPUT_FQ="${DM6_UNMAPPED}"
    else
        echo ""
        echo "[STEP B] dm6 filter SKIPPED (DM6_SPIKEIN_FILTER=false)"
        MM39_INPUT_FQ="${RRNA_UNMAPPED}"
    fi

    # =========================================================================
    # STEP C — Align to mm39
    # =========================================================================
    echo ""
    echo "[STEP C] Aligning to mm39..."

    MM39_RAW_BAM="${BAM_MM39_DIR}/${SAMPLE}_mm39_raw.bam"
    MM39_LOG="${LOG_DIR}/${SAMPLE}_mm39_bowtie2.log"

    if [[ -f "${MM39_RAW_BAM}" && -f "${MM39_RAW_BAM}.bai" ]]; then
        echo "[INFO] mm39 raw BAM exists. Skipping."
    else
        bowtie2 \
            --very-sensitive \
            --no-unal \
            -x "${MM39_BT2_IDX}" \
            -U "${MM39_INPUT_FQ}" \
            -p "${THREADS}" \
            2> "${MM39_LOG}" \
        | samtools view -bS -@ 2 - \
        | samtools sort -@ "${THREADS}" -o "${MM39_RAW_BAM}"

        samtools index "${MM39_RAW_BAM}"
    fi

    MM39_INPUT=$(grep "reads; of these:" "${MM39_LOG}" | awk '{print $1}' | head -1)
    MM39_ALIGNED_1=$(grep "aligned exactly 1 time" "${MM39_LOG}" | awk '{print $1}' | head -1)
    MM39_ALIGNED_MULTI=$(grep "aligned >1 times" "${MM39_LOG}" | awk '{print $1}' | head -1)
    MM39_TOTAL=$(( ${MM39_ALIGNED_1:-0} + ${MM39_ALIGNED_MULTI:-0} ))
    MM39_PCT=$(awk "BEGIN {printf \"%.2f\", ${MM39_TOTAL}/${MM39_INPUT:-1}*100}")

    echo "[SANITY CHECK] mm39 alignment:"
    echo "  Input:        ${MM39_INPUT}"
    echo "  Aligned:      ${MM39_TOTAL} (${MM39_PCT}%)"
    echo "  Unique:       ${MM39_ALIGNED_1:-0}"
    echo "  Multi:        ${MM39_ALIGNED_MULTI:-0}"
    echo "[INTERPRETATION] PRO-seq: >70% good, 50-70% acceptable, <50% investigate"

    # =========================================================================
    # STEP D — Filter: MAPQ + chrM removal
    # =========================================================================
    echo ""
    echo "[STEP D] Filter: MAPQ>=${MAPQ_THRESHOLD}, remove chrM..."

    FINAL_BAM="${BAM_MM39_DIR}/${SAMPLE}_final.bam"

    if [[ -f "${FINAL_BAM}" && -f "${FINAL_BAM}.bai" ]]; then
        echo "[INFO] Final BAM exists. Skipping filter."
    else
        # -F 2820 = exclude unmapped (4) + secondary (256) + QC fail (512) + supplementary (2048)
        # -q THR  = MAPQ threshold
        # awk     = drop chrM (header passes through)
        samtools view -b -q "${MAPQ_THRESHOLD}" -F 2820 "${MM39_RAW_BAM}" \
        | samtools view -h - \
        | awk '$1 ~ /^@/ || $3 != "chrM"' \
        | samtools view -bS - \
        | samtools sort -@ "${THREADS}" -o "${FINAL_BAM}"

        samtools index "${FINAL_BAM}"
    fi

    # =========================================================================
    # Final stats
    # =========================================================================
    FINAL_COUNT=$(samtools view -c "${FINAL_BAM}")
    CHRM_READS=$(samtools view -c "${MM39_RAW_BAM}" chrM 2>/dev/null || echo "0")
    PLUS=$(samtools view -c -F 16 "${FINAL_BAM}")
    MINUS=$(samtools view -c -f 16 "${FINAL_BAM}")
    if (( MINUS > 0 )); then
        STRAND_RATIO=$(awk "BEGIN {printf \"%.3f\", ${PLUS}/${MINUS}}")
    else
        STRAND_RATIO="NA"
    fi

    echo ""
    echo "[SANITY CHECK] Final BAM: ${SAMPLE}"
    printf "  Final reads:    %12d\n" "${FINAL_COUNT}"
    printf "  chrM removed:   %12d\n" "${CHRM_READS}"
    printf "  + strand:       %12d\n" "${PLUS}"
    printf "  - strand:       %12d\n" "${MINUS}"
    echo  "  Strand ratio:   ${STRAND_RATIO} (expect ~1.0 for PRO-seq)"

    # Per-sample QC TSV
    cat > "${QC_DIR}/alignment_stats_${SAMPLE}.tsv" << EOF
sample	trimmed_reads	rRNA_reads	rRNA_pct	dm6_reads	dm6_pct	mm39_input	mm39_aligned	mm39_pct	mm39_unique	mm39_multi	chrM_removed	final_reads	plus_strand	minus_strand	strand_ratio
${SAMPLE}	${TRIMMED_READS}	${RRNA_ALIGNED}	${RRNA_PCT}	${DM6_ALIGNED}	${DM6_PCT}	${MM39_INPUT}	${MM39_TOTAL}	${MM39_PCT}	${MM39_ALIGNED_1:-0}	${MM39_ALIGNED_MULTI:-0}	${CHRM_READS}	${FINAL_COUNT}	${PLUS}	${MINUS}	${STRAND_RATIO}
EOF

    echo "[INFO] Per-sample QC: ${QC_DIR}/alignment_stats_${SAMPLE}.tsv"
done

step_header "Step 03 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output:"
echo "  Final BAMs:     ${BAM_MM39_DIR}/*_final.bam"
echo "  Stats:          ${QC_DIR}/alignment_stats_*.tsv"
echo ""
echo "Next: bash 04_make_bigwigs.sh"
