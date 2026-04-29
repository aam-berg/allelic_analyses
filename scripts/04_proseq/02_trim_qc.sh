#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04_proseq/02_trim_qc.sh — Adapter trimming + post-trim FastQC (SE)
# =============================================================================
#
# WHAT THIS DOES:
#   1. Removes 3' adapter (Illumina TruSeq Small RNA) with cutadapt
#   2. Trims low-quality 3' ends
#   3. Discards reads shorter than MIN_READ_LENGTH after trimming
#   4. Runs FastQC on trimmed reads
#
# WHY:
#   PRO-seq fragments are often shorter than the read length, so the sequencer
#   reads through into the 3' adapter on most reads. Adapter trimming is
#   especially critical for PRO-seq — adapter rates of 30-60% are typical
#   and not a quality concern.
#
# PARALLEL EXECUTION:
#   sbatch --array=0-3 --partition=short --time=4:00:00 --mem=8G \
#       --cpus-per-task=8 -o logs/02_%A_%a.out -e logs/02_%A_%a.err \
#       02_trim_qc.sh
# =============================================================================

source "config.sh"

step_header "PRO-seq Pipeline Step 02: Adapter Trimming + QC"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] cutadapt: $(cutadapt --version)"

create_dirs
resolve_samples "${1:-}"

for SAMPLE in "${RUN_SAMPLES[@]}"; do
    step_header "Trimming ${SAMPLE}"

    INPUT="${FASTQ_DIR}/${SAMPLE}.fastq.gz"
    OUTPUT="${TRIMMED_DIR}/${SAMPLE}_trimmed.fastq.gz"
    REPORT="${LOG_DIR}/${SAMPLE}_cutadapt.log"

    if [[ -f "${OUTPUT}" ]]; then
        echo "[INFO] Trimmed FASTQ exists. Skipping."
        continue
    fi
    if [[ ! -f "${INPUT}" ]]; then
        echo "[ERROR] Input FASTQ missing: ${INPUT}" >&2
        exit 1
    fi

    cutadapt \
        -a "${ADAPTER_3PRIME}" \
        -q "${QUALITY_CUTOFF}" \
        -m "${MIN_READ_LENGTH}" \
        --cores "${THREADS}" \
        -o "${OUTPUT}" \
        "${INPUT}" \
        2>&1 | tee "${REPORT}"

    # Parse stats
    RAW=$(grep "Total reads processed:" "${REPORT}" | awk '{gsub(/,/,""); print $NF}')
    ADAPTED=$(grep "Reads with adapters:" "${REPORT}" | awk '{gsub(/,/,""); print $(NF-1)}')
    TOO_SHORT=$(grep "Reads that were too short:" "${REPORT}" | awk '{gsub(/,/,""); print $(NF-1)}')
    WRITTEN=$(grep "Reads written (passing filters):" "${REPORT}" | awk '{gsub(/,/,""); print $(NF-1)}')

    if [[ -n "${RAW}" && "${RAW}" -gt 0 ]]; then
        PCT_ADAPTER=$(awk "BEGIN {printf \"%.1f\", ${ADAPTED:-0}/${RAW}*100}")
        PCT_SURVIVING=$(awk "BEGIN {printf \"%.1f\", ${WRITTEN:-0}/${RAW}*100}")
    else
        PCT_ADAPTER="NA"; PCT_SURVIVING="NA"
    fi

    echo ""
    echo "[SANITY CHECK] ${SAMPLE}:"
    echo "  Raw reads:           ${RAW}"
    echo "  Reads with adapter:  ${ADAPTED:-NA} (${PCT_ADAPTER}%)"
    echo "  Too short:           ${TOO_SHORT:-NA}"
    echo "  Surviving:           ${WRITTEN:-NA} (${PCT_SURVIVING}%)"
    echo "[INTERPRETATION]"
    echo "  Adapter content 30-60%: NORMAL for PRO-seq"
    echo "  Adapter content >70%:   over-fragmentation, often still OK"
    echo "  Adapter content <10%:   unusual, check adapter sequence"

    echo -e "${SAMPLE}\t${RAW}\t${ADAPTED:-NA}\t${TOO_SHORT:-NA}\t${WRITTEN:-NA}\t${PCT_ADAPTER}\t${PCT_SURVIVING}" \
        > "${QC_DIR}/trimming_stats_${SAMPLE}.tsv"

    echo ""
    echo "[INFO] FastQC on trimmed reads..."
    fastqc \
        --outdir "${FASTQC_TRIM_DIR}" \
        --threads "${THREADS}" \
        --quiet \
        "${OUTPUT}"
done

step_header "Step 02 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output:"
echo "  Trimmed:        ${TRIMMED_DIR}/"
echo "  FastQC trimmed: ${FASTQC_TRIM_DIR}/"
echo "  Trim stats:     ${QC_DIR}/trimming_stats_*.tsv"
echo ""
echo "Next: bash 03_align_filter.sh"
