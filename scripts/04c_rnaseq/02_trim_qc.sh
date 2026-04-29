#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04c_rnaseq/02_trim_qc.sh — PE TruSeq adapter trimming + post-trim FastQC
# =============================================================================
#
# WHAT THIS DOES:
#   1. Removes 3' Illumina TruSeq adapters from both reads with cutadapt
#   2. Trims low-quality 3' bases
#   3. Discards pairs where either read becomes < MIN_READ_LENGTH (20)
#   4. Runs FastQC on trimmed pairs
#
# WHY TRUSEQ ADAPTERS:
#   Library was prepared with Illumina TruSeq Stranded Total RNA kit (per
#   GSE200699 GEO record). The standard 33-mer adapters are anchored on the
#   conserved 13-mer cores, but cutadapt handles partial matches gracefully.
#
# PARALLEL EXECUTION:
#   sbatch --array=0-1 --partition=short --time=2:00:00 --mem=8G \
#       --cpus-per-task=8 -o logs/02_%A_%a.out -e logs/02_%A_%a.err \
#       02_trim_qc.sh
# =============================================================================

source "config.sh"

step_header "RNA-seq Pipeline Step 02: PE Trimming + QC"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] cutadapt: $(cutadapt --version)"

create_dirs
resolve_samples "${1:-}"

for SAMPLE in "${RUN_SAMPLES[@]}"; do
    step_header "Trimming ${SAMPLE}"

    R1_IN="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
    R2_IN="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"
    R1_OUT="${TRIMMED_DIR}/${SAMPLE}_R1_trimmed.fastq.gz"
    R2_OUT="${TRIMMED_DIR}/${SAMPLE}_R2_trimmed.fastq.gz"
    REPORT="${LOG_DIR}/${SAMPLE}_cutadapt.log"

    if [[ -f "${R1_OUT}" && -f "${R2_OUT}" ]]; then
        echo "[INFO] Trimmed FASTQs exist. Skipping."
        continue
    fi
    if [[ ! -f "${R1_IN}" || ! -f "${R2_IN}" ]]; then
        echo "[ERROR] Input FASTQ(s) missing for ${SAMPLE}." >&2
        exit 1
    fi

    cutadapt \
        -a "${ADAPTER_R1}" \
        -A "${ADAPTER_R2}" \
        -q "${QUALITY_CUTOFF}" \
        -m "${MIN_READ_LENGTH}" \
        --cores "${THREADS}" \
        -o "${R1_OUT}" \
        -p "${R2_OUT}" \
        "${R1_IN}" "${R2_IN}" \
        2>&1 | tee "${REPORT}"

    # Parse stats
    RAW_PAIRS=$(grep "Total read pairs processed:" "${REPORT}" \
                | awk '{gsub(/,/,""); print $NF}')
    TOO_SHORT=$(grep "Pairs that were too short:" "${REPORT}" \
                | awk '{gsub(/,/,""); print $(NF-1)}' || echo "0")
    WRITTEN=$(grep "Pairs written (passing filters):" "${REPORT}" \
              | awk '{gsub(/,/,""); print $(NF-1)}')

    if [[ -n "${RAW_PAIRS}" && "${RAW_PAIRS}" -gt 0 ]]; then
        PCT_SURV=$(awk "BEGIN {printf \"%.2f\", ${WRITTEN:-0}/${RAW_PAIRS}*100}")
        PCT_SHORT=$(awk "BEGIN {printf \"%.2f\", ${TOO_SHORT:-0}/${RAW_PAIRS}*100}")
    else
        PCT_SURV="NA"; PCT_SHORT="NA"
    fi

    echo ""
    echo "[SANITY CHECK] ${SAMPLE}:"
    echo "  Raw pairs:        ${RAW_PAIRS}"
    echo "  Too short:        ${TOO_SHORT:-NA} (${PCT_SHORT}%)"
    echo "  Pairs surviving:  ${WRITTEN:-NA} (${PCT_SURV}%)"
    echo "[INTERPRETATION]"
    echo "  Surviving rate >85%: typical for stranded total RNA"
    echo "  Surviving rate <70%: investigate (over-aggressive trim? too-short inserts?)"

    echo -e "${SAMPLE}\t${RAW_PAIRS}\t${TOO_SHORT:-NA}\t${WRITTEN:-NA}\t${PCT_SHORT}\t${PCT_SURV}" \
        > "${QC_DIR}/trimming_stats_${SAMPLE}.tsv"

    echo ""
    echo "[INFO] Running FastQC on trimmed reads..."
    fastqc \
        --outdir "${FASTQC_TRIM_DIR}" \
        --threads "${THREADS}" \
        --quiet \
        "${R1_OUT}" "${R2_OUT}"
done

step_header "Step 02 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output:"
echo "  Trimmed:        ${TRIMMED_DIR}/"
echo "  FastQC trimmed: ${FASTQC_TRIM_DIR}/"
echo ""
echo "Next: bash 03_align_star.sh"
