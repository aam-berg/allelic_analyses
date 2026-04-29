#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 02_trim_qc.sh — Adapter trimming (PE) + post-trim FastQC
# =============================================================================
#
# WHAT THIS DOES:
#   1. Removes 3' Nextera adapters from both reads with cutadapt
#   2. Trims low-quality bases from 3' ends
#   3. Discards read pairs where either read becomes shorter than 20 nt
#   4. Runs FastQC on trimmed reads
#
# WHY NEXTERA ADAPTERS:
#   ATAC-seq libraries are made with Tn5 transposase, which inserts the same
#   Nextera Illumina adapters used in Nextera DNA-seq libraries. The 3' end
#   of any read that runs past the insert is the start of the Nextera adapter
#   on the opposite end of the fragment. The 19-mer "CTGTCTCTTATACACATCT" is
#   the conserved core that's identical between R1 and R2 adapter primers.
#
# WHY -m 20 ON BOTH READS:
#   bowtie2 PE struggles with very short reads. 20 nt is the conventional
#   floor for unique mappability in mouse.
#
# PARALLEL EXECUTION:
#   sbatch --array=0-3 --partition=short --time=4:00:00 --mem=8G \
#       --cpus-per-task=8 -o logs/02_%A_%a.out -e logs/02_%A_%a.err \
#       02_trim_qc.sh
# =============================================================================

source "config.sh"

step_header "ATAC Pipeline Step 02: PE Trimming + QC"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] cutadapt: $(cutadapt --version)"
echo "[INFO] fastqc:   $(fastqc --version 2>&1 | head -1)"

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
        echo "[INFO] Trimmed FASTQs exist for ${SAMPLE}. Skipping."
        continue
    fi

    if [[ ! -f "${R1_IN}" || ! -f "${R2_IN}" ]]; then
        echo "[ERROR] Input FASTQ(s) missing for ${SAMPLE}:" >&2
        [[ -f "${R1_IN}" ]] || echo "  Missing: ${R1_IN}" >&2
        [[ -f "${R2_IN}" ]] || echo "  Missing: ${R2_IN}" >&2
        exit 1
    fi

    # cutadapt PE invocation:
    #   -a / -A         : 3' adapter on R1 / R2 (Nextera core)
    #   -q              : quality trim 3' ends BEFORE adapter removal
    #   -m              : minimum length AFTER trimming; pair dropped if either
    #                     read falls below
    #   --cores 0       : use all available cores (cutadapt picks)
    #   -o / -p         : output R1 / R2
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

    # ---- Parse report ----
    RAW_PAIRS=$(grep "Total read pairs processed:" "${REPORT}" \
                | awk '{gsub(/,/,""); print $NF}')
    TOO_SHORT=$(grep "Pairs that were too short:" "${REPORT}" \
                | awk '{gsub(/,/,""); print $(NF-1)}' || echo "0")
    WRITTEN=$(grep "Pairs written (passing filters):" "${REPORT}" \
              | awk '{gsub(/,/,""); print $(NF-1)}')

    if [[ -n "${RAW_PAIRS}" && "${RAW_PAIRS}" -gt 0 ]]; then
        PCT_SURVIVING=$(awk "BEGIN {printf \"%.2f\", ${WRITTEN:-0}/${RAW_PAIRS}*100}")
        PCT_SHORT=$(awk "BEGIN {printf \"%.2f\", ${TOO_SHORT:-0}/${RAW_PAIRS}*100}")
    else
        PCT_SURVIVING="NA"; PCT_SHORT="NA"
    fi

    echo ""
    echo "[SANITY CHECK] ${SAMPLE}:"
    echo "  Raw pairs:        ${RAW_PAIRS}"
    echo "  Too short:        ${TOO_SHORT:-NA} (${PCT_SHORT}%)"
    echo "  Pairs surviving:  ${WRITTEN:-NA} (${PCT_SURVIVING}%)"
    echo ""
    echo "[INTERPRETATION]"
    echo "  Surviving rate >85%: typical for ATAC"
    echo "  Surviving rate <70%: investigate (over-aggressive Tn5? short fragments?)"

    # Per-sample QC (parallel-safe)
    echo -e "${SAMPLE}\t${RAW_PAIRS}\t${TOO_SHORT:-NA}\t${WRITTEN:-NA}\t${PCT_SHORT}\t${PCT_SURVIVING}" \
        > "${QC_DIR}/trimming_stats_${SAMPLE}.tsv"

    # ---- Post-trim FastQC ----
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
echo "  Trimmed:         ${TRIMMED_DIR}/"
echo "  FastQC trimmed:  ${FASTQC_TRIM_DIR}/"
echo "  Trim stats:      ${QC_DIR}/trimming_stats_*.tsv"
echo ""
echo "Next: bash 03_align_filter.sh"
