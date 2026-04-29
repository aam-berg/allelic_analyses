#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04c_rnaseq/01_download_fastq.sh — PE download with multi-SRR concatenation
# =============================================================================
#
# WHAT THIS DOES (per sample):
#   1. For each SRR in the sample's comma-separated list, downloads PE FASTQ
#   2. Concatenates all R1 files together and all R2 files together,
#      producing one {SAMPLE}_R1.fastq.gz + {SAMPLE}_R2.fastq.gz pair
#   3. Runs FastQC on the concatenated fastqs
#
# WHY MULTI-SRR CONCATENATION:
#   GSE200699 splits each biological replicate across multiple SRR accessions
#   (lane-split technical reps). They should be concatenated at the FASTQ
#   level so STAR sees a single sample with full depth.
#
# IDEMPOTENCY:
#   - Skip download if per-SRR temp files exist
#   - Skip concat if final per-sample fastqs exist
#   - Skip FastQC if HTML report exists
#
# PARALLEL EXECUTION:
#   sbatch --array=0-1 --partition=short --time=4:00:00 --mem=16G \
#       --cpus-per-task=8 -o logs/01_%A_%a.out -e logs/01_%A_%a.err \
#       01_download_fastq.sh
# =============================================================================

source "config.sh"

step_header "RNA-seq Pipeline Step 01: Multi-SRR Download + Concatenation"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] fasterq-dump: $(fasterq-dump --version 2>&1 | head -1 || echo 'unknown')"
echo "[INFO] fastqc:       $(fastqc --version 2>&1 | head -1)"

create_dirs
resolve_samples "${1:-}"

for SAMPLE in "${RUN_SAMPLES[@]}"; do
    get_srrs_for_sample "${SAMPLE}"
    SRR_LIST=("${SRR_ARRAY[@]}")

    step_header "Downloading ${SAMPLE} (${#SRR_LIST[@]} SRRs: ${SRR_LIST[*]})"

    R1_FINAL="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
    R2_FINAL="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"

    if [[ -f "${R1_FINAL}" && -f "${R2_FINAL}" ]]; then
        echo "[INFO] Concatenated FASTQs exist for ${SAMPLE}. Skipping download."
    else
        # Per-SRR temp directory
        SAMPLE_TMP="${FASTQ_DIR}/_tmp_${SAMPLE}"
        mkdir -p "${SAMPLE_TMP}"

        # ---- Download each SRR ----
        for SRR in "${SRR_LIST[@]}"; do
            R1_SRR="${SAMPLE_TMP}/${SRR}_1.fastq"
            R2_SRR="${SAMPLE_TMP}/${SRR}_2.fastq"
            R1_SRR_GZ="${R1_SRR}.gz"
            R2_SRR_GZ="${R2_SRR}.gz"

            if [[ -f "${R1_SRR_GZ}" && -f "${R2_SRR_GZ}" ]]; then
                echo "[INFO]   ${SRR}: per-SRR fastqs exist. Skipping."
                continue
            fi

            echo "[INFO]   Downloading ${SRR}..."
            fasterq-dump "${SRR}" \
                --outdir "${SAMPLE_TMP}" \
                --split-3 \
                --threads "${THREADS}" \
                --temp "${SCRATCH_DIR}/tmp" \
                --progress

            if [[ ! -f "${R1_SRR}" || ! -f "${R2_SRR}" ]]; then
                echo "[ERROR] Expected PE output for ${SRR}; got:" >&2
                ls -la "${SAMPLE_TMP}/${SRR}"* 2>/dev/null >&2 || true
                exit 1
            fi

            # Discard any unpaired remainder
            rm -f "${SAMPLE_TMP}/${SRR}.fastq"

            echo "[INFO]   Compressing ${SRR}..."
            pigz -p "${THREADS}" "${R1_SRR}" "${R2_SRR}"
        done

        # ---- Concatenate ----
        echo ""
        echo "[INFO] Concatenating ${#SRR_LIST[@]} SRRs into ${SAMPLE} R1+R2..."

        # Build ordered file lists
        R1_FILES=()
        R2_FILES=()
        for SRR in "${SRR_LIST[@]}"; do
            R1_FILES+=("${SAMPLE_TMP}/${SRR}_1.fastq.gz")
            R2_FILES+=("${SAMPLE_TMP}/${SRR}_2.fastq.gz")
        done

        cat "${R1_FILES[@]}" > "${R1_FINAL}"
        cat "${R2_FILES[@]}" > "${R2_FINAL}"

        echo "[INFO] Concatenated:"
        echo "   ${R1_FINAL}  ($(du -h "${R1_FINAL}" | cut -f1))"
        echo "   ${R2_FINAL}  ($(du -h "${R2_FINAL}" | cut -f1))"

        # Clean up per-SRR temp files (keep some space)
        rm -rf "${SAMPLE_TMP}"
    fi

    # ---- Sanity check ----
    R1_COUNT=$(($(zcat "${R1_FINAL}" | wc -l) / 4))
    R2_COUNT=$(($(zcat "${R2_FINAL}" | wc -l) / 4))

    echo ""
    echo "[SANITY CHECK] ${SAMPLE}:"
    echo "  R1 reads: ${R1_COUNT}"
    echo "  R2 reads: ${R2_COUNT}"

    if [[ "${R1_COUNT}" -ne "${R2_COUNT}" ]]; then
        echo "[ERROR] R1/R2 read count mismatch for ${SAMPLE}!" >&2
        exit 1
    fi

    # Per-sample QC TSV (overwrite)
    SRR_LIST_STR=$(IFS=,; echo "${SRR_LIST[*]}")
    echo -e "${SAMPLE}\t${SRR_LIST_STR}\t${R1_COUNT}" \
        > "${QC_DIR}/read_counts_raw_${SAMPLE}.tsv"

    # ---- FastQC ----
    echo "[INFO] Running FastQC..."
    fastqc \
        --outdir "${FASTQC_RAW_DIR}" \
        --threads "${THREADS}" \
        --quiet \
        "${R1_FINAL}" "${R2_FINAL}"
done

step_header "Step 01 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output:"
echo "  Concatenated FASTQs: ${FASTQ_DIR}/{SAMPLE}_R{1,2}.fastq.gz"
echo "  FastQC raw:          ${FASTQC_RAW_DIR}/"
echo ""
echo "Next: bash 02_trim_qc.sh"
