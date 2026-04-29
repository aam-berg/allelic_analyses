#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 01_download_fastq.sh — Download paired-end ATAC FASTQ from SRA + raw QC
# =============================================================================
#
# WHAT THIS DOES:
#   1. Downloads R1 and R2 FASTQ for each sample using fasterq-dump
#   2. Renames the SRR-named output files to sample-named files
#   3. Compresses with pigz
#   4. Runs FastQC on raw reads
#
# WHY KEEP BOTH R1 AND R2 (unlike PRO-seq):
#   ATAC inserts span open chromatin regions; both ends are equally
#   informative about cut-site positions. We use both reads for fragment-
#   level analysis (peak calling, fragment-coverage bigWigs).
#
# PARALLEL EXECUTION:
#   sbatch --array=0-3 --partition=short --time=6:00:00 --mem=16G \
#       --cpus-per-task=8 -o logs/01_%A_%a.out -e logs/01_%A_%a.err \
#       01_download_fastq.sh
#
#   # Single sample by index or name:
#   bash 01_download_fastq.sh 0
#   bash 01_download_fastq.sh mESC_ATAC_rep1
# =============================================================================

source "config.sh"

step_header "ATAC Pipeline Step 01: Download FASTQ + Raw QC"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] fasterq-dump: $(fasterq-dump --version 2>&1 | head -1)"
echo "[INFO] fastqc:       $(fastqc --version 2>&1 | head -1)"

create_dirs
resolve_samples "${1:-}"

for SAMPLE in "${RUN_SAMPLES[@]}"; do
    SRR="${SAMPLES[${SAMPLE}]}"
    step_header "Downloading ${SAMPLE} (${SRR})"

    R1="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
    R2="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"

    # ---- Idempotency check ----
    if [[ -f "${R1}" && -f "${R2}" ]]; then
        echo "[INFO] Both FASTQs already exist for ${SAMPLE}. Skipping download."
    else
        echo "[INFO] Downloading ${SRR} with fasterq-dump..."

        # --split-3: produces _1.fastq, _2.fastq for PE; .fastq for SE; and
        # .fastq for unpaired reads (rare). We expect _1 and _2 only.
        fasterq-dump "${SRR}" \
            --outdir "${FASTQ_DIR}" \
            --split-3 \
            --threads "${THREADS}" \
            --temp "${SCRATCH_DIR}/tmp" \
            --progress

        # ---- Detect paired-end output and rename ----
        if [[ -f "${FASTQ_DIR}/${SRR}_1.fastq" && -f "${FASTQ_DIR}/${SRR}_2.fastq" ]]; then
            echo "[INFO] Paired-end output detected (expected)."
            mv "${FASTQ_DIR}/${SRR}_1.fastq" "${FASTQ_DIR}/${SAMPLE}_R1.fastq"
            mv "${FASTQ_DIR}/${SRR}_2.fastq" "${FASTQ_DIR}/${SAMPLE}_R2.fastq"
            # Discard any unpaired remainder file
            rm -f "${FASTQ_DIR}/${SRR}.fastq"
        else
            echo "[ERROR] Expected paired-end output for ${SRR} but did not find" >&2
            echo "        both ${SRR}_1.fastq and ${SRR}_2.fastq." >&2
            ls -la "${FASTQ_DIR}/${SRR}"* 2>/dev/null >&2 || true
            exit 1
        fi

        # ---- Compress ----
        echo "[INFO] Compressing with pigz..."
        pigz -p "${THREADS}" "${FASTQ_DIR}/${SAMPLE}_R1.fastq" \
                              "${FASTQ_DIR}/${SAMPLE}_R2.fastq"
    fi

    # ---- Sanity check: read counts must match ----
    R1_COUNT=$(($(zcat "${R1}" | wc -l) / 4))
    R2_COUNT=$(($(zcat "${R2}" | wc -l) / 4))

    echo "[SANITY CHECK] ${SAMPLE}:"
    echo "  R1: ${R1_COUNT} reads"
    echo "  R2: ${R2_COUNT} reads"

    if [[ "${R1_COUNT}" -ne "${R2_COUNT}" ]]; then
        echo "[ERROR] R1/R2 read count mismatch for ${SAMPLE}!" >&2
        exit 1
    fi

    # Per-sample QC (parallel-safe — overwrite, not append)
    echo -e "${SAMPLE}\t${SRR}\t${R1_COUNT}" > "${QC_DIR}/read_counts_raw_${SAMPLE}.tsv"

    # ---- FastQC ----
    echo "[INFO] Running FastQC on R1 and R2..."
    fastqc \
        --outdir "${FASTQC_RAW_DIR}" \
        --threads "${THREADS}" \
        --quiet \
        "${R1}" "${R2}"
done

step_header "Step 01 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output:"
echo "  Raw FASTQ:  ${FASTQ_DIR}/"
echo "  FastQC raw: ${FASTQC_RAW_DIR}/"
echo ""
echo "Next: bash 02_trim_qc.sh"
