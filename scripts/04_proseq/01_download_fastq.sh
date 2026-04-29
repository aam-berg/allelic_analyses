#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04_proseq/01_download_fastq.sh — Download PRO-seq FASTQ from SRA + raw QC
# =============================================================================
#
# WHAT THIS DOES:
#   1. Downloads raw FASTQ for each sample using fasterq-dump
#   2. If PE data: keep only Read 1 (PRO-seq's informative end)
#   3. Compresses with pigz
#   4. Runs FastQC on the raw reads
#
# WHY KEEP ONLY R1 FOR PE PRO-seq:
#   PRO-seq libraries are sequenced from cDNA (reverse complement of nascent
#   RNA). In PE mode, R1 contains the 3' end of nascent RNA at the read's 5'
#   end — which is the RNAPII active site. R2 reads from the opposite end of
#   the cDNA, which is the 5' end of the original RNA fragment, and is much
#   less informative for pause analysis. So we discard R2.
#
# CHANGES FROM PREVIOUS VERSION:
#   - Per-sample QC TSV uses '>' (overwrite) instead of '>>' (append).
#     The append behavior caused duplicate lines on re-runs.
#   - Sample loading now via samples.tsv (in config.sh).
#   - lib/ utilities are not used here; this is a pipeline-specific download.
#
# PARALLEL EXECUTION:
#   sbatch --array=0-3 --partition=short --time=6:00:00 --mem=16G \
#       --cpus-per-task=8 -o logs/01_%A_%a.out -e logs/01_%A_%a.err \
#       01_download_fastq.sh
# =============================================================================

source "config.sh"

step_header "PRO-seq Pipeline Step 01: Download FASTQ + Raw QC"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] fasterq-dump: $(fasterq-dump --version 2>&1 | head -1 || echo 'unknown')"
echo "[INFO] fastqc:       $(fastqc --version 2>&1 | head -1)"

create_dirs
resolve_samples "${1:-}"

for SAMPLE in "${RUN_SAMPLES[@]}"; do
    SRR="${SAMPLES[${SAMPLE}]}"
    step_header "Downloading ${SAMPLE} (${SRR})"

    OUT_FQ="${FASTQ_DIR}/${SAMPLE}.fastq.gz"

    # ---- Idempotency check ----
    if [[ -f "${OUT_FQ}" ]]; then
        echo "[INFO] FASTQ already exists for ${SAMPLE}. Skipping download."
    else
        echo "[INFO] Downloading ${SRR} with fasterq-dump..."
        fasterq-dump "${SRR}" \
            --outdir "${FASTQ_DIR}" \
            --split-3 \
            --threads "${THREADS}" \
            --temp "${SCRATCH_DIR}/tmp" \
            --progress

        # ---- Detect SE vs PE output ----
        if [[ -f "${FASTQ_DIR}/${SRR}_1.fastq" && -f "${FASTQ_DIR}/${SRR}_2.fastq" ]]; then
            echo "[INFO] Paired-end output detected for ${SRR}."
            echo "[INFO] Keeping R1 only (informative for PRO-seq); discarding R2."
            mv "${FASTQ_DIR}/${SRR}_1.fastq" "${FASTQ_DIR}/${SAMPLE}.fastq"
            rm -f "${FASTQ_DIR}/${SRR}_2.fastq" "${FASTQ_DIR}/${SRR}.fastq"
        elif [[ -f "${FASTQ_DIR}/${SRR}.fastq" ]]; then
            echo "[INFO] Single-end output detected for ${SRR}."
            mv "${FASTQ_DIR}/${SRR}.fastq" "${FASTQ_DIR}/${SAMPLE}.fastq"
        else
            echo "[ERROR] Unexpected output from fasterq-dump for ${SRR}:" >&2
            ls -la "${FASTQ_DIR}/${SRR}"* 2>/dev/null >&2 || true
            exit 1
        fi

        echo "[INFO] Compressing..."
        pigz -p "${THREADS}" "${FASTQ_DIR}/${SAMPLE}.fastq"
    fi

    # ---- Sanity check: read count ----
    READ_COUNT=$(($(zcat "${OUT_FQ}" | wc -l) / 4))
    echo "[SANITY CHECK] ${SAMPLE}: ${READ_COUNT} reads"

    # Per-sample QC: use '>' (overwrite) — fix from old '>>' (append) bug
    echo -e "${SAMPLE}\t${SRR}\t${READ_COUNT}" \
        > "${QC_DIR}/read_counts_raw_${SAMPLE}.tsv"

    # ---- FastQC ----
    echo "[INFO] Running FastQC..."
    fastqc \
        --outdir "${FASTQC_RAW_DIR}" \
        --threads "${THREADS}" \
        --quiet \
        "${OUT_FQ}"

    echo ""
    echo "[INFO] FastQC report: ${FASTQC_RAW_DIR}/${SAMPLE}_fastqc.html"
    echo "[INFO] Expect HIGH adapter contamination (30-60%) for PRO-seq."
done

step_header "Step 01 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output:"
echo "  Raw FASTQ:   ${FASTQ_DIR}/"
echo "  FastQC raw:  ${FASTQC_RAW_DIR}/"
echo ""
echo "Next: bash 02_trim_qc.sh"
