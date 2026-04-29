#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 01_download_fastq.sh — Download raw PRO-seq FASTQ from SRA + initial QC
# =============================================================================
#
# WHAT THIS DOES:
#   1. Downloads raw FASTQ files from SRA using fasterq-dump
#   2. Detects whether data is single-end (SE) or paired-end (PE)
#      - If PE, keeps only Read 1 (see WHY below)
#   3. Runs FastQC on the raw reads for initial quality assessment
#
# WHY:
#   - We start from raw reads so we can align to mm39 and control every step.
#   - PRO-seq libraries are sequenced from cDNA (reverse complement of nascent
#     RNA). In PE mode, Read 1 contains the informative end (3' end of nascent
#     RNA = active transcription site). Read 2 is redundant for standard
#     PRO-seq analysis, so we discard it.
#   - FastQC gives us adapter contamination rates, quality scores, and read
#     length distributions BEFORE trimming, so we can assess raw data quality.
#
# PARALLEL EXECUTION:
#   # As a SLURM array (recommended — runs all 4 samples in parallel):
#   sbatch --array=0-3 --partition=short --time=6:00:00 --mem=16G --cpus-per-task=8 \
#       -o logs/01_%A_%a.out -e logs/01_%A_%a.err 01_download_fastq.sh
#
#   # Single sample by index:
#   bash 01_download_fastq.sh 2        # processes sample index 2 (rep3)
#
#   # Single sample by name:
#   bash 01_download_fastq.sh WT_PROseq_rep1
#
#   # All samples sequentially (no argument):
#   bash 01_download_fastq.sh
# =============================================================================

source "config.sh"

step_header "PRO-seq Pipeline Step 01: Download FASTQ and Raw QC"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] fasterq-dump version: $(fasterq-dump --version 2>&1 | head -1 || echo 'unknown')"
echo "[INFO] fastqc version: $(fastqc --version 2>&1 | head -1)"

create_dirs

# --- Resolve which sample(s) to process ---
resolve_samples "${1:-}"

# =============================================================================
# Process each sample
# =============================================================================
for SAMPLE in "${RUN_SAMPLES[@]}"; do
    SRR="${SAMPLES[$SAMPLE]}"
    step_header "Downloading ${SAMPLE} (${SRR})"

    # --- Check if already downloaded ---
    if [[ -f "${FASTQ_DIR}/${SAMPLE}.fastq.gz" ]]; then
        echo "[INFO] ${SAMPLE}.fastq.gz already exists. Skipping download."
    else
        echo "[INFO] Downloading ${SRR} with fasterq-dump..."
        # fasterq-dump is the modern replacement for fastq-dump; much faster.
        # --split-3: if PE, creates _1.fastq, _2.fastq, and unpaired.fastq
        #            if SE, creates just one .fastq file
        fasterq-dump "${SRR}" \
            --outdir "${FASTQ_DIR}" \
            --split-3 \
            --threads "${THREADS}" \
            --temp "${SCRATCH_DIR}" \
            --progress

        # --- Detect SE vs PE and handle accordingly ---
        if [[ -f "${FASTQ_DIR}/${SRR}_1.fastq" && -f "${FASTQ_DIR}/${SRR}_2.fastq" ]]; then
            echo "[INFO] Paired-end data detected for ${SRR}."
            echo "[INFO] For PRO-seq, we use only Read 1 (contains 3' end of nascent RNA)."
            echo "[INFO] Read 2 is discarded (it reads into the 5' adapter / less informative)."
            mv "${FASTQ_DIR}/${SRR}_1.fastq" "${FASTQ_DIR}/${SAMPLE}.fastq"
            rm -f "${FASTQ_DIR}/${SRR}_2.fastq" "${FASTQ_DIR}/${SRR}.fastq"
        elif [[ -f "${FASTQ_DIR}/${SRR}.fastq" ]]; then
            echo "[INFO] Single-end data detected for ${SRR}."
            mv "${FASTQ_DIR}/${SRR}.fastq" "${FASTQ_DIR}/${SAMPLE}.fastq"
        else
            echo "[ERROR] Unexpected output from fasterq-dump for ${SRR}. Check files:"
            ls -la "${FASTQ_DIR}/${SRR}"* 2>/dev/null || echo "  No files found!"
            exit 1
        fi

        # Compress to save space
        echo "[INFO] Compressing ${SAMPLE}.fastq..."
        pigz -p "${THREADS}" "${FASTQ_DIR}/${SAMPLE}.fastq"
        echo "[INFO] Done: ${FASTQ_DIR}/${SAMPLE}.fastq.gz"
    fi

    # --- Sanity check: read count ---
    READ_COUNT=$(zcat "${FASTQ_DIR}/${SAMPLE}.fastq.gz" | wc -l)
    READ_COUNT=$((READ_COUNT / 4))
    echo "[SANITY CHECK] ${SAMPLE}: ${READ_COUNT} reads"
    # Thread-safe append (one line per sample, flock for array safety)
    echo -e "${SAMPLE}\t${SRR}\t${READ_COUNT}" >> "${QC_DIR}/read_counts_raw_${SAMPLE}.tsv"

    # --- FastQC ---
    echo "[INFO] Running FastQC on ${SAMPLE}..."
    fastqc \
        --outdir "${FASTQC_RAW_DIR}" \
        --threads "${THREADS}" \
        --quiet \
        "${FASTQ_DIR}/${SAMPLE}.fastq.gz"

    echo "[INFO] FastQC report: ${FASTQC_RAW_DIR}/${SAMPLE}_fastqc.html"
    echo ""
    echo "[INFO] Check FastQC for:"
    echo "  - Per-base quality scores (should be >20 across most of the read)"
    echo "  - Adapter content (expect HIGH adapter contamination — normal for PRO-seq"
    echo "    because nascent RNA fragments can be shorter than the read length)"
    echo "  - Sequence length distribution"
done

# =============================================================================
# Done
# =============================================================================
step_header "Step 01 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output:"
echo "  Raw FASTQ:   ${FASTQ_DIR}/"
echo "  FastQC:      ${FASTQC_RAW_DIR}/"
echo ""
echo "Next: Run 02_trim_qc.sh (after all step 01 array tasks finish)"