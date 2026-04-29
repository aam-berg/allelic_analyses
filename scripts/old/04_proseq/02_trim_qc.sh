#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 02_trim_qc.sh — Adapter trimming, quality filtering, and post-trim QC
# =============================================================================
#
# WHAT THIS DOES:
#   1. Removes 3' adapter sequences from reads using cutadapt
#   2. Trims low-quality bases from 3' end
#   3. Discards reads shorter than 20 nt after trimming
#   4. Runs FastQC on trimmed reads
#
# WHY (and how PRO-seq differs from RNA-seq here):
#   - PRO-seq reads come from nascent RNA fragments of variable length. Many
#     fragments are SHORTER than the sequencing read length, so the sequencer
#     reads through the insert and into the 3' adapter. This makes adapter
#     trimming CRITICAL for PRO-seq (even more so than RNA-seq).
#   - We use cutadapt because:
#     (a) It's the field standard for PRO-seq (used in Mahat et al. 2016, and
#         in the original processing of this dataset)
#     (b) It handles partial adapter matches well
#   - The adapter is the Illumina TruSeq Small RNA 3' Adapter
#     (TGGAATTCTCGGGTGCCAAGG), standard for the Mahat PRO-seq protocol.
#   - We do NOT do 5' trimming or UMI extraction because this protocol
#     does not use UMIs (confirmed from GEO metadata).
#   - Min length 20 nt: shorter reads cannot be uniquely mapped to ~2.7 Gb genome.
#
# KEY DIFFERENCE FROM RNA-seq:
#   - RNA-seq: adapter contamination is typically low (<5%) because inserts
#     are usually 200-500 bp.
#   - PRO-seq: adapter contamination is HIGH (often 30-60%) because nascent
#     RNA fragments after base hydrolysis are typically 30-100 nt.
#   - This is expected and NOT a sign of poor library quality.
#
# PARALLEL EXECUTION:
#   sbatch --array=0-3 --partition=short --time=4:00:00 --mem=8G --cpus-per-task=8 \
#       -o logs/02_%A_%a.out -e logs/02_%A_%a.err 02_trim_qc.sh
# =============================================================================

source "config.sh"

step_header "PRO-seq Pipeline Step 02: Adapter Trimming & QC"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] cutadapt version: $(cutadapt --version)"
echo "[INFO] fastqc version: $(fastqc --version 2>&1 | head -1)"

create_dirs
resolve_samples "${1:-}"

for SAMPLE in "${RUN_SAMPLES[@]}"; do
    step_header "Trimming ${SAMPLE}"

    INPUT="${FASTQ_DIR}/${SAMPLE}.fastq.gz"
    OUTPUT="${TRIMMED_DIR}/${SAMPLE}_trimmed.fastq.gz"
    REPORT="${LOG_DIR}/${SAMPLE}_cutadapt.log"

    if [[ -f "${OUTPUT}" ]]; then
        echo "[INFO] Trimmed file already exists: ${OUTPUT}. Skipping."
        continue
    fi

    if [[ ! -f "${INPUT}" ]]; then
        echo "[ERROR] Input file not found: ${INPUT}"
        exit 1
    fi

    # cutadapt parameters explained:
    #   -a ADAPTER : remove 3' adapter (this is where adapter appears for SE reads)
    #   -q 20      : trim bases with Phred < 20 from 3' end BEFORE adapter removal
    #   -m 20      : discard reads shorter than 20 nt after all trimming
    #   --cores    : parallel processing
    #   -o         : output file
    cutadapt \
        -a "${ADAPTER_3PRIME}" \
        -q "${QUALITY_CUTOFF}" \
        -m "${MIN_READ_LENGTH}" \
        --cores "${THREADS}" \
        -o "${OUTPUT}" \
        "${INPUT}" \
        2>&1 | tee "${REPORT}"

    # --- Parse cutadapt report for QC ---
    RAW=$(grep "Total reads processed:" "${REPORT}" | awk '{gsub(/,/,""); print $NF}')
    ADAPTED=$(grep "Reads with adapters:" "${REPORT}" | awk '{gsub(/,/,""); print $(NF-1)}')
    TOO_SHORT=$(grep "Reads that were too short:" "${REPORT}" | awk '{gsub(/,/,""); print $(NF-1)}')
    WRITTEN=$(grep "Reads written (passing filters):" "${REPORT}" | awk '{gsub(/,/,""); print $(NF-1)}')

    if [[ -n "${RAW}" && "${RAW}" -gt 0 ]]; then
        PCT_ADAPTER=$(awk "BEGIN {printf \"%.1f\", ${ADAPTED:-0}/${RAW}*100}")
        PCT_SURVIVING=$(awk "BEGIN {printf \"%.1f\", ${WRITTEN:-0}/${RAW}*100}")
    else
        PCT_ADAPTER="NA"
        PCT_SURVIVING="NA"
    fi

    # Per-sample QC file (safe for parallel writes)
    echo -e "${SAMPLE}\t${RAW}\t${ADAPTED:-NA}\t${TOO_SHORT:-NA}\t${WRITTEN:-NA}\t${PCT_ADAPTER}\t${PCT_SURVIVING}" \
        > "${QC_DIR}/trimming_stats_${SAMPLE}.tsv"

    echo ""
    echo "[SANITY CHECK] ${SAMPLE} trimming summary:"
    echo "  Raw reads:          ${RAW}"
    echo "  Reads with adapter: ${ADAPTED:-NA} (${PCT_ADAPTER}%)"
    echo "  Too short (<${MIN_READ_LENGTH} nt): ${TOO_SHORT:-NA}"
    echo "  Reads surviving:    ${WRITTEN:-NA} (${PCT_SURVIVING}%)"
    echo ""
    echo "[INTERPRETATION]"
    echo "  Adapter content 30-60%: NORMAL for PRO-seq (short nascent RNA fragments)"
    echo "  Adapter content >70%:   Possible over-fragmentation, often still OK"
    echo "  Adapter content <10%:   Unusual — check adapter sequence"
    echo "  Survival rate >60%:     Good"

    # --- FastQC on trimmed reads ---
    echo ""
    echo "[INFO] Running FastQC on trimmed reads..."
    fastqc \
        --outdir "${FASTQC_TRIM_DIR}" \
        --threads "${THREADS}" \
        --quiet \
        "${OUTPUT}"

    echo "[INFO] Post-trim FastQC: ${FASTQC_TRIM_DIR}/${SAMPLE}_trimmed_fastqc.html"
done

step_header "Step 02 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output:"
echo "  Trimmed FASTQ:    ${TRIMMED_DIR}/"
echo "  FastQC (trimmed): ${FASTQC_TRIM_DIR}/"
echo "  Trim stats:       ${QC_DIR}/trimming_stats_*.tsv"
echo ""
echo "Next: Run 03_align_filter.sh (after all step 02 array tasks finish)"