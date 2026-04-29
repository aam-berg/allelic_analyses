#!/usr/bin/env bash
# =============================================================================
# 01_download_and_process.sh — Download, trim, align, filter, and deduplicate
# =============================================================================
#
# OVERVIEW OF WHAT THIS SCRIPT DOES (and why):
#
#   For each PRO-seq replicate:
#
#   STEP A: DOWNLOAD raw FASTQ from SRA
#     Tool: fasterq-dump (SRA Toolkit)
#     Why:  Get the raw sequencing reads.
#
#   STEP B: PRE-TRIM QC
#     Tool: FastQC
#     Why:  See raw read quality, adapter contamination, length distribution
#           before we do anything. Good for debugging.
#
#   STEP C: ADAPTER TRIMMING & QUALITY FILTERING
#     Tool: cutadapt
#     Why:  PRO-seq libraries have a 3' adapter (TruSeq Small RNA adapter:
#           TGGAATTCTCGGGTGCCAAGG) ligated to the 3' end of the nascent RNA.
#           After RT and sequencing, read-through into this adapter will appear
#           at the 3' end of sequencing reads. We must remove it for accurate
#           alignment. We also trim low-quality 3' bases (Q<20) and discard
#           reads shorter than 20 nt (too short for unique mapping).
#
#   STEP D: POST-TRIM QC
#     Tool: FastQC
#     Why:  Verify trimming worked — adapter content should be gone, quality
#           should be uniformly good.
#
#   STEP E: SPIKE-IN FILTERING (align to dm6)
#     Tool: bowtie2
#     Why:  The original experiment spiked in Drosophila cells for normalization.
#           We align to dm6 first and KEEP ONLY THE UNMAPPED READS (i.e., the
#           reads that are NOT Drosophila). This ensures our mouse signal is clean.
#           We also report the spike-in percentage as a QC metric.
#     Settings:
#       --very-sensitive: thorough alignment to catch all spike-in reads
#       --no-mixed --no-discordant: (for PE, if applicable)
#       --un: output unmapped reads (these are our mouse reads)
#
#   STEP F: ALIGN TO mm39
#     Tool: bowtie2
#     Why:  Map the spike-in-filtered reads to the mouse mm39 genome.
#     Settings:
#       --sensitive: good balance of speed and accuracy for typical PRO-seq reads
#       --no-unal: don't output unmapped reads (saves space)
#       --threads: parallelize
#
#   STEP G: POST-ALIGNMENT FILTERING
#     Tool: samtools
#     Why:  Keep only high-quality, uniquely-mapped reads.
#       - Remove unmapped reads (flag 4)
#       - Remove reads with MAPQ < 10 (multimappers/ambiguous alignments)
#       - Sort by coordinate (required for downstream tools)
#       - Index the BAM
#
#   STEP H: PCR DUPLICATE REMOVAL
#     Tool: samtools markdup
#     Why:  PRO-seq libraries undergo PCR amplification (9-11 cycles per the
#           protocol). PCR duplicates inflate read counts at certain positions
#           and must be removed for accurate quantification.
#     Note: For SE reads, samtools markdup identifies duplicates based on
#           alignment position and strand, which is appropriate here.
#
# RUNTIME: ~1-2 hours per sample (depending on download speed and read count)
# MEMORY:  ~8-16 GB peak (bowtie2 alignment)
#
# USAGE:
#   bash 01_download_and_process.sh                  # Process all samples
#   bash 01_download_and_process.sh SRR17640044      # Process one sample
#
# Designed to work with submit_01.sh for SLURM job array parallelism.
# =============================================================================

# ---------------------------------------------------------------------------
# Error handling: we use -eo pipefail but NOT -u.
#   -e          : exit on any command failure
#   -o pipefail : catch failures inside pipes (e.g. bowtie2 | samtools)
#   NOT -u      : unbound variable errors silently kill conda activate,
#                 associative array lookups, and ${1:-} patterns.
# ---------------------------------------------------------------------------
set -eo pipefail

source "config.sh"

echo "============================================================"
echo "PRO-seq Pipeline: Download & Process"
echo "Started: $(date)"
echo "============================================================"

# --- Activate conda environment ---
# conda activate can trip up set -e / set -u depending on conda version,
# so we wrap it:
set +e
source activate "proseq_wasp" 2>/dev/null || conda activate "proseq_wasp"
set -e

echo "[INFO] Conda env: ${CONDA_ENV_NAME}"
echo "[INFO] cutadapt version: $(cutadapt --version)"
echo "[INFO] bowtie2 version: $(bowtie2 --version 2>&1 | head -1)"
echo "[INFO] samtools version: $(samtools --version | head -1)"

# --- Create directories ---
mkdir -p "${FASTQ_DIR}" "${TRIMMED_DIR}" "${ALIGNED_DIR}" "${QC_DIR}"
mkdir -p "${QC_DIR}/fastqc_raw" "${QC_DIR}/fastqc_trimmed" "${QC_DIR}/logs"

# --- Temp directory: use local SSD if available (much faster I/O) ---
SCRATCH_TMPDIR="${TMPDIR:-/tmp}/${USER}_proseq_$$"
mkdir -p "${SCRATCH_TMPDIR}"
trap 'rm -rf "${SCRATCH_TMPDIR}"' EXIT
echo "[INFO] Temp directory: ${SCRATCH_TMPDIR}"

# --- Determine which samples to process ---
# If a specific SRR is given as argument, process only that sample
FILTER_SRR="${1:-}"

# --- Helper: fast read counting ---
# Counting reads via pigz + wc -l is significantly faster than awk on large files.
count_reads_fastq() {
    # Usage: count_reads_fastq file.fastq.gz  OR  count_reads_fastq file.fastq
    local f="$1"
    local lines
    if [[ "${f}" == *.gz ]]; then
        lines=$(pigz -dc "${f}" | wc -l)
    else
        lines=$(wc -l < "${f}")
    fi
    echo $(( lines / 4 ))
}

# =============================================================================
# MAIN PROCESSING LOOP
# =============================================================================
for SAMPLE_NAME in "${!SAMPLES[@]}"; do
    SRR="${SAMPLES[${SAMPLE_NAME}]}"

    # If user specified a specific SRR, skip others
    if [[ -n "${FILTER_SRR}" && "${SRR}" != "${FILTER_SRR}" ]]; then
        continue
    fi

    echo ""
    echo "============================================================"
    echo "PROCESSING: ${SAMPLE_NAME} (${SRR})"
    echo "============================================================"

    # Log file for this sample
    LOG="${QC_DIR}/logs/${SAMPLE_NAME}_processing.log"
    echo "Processing log for ${SAMPLE_NAME} (${SRR})" > "${LOG}"
    echo "Started: $(date)" >> "${LOG}"

    # =========================================================================
    # STEP A: DOWNLOAD from SRA
    # =========================================================================
    echo ""
    echo "--- STEP A: Download FASTQ for ${SRR} ---"

    FASTQ_FILE="${FASTQ_DIR}/${SRR}.fastq"

    if [[ -f "${FASTQ_FILE}" || -f "${FASTQ_FILE}.gz" ]]; then
        echo "[INFO] FASTQ already exists for ${SRR}. Skipping download."
        # Use the gz version if available
        if [[ -f "${FASTQ_FILE}.gz" && ! -f "${FASTQ_FILE}" ]]; then
            echo "[INFO] Decompressing existing .gz file..."
            pigz -dk "${FASTQ_FILE}.gz"
        fi
    else
        echo "[INFO] Downloading ${SRR} with fasterq-dump..."
        echo "[INFO] This may take 10-30 minutes depending on connection speed."

        # fasterq-dump downloads and converts SRA to FASTQ
        # --split-3: auto-detect SE vs PE; SE → single file, PE → _1 and _2 files
        # --outdir: where to save
        # --temp: temp directory — using local scratch for speed
        fasterq-dump "${SRR}" \
            --split-3 \
            --outdir "${FASTQ_DIR}" \
            --temp "${SCRATCH_TMPDIR}" \
            --threads "${THREADS}" \
            --progress

        echo "[INFO] Download complete."
    fi

    # --- Detect SE vs PE ---
    # fasterq-dump with --split-3 produces:
    #   SE: {SRR}.fastq
    #   PE: {SRR}_1.fastq and {SRR}_2.fastq (and possibly {SRR}.fastq for unpaired)
    if [[ -f "${FASTQ_DIR}/${SRR}_1.fastq" && -f "${FASTQ_DIR}/${SRR}_2.fastq" ]]; then
        LAYOUT="PE"
        R1="${FASTQ_DIR}/${SRR}_1.fastq"
        R2="${FASTQ_DIR}/${SRR}_2.fastq"
        echo "[INFO] Detected layout: PAIRED-END"
        echo "[INFO]   R1: ${R1}"
        echo "[INFO]   R2: ${R2}"
        echo "  WARNING: This pipeline is optimized for SE PRO-seq data."
        echo "  For PE data, we will use only R1 (the read containing the insert)."
        echo "  R2 in some PE PRO-seq protocols is the UMI or adapter-only read."
        # For standard PRO-seq PE, R1 is the informative read
        FASTQ_INPUT="${R1}"
    else
        LAYOUT="SE"
        FASTQ_INPUT="${FASTQ_FILE}"
        echo "[INFO] Detected layout: SINGLE-END"
        echo "[INFO]   FASTQ: ${FASTQ_INPUT}"
    fi

    # Sanity check
    RAW_READ_COUNT=$(count_reads_fastq "${FASTQ_INPUT}")
    echo "[SANITY CHECK] Raw read count: ${RAW_READ_COUNT}"
    echo "raw_reads=${RAW_READ_COUNT}" >> "${LOG}"
    echo "[SANITY CHECK] First read:"
    head -4 "${FASTQ_INPUT}" | sed 's/^/    /'

    # =========================================================================
    # STEP B: PRE-TRIM QC with FastQC
    # =========================================================================
    echo ""
    echo "--- STEP B: Pre-trim FastQC ---"

    if [[ -f "${QC_DIR}/fastqc_raw/${SRR}_fastqc.html" || \
          -f "${QC_DIR}/fastqc_raw/$(basename "${FASTQ_INPUT}" .fastq)_fastqc.html" ]]; then
        echo "[INFO] Pre-trim FastQC already exists. Skipping."
    else
        echo "[INFO] Running FastQC on raw reads..."
        fastqc "${FASTQ_INPUT}" \
            --outdir "${QC_DIR}/fastqc_raw" \
            --threads "${THREADS}" \
            --quiet
        echo "[INFO] FastQC complete."
    fi

    # =========================================================================
    # STEP C: ADAPTER TRIMMING & QUALITY FILTERING
    # =========================================================================
    echo ""
    echo "--- STEP C: Adapter trimming & quality filtering (cutadapt) ---"

    TRIMMED_FASTQ="${TRIMMED_DIR}/${SAMPLE_NAME}_trimmed.fastq.gz"
    CUTADAPT_LOG="${QC_DIR}/logs/${SAMPLE_NAME}_cutadapt.log"

    if [[ -f "${TRIMMED_FASTQ}" ]]; then
        echo "[INFO] Trimmed FASTQ already exists. Skipping."
    else
        echo "[INFO] Running cutadapt..."
        echo "[INFO]   Adapter (3'): ${ADAPTER_3PRIME}"
        echo "[INFO]   Quality cutoff: ${QUALITY_CUTOFF}"
        echo "[INFO]   Min length: ${MIN_READ_LENGTH}"

        # cutadapt flags explained:
        #   -a ADAPTER: trim 3' adapter sequence
        #   -q QUALITY: trim low-quality bases from 3' end
        #   -m MIN_LEN: discard reads shorter than this after trimming
        #   --trim-n: trim N bases from ends
        #   -j THREADS: use multiple cores
        #   -o OUTPUT: output file
        cutadapt \
            -a "${ADAPTER_3PRIME}" \
            -q "${QUALITY_CUTOFF}" \
            -m "${MIN_READ_LENGTH}" \
            --trim-n \
            -j "${THREADS}" \
            -o "${TRIMMED_FASTQ}" \
            "${FASTQ_INPUT}" \
            2>&1 | tee "${CUTADAPT_LOG}"

        echo "[INFO] Cutadapt complete."
    fi

    # Parse cutadapt stats
    TRIMMED_READ_COUNT=$(count_reads_fastq "${TRIMMED_FASTQ}")
    echo "[SANITY CHECK] Reads after trimming: ${TRIMMED_READ_COUNT}"
    TRIM_SURVIVAL_PCT=$(echo "scale=1; ${TRIMMED_READ_COUNT} * 100 / ${RAW_READ_COUNT}" | bc)
    echo "[SANITY CHECK] Survival rate: ${TRIM_SURVIVAL_PCT}%"
    echo "trimmed_reads=${TRIMMED_READ_COUNT}" >> "${LOG}"
    echo "trim_survival_pct=${TRIM_SURVIVAL_PCT}" >> "${LOG}"

    # =========================================================================
    # STEP D: POST-TRIM QC
    # =========================================================================
    echo ""
    echo "--- STEP D: Post-trim FastQC ---"

    if [[ -f "${QC_DIR}/fastqc_trimmed/${SAMPLE_NAME}_trimmed_fastqc.html" ]]; then
        echo "[INFO] Post-trim FastQC already exists. Skipping."
    else
        echo "[INFO] Running FastQC on trimmed reads..."
        fastqc "${TRIMMED_FASTQ}" \
            --outdir "${QC_DIR}/fastqc_trimmed" \
            --threads "${THREADS}" \
            --quiet
        echo "[INFO] FastQC complete."
    fi

    # =========================================================================
    # STEP E: SPIKE-IN FILTERING (align to dm6, keep unmapped)
    # =========================================================================
    echo ""
    echo "--- STEP E: Spike-in filtering (align to dm6) ---"

    NOSPIKEIN_FASTQ="${TRIMMED_DIR}/${SAMPLE_NAME}_nospikein.fastq.gz"

    if [[ -f "${NOSPIKEIN_FASTQ}" ]]; then
        echo "[INFO] Spike-in filtered FASTQ already exists. Skipping."
    else
        echo "[INFO] Aligning to dm6 to identify Drosophila spike-in reads..."
        echo "[INFO] Reads that do NOT align to dm6 → mouse reads (kept)"
        echo "[INFO] Reads that DO align to dm6 → spike-in reads (discarded)"

        # We use a named pipe trick to avoid writing a large intermediate file.
        # bowtie2 --un outputs reads that failed to align (= mouse reads).
        # We pipe the SAM output to /dev/null since we don't need it.
        #
        # --very-sensitive: be thorough so we catch all spike-in reads
        # --no-unal: don't include unaligned reads in SAM output (saves time)
        # --un: write unaligned reads (our mouse reads) to this file

        NOSPIKEIN_FASTQ_UNZIPPED="${SCRATCH_TMPDIR}/${SAMPLE_NAME}_nospikein.fastq"

        bowtie2 \
            -x "${DM6_BT2_IDX}" \
            -U "${TRIMMED_FASTQ}" \
            --very-sensitive \
            --no-unal \
            --un "${NOSPIKEIN_FASTQ_UNZIPPED}" \
            --threads "${THREADS}" \
            -S /dev/null \
            2> "${QC_DIR}/logs/${SAMPLE_NAME}_dm6_alignment.log"

        echo "[INFO] Spike-in alignment log:"
        cat "${QC_DIR}/logs/${SAMPLE_NAME}_dm6_alignment.log" | sed 's/^/    /'

        # Compress the output and move to final location
        echo "[INFO] Compressing spike-in-filtered FASTQ..."
        pigz -p "${THREADS}" -c "${NOSPIKEIN_FASTQ_UNZIPPED}" > "${NOSPIKEIN_FASTQ}"
        rm -f "${NOSPIKEIN_FASTQ_UNZIPPED}"
        echo "[INFO] Spike-in filtering complete."
    fi

    # Count reads surviving spike-in filter
    NOSPIKEIN_READ_COUNT=$(count_reads_fastq "${NOSPIKEIN_FASTQ}")
    SPIKEIN_COUNT=$((TRIMMED_READ_COUNT - NOSPIKEIN_READ_COUNT))
    SPIKEIN_PCT=$(echo "scale=2; ${SPIKEIN_COUNT} * 100 / ${TRIMMED_READ_COUNT}" | bc)
    echo "[SANITY CHECK] Reads after spike-in removal: ${NOSPIKEIN_READ_COUNT}"
    echo "[SANITY CHECK] Spike-in reads removed: ${SPIKEIN_COUNT} (${SPIKEIN_PCT}%)"
    echo "spikein_reads=${SPIKEIN_COUNT}" >> "${LOG}"
    echo "spikein_pct=${SPIKEIN_PCT}" >> "${LOG}"
    echo "post_spikein_reads=${NOSPIKEIN_READ_COUNT}" >> "${LOG}"

    # =========================================================================
    # STEP F: ALIGN TO mm39
    # =========================================================================
    echo ""
    echo "--- STEP F: Align to mm39 ---"

    RAW_BAM="${ALIGNED_DIR}/${SAMPLE_NAME}_mm39_raw.bam"

    if [[ -f "${RAW_BAM}" ]]; then
        echo "[INFO] mm39 alignment BAM already exists. Skipping."
    else
        echo "[INFO] Aligning to mm39 with bowtie2..."
        echo "[INFO] Settings: --sensitive, MAPQ filtering done in next step"

        # bowtie2 flags:
        #   -x: index prefix
        #   -U: unpaired reads (SE)
        #   --sensitive: good accuracy/speed tradeoff for typical PRO-seq reads
        #                (--very-sensitive is ~2x slower with marginal gain here)
        #   --no-unal: skip unaligned reads in output (saves space)
        #   --threads: parallelism
        # Pipe to samtools to convert SAM→BAM on the fly (saves disk space)

        bowtie2 \
            -x "${MM39_BT2_IDX}" \
            -U "${NOSPIKEIN_FASTQ}" \
            --sensitive \
            --no-unal \
            --threads "${THREADS}" \
            2> "${QC_DIR}/logs/${SAMPLE_NAME}_mm39_alignment.log" \
        | samtools view -bS -@ "${THREADS}" - \
        > "${RAW_BAM}"

        echo "[INFO] mm39 alignment log:"
        cat "${QC_DIR}/logs/${SAMPLE_NAME}_mm39_alignment.log" | sed 's/^/    /'
        echo "[INFO] Alignment complete."
    fi

    # =========================================================================
    # STEP G: POST-ALIGNMENT FILTERING
    # =========================================================================
    echo ""
    echo "--- STEP G: Filter, sort, index ---"

    FILTERED_BAM="${ALIGNED_DIR}/${SAMPLE_NAME}_mm39_filtered.bam"

    if [[ -f "${FILTERED_BAM}" ]]; then
        echo "[INFO] Filtered BAM already exists. Skipping."
    else
        echo "[INFO] Filtering: remove unmapped (flag 4), MAPQ < ${MAPQ_THRESHOLD}"
        echo "[INFO] Then sorting by coordinate..."

        # samtools view flags:
        #   -b: output BAM
        #   -F 4: exclude unmapped reads (bit flag 4)
        #   -q MAPQ: minimum mapping quality
        #     MAPQ >= 10 effectively removes most multimappers.
        #     bowtie2 assigns MAPQ=0 to reads that map equally well to
        #     multiple locations, and low MAPQ to ambiguous mappings.
        #
        # samtools sort:
        #   sort by coordinate (required for markdup and indexing)
        #   -T: use local scratch for temp files (faster I/O)

        samtools view -b -F 4 -q "${MAPQ_THRESHOLD}" -@ "${THREADS}" "${RAW_BAM}" \
        | samtools sort -@ "${THREADS}" -T "${SCRATCH_TMPDIR}/${SAMPLE_NAME}_sort" \
            -o "${FILTERED_BAM}" -

        samtools index "${FILTERED_BAM}"
        echo "[INFO] Filtering and sorting complete."
    fi

    FILTERED_READ_COUNT=$(samtools view -c "${FILTERED_BAM}")
    echo "[SANITY CHECK] Reads after filtering (unique, mapped, MAPQ>=${MAPQ_THRESHOLD}): ${FILTERED_READ_COUNT}"
    echo "filtered_mapped_reads=${FILTERED_READ_COUNT}" >> "${LOG}"

    # =========================================================================
    # STEP H: PCR DUPLICATE REMOVAL
    # =========================================================================
    echo ""
    echo "--- STEP H: Remove PCR duplicates (samtools markdup) ---"

    DEDUP_BAM="${ALIGNED_DIR}/${SAMPLE_NAME}_mm39_final.bam"
    DEDUP_STATS="${QC_DIR}/logs/${SAMPLE_NAME}_markdup_stats.txt"

    if [[ -f "${DEDUP_BAM}" ]]; then
        echo "[INFO] Deduplicated BAM already exists. Skipping."
    else
        echo "[INFO] Running samtools markdup..."
        echo "[INFO] For SE data, duplicates are defined by same alignment start"
        echo "[INFO] position and strand (which is appropriate for PRO-seq)."

        # samtools markdup:
        #   -r: remove duplicates (not just mark them)
        #   -s: report stats to stderr
        #   -S: mark SE reads (without requiring fixmate). This flag is
        #       critical for SE PRO-seq — without it, markdup may silently
        #       produce an empty or truncated BAM.
        #   -f FILE: write detailed stats to file
        #   Input must be coordinate-sorted (which it is from Step G)
        samtools markdup \
            -r \
            -s \
            -S \
            -f "${DEDUP_STATS}" \
            -@ "${THREADS}" \
            "${FILTERED_BAM}" \
            "${DEDUP_BAM}" \
            2> "${QC_DIR}/logs/${SAMPLE_NAME}_markdup.log"

        # Verify the output BAM was actually created
        if [[ ! -s "${DEDUP_BAM}" ]]; then
            echo "[ERROR] samtools markdup failed — output BAM is missing or empty!"
            echo "[ERROR] Check: ${QC_DIR}/logs/${SAMPLE_NAME}_markdup.log"
            cat "${QC_DIR}/logs/${SAMPLE_NAME}_markdup.log" | sed 's/^/    /'
            exit 1
        fi

        samtools index "${DEDUP_BAM}"
        echo "[INFO] Duplicate removal complete."
    fi

    # Report dedup stats
    FINAL_READ_COUNT=$(samtools view -c "${DEDUP_BAM}")
    DUP_COUNT=$((FILTERED_READ_COUNT - FINAL_READ_COUNT))
    if [[ "${FILTERED_READ_COUNT}" -gt 0 ]]; then
        DUP_PCT=$(echo "scale=2; ${DUP_COUNT} * 100 / ${FILTERED_READ_COUNT}" | bc)
    else
        DUP_PCT="0.00"
    fi
    echo "[SANITY CHECK] Reads after deduplication: ${FINAL_READ_COUNT}"
    echo "[SANITY CHECK] Duplicates removed: ${DUP_COUNT} (${DUP_PCT}%)"
    echo "final_reads=${FINAL_READ_COUNT}" >> "${LOG}"
    echo "duplicates=${DUP_COUNT}" >> "${LOG}"
    echo "duplicate_pct=${DUP_PCT}" >> "${LOG}"

    # Dedup stats summary
    echo "[INFO] samtools markdup stats:"
    cat "${DEDUP_STATS}" 2>/dev/null | sed 's/^/    /' || true

    # =========================================================================
    # Cleanup: compress raw FASTQ to save space
    # =========================================================================
    echo ""
    echo "--- Cleanup ---"
    if [[ -f "${FASTQ_INPUT}" && ! -f "${FASTQ_INPUT}.gz" ]]; then
        echo "[INFO] Compressing raw FASTQ to save space..."
        pigz -p "${THREADS}" "${FASTQ_INPUT}"
    fi
    # Remove the raw (unsorted/unfiltered) BAM — no longer needed
    if [[ -f "${RAW_BAM}" && -f "${DEDUP_BAM}" ]]; then
        echo "[INFO] Removing raw BAM: ${RAW_BAM}"
        rm -f "${RAW_BAM}"
    fi
    # KEEP the filtered (pre-dedup) BAM — it is the input for WASP.
    # Standard positional dedup (Picard) can introduce allelic bias,
    # so the WASP pipeline needs to start from the pre-dedup BAM and
    # apply its own allele-aware deduplication.
    if [[ -f "${FILTERED_BAM}" ]]; then
        echo "[INFO] KEEPING pre-dedup BAM for WASP: ${FILTERED_BAM}"
    fi

    echo ""
    echo "[SUMMARY] ${SAMPLE_NAME} (${SRR}):"
    echo "  Raw reads:            ${RAW_READ_COUNT}"
    echo "  After trimming:       ${TRIMMED_READ_COUNT} (${TRIM_SURVIVAL_PCT}% survived)"
    echo "  Spike-in removed:     ${SPIKEIN_COUNT} (${SPIKEIN_PCT}%)"
    echo "  After mm39 alignment: ${FILTERED_READ_COUNT} (unique, MAPQ>=${MAPQ_THRESHOLD})"
    echo "  After deduplication:  ${FINAL_READ_COUNT} (${DUP_PCT}% duplicates)"
    echo "  Final BAM:            ${DEDUP_BAM}"
    echo ""
    echo "Finished: $(date)" >> "${LOG}"

done

echo ""
echo "============================================================"
echo "All samples processed."
echo "Final BAMs are in: ${ALIGNED_DIR}"
echo "QC logs are in:    ${QC_DIR}/logs"
echo "Finished: $(date)"
echo "============================================================"
echo ""
echo "Next step: Run 02_make_bigwigs.sh"