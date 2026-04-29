#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 03_align_filter.sh — rRNA removal → Spike-in removal → mm39 alignment → Filter
# =============================================================================
#
# WHAT THIS DOES:
#   1. Align trimmed reads to rRNA reference → DISCARD rRNA reads
#   2. Align surviving reads to dm6 (Drosophila spike-in) → DISCARD spike-in reads
#   3. Align surviving reads to mm39 (mouse) genome using bowtie2
#   4. Filter: remove low-MAPQ, secondary/supplementary, and chrM reads
#   5. Generate final filtered BAM + comprehensive QC stats
#
# STEP-BY-STEP RATIONALE:
#
#   STEP 1 — rRNA FILTERING (new):
#     rDNA loci are among the most heavily transcribed regions in the genome.
#     RNA Pol I transcribes the 45S pre-rRNA repeat, and in PRO-seq this
#     produces a substantial number of reads (often 5-20% of total).
#     Problem: some rRNA reads map uniquely to the genome and pass MAPQ
#     filtering, so MAPQ alone doesn't remove them all.
#     Solution: align to an rRNA reference first (28S, 18S, 5.8S, 5S, 45S)
#     and discard aligned reads. This is standard in dedicated PRO-seq
#     pipelines (e.g., proseq2.0 from the Danko lab).
#
#   STEP 2 — SPIKE-IN FILTERING:
#     This PRO-seq experiment includes Drosophila S2 cells as spike-in controls.
#     If we skip this, fly reads would contaminate our mouse signal (especially
#     in repetitive regions). Aligning to dm6 FIRST and removing those reads
#     is standard.
#
#   STEP 3 — BOWTIE2 ALIGNMENT TO mm39:
#     *** CRITICAL DIFFERENCE FROM RNA-seq ***
#     - RNA-seq uses splice-aware aligners (STAR, HISAT2) because mature mRNA
#       has been spliced and reads span exon-exon junctions.
#     - PRO-seq measures NASCENT RNA that has NOT been spliced. Introns are
#       still present.
#     - Therefore we use bowtie2 (non-splice-aware). A splice-aware aligner
#       would INCORRECTLY split reads across non-existent splice junctions.
#
#   STEP 4 — MAPQ FILTERING (>= 10):
#     Removes multi-mapping reads that create ambiguous signal.
#
#   MITOCHONDRIAL READ REMOVAL:
#     chrM genes are highly transcribed and skew normalization.
#
#   NO PCR DEDUPLICATION:
#     PRO-seq measures Pol II positions at single-nucleotide resolution. At
#     strong pause sites, many Pol II molecules sit at the EXACT SAME nucleotide.
#     Position-based dedup would collapse real biological signal into one read.
#     This is standard (Mahat et al. 2016, proseq2.0, original GEO processing).
#
# PARALLEL EXECUTION:
#   sbatch --array=0-3 --partition=short --time=12:00:00 --mem=32G --cpus-per-task=8 \
#       -o logs/03_%A_%a.out -e logs/03_%A_%a.err 03_align_filter.sh
# =============================================================================

source "config.sh"

step_header "PRO-seq Pipeline Step 03: Align & Filter"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] bowtie2 version: $(bowtie2 --version 2>&1 | head -1)"
echo "[INFO] samtools version: $(samtools --version | head -1)"

create_dirs
resolve_samples "${1:-}"

# --- Verify reference indices exist ---
echo ""
echo "[INFO] Checking reference indices..."
for IDX_PREFIX in "${RRNA_BT2_IDX}" "${DM6_BT2_IDX}" "${MM39_BT2_IDX}"; do
    if [[ ! -f "${IDX_PREFIX}.1.bt2" ]]; then
        echo "[ERROR] Missing bowtie2 index: ${IDX_PREFIX}.1.bt2"
        echo "[ERROR] Run 00_setup_references.sh first."
        exit 1
    fi
done
echo "[INFO] All reference indices found."

for SAMPLE in "${RUN_SAMPLES[@]}"; do
    step_header "Processing ${SAMPLE}"

    INPUT="${TRIMMED_DIR}/${SAMPLE}_trimmed.fastq.gz"
    if [[ ! -f "${INPUT}" ]]; then
        echo "[ERROR] Trimmed FASTQ not found: ${INPUT}"
        exit 1
    fi

    # Track read counts through the filtering cascade
    TRIMMED_READS=""
    RRNA_ALIGNED=""
    DM6_ALIGNED=""
    MM39_ALIGNED_1=""
    MM39_ALIGNED_MULTI=""

    # =========================================================================
    # STEP 1: Filter rRNA reads
    # =========================================================================
    echo "[STEP 1] Filtering rRNA reads..."

    RRNA_UNMAPPED="${BAM_RRNA_DIR}/${SAMPLE}_rRNA_unmapped.fastq.gz"
    RRNA_LOG="${LOG_DIR}/${SAMPLE}_rRNA_bowtie2.log"

    if [[ -f "${RRNA_UNMAPPED}" ]]; then
        echo "[INFO] rRNA-filtered reads already exist. Skipping."
    else
        # --very-fast: rRNA sequences are simple; speed > sensitivity here
        # --no-unal:   don't output unmapped to SAM (we want them as FASTQ)
        # --un-gz:     write non-rRNA reads to gzipped FASTQ
        bowtie2 \
            --very-fast \
            --no-unal \
            --un-gz "${RRNA_UNMAPPED}" \
            -x "${RRNA_BT2_IDX}" \
            -U "${INPUT}" \
            -p "${THREADS}" \
            2> "${RRNA_LOG}" \
        | samtools view -bS -@ 2 - > /dev/null   # discard the BAM, we only need unmapped FASTQ
    fi

    # Parse stats
    TRIMMED_READS=$(grep "reads; of these:" "${RRNA_LOG}" | awk '{print $1}')
    RRNA_ALIGNED_1=$(grep "aligned exactly 1 time" "${RRNA_LOG}" | awk '{print $1}')
    RRNA_MULTI=$(grep "aligned >1 times" "${RRNA_LOG}" | awk '{print $1}')
    RRNA_ALIGNED=$((${RRNA_ALIGNED_1:-0} + ${RRNA_MULTI:-0}))
    RRNA_PCT=$(awk "BEGIN {printf \"%.2f\", ${RRNA_ALIGNED}/${TRIMMED_READS:-1}*100}")
    POST_RRNA=$((${TRIMMED_READS:-0} - ${RRNA_ALIGNED}))

    echo "[SANITY CHECK] rRNA filtering:"
    echo "  Input reads:     ${TRIMMED_READS}"
    echo "  rRNA reads:      ${RRNA_ALIGNED} (${RRNA_PCT}%)"
    echo "  Passed to dm6:   ${POST_RRNA}"
    echo ""
    echo "[INTERPRETATION] rRNA rates in PRO-seq:"
    echo "  5-20%:  Typical (rDNA is heavily transcribed by Pol I)"
    echo "  >30%:   High but not unusual for some cell types"
    echo "  <2%:    Surprisingly low — check rRNA index"

    # =========================================================================
    # STEP 2: Filter dm6 spike-in reads
    # =========================================================================
    echo ""
    echo "[STEP 2] Filtering dm6 spike-in reads..."

    DM6_BAM="${BAM_DM6_DIR}/${SAMPLE}_dm6.bam"
    DM6_UNMAPPED="${BAM_DM6_DIR}/${SAMPLE}_dm6_unmapped.fastq.gz"
    DM6_LOG="${LOG_DIR}/${SAMPLE}_dm6_bowtie2.log"

    if [[ -f "${DM6_UNMAPPED}" ]]; then
        echo "[INFO] dm6 unmapped reads already exist. Skipping."
    else
        bowtie2 \
            --very-sensitive \
            --no-unal \
            --un-gz "${DM6_UNMAPPED}" \
            -x "${DM6_BT2_IDX}" \
            -U "${RRNA_UNMAPPED}" \
            -p "${THREADS}" \
            2> "${DM6_LOG}" \
        | samtools view -bS -@ 2 - \
        | samtools sort -@ 2 -o "${DM6_BAM}"

        samtools index "${DM6_BAM}"
    fi

    DM6_INPUT=$(grep "reads; of these:" "${DM6_LOG}" | awk '{print $1}')
    DM6_ALIGNED_1=$(grep "aligned exactly 1 time" "${DM6_LOG}" | awk '{print $1}')
    DM6_MULTI=$(grep "aligned >1 times" "${DM6_LOG}" | awk '{print $1}')
    DM6_ALIGNED=$((${DM6_ALIGNED_1:-0} + ${DM6_MULTI:-0}))
    DM6_PCT=$(awk "BEGIN {printf \"%.2f\", ${DM6_ALIGNED}/${DM6_INPUT:-1}*100}")

    echo "[SANITY CHECK] dm6 spike-in:"
    echo "  Input reads:     ${DM6_INPUT}"
    echo "  Aligned to dm6:  ${DM6_ALIGNED} (${DM6_PCT}%)"
    echo "  Passed to mm39:  $((${DM6_INPUT:-0} - ${DM6_ALIGNED}))"
    echo ""
    echo "[INTERPRETATION] Spike-in rates:"
    echo "  1-10%:  Typical"
    echo "  >20%:   High — spike-in cells may be over-represented"
    echo "  <0.5%:  Very low — spike-in normalization may be unreliable"

    # =========================================================================
    # STEP 3: Align to mm39
    # =========================================================================
    echo ""
    echo "[STEP 3] Aligning to mm39..."

    MM39_RAW_BAM="${BAM_MM39_DIR}/${SAMPLE}_mm39_raw.bam"
    MM39_LOG="${LOG_DIR}/${SAMPLE}_mm39_bowtie2.log"

    if [[ -f "${MM39_RAW_BAM}" ]]; then
        echo "[INFO] mm39 raw BAM already exists. Skipping."
    else
        bowtie2 \
            --very-sensitive \
            --no-unal \
            -x "${MM39_BT2_IDX}" \
            -U "${DM6_UNMAPPED}" \
            -p "${THREADS}" \
            2> "${MM39_LOG}" \
        | samtools view -bS -@ 2 - \
        | samtools sort -@ "${THREADS}" -o "${MM39_RAW_BAM}"

        samtools index "${MM39_RAW_BAM}"
    fi

    MM39_INPUT=$(grep "reads; of these:" "${MM39_LOG}" | awk '{print $1}')
    MM39_ALIGNED_1=$(grep "aligned exactly 1 time" "${MM39_LOG}" | awk '{print $1}')
    MM39_ALIGNED_MULTI=$(grep "aligned >1 times" "${MM39_LOG}" | awk '{print $1}')
    MM39_TOTAL=$(( ${MM39_ALIGNED_1:-0} + ${MM39_ALIGNED_MULTI:-0} ))
    MM39_PCT=$(awk "BEGIN {printf \"%.2f\", ${MM39_TOTAL}/${MM39_INPUT:-1}*100}")

    echo "[SANITY CHECK] mm39 alignment:"
    echo "  Input reads:   ${MM39_INPUT}"
    echo "  Aligned:       ${MM39_TOTAL} (${MM39_PCT}%)"
    echo "  Unique:        ${MM39_ALIGNED_1:-0}"
    echo "  Multi-mapped:  ${MM39_ALIGNED_MULTI:-0}"
    echo ""
    echo "[INTERPRETATION] PRO-seq alignment rates:"
    echo "  >70%:    Good"
    echo "  50-70%:  Acceptable"
    echo "  <50%:    Investigate (contamination? wrong genome?)"

    # =========================================================================
    # STEP 4: Filter — MAPQ, secondary/supplementary, chrM
    # =========================================================================
    echo ""
    echo "[STEP 4] Filtering: MAPQ >= ${MAPQ_THRESHOLD}, remove chrM..."

    FINAL_BAM="${BAM_MM39_DIR}/${SAMPLE}_final.bam"

    if [[ -f "${FINAL_BAM}" ]]; then
        echo "[INFO] Final BAM already exists. Skipping."
    else
        # -F 2820 = exclude unmapped (4) + secondary (256) + QC-fail (512) + supplementary (2048)
        # -q 10   = MAPQ >= 10 (removes multi-mappers)
        # awk      = remove chrM reads
        samtools view -b -q "${MAPQ_THRESHOLD}" -F 2820 "${MM39_RAW_BAM}" \
        | samtools view -h - \
        | awk '$1 ~ /^@/ || $3 != "chrM"' \
        | samtools view -bS - \
        | samtools sort -@ "${THREADS}" -o "${FINAL_BAM}"

        samtools index "${FINAL_BAM}"
    fi


    # =========================================================================
    # OPTIONAL STEP 4b: 5' end deduplication
    # =========================================================================
    if [[ "${DEDUP_5PRIME}" == "true" ]]; then
        echo ""
        echo "[STEP 4b] Optional: 5' end deduplication..."
        echo "[NOTE] This collapses reads sharing the same 5' mapping position + strand."
        echo "[NOTE] Use with caution — at strong Pol II pause sites, identical positions"
        echo "       are EXPECTED biological signal, not PCR artifacts."

        DEDUP_BAM="${BAM_MM39_DIR}/${SAMPLE}_final${DEDUP_SUFFIX}.bam"
        DEDUP_STATS="${LOG_DIR}/${SAMPLE}_markdup_stats.txt"

        if [[ -f "${DEDUP_BAM}" ]]; then
            echo "[INFO] Dedup BAM already exists. Skipping."
        else
            samtools markdup -r -s \
                "${FINAL_BAM}" "${DEDUP_BAM}" \
                2> "${DEDUP_STATS}"

            samtools index "${DEDUP_BAM}"

            # Parse stats
            BEFORE_DEDUP=$(samtools view -c "${FINAL_BAM}")
            AFTER_DEDUP=$(samtools view -c "${DEDUP_BAM}")
            REMOVED=$((BEFORE_DEDUP - AFTER_DEDUP))
            PCT_DUP=$(awk "BEGIN {printf \"%.2f\", ${REMOVED}/${BEFORE_DEDUP:-1}*100}")

            echo "[SANITY CHECK] 5' deduplication: ${SAMPLE}"
            echo "  Before:    ${BEFORE_DEDUP}"
            echo "  After:     ${AFTER_DEDUP}"
            echo "  Removed:   ${REMOVED} (${PCT_DUP}%)"
            echo ""
            echo "[INTERPRETATION]"
            echo "  <10% removed:  Low duplication (typical for good PRO-seq libraries)"
            echo "  10-30%:        Moderate — may indicate some PCR over-amplification"
            echo "  >30%:          High — library complexity may be low"

            echo -e "${SAMPLE}\t${BEFORE_DEDUP}\t${AFTER_DEDUP}\t${REMOVED}\t${PCT_DUP}" \
                > "${QC_DIR}/dedup_stats_${SAMPLE}.tsv"
        fi
    fi



    # --- Collect final QC stats ---
    FINAL_COUNT=$(samtools view -c "${FINAL_BAM}")
    CHRM_READS=$(samtools view -c "${MM39_RAW_BAM}" chrM 2>/dev/null || echo "0")
    PLUS_READS=$(samtools view -c -F 16 "${FINAL_BAM}")
    MINUS_READS=$(samtools view -c -f 16 "${FINAL_BAM}")
    STRAND_RATIO=$(awk "BEGIN {printf \"%.3f\", ${PLUS_READS}/${MINUS_READS:-1}}")

    echo ""
    echo "[SANITY CHECK] Final BAM: ${SAMPLE}"
    echo "  Final read count:  ${FINAL_COUNT}"
    echo "  chrM removed:      ${CHRM_READS}"
    echo "  + strand (BAM):    ${PLUS_READS}"
    echo "  - strand (BAM):    ${MINUS_READS}"
    echo "  Strand ratio:      ${STRAND_RATIO} (expect ~1.0)"

    # --- Write per-sample QC file (safe for parallel writes) ---
    cat > "${QC_DIR}/alignment_stats_${SAMPLE}.tsv" << EOF
sample	trimmed_reads	rRNA_reads	rRNA_pct	dm6_reads	dm6_pct	mm39_input	mm39_aligned	mm39_pct	mm39_unique	mm39_multi	chrM_removed	final_reads	plus_strand	minus_strand	strand_ratio
${SAMPLE}	${TRIMMED_READS}	${RRNA_ALIGNED}	${RRNA_PCT}	${DM6_ALIGNED}	${DM6_PCT}	${MM39_INPUT}	${MM39_TOTAL}	${MM39_PCT}	${MM39_ALIGNED_1:-0}	${MM39_ALIGNED_MULTI:-0}	${CHRM_READS}	${FINAL_COUNT}	${PLUS_READS}	${MINUS_READS}	${STRAND_RATIO}
EOF

    echo ""
    echo "[INFO] Per-sample stats: ${QC_DIR}/alignment_stats_${SAMPLE}.tsv"

done

step_header "Step 03 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output:"
echo "  Final BAMs:  ${BAM_MM39_DIR}/*_final.bam"
echo "  Stats:       ${QC_DIR}/alignment_stats_*.tsv"
echo ""
echo "Next: Run 04_make_bigwigs.sh (after all step 03 array tasks finish)"