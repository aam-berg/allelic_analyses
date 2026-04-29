#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04_proseq/06_wasp_filter.sh — WASP mapping bias correction (single-end)
# =============================================================================
#
# WHAT THIS DOES (per sample):
#   1. find_intersecting_snps.py (SE mode):
#        Identifies reads overlapping het SNPs. Outputs:
#          {sample}.keep.bam       : reads NOT overlapping SNPs
#          {sample}.to.remap.bam   : reads overlapping SNPs (the originals)
#          {sample}.remap.fq.gz    : allele-swapped versions for remapping
#
#   2. Re-map allele-swapped reads with bowtie2 SE, IDENTICAL to step 03's
#      bowtie2 invocation (--very-sensitive). Drift would invalidate WASP.
#
#   3. filter_remapped_reads.py: keeps only reads that mapped to the SAME
#      position regardless of allele swap.
#
#   4. Merge keep + remap.keep -> {SAMPLE}_wasp.bam
#
# NOTE — DIFFERENCE FROM 03_atac/06_wasp_filter.sh:
#   - Single-end here vs paired-end in ATAC. No --is_paired_end flag.
#   - Single FASTQ (.remap.fq.gz) rather than a pair (.remap.fq1.gz /
#     .remap.fq2.gz).
#   - The remap bowtie2 invocation matches step 03 here (--very-sensitive
#     for SE; ATAC uses --very-sensitive -X 1000 --no-mixed --no-discordant
#     for PE).
#
# PARALLEL EXECUTION:
#   sbatch --array=0-3 --partition=medium --time=8:00:00 --mem=32G \
#       --cpus-per-task=8 -o logs/06_%A_%a.out -e logs/06_%A_%a.err \
#       06_wasp_filter.sh
# =============================================================================

source "config.sh"

step_header "PRO-seq Pipeline Step 06: WASP Mapping Bias Correction (SE)"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] python:   $(python --version 2>&1)"
echo "[INFO] pysam:    $(python -c 'import pysam; print(pysam.__version__)')"
echo "[INFO] bowtie2:  $(bowtie2 --version 2>&1 | head -1)"

create_dirs
resolve_samples "${1:-}"

for SCRIPT in "${WASP_FIND}" "${WASP_FILTER}"; do
    [[ -f "${SCRIPT}" ]] || \
        { echo "[ERROR] WASP missing: ${SCRIPT}. Run 05_wasp_setup.sh." >&2; exit 1; }
done
ls "${WASP_SNP_DIR}"/chr*.snps.txt.gz 1>/dev/null 2>&1 || \
    { echo "[ERROR] No WASP SNP files. Run 05_wasp_setup.sh." >&2; exit 1; }

for SAMPLE in "${RUN_SAMPLES[@]}"; do
    step_header "WASP filtering: ${SAMPLE}"

    INPUT_BAM="${BAM_MM39_DIR}/${SAMPLE}_final.bam"
    WORKDIR="${WASP_INTERMEDIATE_DIR}/${SAMPLE}"
    FINAL_WASP_BAM="${WASP_FILTERED_DIR}/${SAMPLE}_wasp.bam"

    if [[ -f "${FINAL_WASP_BAM}" && -f "${FINAL_WASP_BAM}.bai" ]]; then
        echo "[INFO] WASP-filtered BAM exists. Skipping."
        continue
    fi
    [[ -f "${INPUT_BAM}" ]] || \
        { echo "[ERROR] Input BAM missing: ${INPUT_BAM}" >&2; exit 1; }

    mkdir -p "${WORKDIR}"

    INPUT_COUNT=$(samtools view -c "${INPUT_BAM}")
    echo "[INFO] ${SAMPLE} input: ${INPUT_COUNT} reads"

    # ---- WASP expects outputs in same dir as input; symlink in WORKDIR ----
    WASP_INPUT="${WORKDIR}/${SAMPLE}.bam"
    [[ -e "${WASP_INPUT}"     ]] || ln -sf "${INPUT_BAM}"     "${WASP_INPUT}"
    [[ -e "${WASP_INPUT}.bai" ]] || ln -sf "${INPUT_BAM}.bai" "${WASP_INPUT}.bai"

    KEEP_BAM="${WORKDIR}/${SAMPLE}.keep.bam"
    TO_REMAP_BAM="${WORKDIR}/${SAMPLE}.to.remap.bam"
    REMAP_FQ="${WORKDIR}/${SAMPLE}.remap.fq.gz"

    # =========================================================================
    # STEP 1 — find_intersecting_snps (SE mode; no --is_paired_end)
    # =========================================================================
    if [[ -f "${KEEP_BAM}" && -f "${REMAP_FQ}" ]]; then
        echo "[STEP 1] find_intersecting_snps output exists. Skipping."
    else
        echo "[STEP 1] Running find_intersecting_snps.py (SE)..."
        python "${WASP_FIND}" \
            --is_sorted \
            --output_dir "${WORKDIR}" \
            --snp_dir "${WASP_SNP_DIR}" \
            "${WASP_INPUT}"
    fi

    KEEP_COUNT=$(samtools view -c "${KEEP_BAM}")
    REMAP_COUNT=$(samtools view -c "${TO_REMAP_BAM}")
    REMAP_FQ_ENTRIES=$(zcat "${REMAP_FQ}" | awk 'NR%4==1' | wc -l)

    echo ""
    echo "[SANITY CHECK] After find_intersecting_snps:"
    echo "  Reads NOT overlapping SNPs:  ${KEEP_COUNT}"
    echo "  Reads overlapping SNPs:      ${REMAP_COUNT}"
    echo "  FASTQ entries to remap:      ${REMAP_FQ_ENTRIES}"
    echo "    (>= read count: reads at multiple SNPs produce multiple"
    echo "     allele-swapped entries)"
    if (( KEEP_COUNT + REMAP_COUNT > 0 )); then
        FRAC=$(awk "BEGIN {printf \"%.2f\", \
            ${REMAP_COUNT}/(${KEEP_COUNT}+${REMAP_COUNT})*100}")
        echo "  Fraction overlapping SNPs:   ${FRAC}%"
    fi

    # =========================================================================
    # STEP 2 — Remap allele-swapped reads (IDENTICAL bowtie2 args to step 03)
    # =========================================================================
    REMAP_BAM="${WORKDIR}/${SAMPLE}.remap.bam"
    REMAP_LOG="${LOG_DIR}/${SAMPLE}_remap_bowtie2.log"

    if [[ -f "${REMAP_BAM}" && -f "${REMAP_BAM}.bai" ]]; then
        echo ""
        echo "[STEP 2] Remap BAM exists. Skipping."
    else
        echo ""
        echo "[STEP 2] Re-mapping allele-swapped reads..."
        # CRITICAL: identical to step 03's mm39 bowtie2 args
        bowtie2 \
            --very-sensitive \
            --no-unal \
            -x "${MM39_BT2_IDX}" \
            -U "${REMAP_FQ}" \
            -p "${THREADS}" \
            2> "${REMAP_LOG}" \
        | samtools view -bS -@ 2 - \
        | samtools sort -@ "${THREADS}" -o "${REMAP_BAM}"
        samtools index "${REMAP_BAM}"
    fi

    REMAP_ALIGNED=$(samtools view -c "${REMAP_BAM}")
    echo "[SANITY CHECK] Remapping:"
    echo "  Input to remap:    ${REMAP_COUNT}"
    echo "  Re-mapped:         ${REMAP_ALIGNED}"

    # =========================================================================
    # STEP 3 — filter_remapped_reads
    # =========================================================================
    REMAP_KEEP="${WORKDIR}/${SAMPLE}.remap.keep.bam"

    if [[ -f "${REMAP_KEEP}" ]]; then
        echo ""
        echo "[STEP 3] Filtered remap BAM exists. Skipping."
    else
        echo ""
        echo "[STEP 3] Filtering remapped reads (must map to same position)..."
        python "${WASP_FILTER}" \
            "${TO_REMAP_BAM}" \
            "${REMAP_BAM}" \
            "${REMAP_KEEP}"
    fi

    REMAP_KEEP_COUNT=$(samtools view -c "${REMAP_KEEP}")
    REMAP_REMOVED=$((REMAP_COUNT - REMAP_KEEP_COUNT))
    if (( REMAP_COUNT > 0 )); then
        REMAP_REMOVED_PCT=$(awk "BEGIN {printf \"%.2f\", ${REMAP_REMOVED}/${REMAP_COUNT}*100}")
    else
        REMAP_REMOVED_PCT="NA"
    fi

    echo "[SANITY CHECK] WASP filtering:"
    echo "  SNP-overlapping reads input:   ${REMAP_COUNT}"
    echo "  Passed WASP (same position):   ${REMAP_KEEP_COUNT}"
    echo "  Removed (mapping bias):        ${REMAP_REMOVED} (${REMAP_REMOVED_PCT}%)"
    echo "[INTERPRETATION] Mouse F1 hybrid: 5-15% typical; >25% high"

    # =========================================================================
    # STEP 4 — Merge keep + remap.keep
    # =========================================================================
    echo ""
    echo "[STEP 4] Merging into final WASP BAM..."

    REMAP_KEEP_SORTED="${WORKDIR}/${SAMPLE}.remap.keep.sorted.bam"
    samtools sort -@ "${THREADS}" -o "${REMAP_KEEP_SORTED}" "${REMAP_KEEP}"
    samtools index "${REMAP_KEEP_SORTED}"

    samtools merge -@ "${THREADS}" -f "${WORKDIR}/${SAMPLE}_merged.bam" \
        "${KEEP_BAM}" "${REMAP_KEEP_SORTED}"
    samtools sort -@ "${THREADS}" -o "${FINAL_WASP_BAM}" \
        "${WORKDIR}/${SAMPLE}_merged.bam"
    samtools index "${FINAL_WASP_BAM}"

    rm -f "${WORKDIR}/${SAMPLE}_merged.bam"

    # ---- Final stats ----
    FINAL_COUNT=$(samtools view -c "${FINAL_WASP_BAM}")
    TOTAL_REMOVED=$((INPUT_COUNT - FINAL_COUNT))
    if (( INPUT_COUNT > 0 )); then
        TOTAL_REMOVED_PCT=$(awk "BEGIN {printf \"%.2f\", ${TOTAL_REMOVED}/${INPUT_COUNT}*100}")
    else
        TOTAL_REMOVED_PCT="NA"
    fi

    PLUS=$(samtools view -c -F 16 "${FINAL_WASP_BAM}")
    MINUS=$(samtools view -c -f 16 "${FINAL_WASP_BAM}")
    if (( MINUS > 0 )); then
        RATIO=$(awk "BEGIN {printf \"%.3f\", ${PLUS}/${MINUS}}")
    else
        RATIO="NA"
    fi

    echo ""
    echo "[SANITY CHECK] Final WASP BAM: ${SAMPLE}"
    echo "  Input (pre-WASP):    ${INPUT_COUNT}"
    echo "  Output (post-WASP):  ${FINAL_COUNT}"
    echo "  Total removed:       ${TOTAL_REMOVED} (${TOTAL_REMOVED_PCT}%)"
    echo "  Strand ratio:        ${RATIO} (still ~1.0 expected)"

    cat > "${QC_DIR}/wasp_stats_${SAMPLE}.tsv" << EOF
sample	input	no_snp	snp_overlap	snp_remapped	snp_kept	snp_removed	snp_removed_pct	final	total_removed_pct	strand_ratio
${SAMPLE}	${INPUT_COUNT}	${KEEP_COUNT}	${REMAP_COUNT}	${REMAP_ALIGNED}	${REMAP_KEEP_COUNT}	${REMAP_REMOVED}	${REMAP_REMOVED_PCT}	${FINAL_COUNT}	${TOTAL_REMOVED_PCT}	${RATIO}
EOF
done

step_header "Step 06 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output:"
echo "  WASP BAMs:  ${WASP_FILTERED_DIR}/"
echo "  QC stats:   ${QC_DIR}/wasp_stats_*.tsv"
echo ""
echo "Next: bash 07_allele_specific.sh"
