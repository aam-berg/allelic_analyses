#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 06_wasp_filter.sh — WASP mapping bias correction (paired-end)
# =============================================================================
#
# WHAT THIS DOES (per sample):
#   1. find_intersecting_snps.py (PE mode):
#        Identifies fragments where at least one mate covers a het SNP.
#        Outputs:
#          - {sample}.keep.bam       : fragments with no SNP overlap
#          - {sample}.to.remap.bam   : fragments overlapping SNPs (the originals)
#          - {sample}.remap.fq1.gz   : R1 reads with allele swapped (one fastq
#          - {sample}.remap.fq2.gz     entry per allele combination)
#
#   2. Re-map the allele-swapped reads with bowtie2 PE, using the SAME
#      parameters as in step 03. This is critical: any alignment-parameter
#      drift would invalidate WASP's logic.
#
#   3. filter_remapped_reads.py:
#        Compares re-mapped positions against original positions; keeps
#        only fragments that mapped to the SAME position regardless of allele.
#
#   4. Merge {sample}.keep.bam + filtered remap BAM -> {sample}_wasp.bam
#
# WHY THIS ORDER:
#   WASP runs on the FINAL filtered BAM (post MAPQ + chrM + dedup), but BEFORE
#   peak calling has already happened in step 04. Peaks are called on the
#   pre-WASP BAM because peak boundaries shouldn't depend on allele balance.
#   Allele-specific quantification within those peaks happens later, in
#   06_allele_pairs/.
#
# WHY --no-mixed --no-discordant ON REMAP:
#   Identical to step 03's bowtie2 invocation. Any difference would break
#   WASP's logic.
#
# PARALLEL EXECUTION:
#   sbatch --array=0-3 --partition=medium --time=8:00:00 --mem=32G \
#       --cpus-per-task=8 -o logs/06_%A_%a.out -e logs/06_%A_%a.err \
#       06_wasp_filter.sh
# =============================================================================

source "config.sh"

step_header "ATAC Pipeline Step 06: WASP Mapping Bias Correction (PE)"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] python:   $(python --version 2>&1)"
echo "[INFO] pysam:    $(python -c 'import pysam; print(pysam.__version__)')"
echo "[INFO] bowtie2:  $(bowtie2 --version 2>&1 | head -1)"

create_dirs
resolve_samples "${1:-}"

# Verify WASP scripts exist
for SCRIPT in "${WASP_FIND}" "${WASP_FILTER}"; do
    if [[ ! -f "${SCRIPT}" ]]; then
        echo "[ERROR] WASP script not found: ${SCRIPT}" >&2
        echo "[ERROR] Run 05_wasp_setup.sh first." >&2
        exit 1
    fi
done
if ! ls "${WASP_SNP_DIR}"/chr*.snps.txt.gz 1>/dev/null 2>&1; then
    echo "[ERROR] No WASP SNP files in ${WASP_SNP_DIR}" >&2
    echo "[ERROR] Run 05_wasp_setup.sh first." >&2
    exit 1
fi

for SAMPLE in "${RUN_SAMPLES[@]}"; do
    step_header "WASP filtering: ${SAMPLE}"

    INPUT_BAM="${BAM_FINAL_DIR}/${SAMPLE}_final.bam"
    WORKDIR="${WASP_INTERMEDIATE_DIR}/${SAMPLE}"
    FINAL_WASP_BAM="${WASP_FILTERED_DIR}/${SAMPLE}_wasp.bam"

    if [[ -f "${FINAL_WASP_BAM}" ]]; then
        echo "[INFO] WASP-filtered BAM exists for ${SAMPLE}: ${FINAL_WASP_BAM}"
        echo "[INFO] Skipping. Delete to re-run."
        continue
    fi

    if [[ ! -f "${INPUT_BAM}" ]]; then
        echo "[ERROR] Input BAM not found: ${INPUT_BAM}" >&2
        exit 1
    fi

    mkdir -p "${WORKDIR}"

    INPUT_COUNT=$(samtools view -c "${INPUT_BAM}")
    echo "[INFO] ${SAMPLE} input: ${INPUT_COUNT} reads"

    # ===========================================================================
    # STEP 1 — find_intersecting_snps (PE mode)
    # ===========================================================================
    # WASP names outputs after the input filename. We symlink into WORKDIR
    # so all outputs land there.
    WASP_INPUT="${WORKDIR}/${SAMPLE}.bam"
    [[ -e "${WASP_INPUT}"    ]] || ln -sf "${INPUT_BAM}"     "${WASP_INPUT}"
    [[ -e "${WASP_INPUT}.bai" ]] || ln -sf "${INPUT_BAM}.bai" "${WASP_INPUT}.bai"

    KEEP_BAM="${WORKDIR}/${SAMPLE}.keep.bam"
    TO_REMAP_BAM="${WORKDIR}/${SAMPLE}.to.remap.bam"
    REMAP_FQ1="${WORKDIR}/${SAMPLE}.remap.fq1.gz"
    REMAP_FQ2="${WORKDIR}/${SAMPLE}.remap.fq2.gz"

    if [[ -f "${KEEP_BAM}" && -f "${REMAP_FQ1}" && -f "${REMAP_FQ2}" ]]; then
        echo "[STEP 1] find_intersecting_snps output exists. Skipping."
    else
        echo "[STEP 1] Running find_intersecting_snps.py --is_paired_end..."
        python "${WASP_FIND}" \
            --is_paired_end \
            --is_sorted \
            --output_dir "${WORKDIR}" \
            --snp_dir "${WASP_SNP_DIR}" \
            "${WASP_INPUT}"
    fi

    KEEP_COUNT=$(samtools view -c "${KEEP_BAM}")
    REMAP_COUNT=$(samtools view -c "${TO_REMAP_BAM}")

    echo ""
    echo "[SANITY CHECK] After find_intersecting_snps:"
    echo "  Reads NOT overlapping SNPs:   ${KEEP_COUNT}"
    echo "  Reads overlapping SNPs:       ${REMAP_COUNT}"
    if [[ $((KEEP_COUNT + REMAP_COUNT)) -gt 0 ]]; then
        FRAC_SNP=$(awk "BEGIN {printf \"%.2f\", \
            ${REMAP_COUNT}/(${KEEP_COUNT}+${REMAP_COUNT})*100}")
        echo "  Fraction overlapping SNPs:    ${FRAC_SNP}%"
    fi

    # ===========================================================================
    # STEP 2 — Remap with IDENTICAL bowtie2 parameters
    # ===========================================================================
    REMAP_BAM="${WORKDIR}/${SAMPLE}.remap.bam"
    REMAP_LOG="${LOG_DIR}/${SAMPLE}_remap_bowtie2.log"

    if [[ -f "${REMAP_BAM}" && -f "${REMAP_BAM}.bai" ]]; then
        echo ""
        echo "[STEP 2] Remap BAM exists. Skipping."
    else
        echo ""
        echo "[STEP 2] Re-mapping allele-swapped reads with bowtie2 PE..."
        # CRITICAL: same parameters as step 03's alignment
        bowtie2 \
            --very-sensitive \
            -X "${MAX_INSERT}" \
            --no-mixed \
            --no-discordant \
            --no-unal \
            -x "${MM39_BT2_IDX}" \
            -1 "${REMAP_FQ1}" -2 "${REMAP_FQ2}" \
            -p "${THREADS}" \
            2> "${REMAP_LOG}" \
        | samtools view -bS -@ 2 - \
        | samtools sort -@ "${THREADS}" -o "${REMAP_BAM}"

        samtools index "${REMAP_BAM}"
    fi

    REMAP_ALIGNED=$(samtools view -c "${REMAP_BAM}")
    echo "[SANITY CHECK] Remapping:"
    echo "  Input reads to remap:    ${REMAP_COUNT}"
    echo "  Re-mapped successfully:  ${REMAP_ALIGNED}"

    # ===========================================================================
    # STEP 3 — filter_remapped_reads
    # ===========================================================================
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

    echo "[SANITY CHECK] WASP filtering of SNP-overlapping reads:"
    echo "  SNP-overlapping reads input:   ${REMAP_COUNT}"
    echo "  Passed WASP (same position):   ${REMAP_KEEP_COUNT}"
    echo "  Removed (mapping bias):        ${REMAP_REMOVED} (${REMAP_REMOVED_PCT}%)"
    echo ""
    echo "[INTERPRETATION] Mapping bias removal rates:"
    echo "  5-15%: typical for mouse F1 hybrids"
    echo "  >25%:  high; may indicate divergent regions or alignment issues"
    echo "  <5%:   very low; mapping was already mostly bias-free"

    # ===========================================================================
    # STEP 4 — Merge keep + remap.keep into final WASP BAM
    # ===========================================================================
    echo ""
    echo "[STEP 4] Merging keep + remap.keep into final WASP BAM..."

    REMAP_KEEP_SORTED="${WORKDIR}/${SAMPLE}.remap.keep.sorted.bam"
    samtools sort -@ "${THREADS}" -o "${REMAP_KEEP_SORTED}" "${REMAP_KEEP}"
    samtools index "${REMAP_KEEP_SORTED}"

    samtools merge -@ "${THREADS}" -f "${WORKDIR}/${SAMPLE}_merged.bam" \
        "${KEEP_BAM}" "${REMAP_KEEP_SORTED}"

    samtools sort -@ "${THREADS}" -o "${FINAL_WASP_BAM}" \
        "${WORKDIR}/${SAMPLE}_merged.bam"
    samtools index "${FINAL_WASP_BAM}"

    rm -f "${WORKDIR}/${SAMPLE}_merged.bam"

    # ===========================================================================
    # Final stats
    # ===========================================================================
    FINAL_COUNT=$(samtools view -c "${FINAL_WASP_BAM}")
    TOTAL_REMOVED=$((INPUT_COUNT - FINAL_COUNT))
    if (( INPUT_COUNT > 0 )); then
        TOTAL_REMOVED_PCT=$(awk "BEGIN {printf \"%.2f\", ${TOTAL_REMOVED}/${INPUT_COUNT}*100}")
    else
        TOTAL_REMOVED_PCT="NA"
    fi

    echo ""
    echo "[SANITY CHECK] Final WASP BAM: ${SAMPLE}"
    echo "  Input (pre-WASP):    ${INPUT_COUNT}"
    echo "  Output (post-WASP):  ${FINAL_COUNT}"
    echo "  Total removed:       ${TOTAL_REMOVED} (${TOTAL_REMOVED_PCT}%)"

    # Strand balance still ~1.0 (ATAC; nothing in WASP biases strand)
    PLUS=$(samtools view -c -F 16 "${FINAL_WASP_BAM}")
    MINUS=$(samtools view -c -f 16 "${FINAL_WASP_BAM}")
    if (( MINUS > 0 )); then
        RATIO=$(awk "BEGIN {printf \"%.3f\", ${PLUS}/${MINUS}}")
    else
        RATIO="NA"
    fi
    echo "  Strand ratio (post-WASP): ${RATIO}"

    cat > "${QC_DIR}/wasp_stats_${SAMPLE}.tsv" << EOF
sample	input	no_snp	snp_overlap	snp_remapped	snp_kept	snp_removed	snp_removed_pct	final	total_removed_pct	strand_ratio
${SAMPLE}	${INPUT_COUNT}	${KEEP_COUNT}	${REMAP_COUNT}	${REMAP_ALIGNED}	${REMAP_KEEP_COUNT}	${REMAP_REMOVED}	${REMAP_REMOVED_PCT}	${FINAL_COUNT}	${TOTAL_REMOVED_PCT}	${RATIO}
EOF
done

step_header "Step 06 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output:"
echo "  WASP-filtered BAMs:  ${WASP_FILTERED_DIR}/"
echo "  QC stats:            ${QC_DIR}/wasp_stats_*.tsv"
echo ""
echo "Next: bash 07_allele_specific.sh"
