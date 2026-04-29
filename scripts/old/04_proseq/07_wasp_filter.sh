#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 07_wasp_filter.sh — WASP mapping bias correction
# =============================================================================
#
# WHAT THIS DOES (per sample):
#   1. find_intersecting_snps.py:
#      - Identifies reads overlapping het SNPs
#      - For each such read, creates a version with the allele swapped
#      - Outputs: reads NOT overlapping SNPs (.keep.bam),
#                 reads to remap (.remap.fq.gz), and the originals (.to.remap.bam)
#
#   2. Bowtie2 remapping:
#      - Re-aligns the allele-swapped reads using the SAME parameters
#        as the original alignment (--very-sensitive)
#      - This tests whether the read would map to the same position
#        regardless of which allele is present
#
#   3. filter_remapped_reads.py:
#      - Compares the original and remapped positions
#      - Keeps only reads that mapped to the SAME position after allele swap
#      - Removes reads where mapping was allele-dependent (= mapping bias)
#
#   4. Merge:
#      - Combines non-SNP-overlapping reads (.keep.bam) with the
#        bias-filtered SNP-overlapping reads (.remap.keep.bam)
#      - Result: a WASP-filtered BAM with no mapping bias at het SNPs
#
# WHY THIS ORDER MATTERS:
#   We must remap with IDENTICAL bowtie2 parameters. If the aligner settings
#   differ between original and remapped alignment, reads might map differently
#   for reasons unrelated to the allele — invalidating the WASP test.
#
# PARALLEL EXECUTION:
#   sbatch --array=0-3 --partition=medium --time=8:00:00 --mem=32G --cpus-per-task=8 \
#       -o logs/07_%A_%a.out -e logs/07_%A_%a.err 07_wasp_filter.sh
# =============================================================================

source "config_wasp.sh"

step_header "PRO-seq Pipeline Step 07: WASP Mapping Bias Correction"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] Python: $(python --version)"
echo "[INFO] pysam:  $(python -c 'import pysam; print(pysam.__version__)')"
echo "[INFO] bowtie2: $(bowtie2 --version 2>&1 | head -1)"

create_wasp_dirs
resolve_samples "${1:-}"

# Verify WASP scripts
for F in "${WASP_FIND}" "${WASP_FILTER}"; do
    if [[ ! -f "${F}" ]]; then
        echo "[ERROR] WASP script not found: ${F}"
        echo "[ERROR] Run 06_wasp_setup.sh first."
        exit 1
    fi
done

# Verify SNP directory
if ! ls "${WASP_SNP_DIR}"/chr*.snps.txt.gz 1>/dev/null 2>&1; then
    echo "[ERROR] No SNP files in ${WASP_SNP_DIR}"
    echo "[ERROR] Run 06_wasp_setup.sh first."
    exit 1
fi

for SAMPLE in "${RUN_SAMPLES[@]}"; do
    step_header "WASP filtering: ${SAMPLE}"

    INPUT_BAM="${BAM_MM39_DIR}/${SAMPLE}_final.bam"
    WORKDIR="${WASP_INTERMEDIATE_DIR}/${SAMPLE}"
    FINAL_WASP_BAM="${WASP_FILTERED_DIR}/${SAMPLE}_wasp.bam"

    if [[ -f "${FINAL_WASP_BAM}" ]]; then
        echo "[INFO] WASP-filtered BAM already exists: ${FINAL_WASP_BAM}"
        echo "[INFO] Skipping. Delete to re-run."
        continue
    fi

    if [[ ! -f "${INPUT_BAM}" ]]; then
        echo "[ERROR] Input BAM not found: ${INPUT_BAM}"
        exit 1
    fi

    mkdir -p "${WORKDIR}"

    INPUT_COUNT=$(samtools view -c "${INPUT_BAM}")
    echo "[INFO] Input BAM: ${INPUT_COUNT} reads"

    # =========================================================================
    # STEP 1: find_intersecting_snps.py
    # =========================================================================
    echo ""
    echo "[STEP 1] Finding reads that intersect het SNPs..."

    # WASP prefixes output with the input filename.
    # We'll work in WORKDIR with a clean copy/link.
    WASP_INPUT="${WORKDIR}/${SAMPLE}.bam"
    if [[ ! -f "${WASP_INPUT}" ]]; then
        ln -sf "${INPUT_BAM}" "${WASP_INPUT}"
        ln -sf "${INPUT_BAM}.bai" "${WASP_INPUT}.bai"
    fi

    KEEP_BAM="${WORKDIR}/${SAMPLE}.keep.bam"
    TO_REMAP_BAM="${WORKDIR}/${SAMPLE}.to.remap.bam"
    REMAP_FQ="${WORKDIR}/${SAMPLE}.remap.fq.gz"
    REMAP_NUM="${WORKDIR}/${SAMPLE}.to.remap.num.gz"

    if [[ -f "${KEEP_BAM}" && -f "${REMAP_FQ}" ]]; then
        echo "[INFO] find_intersecting_snps output exists. Skipping."
    else
        python "${WASP_FIND}" \
            --is_sorted \
            --output_dir "${WORKDIR}" \
            --snp_dir "${WASP_SNP_DIR}" \
            "${WASP_INPUT}"

        echo "[INFO] find_intersecting_snps.py complete."
    fi

    # QC: how many reads overlap SNPs?
    KEEP_COUNT=$(samtools view -c "${KEEP_BAM}")
    # Use the .to.remap.bam for the true read count (not FASTQ, which has
    # multiple entries per read when a read overlaps >1 SNP due to
    # combinatorial allele swapping).
    REMAP_COUNT=$(samtools view -c "${TO_REMAP_BAM}")
    REMAP_FQ_ENTRIES=$(zcat "${REMAP_FQ}" | awk 'NR%4==1' | wc -l)
    echo "[SANITY CHECK] find_intersecting_snps:"
    echo "  Reads NOT overlapping SNPs (keep):    ${KEEP_COUNT}"
    echo "  Reads overlapping SNPs (to remap):    ${REMAP_COUNT}"
    echo "  FASTQ entries for remapping:          ${REMAP_FQ_ENTRIES}"
    echo "    (higher than read count because reads at multiple SNPs"
    echo "     generate multiple allele-swapped versions)"
    echo "  Fraction overlapping SNPs:            $(awk "BEGIN {printf \"%.2f%%\", ${REMAP_COUNT}/(${KEEP_COUNT}+${REMAP_COUNT})*100}")"

    # =========================================================================
    # STEP 2: Re-map allele-swapped reads with IDENTICAL bowtie2 parameters
    # =========================================================================
    echo ""
    echo "[STEP 2] Re-mapping allele-swapped reads..."

    REMAP_BAM="${WORKDIR}/${SAMPLE}.remap.bam"

    if [[ -f "${REMAP_BAM}" ]]; then
        echo "[INFO] Remapped BAM exists. Skipping."
    else
        # CRITICAL: Use the EXACT SAME alignment parameters as step 03.
        # --very-sensitive, single-end (-U), same genome index.
        # Any difference could cause reads to map differently for reasons
        # OTHER than the allele, which would break WASP's logic.
        bowtie2 \
            --very-sensitive \
            --no-unal \
            -x "${MM39_BT2_IDX}" \
            -U "${REMAP_FQ}" \
            -p "${THREADS}" \
            2> "${WASP_QC_DIR}/${SAMPLE}_remap_bowtie2.log" \
        | samtools view -bS -@ 2 - \
        | samtools sort -@ "${THREADS}" -o "${REMAP_BAM}"

        samtools index "${REMAP_BAM}"
        echo "[INFO] Remapping complete."
    fi

    REMAP_ALIGNED=$(samtools view -c "${REMAP_BAM}")
    echo "[SANITY CHECK] Remapping:"
    echo "  Input reads to remap:   ${REMAP_COUNT}"
    echo "  Remapped successfully:  ${REMAP_ALIGNED}"

    # =========================================================================
    # STEP 3: filter_remapped_reads.py
    # =========================================================================
    echo ""
    echo "[STEP 3] Filtering remapped reads (keep only same-position maps)..."

    REMAP_KEEP="${WORKDIR}/${SAMPLE}.remap.keep.bam"

    if [[ -f "${REMAP_KEEP}" ]]; then
        echo "[INFO] Filtered remap BAM exists. Skipping."
    else
        # filter_remapped_reads.py takes exactly 3 positional args:
        #   to_remap_bam, remap_bam, keep_bam
        # It finds the .to.remap.num.gz file automatically by naming convention.
        python "${WASP_FILTER}" \
            "${TO_REMAP_BAM}" \
            "${REMAP_BAM}" \
            "${REMAP_KEEP}"

        echo "[INFO] filter_remapped_reads.py complete."
    fi

    REMAP_KEEP_COUNT=$(samtools view -c "${REMAP_KEEP}")
    REMAP_REMOVED=$((REMAP_COUNT - REMAP_KEEP_COUNT))
    REMAP_REMOVED_PCT=$(awk "BEGIN {printf \"%.2f\", ${REMAP_REMOVED}/${REMAP_COUNT:-1}*100}")

    echo "[SANITY CHECK] WASP filtering at SNP-overlapping reads:"
    echo "  SNP-overlapping reads input:  ${REMAP_COUNT}"
    echo "  Passed WASP filter (same pos): ${REMAP_KEEP_COUNT}"
    echo "  Removed (mapping bias):        ${REMAP_REMOVED} (${REMAP_REMOVED_PCT}%)"
    echo ""
    echo "[INTERPRETATION] Mapping bias removal rates:"
    echo "  5-15%:  Typical for mouse F1 hybrids"
    echo "  >25%:   High — may indicate highly divergent regions"
    echo "  <5%:    Very low — good, little mapping bias"

    # =========================================================================
    # STEP 4: Merge keep + remap.keep → final WASP-filtered BAM
    # =========================================================================
    echo ""
    echo "[STEP 4] Merging into final WASP-filtered BAM..."

    # Sort remap.keep.bam first (may not be coordinate-sorted)
    REMAP_KEEP_SORTED="${WORKDIR}/${SAMPLE}.remap.keep.sorted.bam"
    samtools sort -@ "${THREADS}" -o "${REMAP_KEEP_SORTED}" "${REMAP_KEEP}"
    samtools index "${REMAP_KEEP_SORTED}"

    # Merge the two BAMs
    samtools merge -@ "${THREADS}" -f "${WORKDIR}/${SAMPLE}_merged_unsorted.bam" \
        "${KEEP_BAM}" "${REMAP_KEEP_SORTED}"

    # Sort and index
    samtools sort -@ "${THREADS}" -o "${FINAL_WASP_BAM}" \
        "${WORKDIR}/${SAMPLE}_merged_unsorted.bam"
    samtools index "${FINAL_WASP_BAM}"

    # Clean up large intermediate
    rm -f "${WORKDIR}/${SAMPLE}_merged_unsorted.bam"

    FINAL_COUNT=$(samtools view -c "${FINAL_WASP_BAM}")
    TOTAL_REMOVED=$((INPUT_COUNT - FINAL_COUNT))
    TOTAL_REMOVED_PCT=$(awk "BEGIN {printf \"%.2f\", ${TOTAL_REMOVED}/${INPUT_COUNT}*100}")

    echo ""
    echo "[SANITY CHECK] Final WASP-filtered BAM: ${SAMPLE}"
    echo "  Input (pre-WASP):     ${INPUT_COUNT}"
    echo "  Output (post-WASP):   ${FINAL_COUNT}"
    echo "  Total removed:        ${TOTAL_REMOVED} (${TOTAL_REMOVED_PCT}%)"
    echo "  Reads not at SNPs:    ${KEEP_COUNT}"
    echo "  Reads at SNPs kept:   ${REMAP_KEEP_COUNT}"
    echo "  Reads at SNPs removed: ${REMAP_REMOVED}"

    # Strand check
    PLUS=$(samtools view -c -F 16 "${FINAL_WASP_BAM}")
    MINUS=$(samtools view -c -f 16 "${FINAL_WASP_BAM}")
    RATIO=$(awk "BEGIN {printf \"%.3f\", ${PLUS}/${MINUS:-1}}")
    echo "  Strand ratio: ${RATIO} (should still be ~1.0)"

    # Write per-sample QC
    cat > "${WASP_QC_DIR}/wasp_stats_${SAMPLE}.tsv" << EOF
sample	input_reads	no_snp_reads	snp_reads	snp_passed_wasp	snp_removed	snp_removed_pct	final_reads	total_removed_pct	strand_ratio
${SAMPLE}	${INPUT_COUNT}	${KEEP_COUNT}	${REMAP_COUNT}	${REMAP_KEEP_COUNT}	${REMAP_REMOVED}	${REMAP_REMOVED_PCT}	${FINAL_COUNT}	${TOTAL_REMOVED_PCT}	${RATIO}
EOF

done

step_header "Step 07 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output:"
echo "  WASP-filtered BAMs: ${WASP_FILTERED_DIR}/"
echo "  QC stats:           ${WASP_QC_DIR}/"
echo ""
echo "Next: Run 08_allele_specific.sh (after all step 07 array tasks finish)"