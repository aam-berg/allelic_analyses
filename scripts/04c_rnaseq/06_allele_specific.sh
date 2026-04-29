#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04c_rnaseq/06_allele_specific.sh — WASP-tag filter + PE allele split
# =============================================================================
#
# WHAT THIS DOES (per sample):
#   1. Filters the STAR BAM to keep only reads passing the WASP test:
#        - reads with vW:i:1 (passed WASP)
#        - reads with NO vW tag (don't overlap any SNP, kept by default)
#      Reads with vW:i:2-7 are removed (mapping-bias-suspect).
#   2. Calls lib/split_by_allele.py --paired-end on the filtered BAM, which
#      produces:
#          {SAMPLE}_ref.bam        fragments supporting REF (129S1) allele
#          {SAMPLE}_alt.bam        fragments supporting ALT (CAST) allele
#          {SAMPLE}_nosnp.bam      fragments not overlapping any het SNP
#          {SAMPLE}_ambiguous.bam  fragments with conflicting allele evidence
#          {SAMPLE}_allele_counts.tsv   per-SNP allele count summary
#   3. Sorts and indexes each output BAM.
#
# THIS REPLACES THE EXPLICIT find/remap/filter CYCLE used in 03_atac/ and
# 04_proseq/ — STAR did the WASP work during alignment, we just need to
# select the passing reads.
#
# WHY KEEP THE WASP-PASSED BAM AS A SEPARATE FILE:
#   The full STAR BAM is the right input for gene-level expression (step 04)
#   and strand bigWigs (step 05). The WASP-filtered BAM is the right input
#   for allele-specific work and is a direct input to lib/split_by_allele.py.
#
# WHY split_by_allele.py --paired-end:
#   PE-aware classification: both mates of a fragment must agree on allele
#   evidence, and both go to the same output BAM. This is identical to how
#   03_atac/ uses it.
#
# PARALLEL EXECUTION:
#   sbatch --array=0-1 --partition=short --time=2:00:00 --mem=16G \
#       --cpus-per-task=8 -o logs/06_%A_%a.out -e logs/06_%A_%a.err \
#       06_allele_specific.sh
# =============================================================================

source "config.sh"

step_header "RNA-seq Pipeline Step 06: WASP Filter + Allele Split (PE)"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] python:   $(python --version 2>&1)"
echo "[INFO] pysam:    $(python -c 'import pysam; print(pysam.__version__)')"
echo "[INFO] samtools: $(samtools --version | head -1)"

create_dirs
resolve_samples "${1:-}"

[[ -f "${SPLIT_BY_ALLELE_PY}" ]] || \
    { echo "[ERROR] split_by_allele.py missing: ${SPLIT_BY_ALLELE_PY}" >&2; exit 1; }

ls "${WASP_SNP_DIR}"/chr*.snps.txt.gz 1>/dev/null 2>&1 || \
    { echo "[ERROR] WASP SNP files missing. Run 00_setup_references.sh first." >&2; exit 1; }

for SAMPLE in "${RUN_SAMPLES[@]}"; do
    step_header "WASP filter + allele split: ${SAMPLE}"

    INPUT_BAM="${BAM_STAR_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    [[ -f "${INPUT_BAM}" ]] || { echo "[ERROR] STAR BAM missing: ${INPUT_BAM}" >&2; exit 1; }

    FILTERED_BAM="${BAM_FINAL_DIR}/${SAMPLE}_final.bam"

    # =========================================================================
    # STEP A — WASP filter via vW tag
    # =========================================================================
    if [[ -f "${FILTERED_BAM}" && -f "${FILTERED_BAM}.bai" ]]; then
        echo "[STEP A] WASP-filtered BAM exists. Skipping."
    else
        echo "[STEP A] Filtering by WASP tag..."

        INPUT_COUNT=$(samtools view -c "${INPUT_BAM}")

        # Keep reads where:
        #   - vW tag is absent, OR
        #   - vW:i:1 (passed WASP)
        # awk processes header (^@) unchanged, then for each alignment line
        # checks for any vW:i:N tag. Keep if absent or value is 1.
        samtools view -h "${INPUT_BAM}" \
        | awk '
            /^@/ { print; next }
            {
                wasp_pass = 1
                for (i=12; i<=NF; i++) {
                    if ($i ~ /^vW:i:/) {
                        val = substr($i, 6)
                        if (val != "1") { wasp_pass = 0 }
                        break
                    }
                }
                if (wasp_pass) print
            }' \
        | samtools view -bS -@ 2 - \
        | samtools sort -@ "${THREADS}" -o "${FILTERED_BAM}"

        samtools index "${FILTERED_BAM}"

        FILTERED_COUNT=$(samtools view -c "${FILTERED_BAM}")
        REMOVED=$((INPUT_COUNT - FILTERED_COUNT))
        if (( INPUT_COUNT > 0 )); then
            REMOVED_PCT=$(awk "BEGIN {printf \"%.2f\", ${REMOVED}/${INPUT_COUNT}*100}")
        else
            REMOVED_PCT="NA"
        fi

        echo ""
        echo "[SANITY CHECK] WASP filter:"
        printf "  Input reads:     %12d\n" "${INPUT_COUNT}"
        printf "  Kept (WASP-pass + no-tag): %12d\n" "${FILTERED_COUNT}"
        printf "  Removed:         %12d  (%s%%)\n" "${REMOVED}" "${REMOVED_PCT}"
        echo "[INTERPRETATION]"
        echo "  Removal 2-10%: typical for mouse F1 hybrids"
        echo "  Removal >15%:  high; investigate (highly divergent regions?)"
        echo "  Removal <1%:   low; mapping was largely bias-free"

        cat > "${QC_DIR}/wasp_filter_${SAMPLE}.tsv" << EOF
sample	input	kept	removed	removed_pct
${SAMPLE}	${INPUT_COUNT}	${FILTERED_COUNT}	${REMOVED}	${REMOVED_PCT}
EOF
    fi

    # =========================================================================
    # STEP B — Allele split via lib/split_by_allele.py --paired-end
    # =========================================================================
    OUT_PREFIX="${ALLELE_DIR}/${SAMPLE}"
    REF_BAM="${OUT_PREFIX}_ref.bam"
    ALT_BAM="${OUT_PREFIX}_alt.bam"

    if [[ -f "${REF_BAM}.bai" && -f "${ALT_BAM}.bai" ]]; then
        echo ""
        echo "[STEP B] Allele-split BAMs (sorted+indexed) exist. Skipping."
    else
        echo ""
        echo "[STEP B] Splitting WASP-passed BAM by allele (PE mode)..."

        python "${SPLIT_BY_ALLELE_PY}" \
            --bam "${FILTERED_BAM}" \
            --snp_dir "${WASP_SNP_DIR}" \
            --output_prefix "${OUT_PREFIX}" \
            --paired-end

        # Sort + index each output BAM
        echo ""
        echo "[INFO] Sorting and indexing allele BAMs..."
        for ALLELE in ref alt nosnp ambiguous; do
            BAM="${OUT_PREFIX}_${ALLELE}.bam"
            if [[ ! -f "${BAM}" ]]; then
                echo "  [WARN] ${ALLELE}.bam not produced (zero reads of this class)"
                continue
            fi
            TMP="${OUT_PREFIX}_${ALLELE}.sorted.bam"
            samtools sort -@ "${THREADS}" -o "${TMP}" "${BAM}"
            mv "${TMP}" "${BAM}"
            samtools index "${BAM}"
            COUNT=$(samtools view -c "${BAM}")
            printf "  %-10s  %12d reads\n" "${ALLELE}:" "${COUNT}"
        done

        # ---- Sanity check ----
        REF_C=$(samtools view -c "${REF_BAM}" 2>/dev/null || echo 0)
        ALT_C=$(samtools view -c "${ALT_BAM}" 2>/dev/null || echo 0)
        ALLELIC=$((REF_C + ALT_C))

        if (( ALLELIC > 0 )); then
            REF_PCT=$(awk "BEGIN {printf \"%.2f\", ${REF_C}/${ALLELIC}*100}")
            echo ""
            echo "  Among allelic reads (ref+alt only):"
            echo "    ref fraction: ${REF_PCT}%   (expect ~50% after WASP)"

            DEV=$(awk "BEGIN {x=${REF_PCT}-50; if(x<0)x=-x; print x}")
            if (( $(awk "BEGIN {print (${DEV}>5)?1:0}") )); then
                echo "  [WARNING] Ref fraction deviates from 50% by ${DEV} pp."
            else
                echo "  [OK] Ref fraction within 5pp of 50%."
            fi
        fi

        # Per-sample QC TSV with coverage thresholds
        COUNTS_FILE="${OUT_PREFIX}_allele_counts.tsv"
        if [[ -f "${COUNTS_FILE}" ]]; then
            SNPS_GE5=$(awk 'NR>1 && ($5+$6) >= 5' "${COUNTS_FILE}" | wc -l)
            SNPS_GE10=$(awk 'NR>1 && ($5+$6) >= 10' "${COUNTS_FILE}" | wc -l)
            SNPS_GE20=$(awk 'NR>1 && ($5+$6) >= 20' "${COUNTS_FILE}" | wc -l)
        else
            SNPS_GE5=0; SNPS_GE10=0; SNPS_GE20=0
        fi

        cat > "${QC_DIR}/allele_split_${SAMPLE}.tsv" << EOF
sample	ref	alt	allelic	ref_pct	snps_ge5	snps_ge10	snps_ge20
${SAMPLE}	${REF_C}	${ALT_C}	${ALLELIC}	${REF_PCT:-NA}	${SNPS_GE5}	${SNPS_GE10}	${SNPS_GE20}
EOF
    fi
done

step_header "Step 06 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output:"
echo "  WASP-filtered BAMs: ${BAM_FINAL_DIR}/{SAMPLE}_final.bam"
echo "  Allele-split BAMs:  ${ALLELE_DIR}/{SAMPLE}_{ref,alt,nosnp,ambiguous}.bam"
echo "  Per-SNP counts:     ${ALLELE_DIR}/{SAMPLE}_allele_counts.tsv"
echo ""
echo "Next: bash 07_make_bigwigs_allele.sh"
