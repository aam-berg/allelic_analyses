#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 07_allele_specific.sh — Split WASP-filtered BAMs by allele
# =============================================================================
#
# WHAT THIS DOES:
#   For each sample, calls lib/split_by_allele.py in paired-end mode on the
#   WASP-filtered BAM. Produces:
#       {SAMPLE}_ref.bam        fragments where both mates support REF allele(s)
#       {SAMPLE}_alt.bam        fragments where both mates support ALT allele(s)
#       {SAMPLE}_nosnp.bam      fragments not overlapping any het SNP
#       {SAMPLE}_ambiguous.bam  fragments with conflicting allele evidence
#       {SAMPLE}_allele_counts.tsv   per-SNP allele count summary
#
# WHY PAIRED-END MODE:
#   In split_by_allele.py PE mode, both mates of a fragment are buffered and
#   the fragment is classified based on the union of their SNP evidence. Then
#   BOTH mates of the fragment are written to the same output BAM. This:
#     - Avoids double-counting (one fragment = one observation, not two reads)
#     - Keeps output BAMs properly paired (downstream tools won't choke)
#     - Resolves the case where R1 covers SNP A and R2 covers SNP B (the
#       fragment contributes evidence at both SNPs, classified as ref/alt
#       only if all of that evidence is consistent)
#
# WHY THIS RUNS AFTER WASP:
#   WASP removes reads that show mapping bias to one allele. After WASP, the
#   remaining reads at SNP-overlapping positions reflect true allelic balance
#   (at least, with no mapping confound). Splitting by allele after WASP gives
#   us a clean ref-only and alt-only set for downstream quantification.
#
# EXPECTED OUTPUT BALANCE:
#   For F1 hybrid het SNPs after WASP, we expect roughly 50/50 ref:alt
#   genome-wide. The script reports the global ref fraction; warnings appear
#   if it deviates from 50% by more than 5 percentage points (as a
#   sanity check on WASP behavior).
#
# PARALLEL EXECUTION:
#   sbatch --array=0-3 --partition=short --time=4:00:00 --mem=16G \
#       --cpus-per-task=8 -o logs/07_%A_%a.out -e logs/07_%A_%a.err \
#       07_allele_specific.sh
# =============================================================================

source "config.sh"

step_header "ATAC Pipeline Step 07: Allele-Specific Splitting"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] python:   $(python --version 2>&1)"
echo "[INFO] pysam:    $(python -c 'import pysam; print(pysam.__version__)')"
echo "[INFO] samtools: $(samtools --version | head -1)"
echo "[INFO] script:   ${SPLIT_BY_ALLELE_PY}"

create_dirs
resolve_samples "${1:-}"

# Verify the shared script exists
if [[ ! -f "${SPLIT_BY_ALLELE_PY}" ]]; then
    echo "[ERROR] split_by_allele.py not found at: ${SPLIT_BY_ALLELE_PY}" >&2
    echo "[ERROR] The lib/ directory should be at: ${LIB_DIR}" >&2
    exit 1
fi

for SAMPLE in "${RUN_SAMPLES[@]}"; do
    step_header "Splitting ${SAMPLE} by allele"

    WASP_BAM="${WASP_FILTERED_DIR}/${SAMPLE}_wasp.bam"
    OUT_PREFIX="${ALLELE_DIR}/${SAMPLE}"

    REF_BAM="${OUT_PREFIX}_ref.bam"
    ALT_BAM="${OUT_PREFIX}_alt.bam"

    if [[ -f "${REF_BAM}" && -f "${ALT_BAM}" \
       && -f "${REF_BAM}.bai" && -f "${ALT_BAM}.bai" ]]; then
        echo "[INFO] Allele-split BAMs (sorted+indexed) exist for ${SAMPLE}. Skipping."
        continue
    fi

    if [[ ! -f "${WASP_BAM}" ]]; then
        echo "[ERROR] WASP BAM missing for ${SAMPLE}: ${WASP_BAM}" >&2
        echo "[ERROR] Run 06_wasp_filter.sh first." >&2
        exit 1
    fi

    # ---- Run the splitter ----
    echo "[INFO] Running split_by_allele.py --paired-end..."
    python "${SPLIT_BY_ALLELE_PY}" \
        --bam "${WASP_BAM}" \
        --snp_dir "${WASP_SNP_DIR}" \
        --output_prefix "${OUT_PREFIX}" \
        --paired-end

    # ---- Sort + index each output BAM ----
    # split_by_allele writes BAMs in input order (which was coord-sorted, so
    # output is already roughly sorted), but we re-sort to guarantee correct
    # ordering for samtools index.
    echo ""
    echo "[INFO] Sorting and indexing allele BAMs..."
    for ALLELE in ref alt nosnp ambiguous; do
        BAM="${OUT_PREFIX}_${ALLELE}.bam"
        if [[ ! -f "${BAM}" ]]; then
            echo "  [WARN] ${ALLELE}.bam not produced (this is normal if zero reads of this class)"
            continue
        fi
        TMP_SORTED="${OUT_PREFIX}_${ALLELE}.sorted.bam"
        samtools sort -@ "${THREADS}" -o "${TMP_SORTED}" "${BAM}"
        mv "${TMP_SORTED}" "${BAM}"
        samtools index "${BAM}"
        COUNT=$(samtools view -c "${BAM}")
        printf "  %-10s  %12d reads\n" "${ALLELE}:" "${COUNT}"
    done

    # ---- Sanity check: ref/alt balance after WASP ----
    REF_COUNT=$(samtools view -c "${REF_BAM}" 2>/dev/null || echo 0)
    ALT_COUNT=$(samtools view -c "${ALT_BAM}" 2>/dev/null || echo 0)
    NOSNP_COUNT=$(samtools view -c "${OUT_PREFIX}_nosnp.bam" 2>/dev/null || echo 0)
    AMB_COUNT=$(samtools view -c "${OUT_PREFIX}_ambiguous.bam" 2>/dev/null || echo 0)

    TOTAL=$((REF_COUNT + ALT_COUNT + NOSNP_COUNT + AMB_COUNT))
    ALLELIC=$((REF_COUNT + ALT_COUNT))

    echo ""
    echo "[SANITY CHECK] ${SAMPLE}:"
    if (( TOTAL > 0 )); then
        REF_FRAC_TOT=$(awk "BEGIN {printf \"%.2f\", ${REF_COUNT}/${TOTAL}*100}")
        ALT_FRAC_TOT=$(awk "BEGIN {printf \"%.2f\", ${ALT_COUNT}/${TOTAL}*100}")
        NOSNP_FRAC=$(awk "BEGIN {printf \"%.2f\", ${NOSNP_COUNT}/${TOTAL}*100}")
        AMB_FRAC=$(awk "BEGIN {printf \"%.2f\", ${AMB_COUNT}/${TOTAL}*100}")
        echo "  Class fractions:"
        printf "    ref:        %12d  (%s%%)\n" "${REF_COUNT}"   "${REF_FRAC_TOT}"
        printf "    alt:        %12d  (%s%%)\n" "${ALT_COUNT}"   "${ALT_FRAC_TOT}"
        printf "    nosnp:      %12d  (%s%%)\n" "${NOSNP_COUNT}" "${NOSNP_FRAC}"
        printf "    ambiguous:  %12d  (%s%%)\n" "${AMB_COUNT}"   "${AMB_FRAC}"
    fi

    REF_FRAC_OF_ALLELIC="NA"
    if (( ALLELIC > 0 )); then
        REF_FRAC_OF_ALLELIC=$(awk "BEGIN {printf \"%.2f\", ${REF_COUNT}/${ALLELIC}*100}")
        echo ""
        echo "  Among allelic reads (ref+alt only):"
        echo "    ref fraction: ${REF_FRAC_OF_ALLELIC}%   (expect ~50% after WASP)"

        DEV=$(awk "BEGIN {x=${REF_FRAC_OF_ALLELIC}-50; if(x<0)x=-x; print x}")
        if (( $(awk "BEGIN {print (${DEV}>5)?1:0}") )); then
            echo "  [WARNING] Ref fraction deviates from 50% by ${DEV} percentage points."
            echo "             Investigate: WASP working correctly? heterozygous SNPs verified?"
        else
            echo "  [OK] Ref fraction within 5pp of 50%."
        fi
    fi

    # ---- Per-sample QC ----
    cat > "${QC_DIR}/allele_split_${SAMPLE}.tsv" << EOF
sample	total	ref	alt	nosnp	ambiguous	ref_fraction_of_allelic
${SAMPLE}	${TOTAL}	${REF_COUNT}	${ALT_COUNT}	${NOSNP_COUNT}	${AMB_COUNT}	${REF_FRAC_OF_ALLELIC}
EOF

    echo ""
    echo "  Per-SNP allele counts: ${OUT_PREFIX}_allele_counts.tsv"
done

step_header "Step 07 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output:"
echo "  Allele-split BAMs:    ${ALLELE_DIR}/"
echo "  Per-SNP allele counts: ${ALLELE_DIR}/{SAMPLE}_allele_counts.tsv"
echo "  QC stats:              ${QC_DIR}/allele_split_*.tsv"
echo ""
echo "Next: bash 08_make_bigwigs.sh"
