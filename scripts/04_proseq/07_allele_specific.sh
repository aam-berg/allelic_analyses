#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04_proseq/07_allele_specific.sh — Split WASP BAMs by allele (SE)
# =============================================================================
#
# WHAT THIS DOES:
#   For each sample, calls lib/split_by_allele.py in single-end mode on the
#   WASP-filtered BAM. Produces:
#       {SAMPLE}_ref.bam        reads supporting REF (129S1) allele
#       {SAMPLE}_alt.bam        reads supporting ALT (CAST) allele
#       {SAMPLE}_nosnp.bam      reads not overlapping any het SNP
#       {SAMPLE}_ambiguous.bam  reads with conflicting allele evidence
#       {SAMPLE}_allele_counts.tsv   per-SNP allele count summary
#
#   Sorts and indexes the ref/alt BAMs (used by step 08 for bigWigs).
#   Bigwig generation is in step 08 (separate to keep this stage cheap).
#
# CHANGES FROM PREVIOUS VERSION (was 08_allele_specific.sh):
#   - split_by_allele.py is now in lib/ (shared with 03_atac/) and called
#     with --single-end.
#   - BIGWIG GENERATION MOVED TO STEP 08. Previously this step did both;
#     splitting it apart makes the step lighter, parallelizable, and matches
#     03_atac/ structure.
#   - Fixed bug: previous version had a printf line that always printed
#     "SNPs with >=5 reads: 0" (literal 0, ignoring data), then printed the
#     same statistic correctly. Removed the literal 0 line.
#
# PARALLEL EXECUTION:
#   sbatch --array=0-3 --partition=short --time=4:00:00 --mem=16G \
#       --cpus-per-task=8 -o logs/07_%A_%a.out -e logs/07_%A_%a.err \
#       07_allele_specific.sh
# =============================================================================

source "config.sh"

step_header "PRO-seq Pipeline Step 07: Allele-Specific Splitting (SE)"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] python:   $(python --version 2>&1)"
echo "[INFO] pysam:    $(python -c 'import pysam; print(pysam.__version__)')"
echo "[INFO] script:   ${SPLIT_BY_ALLELE_PY}"

create_dirs
resolve_samples "${1:-}"

[[ -f "${SPLIT_BY_ALLELE_PY}" ]] || \
    { echo "[ERROR] split_by_allele.py not found at: ${SPLIT_BY_ALLELE_PY}" >&2; exit 1; }

for SAMPLE in "${RUN_SAMPLES[@]}"; do
    step_header "Splitting ${SAMPLE} by allele"

    WASP_BAM="${WASP_FILTERED_DIR}/${SAMPLE}_wasp.bam"
    OUT_PREFIX="${ALLELE_DIR}/${SAMPLE}"
    REF_BAM="${OUT_PREFIX}_ref.bam"
    ALT_BAM="${OUT_PREFIX}_alt.bam"

    if [[ -f "${REF_BAM}.bai" && -f "${ALT_BAM}.bai" ]]; then
        echo "[INFO] Allele-split BAMs (sorted + indexed) exist. Skipping."
        continue
    fi
    [[ -f "${WASP_BAM}" ]] || \
        { echo "[ERROR] WASP BAM missing: ${WASP_BAM}" >&2; exit 1; }

    # ---- Run the splitter ----
    echo "[INFO] Running split_by_allele.py --single-end..."
    python "${SPLIT_BY_ALLELE_PY}" \
        --bam "${WASP_BAM}" \
        --snp_dir "${WASP_SNP_DIR}" \
        --output_prefix "${OUT_PREFIX}" \
        --single-end

    # ---- Sort and index each output BAM ----
    echo ""
    echo "[INFO] Sorting and indexing allele BAMs..."
    for ALLELE in ref alt nosnp ambiguous; do
        BAM="${OUT_PREFIX}_${ALLELE}.bam"
        if [[ ! -f "${BAM}" ]]; then
            echo "  [WARN] ${ALLELE}.bam not produced (probably zero reads of this class)"
            continue
        fi
        TMP="${OUT_PREFIX}_${ALLELE}.sorted.bam"
        samtools sort -@ "${THREADS}" -o "${TMP}" "${BAM}"
        mv "${TMP}" "${BAM}"
        samtools index "${BAM}"
        COUNT=$(samtools view -c "${BAM}")
        printf "  %-10s  %12d reads\n" "${ALLELE}:" "${COUNT}"
    done

    # ---- Sanity checks ----
    REF_C=$(samtools view -c "${REF_BAM}" 2>/dev/null || echo 0)
    ALT_C=$(samtools view -c "${ALT_BAM}" 2>/dev/null || echo 0)
    NOSNP_C=$(samtools view -c "${OUT_PREFIX}_nosnp.bam" 2>/dev/null || echo 0)
    AMB_C=$(samtools view -c "${OUT_PREFIX}_ambiguous.bam" 2>/dev/null || echo 0)
    TOTAL=$((REF_C + ALT_C + NOSNP_C + AMB_C))
    ALLELIC=$((REF_C + ALT_C))

    echo ""
    echo "[SANITY CHECK] ${SAMPLE}:"
    if (( TOTAL > 0 )); then
        echo "  Class fractions:"
        for pair in "ref:${REF_C}" "alt:${ALT_C}" "nosnp:${NOSNP_C}" "ambiguous:${AMB_C}"; do
            label=${pair%:*}; val=${pair#*:}
            pct=$(awk "BEGIN {printf \"%.2f\", ${val}/${TOTAL}*100}")
            printf "    %-10s  %12d  (%s%%)\n" "${label}" "${val}" "${pct}"
        done
    fi

    REF_FRAC_OF_ALLELIC="NA"
    if (( ALLELIC > 0 )); then
        REF_FRAC_OF_ALLELIC=$(awk "BEGIN {printf \"%.2f\", ${REF_C}/${ALLELIC}*100}")
        echo ""
        echo "  Among allelic reads (ref+alt only):"
        echo "    ref fraction: ${REF_FRAC_OF_ALLELIC}%   (expect ~50% after WASP)"

        DEV=$(awk "BEGIN {x=${REF_FRAC_OF_ALLELIC}-50; if(x<0)x=-x; print x}")
        if (( $(awk "BEGIN {print (${DEV}>5)?1:0}") )); then
            echo "  [WARNING] Ref fraction deviates from 50% by ${DEV} pp."
        else
            echo "  [OK] Ref fraction within 5pp of 50%."
        fi
    fi

    # ---- SNP coverage distribution (correct version; no literal 0 bug) ----
    COUNTS_FILE="${OUT_PREFIX}_allele_counts.tsv"
    if [[ -f "${COUNTS_FILE}" ]]; then
        TOTAL_SNPS_COVERED=$(awk 'NR>1 && ($5+$6) > 0' "${COUNTS_FILE}" | wc -l)
        SNPS_GE5=$(awk 'NR>1 && ($5+$6) >= 5' "${COUNTS_FILE}" | wc -l)
        SNPS_GE10=$(awk 'NR>1 && ($5+$6) >= 10' "${COUNTS_FILE}" | wc -l)
        SNPS_GE20=$(awk 'NR>1 && ($5+$6) >= 20' "${COUNTS_FILE}" | wc -l)

        echo ""
        echo "  Coverage at het SNPs:"
        awk 'NR>1 {total=$5+$6; if(total>=1) print total}' "${COUNTS_FILE}" \
        | sort -n \
        | awk '
            BEGIN {n=0}
            {vals[n++]=$1}
            END {
                if (n==0) {print "    No covered SNPs"; exit}
                printf "    Min:    %d\n", vals[0]
                printf "    Q1:     %d\n", vals[int(n*0.25)]
                printf "    Median: %d\n", vals[int(n*0.5)]
                printf "    Q3:     %d\n", vals[int(n*0.75)]
                printf "    Max:    %d\n", vals[n-1]
            }
        '
        echo "    SNPs with coverage:    ${TOTAL_SNPS_COVERED}"
        echo "    SNPs with >=5 reads:   ${SNPS_GE5}"
        echo "    SNPs with >=10 reads:  ${SNPS_GE10}"
        echo "    SNPs with >=20 reads:  ${SNPS_GE20}"
    fi

    # ---- Per-sample QC ----
    cat > "${QC_DIR}/allele_split_${SAMPLE}.tsv" << EOF
sample	total	ref	alt	nosnp	ambiguous	ref_fraction_of_allelic	snps_ge5	snps_ge10	snps_ge20
${SAMPLE}	${TOTAL}	${REF_C}	${ALT_C}	${NOSNP_C}	${AMB_C}	${REF_FRAC_OF_ALLELIC}	${SNPS_GE5:-0}	${SNPS_GE10:-0}	${SNPS_GE20:-0}
EOF

done

step_header "Step 07 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output:"
echo "  Allele BAMs:           ${ALLELE_DIR}/{SAMPLE}_{ref,alt,nosnp,ambiguous}.bam"
echo "  Per-SNP allele counts: ${ALLELE_DIR}/{SAMPLE}_allele_counts.tsv"
echo "  QC stats:              ${QC_DIR}/allele_split_*.tsv"
echo ""
echo "Next: bash 08_make_bigwigs_allele.sh"
