#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 09_merge_qc.sh — Merge across replicates + compile final QC report
# =============================================================================
#
# WHAT THIS DOES (single non-array job, runs once at end):
#   1. Merges per-replicate final BAMs into a combined merged_final.bam
#      and a corresponding merged.bw bigWig.
#   2. Merges per-replicate ref BAMs and alt BAMs into merged_ref.bam +
#      merged_alt.bam, with corresponding bigWigs.
#   3. Compiles all per-sample QC TSVs (read counts, trimming, alignment,
#      WASP, allele-split, peak counts) into a single human-readable QC
#      report at ${QC_DIR}/pipeline_qc_report.txt.
#
# WHY THIS IS A SEPARATE STEP, NOT PART OF 08:
#   Steps 01-08 are per-sample (parallel-safe SLURM arrays). Merging across
#   samples is inherently a single-job operation that must wait for all
#   per-sample work to complete. Splitting it out cleanly avoids race
#   conditions and makes the dependency graph explicit:
#     08 (array, depends on 07) -> 09 (single, depends on 08).
#
# WHEN TO USE THE MERGED BIGWIGS:
#   The merged signal has 4x the depth of any individual replicate, so it's
#   the cleanest track for browser inspection and motif-aggregate plots in
#   06_allele_pairs/ and 07_analysis/. Per-replicate tracks remain useful
#   for rep-level reproducibility checks.
#
# SBATCH RESOURCES:
#   sbatch --partition=short --time=2:00:00 --mem=16G --cpus-per-task=4 \
#       -o logs/09_%j.out -e logs/09_%j.err 09_merge_qc.sh
# =============================================================================

source "config.sh"

step_header "ATAC Pipeline Step 09: Merge + Final QC"

source activate "${CONDA_ENV_NAME}"
create_dirs

# -----------------------------------------------------------------------------
# Helper: BAM -> fragment-coverage bigWig (same as step 08)
# -----------------------------------------------------------------------------
make_atac_bigwig() {
    local BAM="$1"
    local OUTBW="$2"
    local LABEL="$3"

    if [[ -f "${OUTBW}" ]]; then
        echo "  [SKIP] ${LABEL}: bigWig already exists."
        return 0
    fi
    if [[ ! -f "${BAM}" ]]; then
        echo "  [SKIP] ${LABEL}: input BAM missing (${BAM})"
        return 0
    fi

    local TMPBG="${OUTBW%.bw}.bedGraph"
    bedtools genomecov -ibam "${BAM}" -bg -pc \
    | sort -k1,1 -k2,2n \
    | awk 'NR==FNR {chr[$1]=1; next} ($1 in chr)' "${MM39_CHROM_SIZES}" - \
    > "${TMPBG}"
    bedGraphToBigWig "${TMPBG}" "${MM39_CHROM_SIZES}" "${OUTBW}"
    rm -f "${TMPBG}"
    echo "  [OK]   ${LABEL}: ${OUTBW}"
}

# -----------------------------------------------------------------------------
# 1. Merged final BAM + bigWig (across all replicates)
# -----------------------------------------------------------------------------
step_header "Merging final BAMs across replicates"

MERGED_FINAL_BAM="${BAM_FINAL_DIR}/merged_final.bam"
MERGED_FINAL_BW="${BIGWIG_MERGED_DIR}/merged.bw"

if [[ -f "${MERGED_FINAL_BAM}" && -f "${MERGED_FINAL_BAM}.bai" ]]; then
    echo "[INFO] merged_final.bam exists. Skipping merge."
else
    INPUT_BAMS=()
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        BAM="${BAM_FINAL_DIR}/${SAMPLE}_final.bam"
        if [[ ! -f "${BAM}" ]]; then
            echo "[ERROR] Missing per-rep final BAM: ${BAM}" >&2
            exit 1
        fi
        INPUT_BAMS+=("${BAM}")
    done

    echo "[INFO] Merging ${#INPUT_BAMS[@]} BAMs..."
    samtools merge -@ "${THREADS}" -f "${MERGED_FINAL_BAM}" "${INPUT_BAMS[@]}"
    samtools index "${MERGED_FINAL_BAM}"
    echo "[INFO] Indexed: ${MERGED_FINAL_BAM}"
fi

make_atac_bigwig "${MERGED_FINAL_BAM}" "${MERGED_FINAL_BW}" "merged (all reps)"

# -----------------------------------------------------------------------------
# 2. Merged allele-specific BAMs + bigWigs
# -----------------------------------------------------------------------------
step_header "Merging allele-specific BAMs"

for ALLELE in ref alt; do
    INPUT_BAMS=()
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        BAM="${ALLELE_DIR}/${SAMPLE}_${ALLELE}.bam"
        [[ -f "${BAM}" ]] && INPUT_BAMS+=("${BAM}")
    done

    if [[ ${#INPUT_BAMS[@]} -eq 0 ]]; then
        echo "[WARN] No per-rep ${ALLELE}.bam files found; skipping merge for ${ALLELE}."
        continue
    fi

    MERGED_ALLELE_BAM="${ALLELE_DIR}/merged_${ALLELE}.bam"
    MERGED_ALLELE_BW="${BIGWIG_MERGED_DIR}/merged_${ALLELE}.bw"

    if [[ -f "${MERGED_ALLELE_BAM}" && -f "${MERGED_ALLELE_BAM}.bai" ]]; then
        echo "[INFO] merged_${ALLELE}.bam exists. Skipping merge."
    else
        echo "[INFO] Merging ${#INPUT_BAMS[@]} ${ALLELE} BAMs..."
        samtools merge -@ "${THREADS}" -f "${MERGED_ALLELE_BAM}" "${INPUT_BAMS[@]}"
        samtools index "${MERGED_ALLELE_BAM}"
    fi

    make_atac_bigwig "${MERGED_ALLELE_BAM}" "${MERGED_ALLELE_BW}" "merged ${ALLELE}"
done

# Sanity check on merged ref/alt balance
if [[ -f "${ALLELE_DIR}/merged_ref.bam" && -f "${ALLELE_DIR}/merged_alt.bam" ]]; then
    REF_TOTAL=$(samtools view -c "${ALLELE_DIR}/merged_ref.bam")
    ALT_TOTAL=$(samtools view -c "${ALLELE_DIR}/merged_alt.bam")
    ALLELIC=$((REF_TOTAL + ALT_TOTAL))
    if (( ALLELIC > 0 )); then
        REF_PCT=$(awk "BEGIN {printf \"%.2f\", ${REF_TOTAL}/${ALLELIC}*100}")
        echo ""
        echo "[SANITY CHECK] Merged allelic balance:"
        printf "  ref reads:    %12d\n" "${REF_TOTAL}"
        printf "  alt reads:    %12d\n" "${ALT_TOTAL}"
        printf "  ref fraction: %s%%   (expect ~50%% after WASP)\n" "${REF_PCT}"
    fi
fi

# -----------------------------------------------------------------------------
# 3. Compile final QC report
# -----------------------------------------------------------------------------
step_header "Compiling pipeline QC report"

QC_REPORT="${QC_DIR}/pipeline_qc_report.txt"
{
    echo "============================================================"
    echo "F121-9 ATAC-seq Pipeline — QC Report"
    echo "Generated: $(date)"
    echo "Pipeline:  ${SCRIPT_DIR}"
    echo "Output:    ${SCRATCH_DIR}"
    echo "Samples:   ${SAMPLE_ORDER[*]}"
    echo "============================================================"
    echo ""

    # --- Read counts (raw) ---
    echo "--- 1. Raw read counts (from step 01) ---"
    echo -e "sample\tsrr\traw_reads"
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        F="${QC_DIR}/read_counts_raw_${SAMPLE}.tsv"
        [[ -f "${F}" ]] && cat "${F}" || echo -e "${SAMPLE}\t?\t?"
    done
    echo ""

    # --- Trimming ---
    echo "--- 2. Trimming (from step 02) ---"
    echo -e "sample\traw_pairs\ttoo_short\twritten\tpct_short\tpct_surviving"
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        F="${QC_DIR}/trimming_stats_${SAMPLE}.tsv"
        [[ -f "${F}" ]] && cat "${F}" || echo -e "${SAMPLE}\t?\t?\t?\t?\t?"
    done
    echo ""

    # --- Alignment ---
    echo "--- 3. Alignment + filtering (from step 03) ---"
    HEADER_PRINTED=0
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        F="${QC_DIR}/alignment_stats_${SAMPLE}.tsv"
        if [[ -f "${F}" ]]; then
            if (( HEADER_PRINTED == 0 )); then
                head -1 "${F}"
                HEADER_PRINTED=1
            fi
            tail -n +2 "${F}"
        fi
    done
    echo ""

    # --- Peak counts (per-rep + consensus) ---
    echo "--- 4. Peak calling (from step 04) ---"
    echo -e "sample\tn_peaks\tmedian_width_bp"
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        F="${QC_DIR}/peak_counts_${SAMPLE}.tsv"
        [[ -f "${F}" ]] && cat "${F}" || echo -e "${SAMPLE}\t?\t?"
    done
    if [[ -f "${QC_DIR}/consensus_peak_count.tsv" ]]; then
        echo ""
        echo "Consensus peak set:"
        cat "${QC_DIR}/consensus_peak_count.tsv" \
            | awk 'BEGIN{print "  set\tn_consensus\tn_per_rep_combined"}{print "  "$0}'
    fi
    echo ""

    # --- WASP ---
    echo "--- 5. WASP mapping bias correction (from step 06) ---"
    HEADER_PRINTED=0
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        F="${QC_DIR}/wasp_stats_${SAMPLE}.tsv"
        if [[ -f "${F}" ]]; then
            if (( HEADER_PRINTED == 0 )); then
                head -1 "${F}"
                HEADER_PRINTED=1
            fi
            tail -n +2 "${F}"
        fi
    done
    echo ""

    # --- Allele split ---
    echo "--- 6. Allele-specific splitting (from step 07) ---"
    HEADER_PRINTED=0
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        F="${QC_DIR}/allele_split_${SAMPLE}.tsv"
        if [[ -f "${F}" ]]; then
            if (( HEADER_PRINTED == 0 )); then
                head -1 "${F}"
                HEADER_PRINTED=1
            fi
            tail -n +2 "${F}"
        fi
    done
    echo ""

    echo "============================================================"
    echo "End of report. Source TSVs: ${QC_DIR}/"
    echo "============================================================"
} > "${QC_REPORT}"

echo "[INFO] Wrote QC report: ${QC_REPORT}"
echo ""
echo "----- Preview (first 60 lines) -----"
head -60 "${QC_REPORT}"
echo "----- ... -----"

step_header "Step 09 COMPLETE — Pipeline finished"
echo ""
echo "Final outputs:"
echo "  Merged BAMs:         ${BAM_FINAL_DIR}/merged_final.bam"
echo "                       ${ALLELE_DIR}/merged_{ref,alt}.bam"
echo "  Merged bigWigs:      ${BIGWIG_MERGED_DIR}/merged{,_ref,_alt}.bw"
echo "  Per-rep bigWigs:     ${BIGWIG_IND_DIR}/"
echo "  Allele bigWigs:      ${BIGWIG_ALLELE_DIR}/"
echo "  Consensus peaks:     ${PEAKS_CONSENSUS_DIR}/consensus_peaks.bed"
echo "  Final QC report:     ${QC_REPORT}"
echo ""
echo "Downstream consumers:"
echo "  05_motif_annot/ ← peaks/consensus/consensus_peaks.bed"
echo "  06_allele_pairs/ ← allele_specific/{SAMPLE}_{ref,alt}.bam, bigwig/allele_specific/*.bw"
