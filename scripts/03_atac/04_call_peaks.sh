#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04_call_peaks.sh — MACS2 per-replicate + consensus peak calling
# =============================================================================
#
# This script has TWO modes, selected by the first argument:
#
#   per-rep    Calls MACS2 on one sample (array-friendly).
#   consensus  Builds the consensus set from all per-rep results (single job).
#
# WHY MACS2 PARAMETERS:
#   --nomodel             ATAC fragments come from Tn5 cut events; no need to
#                         build a fragment-size model from peaks
#   --shift -75           Shift each cut site 75 bp upstream...
#   --extsize 150         ...and extend 150 bp, simulating a fragment centered
#                         on the cut site. Standard ATAC narrow-peak setup.
#   -f BAMPE              Use insert (R1+R2) information; lets MACS2 ignore
#                         the shift/extsize and use real fragments instead.
#                         (When -f BAMPE is given, --shift and --extsize are
#                         ignored by MACS2; we leave them for documentation.)
#   -q 0.01               narrow-peak FDR threshold
#   --keep-dup all        we already removed dups in step 03; using "all"
#                         tells MACS2 not to dedup again
#
# HOW THE CONSENSUS IS BUILT:
#   1. Concatenate all per-rep narrowPeak files, tagging each interval with
#      its replicate name.
#   2. Sort and merge overlapping intervals across all replicates with
#      bedtools merge -c 4 -o distinct (this gives a comma-separated list
#      of contributing replicates per merged interval).
#   3. Keep merged intervals where at least MIN_REPS_FOR_CONSENSUS distinct
#      replicates contributed.
#
# This is a simpler and more reproducible alternative to IDR for our
# downstream needs (we want a stable peak set to overlap with motifs).
#
# PARALLEL EXECUTION:
#   # Per-rep (array of 4):
#   sbatch --array=0-3 --partition=short --time=2:00:00 --mem=8G \
#       --cpus-per-task=4 -o logs/04_perrep_%A_%a.out -e logs/04_perrep_%A_%a.err \
#       --wrap="bash 04_call_peaks.sh per-rep"
#
#   # Consensus (single, depends on per-rep array):
#   sbatch --partition=short --time=00:30:00 --mem=8G --cpus-per-task=2 \
#       --dependency=afterok:${PERREP_JOBID} \
#       -o logs/04_consensus_%j.out -e logs/04_consensus_%j.err \
#       --wrap="bash 04_call_peaks.sh consensus"
# =============================================================================

source "config.sh"

MODE="${1:-}"
shift || true

if [[ -z "${MODE}" ]]; then
    echo "Usage: $0 {per-rep|consensus} [sample_index_or_name]" >&2
    exit 1
fi

step_header "ATAC Pipeline Step 04: Peak Calling (mode=${MODE})"

source activate "${CONDA_ENV_NAME}"
create_dirs

# Verify MACS3 is available
if ! command -v macs3 &>/dev/null; then
    echo "[ERROR] macs3 not found. Install with:" >&2
    echo "        pip install macs3" >&2
    exit 1
fi
echo "[INFO] macs3: $(macs3 --version 2>&1)"

# =============================================================================
# MODE: per-rep
# =============================================================================
case "${MODE}" in
    per-rep)
        resolve_samples "${1:-}"

        for SAMPLE in "${RUN_SAMPLES[@]}"; do
            FINAL_BAM="${BAM_FINAL_DIR}/${SAMPLE}_final.bam"
            OUT_DIR="${PEAKS_PER_REP_DIR}/${SAMPLE}"
            PEAKS_FILE="${OUT_DIR}/${SAMPLE}_peaks.narrowPeak"

            if [[ ! -f "${FINAL_BAM}" ]]; then
                echo "[ERROR] Final BAM missing for ${SAMPLE}: ${FINAL_BAM}" >&2
                exit 1
            fi

            if [[ -f "${PEAKS_FILE}" ]]; then
                echo "[INFO] Peaks already exist for ${SAMPLE}: ${PEAKS_FILE}"
                echo "[INFO] Skipping. Delete to re-run."
                continue
            fi

            mkdir -p "${OUT_DIR}"

            echo ""
            echo "[INFO] Calling MACS2 peaks for ${SAMPLE}..."
            macs3 callpeak \
                -t "${FINAL_BAM}" \
                -f BAMPE \
                -g "${MACS2_GSIZE}" \
                -n "${SAMPLE}" \
                --outdir "${OUT_DIR}" \
                --nomodel \
                --shift -75 \
                --extsize 150 \
                -q "${MACS2_QVALUE}" \
                --keep-dup all \
                2>&1 | tee "${LOG_DIR}/${SAMPLE}_macs3.log"

            N_PEAKS=$(wc -l < "${PEAKS_FILE}")
            MEDIAN_W=$(awk '{print $3-$2}' "${PEAKS_FILE}" \
                       | sort -n \
                       | awk 'BEGIN{c=0} {a[c++]=$1} END{
                           if (c==0) {print "NA"; exit}
                           if (c%2) print a[int(c/2)]
                           else printf "%.0f", (a[c/2-1]+a[c/2])/2
                         }')

            echo ""
            echo "[SANITY CHECK] ${SAMPLE}:"
            echo "  Peaks called:    ${N_PEAKS}"
            echo "  Median width:    ${MEDIAN_W} bp"
            echo "[INTERPRETATION]"
            echo "  ATAC mESC typically: 30,000-80,000 peaks at q=0.01"
            echo "  Median width:        300-700 bp common"

            # Per-sample QC
            echo -e "${SAMPLE}\t${N_PEAKS}\t${MEDIAN_W}" \
                > "${QC_DIR}/peak_counts_${SAMPLE}.tsv"
        done
        ;;

    # =============================================================================
    # MODE: consensus
    # =============================================================================
    consensus)
        echo ""
        echo "[INFO] Building consensus peak set."
        echo "[INFO] Threshold: peak must be in at least ${MIN_REPS_FOR_CONSENSUS} of ${#SAMPLE_ORDER[@]} replicates."

        OUT_DIR="${PEAKS_CONSENSUS_DIR}"
        mkdir -p "${OUT_DIR}"

        ALL_PEAKS="${OUT_DIR}/all_peaks_concatenated.bed"
        CONSENSUS="${OUT_DIR}/consensus_peaks.bed"

        # ---- Step 1: concat all per-rep peaks, tagged with sample name ----
        : > "${ALL_PEAKS}"
        for SAMPLE in "${SAMPLE_ORDER[@]}"; do
            PEAK_FILE="${PEAKS_PER_REP_DIR}/${SAMPLE}/${SAMPLE}_peaks.narrowPeak"
            if [[ ! -f "${PEAK_FILE}" ]]; then
                echo "[ERROR] Missing per-rep peaks: ${PEAK_FILE}" >&2
                echo "[ERROR] Run 'bash 04_call_peaks.sh per-rep' first." >&2
                exit 1
            fi
            awk -v s="${SAMPLE}" 'BEGIN{OFS="\t"} {print $1, $2, $3, s}' "${PEAK_FILE}" \
                >> "${ALL_PEAKS}"
        done

        TOTAL_PER_REP=$(wc -l < "${ALL_PEAKS}")
        echo "[INFO] Combined per-rep peaks: ${TOTAL_PER_REP}"

        # ---- Step 2: sort, merge, count distinct contributing replicates ----
        sort -k1,1 -k2,2n -o "${ALL_PEAKS}" "${ALL_PEAKS}"

        # bedtools merge -c 4 -o distinct gives merged intervals with a
        # comma-separated list of contributing rep names
        # awk then counts those names and keeps intervals meeting threshold
        bedtools merge -i "${ALL_PEAKS}" -c 4 -o distinct \
        | awk -v min="${MIN_REPS_FOR_CONSENSUS}" 'BEGIN{OFS="\t"} {
            n = split($4, a, ",")
            if (n >= min) print $1, $2, $3, $4, n
          }' > "${CONSENSUS}"

        N_CONSENSUS=$(wc -l < "${CONSENSUS}")
        echo ""
        echo "[SANITY CHECK]"
        echo "  Consensus peaks (>=${MIN_REPS_FOR_CONSENSUS} reps):  ${N_CONSENSUS}"
        echo ""
        echo "  Distribution of replicate-support counts:"
        awk '{print $5}' "${CONSENSUS}" | sort | uniq -c | sed 's/^/    /'
        echo ""
        echo "  Median consensus peak width:"
        awk '{print $3-$2}' "${CONSENSUS}" \
            | sort -n \
            | awk 'BEGIN{c=0} {a[c++]=$1} END{
                if (c==0) print "  (no peaks)"
                else if (c%2) printf "    %d bp\n", a[int(c/2)]
                else printf "    %.0f bp\n", (a[c/2-1]+a[c/2])/2
              }'

        # Final QC line
        echo -e "consensus\t${N_CONSENSUS}\t${TOTAL_PER_REP}" \
            > "${QC_DIR}/consensus_peak_count.tsv"

        # Cleanup intermediate
        rm -f "${ALL_PEAKS}"
        ;;

    *)
        echo "[ERROR] Unknown mode: ${MODE}. Use 'per-rep' or 'consensus'." >&2
        exit 1
        ;;
esac

step_header "Step 04 COMPLETE (mode=${MODE})"
case "${MODE}" in
    per-rep)
        echo "Output: ${PEAKS_PER_REP_DIR}/*"
        echo "Next: bash 04_call_peaks.sh consensus  (after all per-rep tasks complete)"
        ;;
    consensus)
        echo "Output: ${PEAKS_CONSENSUS_DIR}/consensus_peaks.bed"
        echo "Next: bash 05_wasp_setup.sh"
        ;;
esac
