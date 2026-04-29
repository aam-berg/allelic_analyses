#!/bin/bash
# =============================================================================
# run_pipeline.sh — Master orchestrator for MOODS genome-wide motif scanning
#
# Pipeline steps:
#   Step 0: Compute genome background frequencies, list chromosomes
#   Step 1: Scan each chromosome against all motifs (SLURM array job)
#   Step 2: Merge results, palindrome dedup (optional), QC report
#   Step 3: Comprehensive QC report (genomic context, TSS dist, plots)
#
# All paths come from config.sh — no hardcoded paths in this script.
#
# Usage:
#   bash run_pipeline.sh                         # dry run
#   bash run_pipeline.sh --submit                # submit
#   bash run_pipeline.sh --submit --no-palindrome-dedup
# =============================================================================

set -euo pipefail

SCRIPT_DIR="."
source "${SCRIPT_DIR}/config.sh"

# ---- Args ----
DRY_RUN=true
USE_PALINDROME_DEDUP=true
for arg in "$@"; do
    case "$arg" in
        --submit) DRY_RUN=false ;;
        --no-palindrome-dedup) USE_PALINDROME_DEDUP=false ;;
    esac
done

mkdir -p "${PREPDIR}" "${CHROMDIR}" "${MERGEDIR}" "${QCDIR}" "${LOGDIR}"

echo "============================================"
echo "MOODS Genome-Wide Motif Scanning Pipeline"
echo "============================================"
echo "  Genome:        ${GENOME}"
echo "  Chrom sizes:   ${CHROM_SIZES}"
echo "  MEME:          ${MEME_FILE}"
echo "  GTF:           ${GTF_FILE}"
echo "  P-value:       ${PVALUE}"
echo "  Output root:   ${BASE_OUTDIR}"
echo "  Promoter win:  [-${PROMOTER_UPSTREAM}, +${PROMOTER_DOWNSTREAM}]"
echo "  Dry run:       ${DRY_RUN}"
echo "  Dedup palins:  ${USE_PALINDROME_DEDUP}"
echo ""

# ============================================================================
# Step 0
# ============================================================================
echo "Step 0: Background frequencies & chromosome manifest"

if [[ -f "${PREPDIR}/background.txt" && -f "${PREPDIR}/chromosomes.txt" ]]; then
    echo "  Already computed in ${PREPDIR}. (Delete to force redo.)"
else
    if [[ "${DRY_RUN}" == true ]]; then
        echo "  [DRY] python ${SCRIPT_DIR}/00_prepare.py ${GENOME} ${PREPDIR}"
    else
        python "${SCRIPT_DIR}/00_prepare.py" "${GENOME}" "${PREPDIR}"
    fi
fi

if [[ -f "${PREPDIR}/chromosomes.txt" ]]; then
    N_CHROMS=$(wc -l < "${PREPDIR}/chromosomes.txt")
else
    N_CHROMS=21  # placeholder for dry-run before step 0 has produced output
fi
echo "  Chromosomes to scan: ${N_CHROMS}"
echo ""

# ============================================================================
# Pre-flight: motif sanity check (produces palindrome table)
# ============================================================================
echo "Pre-flight: motif sanity check"

# Track whether dedup will end up being applied so the dry-run accurately
# predicts the --submit behavior (the sanity check produces the palindrome
# table, but only actually runs in --submit mode).
DEDUP_WILL_APPLY=false

if [[ -f "${PALINDROME_TABLE}" ]]; then
    echo "  Palindrome table exists: ${PALINDROME_TABLE}"
    [[ "${USE_PALINDROME_DEDUP}" == true ]] && DEDUP_WILL_APPLY=true
elif [[ "${USE_PALINDROME_DEDUP}" == true ]]; then
    if [[ "${DRY_RUN}" == true ]]; then
        echo "  [DRY] bash ${SCRIPT_DIR}/0x_run_sanity_check.sh"
        echo "        (produces ${PALINDROME_TABLE} synchronously in --submit mode)"
        DEDUP_WILL_APPLY=true
    else
        bash "${SCRIPT_DIR}/0x_run_sanity_check.sh"
        [[ -f "${PALINDROME_TABLE}" ]] && DEDUP_WILL_APPLY=true
    fi
fi

if [[ "${DEDUP_WILL_APPLY}" == true ]]; then
    PALINDROME_ARG="--palindrome-table ${PALINDROME_TABLE}"
    echo "  Palindrome dedup will be applied at step 2."
else
    PALINDROME_ARG=""
    echo "  Palindrome dedup will NOT be applied."
fi
echo ""

# ============================================================================
# Step 1: MOODS scan array job
# ============================================================================
echo "Step 1: MOODS scan (array, one task per chromosome)"

if [[ "${DRY_RUN}" == true ]]; then
    echo "  [DRY] sbatch --array=1-${N_CHROMS} ${SCRIPT_DIR}/01_scan_moods.sbatch"
    SCAN_JOBID="DRYRUN"
else
    SCAN_JOBID=$(sbatch \
        --array=1-${N_CHROMS} \
        --parsable \
        -o "${LOGDIR}/moods_scan_%A_%a.out" \
        -e "${LOGDIR}/moods_scan_%A_%a.err" \
        "${SCRIPT_DIR}/01_scan_moods.sbatch")
    echo "  Submitted: ${SCAN_JOBID}"
fi
echo ""

# ============================================================================
# Step 2: Post-processing (depends on scan)
# ============================================================================
echo "Step 2: Post-processing (depends on step 1)"

POSTPROCESS_CMD="python ${SCRIPT_DIR}/02_postprocess.py \
    --input-dir ${CHROMDIR} \
    --output-dir ${MERGEDIR} \
    --background ${PREPDIR}/background.txt \
    --meme ${MEME_FILE} \
    --pvalue ${PVALUE} \
    --chrom-sizes ${CHROM_SIZES} \
    ${PALINDROME_ARG}"

if [[ "${DRY_RUN}" == true ]]; then
    echo "  [DRY] ${POSTPROCESS_CMD}"
    POST_JOBID="DRYRUN"
else
    POST_JOBID=$(sbatch \
        --parsable \
        --dependency=afterok:${SCAN_JOBID} \
        --job-name=moods_postprocess \
        --partition="${POST_PARTITION}" \
        --time="${POST_TIME}" \
        --mem="${POST_MEM}" \
        -o "${LOGDIR}/moods_postprocess_%j.out" \
        -e "${LOGDIR}/moods_postprocess_%j.err" \
        --wrap="set -euo pipefail; source activate ${CONDA_ENV}; ${POSTPROCESS_CMD}")
    echo "  Submitted: ${POST_JOBID} (depends on ${SCAN_JOBID})"
fi
echo ""

# ============================================================================
# Step 3: Comprehensive QC (depends on postprocess)
# ============================================================================
echo "Step 3: Comprehensive QC report (depends on step 2)"

QC_CMD="python ${SCRIPT_DIR}/03_qc_report.py \
    --motif-dir ${MERGEDIR}/per_motif \
    --merged-bed ${MERGEDIR}/all_motifs_merged.bed.gz \
    --summary-tsv ${MERGEDIR}/motif_hit_summary.tsv \
    --gtf ${GTF_FILE} \
    --chrom-sizes ${CHROM_SIZES} \
    --meme ${MEME_FILE} \
    --background ${PREPDIR}/background.txt \
    --pvalue ${PVALUE} \
    --promoter-upstream ${PROMOTER_UPSTREAM} \
    --promoter-downstream ${PROMOTER_DOWNSTREAM} \
    --output-dir ${QCDIR}"

if [[ "${DRY_RUN}" == true ]]; then
    echo "  [DRY] ${QC_CMD}"
else
    QC_JOBID=$(sbatch \
        --parsable \
        --dependency=afterok:${POST_JOBID} \
        --job-name=moods_qc \
        --partition="${QC_PARTITION}" \
        --time="${QC_TIME}" \
        --mem="${QC_MEM}" \
        -o "${LOGDIR}/moods_qc_%j.out" \
        -e "${LOGDIR}/moods_qc_%j.err" \
        --wrap="set -euo pipefail; source activate ${CONDA_ENV}; ${QC_CMD}")
    echo "  Submitted: ${QC_JOBID} (depends on ${POST_JOBID})"
fi
echo ""

echo "============================================"
if [[ "${DRY_RUN}" == true ]]; then
    echo "Dry run complete. Re-run with --submit to launch jobs."
else
    echo "All jobs submitted."
fi
echo "============================================"
