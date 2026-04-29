#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# run_pipeline.sh — Orchestrator for the F121-9 ATAC pipeline
# =============================================================================
#
# WHAT THIS DOES:
#   Submits the full ATAC pipeline (steps 00-09) as SLURM jobs with proper
#   dependency chains, so each stage waits for its prerequisites.
#
# MODES:
#   --dry-run (default)   Print what WOULD be submitted, then stop.
#   --submit              Actually submit jobs.
#
# DEPENDENCY GRAPH (linearized):
#   00 setup
#   ↓
#   01 download (array of 4)
#   ↓
#   02 trim (array of 4)
#   ↓
#   03 align (array of 4)
#   ├──┬──→ 04a peaks per-rep (array of 4) → 04b consensus (single)
#   │  │
#   │  └──→ 05 wasp setup (single, parallel to 04)
#   │           ↓
#   │       06 wasp filter (array of 4)
#   │           ↓
#   │       07 allele split (array of 4)
#   │           ↓
#   └──────→ 08 bigwigs (array of 4, depends on both 03 and 07)
#                ↓
#            09 merge + qc (single)
#
# Steps 04 (peaks) and 05+06+07+08 (allele track) run in parallel after 03.
# Step 09 waits for everything.
#
# USAGE:
#   bash run_pipeline.sh                 # dry-run (default)
#   bash run_pipeline.sh --submit        # actually submit
#   bash run_pipeline.sh --submit --start-from 03   # skip 00-02 (already done)
#
# Per-step resource specs are at the top of each section below; tune them if
# your data or cluster differ.
# =============================================================================

source "config.sh"

# ----- CLI -----
SUBMIT=false
START_FROM=00
while [[ $# -gt 0 ]]; do
    case "$1" in
        --submit)        SUBMIT=true; shift ;;
        --dry-run)       SUBMIT=false; shift ;;
        --start-from)    START_FROM="$2"; shift 2 ;;
        -h|--help)
            grep -E '^#( |$)' "$0" | sed -E 's/^# ?//' | head -50
            exit 0 ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            exit 1 ;;
    esac
done

NSAMPLES=${#SAMPLE_ORDER[@]}
ARRAY_SPEC="0-$((NSAMPLES - 1))"

# Compute partition / log dir
mkdir -p "${LOG_DIR}"

# Pretty-print plan
echo ""
echo "============================================================"
echo "F121-9 ATAC Pipeline — Submission Plan"
echo "============================================================"
echo "  Mode:         $( ${SUBMIT} && echo SUBMIT || echo "DRY RUN" )"
echo "  Start from:   ${START_FROM}"
echo "  Samples:      ${NSAMPLES} (${SAMPLE_ORDER[*]})"
echo "  Array spec:   ${ARRAY_SPEC}"
echo "  Log dir:      ${LOG_DIR}"
echo "  Scratch:      ${SCRATCH_DIR}"
echo "============================================================"

# -----------------------------------------------------------------------------
# Submission helper
# -----------------------------------------------------------------------------
# Usage: submit_step STEP_NUM SBATCH_ARGS_STRING SCRIPT_AND_ARGS
#   Returns the JOB_ID via SUBMITTED_JOBID global, or NA in dry-run.
SUBMITTED_JOBID="NA"
submit_step() {
    local step_label="$1"
    local sbatch_args="$2"
    local script_args="$3"

    echo ""
    echo "--- ${step_label} ---"
    echo "  sbatch ${sbatch_args} \\"
    echo "         ${script_args}"

    if ${SUBMIT}; then
        # Inject standard log paths
        local log_args="-o ${LOG_DIR}/${step_label}_%A_%a.out -e ${LOG_DIR}/${step_label}_%A_%a.err"
        # If not array, use single-job log naming
        if [[ "${sbatch_args}" != *"--array"* ]]; then
            log_args="-o ${LOG_DIR}/${step_label}_%j.out -e ${LOG_DIR}/${step_label}_%j.err"
        fi

        local out
        out=$(sbatch ${sbatch_args} ${log_args} --parsable ${script_args})
        SUBMITTED_JOBID="${out%%;*}"   # in case --parsable returns "JOBID;cluster"
        echo "  -> Submitted as job ${SUBMITTED_JOBID}"
    else
        SUBMITTED_JOBID="DRY_$(date +%N)_${step_label}"
        echo "  (dry run; would set ${step_label}=${SUBMITTED_JOBID})"
    fi
}

# Boolean: should we run a step at >= START_FROM?
should_run() {
    local step="$1"
    [[ "${step}" > "${START_FROM}" || "${step}" == "${START_FROM}" ]]
}

# -----------------------------------------------------------------------------
# 00 — Reference setup (single)
# -----------------------------------------------------------------------------
JOB_00=""
if should_run 00; then
    submit_step "00_setup" \
        "--partition=short --time=2:00:00 --mem=16G --cpus-per-task=8 --job-name=atac_00" \
        "${SCRIPT_DIR}/00_setup_references.sh"
    JOB_00="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 01 — Download (array of NSAMPLES)
# -----------------------------------------------------------------------------
JOB_01=""
if should_run 01; then
    DEPS=""
    [[ -n "${JOB_00}" ]] && DEPS="--dependency=afterok:${JOB_00}"
    submit_step "01_download" \
        "--array=${ARRAY_SPEC} --partition=short --time=6:00:00 --mem=16G --cpus-per-task=8 --job-name=atac_01 ${DEPS}" \
        "${SCRIPT_DIR}/01_download_fastq.sh"
    JOB_01="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 02 — Trimming (array of NSAMPLES)
# -----------------------------------------------------------------------------
JOB_02=""
if should_run 02; then
    DEPS=""
    [[ -n "${JOB_01}" ]] && DEPS="--dependency=afterok:${JOB_01}"
    submit_step "02_trim" \
        "--array=${ARRAY_SPEC} --partition=short --time=4:00:00 --mem=8G --cpus-per-task=8 --job-name=atac_02 ${DEPS}" \
        "${SCRIPT_DIR}/02_trim_qc.sh"
    JOB_02="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 03 — Align + filter + dedup (array of NSAMPLES)
# -----------------------------------------------------------------------------
JOB_03=""
if should_run 03; then
    DEPS=""
    [[ -n "${JOB_02}" ]] && DEPS="--dependency=afterok:${JOB_02}"
    submit_step "03_align" \
        "--array=${ARRAY_SPEC} --partition=short --time=12:00:00 --mem=32G --cpus-per-task=8 --job-name=atac_03 ${DEPS}" \
        "${SCRIPT_DIR}/03_align_filter.sh"
    JOB_03="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 04a — Peaks per replicate (array of NSAMPLES)
# 04b — Consensus peaks (single, depends on 04a)
# These run in PARALLEL with 05+06+07 after 03.
# -----------------------------------------------------------------------------
JOB_04A=""
JOB_04B=""
if should_run 04; then
    DEPS=""
    [[ -n "${JOB_03}" ]] && DEPS="--dependency=afterok:${JOB_03}"
    submit_step "04a_peaks_perrep" \
        "--array=${ARRAY_SPEC} --partition=short --time=2:00:00 --mem=8G --cpus-per-task=4 --job-name=atac_04a ${DEPS}" \
        "${SCRIPT_DIR}/04_call_peaks.sh per-rep"
    JOB_04A="${SUBMITTED_JOBID}"

    DEPS_B=""
    [[ -n "${JOB_04A}" ]] && DEPS_B="--dependency=afterok:${JOB_04A}"
    submit_step "04b_peaks_consensus" \
        "--partition=short --time=00:30:00 --mem=8G --cpus-per-task=2 --job-name=atac_04b ${DEPS_B}" \
        "${SCRIPT_DIR}/04_call_peaks.sh consensus"
    JOB_04B="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 05 — WASP setup (single, parallel to 04)
# -----------------------------------------------------------------------------
JOB_05=""
if should_run 05; then
    DEPS=""
    [[ -n "${JOB_03}" ]] && DEPS="--dependency=afterok:${JOB_03}"
    submit_step "05_wasp_setup" \
        "--partition=short --time=1:00:00 --mem=8G --cpus-per-task=4 --job-name=atac_05 ${DEPS}" \
        "${SCRIPT_DIR}/05_wasp_setup.sh"
    JOB_05="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 06 — WASP filter (array of NSAMPLES, depends on 05)
# -----------------------------------------------------------------------------
JOB_06=""
if should_run 06; then
    DEPS=""
    [[ -n "${JOB_05}" ]] && DEPS="--dependency=afterok:${JOB_05}"
    submit_step "06_wasp_filter" \
        "--array=${ARRAY_SPEC} --partition=short --time=8:00:00 --mem=32G --cpus-per-task=8 --job-name=atac_06 ${DEPS}" \
        "${SCRIPT_DIR}/06_wasp_filter.sh"
    JOB_06="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 07 — Allele-specific splitting (array of NSAMPLES, depends on 06)
# -----------------------------------------------------------------------------
JOB_07=""
if should_run 07; then
    DEPS=""
    [[ -n "${JOB_06}" ]] && DEPS="--dependency=afterok:${JOB_06}"
    submit_step "07_allele_split" \
        "--array=${ARRAY_SPEC} --partition=short --time=4:00:00 --mem=16G --cpus-per-task=8 --job-name=atac_07 ${DEPS}" \
        "${SCRIPT_DIR}/07_allele_specific.sh"
    JOB_07="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 08 — Per-sample bigWigs (array; depends on 07 for allele tracks AND 03 for
#                          per-rep tracks. We use 07 since 07 depends on 03.)
# -----------------------------------------------------------------------------
JOB_08=""
if should_run 08; then
    DEPS=""
    [[ -n "${JOB_07}" ]] && DEPS="--dependency=afterok:${JOB_07}"
    submit_step "08_bigwigs" \
        "--array=${ARRAY_SPEC} --partition=short --time=2:00:00 --mem=16G --cpus-per-task=4 --job-name=atac_08 ${DEPS}" \
        "${SCRIPT_DIR}/08_make_bigwigs.sh"
    JOB_08="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 09 — Merge + final QC (single, depends on 08 AND 04b)
# -----------------------------------------------------------------------------
JOB_09=""
if should_run 09; then
    DEPS_LIST=()
    [[ -n "${JOB_08}" ]] && DEPS_LIST+=("${JOB_08}")
    [[ -n "${JOB_04B}" ]] && DEPS_LIST+=("${JOB_04B}")
    DEPS=""
    if (( ${#DEPS_LIST[@]} > 0 )); then
        DEPS="--dependency=afterok:$(IFS=:; echo "${DEPS_LIST[*]}")"
    fi
    submit_step "09_merge_qc" \
        "--partition=short --time=2:00:00 --mem=16G --cpus-per-task=4 --job-name=atac_09 ${DEPS}" \
        "${SCRIPT_DIR}/09_merge_qc.sh"
    JOB_09="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "Submission summary"
echo "============================================================"
printf "  %-25s %s\n" "00 setup:"            "${JOB_00:-(skipped)}"
printf "  %-25s %s\n" "01 download:"         "${JOB_01:-(skipped)}"
printf "  %-25s %s\n" "02 trim:"             "${JOB_02:-(skipped)}"
printf "  %-25s %s\n" "03 align:"            "${JOB_03:-(skipped)}"
printf "  %-25s %s\n" "04a peaks per-rep:"   "${JOB_04A:-(skipped)}"
printf "  %-25s %s\n" "04b consensus:"       "${JOB_04B:-(skipped)}"
printf "  %-25s %s\n" "05 wasp setup:"       "${JOB_05:-(skipped)}"
printf "  %-25s %s\n" "06 wasp filter:"      "${JOB_06:-(skipped)}"
printf "  %-25s %s\n" "07 allele split:"     "${JOB_07:-(skipped)}"
printf "  %-25s %s\n" "08 bigwigs:"          "${JOB_08:-(skipped)}"
printf "  %-25s %s\n" "09 merge + qc:"       "${JOB_09:-(skipped)}"
echo ""

if ${SUBMIT}; then
    echo "Monitor: squeue -u \$(whoami)"
    echo "Logs:    ${LOG_DIR}/"
else
    echo "Dry run complete. Re-run with --submit to actually launch."
fi
