#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04_proseq/run_pipeline.sh — Orchestrator for the F121-9 PRO-seq pipeline
# =============================================================================
#
# WHAT THIS DOES:
#   Submits the full PRO-seq pipeline (steps 00-09) as SLURM jobs with
#   proper dependency chains. Each stage waits for its prerequisites.
#
# MODES:
#   --dry-run (default)   Print what WOULD be submitted, then stop.
#   --submit              Actually submit the jobs.
#
# DEPENDENCY GRAPH:
#   00 setup
#   ↓
#   01 download (array of 4)
#   ↓
#   02 trim (array of 4)
#   ↓
#   03 align (array of 4)
#   ├─→ 04 per-rep bigwigs (array of 4) ─────────┐
#   │                                            │
#   └─→ 05 wasp setup (single)                   │
#         ↓                                      │
#       06 wasp filter (array of 4)              │
#         ↓                                      │
#       07 allele split (array of 4)             │
#         ↓                                      │
#       08 allele bigwigs (array of 4) ──────────┤
#                                                ↓
#                                           09 merge + qc (single)
#
# Steps 04 (per-rep bigwigs) and 05→06→07→08 (allele track) run in parallel
# after step 03 completes. Step 09 waits for both branches.
#
# USAGE:
#   bash run_pipeline.sh                          # dry-run (default)
#   bash run_pipeline.sh --submit                 # actually submit
#   bash run_pipeline.sh --submit --start-from 05 # skip 00-04 (already done)
#
# Per-step resource specs are below; tune them if your data or cluster differ.
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

mkdir -p "${LOG_DIR}"

# Pretty-print plan
echo ""
echo "============================================================"
echo "F121-9 PRO-seq Pipeline — Submission Plan"
echo "============================================================"
echo "  Mode:         $( ${SUBMIT} && echo SUBMIT || echo "DRY RUN" )"
echo "  Start from:   ${START_FROM}"
echo "  Samples:      ${NSAMPLES} (${SAMPLE_ORDER[*]})"
echo "  Array spec:   ${ARRAY_SPEC}"
echo "  Log dir:      ${LOG_DIR}"
echo "  Scratch:      ${SCRATCH_DIR}"
echo "  dm6 spike-in filter: ${DM6_SPIKEIN_FILTER}"
echo "============================================================"

# -----------------------------------------------------------------------------
# Submission helper (same pattern as 03_atac/run_pipeline.sh)
# -----------------------------------------------------------------------------
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
        local log_args="-o ${LOG_DIR}/${step_label}_%A_%a.out -e ${LOG_DIR}/${step_label}_%A_%a.err"
        if [[ "${sbatch_args}" != *"--array"* ]]; then
            log_args="-o ${LOG_DIR}/${step_label}_%j.out -e ${LOG_DIR}/${step_label}_%j.err"
        fi

        local out
        out=$(sbatch ${sbatch_args} ${log_args} --parsable ${script_args})
        SUBMITTED_JOBID="${out%%;*}"
        echo "  -> Submitted as job ${SUBMITTED_JOBID}"
    else
        SUBMITTED_JOBID="DRY_$(date +%N)_${step_label}"
        echo "  (dry run; would set ${step_label}=${SUBMITTED_JOBID})"
    fi
}

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
        "--partition=short --time=2:00:00 --mem=16G --cpus-per-task=8 --job-name=proseq_00" \
        "${SCRIPT_DIR}/00_setup_references.sh"
    JOB_00="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 01 — Download (array)
# -----------------------------------------------------------------------------
JOB_01=""
if should_run 01; then
    DEPS=""
    [[ -n "${JOB_00}" ]] && DEPS="--dependency=afterok:${JOB_00}"
    submit_step "01_download" \
        "--array=${ARRAY_SPEC} --partition=short --time=6:00:00 --mem=16G --cpus-per-task=8 --job-name=proseq_01 ${DEPS}" \
        "${SCRIPT_DIR}/01_download_fastq.sh"
    JOB_01="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 02 — Trimming (array)
# -----------------------------------------------------------------------------
JOB_02=""
if should_run 02; then
    DEPS=""
    [[ -n "${JOB_01}" ]] && DEPS="--dependency=afterok:${JOB_01}"
    submit_step "02_trim" \
        "--array=${ARRAY_SPEC} --partition=short --time=4:00:00 --mem=8G --cpus-per-task=8 --job-name=proseq_02 ${DEPS}" \
        "${SCRIPT_DIR}/02_trim_qc.sh"
    JOB_02="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 03 — Align + filter (array). Three sequential bowtie2 runs (rRNA, dm6?, mm39)
#       so this gets the longest single-task wall time.
# -----------------------------------------------------------------------------
JOB_03=""
if should_run 03; then
    DEPS=""
    [[ -n "${JOB_02}" ]] && DEPS="--dependency=afterok:${JOB_02}"
    submit_step "03_align" \
        "--array=${ARRAY_SPEC} --partition=medium --time=12:00:00 --mem=32G --cpus-per-task=8 --job-name=proseq_03 ${DEPS}" \
        "${SCRIPT_DIR}/03_align_filter.sh"
    JOB_03="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 04 — Per-rep bigWigs (array, depends on 03; PARALLEL with 05->06->07->08)
# -----------------------------------------------------------------------------
JOB_04=""
if should_run 04; then
    DEPS=""
    [[ -n "${JOB_03}" ]] && DEPS="--dependency=afterok:${JOB_03}"
    submit_step "04_perrep_bw" \
        "--array=${ARRAY_SPEC} --partition=short --time=2:00:00 --mem=8G --cpus-per-task=4 --job-name=proseq_04 ${DEPS}" \
        "${SCRIPT_DIR}/04_make_bigwigs.sh"
    JOB_04="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 05 — WASP setup (single, depends on 03; parallel to 04)
# -----------------------------------------------------------------------------
JOB_05=""
if should_run 05; then
    DEPS=""
    [[ -n "${JOB_03}" ]] && DEPS="--dependency=afterok:${JOB_03}"
    submit_step "05_wasp_setup" \
        "--partition=short --time=1:00:00 --mem=8G --cpus-per-task=4 --job-name=proseq_05 ${DEPS}" \
        "${SCRIPT_DIR}/05_wasp_setup.sh"
    JOB_05="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 06 — WASP filter (array, depends on 05)
# -----------------------------------------------------------------------------
JOB_06=""
if should_run 06; then
    DEPS=""
    [[ -n "${JOB_05}" ]] && DEPS="--dependency=afterok:${JOB_05}"
    submit_step "06_wasp_filter" \
        "--array=${ARRAY_SPEC} --partition=short --time=8:00:00 --mem=32G --cpus-per-task=8 --job-name=proseq_06 ${DEPS}" \
        "${SCRIPT_DIR}/06_wasp_filter.sh"
    JOB_06="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 07 — Allele split (array, depends on 06)
# -----------------------------------------------------------------------------
JOB_07=""
if should_run 07; then
    DEPS=""
    [[ -n "${JOB_06}" ]] && DEPS="--dependency=afterok:${JOB_06}"
    submit_step "07_allele_split" \
        "--array=${ARRAY_SPEC} --partition=short --time=4:00:00 --mem=16G --cpus-per-task=8 --job-name=proseq_07 ${DEPS}" \
        "${SCRIPT_DIR}/07_allele_specific.sh"
    JOB_07="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 08 — Allele bigWigs (array, depends on 07)
# -----------------------------------------------------------------------------
JOB_08=""
if should_run 08; then
    DEPS=""
    [[ -n "${JOB_07}" ]] && DEPS="--dependency=afterok:${JOB_07}"
    submit_step "08_allele_bw" \
        "--array=${ARRAY_SPEC} --partition=short --time=2:00:00 --mem=8G --cpus-per-task=4 --job-name=proseq_08 ${DEPS}" \
        "${SCRIPT_DIR}/08_make_bigwigs_allele.sh"
    JOB_08="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 09 — Merge + final QC (single, depends on BOTH 04 AND 08)
# -----------------------------------------------------------------------------
JOB_09=""
if should_run 09; then
    DEPS_LIST=()
    [[ -n "${JOB_04}" ]] && DEPS_LIST+=("${JOB_04}")
    [[ -n "${JOB_08}" ]] && DEPS_LIST+=("${JOB_08}")
    DEPS=""
    if (( ${#DEPS_LIST[@]} > 0 )); then
        DEPS="--dependency=afterok:$(IFS=:; echo "${DEPS_LIST[*]}")"
    fi
    submit_step "09_merge_qc" \
        "--partition=short --time=2:00:00 --mem=16G --cpus-per-task=4 --job-name=proseq_09 ${DEPS}" \
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
printf "  %-25s %s\n" "04 per-rep bigwigs:"  "${JOB_04:-(skipped)}"
printf "  %-25s %s\n" "05 wasp setup:"       "${JOB_05:-(skipped)}"
printf "  %-25s %s\n" "06 wasp filter:"      "${JOB_06:-(skipped)}"
printf "  %-25s %s\n" "07 allele split:"     "${JOB_07:-(skipped)}"
printf "  %-25s %s\n" "08 allele bigwigs:"   "${JOB_08:-(skipped)}"
printf "  %-25s %s\n" "09 merge + qc:"       "${JOB_09:-(skipped)}"
echo ""

if ${SUBMIT}; then
    echo "Monitor: squeue -u \$(whoami)"
    echo "Logs:    ${LOG_DIR}/"
else
    echo "Dry run complete. Re-run with --submit to actually launch."
fi
