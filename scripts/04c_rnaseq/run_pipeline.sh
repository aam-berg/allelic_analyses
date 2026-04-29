#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04c_rnaseq/run_pipeline.sh — Orchestrator for the F121-9 RNA-seq pipeline
# =============================================================================
#
# WHAT THIS DOES:
#   Submits the full RNA-seq pipeline (steps 00-08) as SLURM jobs with
#   proper dependency chains.
#
# MODES:
#   --dry-run (default)   Print what WOULD be submitted, then stop.
#   --submit              Actually submit the jobs.
#
# DEPENDENCY GRAPH:
#   00 setup
#   ↓
#   01 download (array of 2)
#   ↓
#   02 trim (array of 2)
#   ↓
#   03 STAR align + WASP (array of 2; longest single step)
#   ├─→ 04 quantify (array of 2) ─────────────────────┐
#   ├─→ 05 stranded bigwigs (array of 2) ─────────────┤
#   └─→ 06 WASP filter + allele split (array of 2) ───┤
#         ↓                                            │
#       07 allele bigwigs (array of 2) ────────────────┤
#                                                      ↓
#                                                 08 merge + qc (single)
#
# Steps 04, 05, 06 all run in parallel after step 03 completes.
# Step 07 depends on step 06 (needs allele BAMs).
# Step 08 waits for 04, 05, 07 (all the per-sample work).
#
# USAGE:
#   bash run_pipeline.sh                          # dry-run (default)
#   bash run_pipeline.sh --submit                 # actually submit
#   bash run_pipeline.sh --submit --start-from 03 # skip 00-02 (already done)
#
# Per-step resource specs are below; tune them if your cluster differs.
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

echo ""
echo "============================================================"
echo "F121-9 RNA-seq Pipeline — Submission Plan"
echo "============================================================"
echo "  Mode:         $( ${SUBMIT} && echo SUBMIT || echo "DRY RUN" )"
echo "  Start from:   ${START_FROM}"
echo "  Samples:      ${NSAMPLES} (${SAMPLE_ORDER[*]})"
echo "  Array spec:   ${ARRAY_SPEC}"
echo "  Log dir:      ${LOG_DIR}"
echo "  Scratch:      ${SCRATCH_DIR}"
echo "  Library:      TruSeq Stranded Total RNA (reverse-stranded), PE"
echo "============================================================"

# -----------------------------------------------------------------------------
# Submission helper
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
        "--partition=short --time=4:00:00 --mem=64G --cpus-per-task=8 --job-name=rnaseq_00" \
        "${SCRIPT_DIR}/00_setup_references.sh"
    JOB_00="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 01 — Multi-SRR download (array)
# -----------------------------------------------------------------------------
JOB_01=""
if should_run 01; then
    DEPS=""
    [[ -n "${JOB_00}" ]] && DEPS="--dependency=afterok:${JOB_00}"
    submit_step "01_download" \
        "--array=${ARRAY_SPEC} --partition=short --time=4:00:00 --mem=16G --cpus-per-task=8 --job-name=rnaseq_01 ${DEPS}" \
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
        "--array=${ARRAY_SPEC} --partition=short --time=2:00:00 --mem=8G --cpus-per-task=8 --job-name=rnaseq_02 ${DEPS}" \
        "${SCRIPT_DIR}/02_trim_qc.sh"
    JOB_02="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 03 — STAR align + WASP (array). Longest single step; needs ample memory.
# -----------------------------------------------------------------------------
JOB_03=""
if should_run 03; then
    DEPS=""
    [[ -n "${JOB_02}" ]] && DEPS="--dependency=afterok:${JOB_02}"
    submit_step "03_star_align" \
        "--array=${ARRAY_SPEC} --partition=medium --time=24:00:00 --mem=64G --cpus-per-task=8 --job-name=rnaseq_03 ${DEPS}" \
        "${SCRIPT_DIR}/03_align_star.sh"
    JOB_03="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 04 — Quantification (array, depends on 03; PARALLEL with 05 and 06)
# -----------------------------------------------------------------------------
JOB_04=""
if should_run 04; then
    DEPS=""
    [[ -n "${JOB_03}" ]] && DEPS="--dependency=afterok:${JOB_03}"
    submit_step "04_quantify" \
        "--array=${ARRAY_SPEC} --partition=short --time=1:00:00 --mem=8G --cpus-per-task=8 --job-name=rnaseq_04 ${DEPS}" \
        "${SCRIPT_DIR}/04_quantify.sh"
    JOB_04="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 05 — Stranded bigwigs (array, depends on 03; PARALLEL with 04 and 06)
# -----------------------------------------------------------------------------
JOB_05=""
if should_run 05; then
    DEPS=""
    [[ -n "${JOB_03}" ]] && DEPS="--dependency=afterok:${JOB_03}"
    submit_step "05_bigwigs" \
        "--array=${ARRAY_SPEC} --partition=short --time=2:00:00 --mem=16G --cpus-per-task=4 --job-name=rnaseq_05 ${DEPS}" \
        "${SCRIPT_DIR}/05_make_bigwigs.sh"
    JOB_05="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 06 — WASP tag filter + allele split (array, depends on 03; PARALLEL with 04, 05)
# -----------------------------------------------------------------------------
JOB_06=""
if should_run 06; then
    DEPS=""
    [[ -n "${JOB_03}" ]] && DEPS="--dependency=afterok:${JOB_03}"
    submit_step "06_allele_split" \
        "--array=${ARRAY_SPEC} --partition=short --time=2:00:00 --mem=16G --cpus-per-task=8 --job-name=rnaseq_06 ${DEPS}" \
        "${SCRIPT_DIR}/06_allele_specific.sh"
    JOB_06="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 07 — Per-allele bigwigs (array, depends on 06)
# -----------------------------------------------------------------------------
JOB_07=""
if should_run 07; then
    DEPS=""
    [[ -n "${JOB_06}" ]] && DEPS="--dependency=afterok:${JOB_06}"
    submit_step "07_allele_bigwigs" \
        "--array=${ARRAY_SPEC} --partition=short --time=2:00:00 --mem=16G --cpus-per-task=4 --job-name=rnaseq_07 ${DEPS}" \
        "${SCRIPT_DIR}/07_make_bigwigs_allele.sh"
    JOB_07="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# 08 — Final merge + QC + expression summary (single, depends on 04, 05, 07)
# -----------------------------------------------------------------------------
JOB_08=""
if should_run 08; then
    DEPS_LIST=()
    [[ -n "${JOB_04}" ]] && DEPS_LIST+=("${JOB_04}")
    [[ -n "${JOB_05}" ]] && DEPS_LIST+=("${JOB_05}")
    [[ -n "${JOB_07}" ]] && DEPS_LIST+=("${JOB_07}")
    DEPS=""
    if (( ${#DEPS_LIST[@]} > 0 )); then
        DEPS="--dependency=afterok:$(IFS=:; echo "${DEPS_LIST[*]}")"
    fi
    submit_step "08_merge_qc" \
        "--partition=short --time=2:00:00 --mem=16G --cpus-per-task=4 --job-name=rnaseq_08 ${DEPS}" \
        "${SCRIPT_DIR}/08_merge_qc.sh"
    JOB_08="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "Submission summary"
echo "============================================================"
printf "  %-25s %s\n" "00 setup:"               "${JOB_00:-(skipped)}"
printf "  %-25s %s\n" "01 download:"            "${JOB_01:-(skipped)}"
printf "  %-25s %s\n" "02 trim:"                "${JOB_02:-(skipped)}"
printf "  %-25s %s\n" "03 STAR align + WASP:"   "${JOB_03:-(skipped)}"
printf "  %-25s %s\n" "04 quantify:"            "${JOB_04:-(skipped)}"
printf "  %-25s %s\n" "05 stranded bigwigs:"    "${JOB_05:-(skipped)}"
printf "  %-25s %s\n" "06 WASP filter + split:" "${JOB_06:-(skipped)}"
printf "  %-25s %s\n" "07 allele bigwigs:"      "${JOB_07:-(skipped)}"
printf "  %-25s %s\n" "08 merge + QC:"          "${JOB_08:-(skipped)}"
echo ""

if ${SUBMIT}; then
    echo "Monitor: squeue -u \$(whoami)"
    echo "Logs:    ${LOG_DIR}/"
else
    echo "Dry run complete. Re-run with --submit to actually launch."
fi
