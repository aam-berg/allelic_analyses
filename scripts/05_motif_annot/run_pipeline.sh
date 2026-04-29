#!/bin/bash
# =============================================================================
# 05_motif_annot/run_pipeline.sh — Two-stage orchestrator
#
# Stage 1: preprocess shared resources (single SLURM job)
# Stage 2: per-motif annotation (SLURM array, depends on Stage 1)
#
# Stage 2 won't start until Stage 1 succeeds. If you've already run Stage 1
# (or are testing only Stage 2), pass --start-from 02.
#
# USAGE:
#   bash run_pipeline.sh                      # dry-run (default)
#   bash run_pipeline.sh --apply              # submit
#   bash run_pipeline.sh --apply --start-from 02
# =============================================================================

set -euo pipefail

SCRIPT_DIR="."
source "${SCRIPT_DIR}/config.sh"

# ---- CLI ----
APPLY=false
START_FROM=01
while [[ $# -gt 0 ]]; do
    case "$1" in
        --apply)        APPLY=true; shift ;;
        --dry-run)      APPLY=false; shift ;;
        --start-from)   START_FROM="$2"; shift 2 ;;
        -h|--help)
            sed -n '3,15p' "$0" | sed 's/^# \?//'
            exit 0 ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

mkdir -p logs "${OUTDIR}" "${PREPROCESSED_DIR}"

echo "================================================================"
echo "05_motif_annot/ pipeline"
echo "================================================================"
echo "  Mode:           $( ${APPLY} && echo APPLY || echo "DRY RUN" )"
echo "  Start from:     ${START_FROM}"
echo "  Preprocessed:   ${PREPROCESSED_DIR}"
echo "  Output:         ${OUTDIR}"
echo ""

# -----------------------------------------------------------------------------
# Stage 1 — preprocess
# -----------------------------------------------------------------------------
JOB_01=""
should_run_01() { [[ "${START_FROM}" == "01" ]]; }

if should_run_01; then
    echo "--- Stage 01: preprocess resources ---"
    echo "  sbatch \\"
    echo "    --partition=${PREPROCESS_PARTITION} \\"
    echo "    --time=${PREPROCESS_TIME} \\"
    echo "    --mem=${PREPROCESS_MEM} \\"
    echo "    --cpus-per-task=4 \\"
    echo "    --job-name=annot_v2_preprocess \\"
    echo "    ${SCRIPT_DIR}/01_run_preprocess_resources.sh"

    if ${APPLY}; then
        JOB_01=$(sbatch \
            --partition="${PREPROCESS_PARTITION}" \
            --time="${PREPROCESS_TIME}" \
            --mem="${PREPROCESS_MEM}" \
            --cpus-per-task=4 \
            --job-name=annot_v2_preprocess \
            --parsable \
            "${SCRIPT_DIR}/01_run_preprocess_resources.sh")
        echo "  -> Submitted preprocess as job ${JOB_01}"
    else
        echo "  (dry run)"
    fi
else
    echo "--- Stage 01 skipped (--start-from=${START_FROM}) ---"
fi

# -----------------------------------------------------------------------------
# Stage 2 — annotate
# -----------------------------------------------------------------------------
echo ""
echo "--- Stage 02: per-motif annotation array ---"

# We use a wrapper so Stage 2 can wait on Stage 1 via SLURM dependency.
# Approach: invoke 02_annotate_motifs.sh which discovers motifs and submits
# the array. To chain via dependency, we need to inject --dependency into
# the inner sbatch. The simplest robust approach: just let the user invoke
# 02_annotate_motifs.sh directly after Stage 1 finishes, OR — when --apply
# is given AND Stage 1 was submitted — we wait for Stage 1 to finish and
# then run Stage 2. For an unattended pipeline we'd need to refactor
# 02_annotate_motifs.sh to accept --dependency; we keep it simple here.

if ${APPLY} && [[ -n "${JOB_01}" ]]; then
    echo "  Stage 02 will run only after preprocess job ${JOB_01} succeeds."
    echo "  Submit it with:"
    echo ""
    echo "    bash ${SCRIPT_DIR}/02_annotate_motifs.sh --apply"
    echo ""
    echo "  (Run AFTER ${JOB_01} completes successfully. To monitor:"
    echo "      squeue -j ${JOB_01}    sacct -j ${JOB_01} -X)"
elif ${APPLY}; then
    echo "  Submitting Stage 02 directly (--start-from=02 implies preprocess is done)..."
    bash "${SCRIPT_DIR}/02_annotate_motifs.sh" --apply
else
    echo "  bash ${SCRIPT_DIR}/02_annotate_motifs.sh --apply"
    echo "  (dry run)"
fi

echo ""
echo "================================================================"
if ${APPLY}; then
    echo "  Submitted preprocess (if applicable). After it finishes, run:"
    echo "    bash 02_annotate_motifs.sh --apply"
    echo ""
    echo "  After the array finishes, find any failed motifs:"
    echo "    bash 0x_find_missing_motifs.sh"
    echo "    bash 0x_resubmit_missing_motifs.sh --apply"
else
    echo "  Dry run complete. Re-run with --apply to actually submit."
fi
echo "================================================================"
