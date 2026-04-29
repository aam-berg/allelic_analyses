#!/bin/bash
# =============================================================================
# 0x_resubmit_missing_motifs.sh — Re-submit SLURM array for missing motifs only
#
# Reads missing list produced by 0x_find_missing_motifs.sh and submits a new
# array job with bumped time/memory limits (see RETRY_* in config.sh).
#
# USAGE:
#   bash 0x_resubmit_missing_motifs.sh                # submit
#   bash 0x_resubmit_missing_motifs.sh --dry-run      # preview
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

MISSING_LIST="${OUTDIR}/missing_motifs.txt"

# ---- Args ----
DRY_RUN=false
for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=true ;;
    esac
done

# ---- Validate ----
if [[ ! -f "${MISSING_LIST}" ]]; then
    echo "ERROR: ${MISSING_LIST} not found." >&2
    echo "       Run 0x_find_missing_motifs.sh first." >&2
    exit 1
fi

N_MISSING=$(wc -l < "${MISSING_LIST}")

echo "=============================================="
echo "Re-submit missing motifs"
echo "Date: $(date)"
echo "=============================================="
echo "  Missing list: ${MISSING_LIST}"
echo "  Count:        ${N_MISSING}"
echo "  Partition:    ${RETRY_PARTITION}"
echo "  Time limit:   ${RETRY_TIME}"
echo "  Memory limit: ${RETRY_MEM}"
echo ""

if [[ ${N_MISSING} -eq 0 ]]; then
    echo "Nothing to resubmit — all motifs completed."
    exit 0
fi

# ---- Generate the retry array job script ----
mkdir -p "${OUTDIR}/logs"
JOB_SCRIPT="${OUTDIR}/run_annotate_missing.sh"

cat > "${JOB_SCRIPT}" << SLURM_EOF
#!/bin/bash
#SBATCH --job-name=annot_retry
#SBATCH --partition=${RETRY_PARTITION}
#SBATCH --time=${RETRY_TIME}
#SBATCH --mem=${RETRY_MEM}
#SBATCH --cpus-per-task=1
#SBATCH --output=${OUTDIR}/logs/annot_retry_%A_%a.out
#SBATCH --error=${OUTDIR}/logs/annot_retry_%A_%a.err

set -euo pipefail

echo "================================================================"
echo "Motif Annotation RETRY — Array Task \${SLURM_ARRAY_TASK_ID}"
echo "Job ID:    \${SLURM_ARRAY_JOB_ID}_\${SLURM_ARRAY_TASK_ID}"
echo "Node:      \$(hostname)"
echo "Date:      \$(date)"
echo "================================================================"

MOTIF_BED=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${MISSING_LIST}")
if [[ -z "\${MOTIF_BED}" ]]; then
    echo "ERROR: No file for array task \${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

MOTIF_ID=\$(basename "\${MOTIF_BED}" .bed)
echo "Motif: \${MOTIF_ID}"

# Remove any partial output so the R script doesn't skip it
rm -f "${OUTDIR}/\${MOTIF_ID}_annotated.tsv.gz"

set +u
source activate ${CONDA_ENV}
set -u

Rscript ${SCRIPT_DIR}/02_annotate_motifs.R \\
    --motif_bed "\${MOTIF_BED}" \\
    --genome_fa ${GENOME_FA} \\
    --preprocessed_dir ${PREPROCESSED_DIR} \\
    --vcf_dir ${VCF_DIR} \\
    --hybrids "${HYBRIDS}" \\
    --outdir ${OUTDIR} \\
    --flank_distances "${FLANK_DISTANCES}" \\
    --promoter_upstream ${PROMOTER_UPSTREAM} \\
    --promoter_downstream ${PROMOTER_DOWNSTREAM}

echo "Finished: \$(date)"
SLURM_EOF

chmod +x "${JOB_SCRIPT}"
echo "Generated retry script: ${JOB_SCRIPT}"
echo ""

# ---- Submit ----
if [[ "${DRY_RUN}" == true ]]; then
    echo "=== DRY RUN ==="
    echo "Would run: sbatch --array=1-${N_MISSING}%${RETRY_MAX_CONCURRENT} ${JOB_SCRIPT}"
else
    JOB_ID=$(sbatch --array=1-${N_MISSING}%${RETRY_MAX_CONCURRENT} \
        "${JOB_SCRIPT}" | awk '{print $NF}')
    echo "=============================================="
    echo "Submitted retry job: ${JOB_ID}"
    echo "Array:               1-${N_MISSING}"
    echo "=============================================="
    echo ""
    echo "Monitor: squeue -u \$(whoami) -j ${JOB_ID}"
    echo "Logs:    ${OUTDIR}/logs/"
fi
