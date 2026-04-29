#!/bin/bash
# =============================================================================
# 0x_resubmit_missing_motifs.sh — Resubmit failed/missing motifs with extra resources
#
# Reads ${OUTDIR}/_missing_motifs.txt (built by 0x_find_missing_motifs.sh)
# and submits a SLURM array using RETRY_PARTITION / RETRY_TIME / RETRY_MEM
# from config.sh — a higher-resource queue meant to absorb timeouts and OOM
# failures from the main run.
#
# USAGE:
#   bash 0x_resubmit_missing_motifs.sh                # dry-run
#   bash 0x_resubmit_missing_motifs.sh --apply
# =============================================================================

set -euo pipefail

SCRIPT_DIR="."
source "${SCRIPT_DIR}/config.sh"

# ---- CLI ----
APPLY=false
PARTITION="${RETRY_PARTITION}"
TIME="${RETRY_TIME}"
MEM="${RETRY_MEM}"
MAX_CONCURRENT="${RETRY_MAX_CONCURRENT}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --apply)            APPLY=true; shift ;;
        --dry-run)          APPLY=false; shift ;;
        --partition)        PARTITION="$2"; shift 2 ;;
        --time)             TIME="$2"; shift 2 ;;
        --mem)              MEM="$2"; shift 2 ;;
        --max-concurrent)   MAX_CONCURRENT="$2"; shift 2 ;;
        -h|--help)
            sed -n '3,15p' "$0" | sed 's/^# \?//'
            exit 0 ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

MISSING_LIST="${OUTDIR}/_missing_motifs.txt"
RETRY_MANIFEST="${OUTDIR}/_retry_manifest.txt"

if [[ ! -f "${MISSING_LIST}" ]]; then
    echo "[ERROR] Missing list not found: ${MISSING_LIST}" >&2
    echo "        Run 0x_find_missing_motifs.sh first." >&2
    exit 1
fi

# ---- Build retry manifest with re-indexed task IDs (0..N-1) ----
> "${RETRY_MANIFEST}"
N_RETRY=0
while IFS=$'\t' read -r _orig_id motif_id bed_path; do
    echo -e "${N_RETRY}\t${motif_id}\t${bed_path}" >> "${RETRY_MANIFEST}"
    N_RETRY=$((N_RETRY + 1))
done < "${MISSING_LIST}"

if (( N_RETRY == 0 )); then
    echo "[OK] No motifs to retry."
    rm -f "${RETRY_MANIFEST}"
    exit 0
fi

echo "================================================================"
echo "Resubmit plan"
echo "================================================================"
echo "  Mode:                $( ${APPLY} && echo APPLY || echo "DRY RUN" )"
echo "  Motifs to retry:     ${N_RETRY}"
echo "  Retry manifest:      ${RETRY_MANIFEST}"
echo "  Partition:           ${PARTITION}"
echo "  Time:                ${TIME}"
echo "  Memory:              ${MEM}"
echo "  Max concurrent:      ${MAX_CONCURRENT}"
echo ""

# ---- Per-task wrapper (separate from primary so logs don't collide) ----
WRAPPER="${OUTDIR}/_retry_one_motif.sh"
cat > "${WRAPPER}" << EOF
#!/bin/bash
#SBATCH --output=logs/retry_%A_%a.out
#SBATCH --error=logs/retry_%A_%a.err
set -euo pipefail

SCRIPT_DIR="${SCRIPT_DIR}"
MANIFEST="${RETRY_MANIFEST}"

LINE=\$(awk -v idx="\${SLURM_ARRAY_TASK_ID}" '\$1 == idx' "\${MANIFEST}")
if [[ -z "\${LINE}" ]]; then
    echo "[ERROR] No retry manifest entry for task \${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi
MOTIF_ID=\$(echo "\${LINE}" | cut -f2)
MOTIF_BED=\$(echo "\${LINE}" | cut -f3)

echo "RETRY task \${SLURM_ARRAY_TASK_ID}: motif=\${MOTIF_ID}"
echo "BED:  \${MOTIF_BED}"
echo "Date: \$(date)"

set +u
source activate "${CONDA_ENV}"
set -u

OUT="${OUTDIR}/\${MOTIF_ID}_annotated.tsv.gz"
if [[ -f "\${OUT}" ]]; then
    echo "[INFO] Already annotated: \${OUT}. Skipping."
    exit 0
fi

Rscript "\${SCRIPT_DIR}/02_annotate_motifs.R" \\
    --motif_bed "\${MOTIF_BED}" \\
    --genome_fa "${GENOME_FA}" \\
    --preprocessed_dir "${PREPROCESSED_DIR}" \\
    --vcf_dir "${VCF_DIR}" \\
    --hybrids "${HYBRIDS}" \\
    --flank_distances "${FLANK_DISTANCES}" \\
    --promoter_upstream ${PROMOTER_UPSTREAM} \\
    --promoter_downstream ${PROMOTER_DOWNSTREAM} \\
    --outdir "${OUTDIR}"

echo ""
echo "RETRY task \${SLURM_ARRAY_TASK_ID} done: \$(date)"
EOF
chmod +x "${WRAPPER}"

ARRAY_SPEC="0-$((N_RETRY - 1))%${MAX_CONCURRENT}"

echo "Would submit:"
echo "  sbatch \\"
echo "    --array=${ARRAY_SPEC} \\"
echo "    --partition=${PARTITION} \\"
echo "    --time=${TIME} \\"
echo "    --mem=${MEM} \\"
echo "    --cpus-per-task=2 \\"
echo "    --job-name=motif_annot_retry \\"
echo "    ${WRAPPER}"

if ${APPLY}; then
    echo ""
    JOB_ID=$(sbatch \
        --array="${ARRAY_SPEC}" \
        --partition="${PARTITION}" \
        --time="${TIME}" \
        --mem="${MEM}" \
        --cpus-per-task=2 \
        --job-name=motif_annot_retry \
        --parsable \
        "${WRAPPER}")
    echo "[OK] Submitted retry as job ${JOB_ID}"
    echo "Monitor: squeue -j ${JOB_ID}"
else
    echo ""
    echo "Dry run; re-run with --apply to actually submit."
fi
