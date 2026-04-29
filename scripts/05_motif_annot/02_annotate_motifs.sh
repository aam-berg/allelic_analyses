#!/bin/bash
# =============================================================================
# 02_annotate_motifs.sh — Submit per-motif annotation as a SLURM array
#
# Discovers all .bed files in MOTIF_DIR (one per motif archetype) and submits
# a SLURM array job. Each task in the array runs 02_annotate_motifs.R on one
# motif BED and produces ${motif_id}_annotated.tsv.gz in OUTDIR.
#
# USAGE:
#   bash 02_annotate_motifs.sh                # dry-run (default)
#   bash 02_annotate_motifs.sh --apply        # actually submit
#   bash 02_annotate_motifs.sh --apply --max-concurrent 200
# =============================================================================

set -euo pipefail

SCRIPT_DIR="."
source "${SCRIPT_DIR}/config.sh"

# ---- Parse CLI ----
APPLY=false
PARTITION="${DEFAULT_PARTITION}"
TIME="${DEFAULT_TIME}"
MEM="${DEFAULT_MEM}"
MAX_CONCURRENT="${DEFAULT_MAX_CONCURRENT}"

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

# ---- Discover motif BEDs ----
mkdir -p "${OUTDIR}" logs

if [[ ! -d "${MOTIF_DIR}" ]]; then
    echo "[ERROR] MOTIF_DIR not found: ${MOTIF_DIR}" >&2
    exit 1
fi

mapfile -t BED_FILES < <(find "${MOTIF_DIR}" -maxdepth 1 -name "*.bed" -type f | sort)
N_MOTIFS=${#BED_FILES[@]}

if (( N_MOTIFS == 0 )); then
    echo "[ERROR] No .bed files found in ${MOTIF_DIR}" >&2
    exit 1
fi

echo "================================================================"
echo "05_motif_annot/02_annotate_motifs.sh — submission plan"
echo "================================================================"
echo "  Mode:                $( ${APPLY} && echo APPLY || echo "DRY RUN" )"
echo "  Motifs found:        ${N_MOTIFS}"
echo "  Motif dir:           ${MOTIF_DIR}"
echo "  Output dir:          ${OUTDIR}"
echo "  Preprocessed dir:    ${PREPROCESSED_DIR}"
echo "  Hybrids:             ${HYBRIDS}"
echo "  Flank distances:     ${FLANK_DISTANCES}"
echo "  Promoter window:     -${PROMOTER_UPSTREAM} / +${PROMOTER_DOWNSTREAM}"
echo "  Partition / time:    ${PARTITION} / ${TIME}"
echo "  Memory:              ${MEM}"
echo "  Max concurrent:      ${MAX_CONCURRENT}"
echo "================================================================"

# ---- Pre-flight: is preprocessing done? ----
REQUIRED_FILES=(
    "tx_metadata.rds"
    "gene_metadata.rds"
    "transcripts_gr.rds"
    "tss_per_transcript_gr.rds"
    "genic_5UTR_gr.rds"
    "genic_3UTR_gr.rds"
    "genic_CDS_gr.rds"
    "genic_exon_gr.rds"
    "genic_intron_gr.rds"
    "splice_donors_gr.rds"
    "splice_acceptors_gr.rds"
)

MISSING=()
for f in "${REQUIRED_FILES[@]}"; do
    [[ -f "${PREPROCESSED_DIR}/${f}" ]] || MISSING+=("${f}")
done

if (( ${#MISSING[@]} > 0 )); then
    echo ""
    echo "[ERROR] Required preprocessed files missing:" >&2
    for f in "${MISSING[@]}"; do echo "          ${PREPROCESSED_DIR}/${f}" >&2; done
    echo "" >&2
    echo "Run: sbatch 01_run_preprocess_resources.sh" >&2
    exit 1
fi

# Optional resources just warn
for f in atac_consensus_gr.rds rnaseq_expression.rds; do
    if [[ ! -f "${PREPROCESSED_DIR}/${f}" ]]; then
        echo ""
        echo "[WARN] Optional resource missing: ${f}"
        echo "       Annotations from this resource will be NA."
    fi
done

# ---- Build manifest of motifs with task indices ----
MANIFEST="${OUTDIR}/_motif_manifest.txt"
mkdir -p "${OUTDIR}"
> "${MANIFEST}"
for i in "${!BED_FILES[@]}"; do
    bed="${BED_FILES[$i]}"
    motif_id=$(basename "${bed}" .bed)
    echo -e "${i}\t${motif_id}\t${bed}" >> "${MANIFEST}"
done

echo ""
echo "Wrote manifest with ${N_MOTIFS} entries: ${MANIFEST}"
echo "  First few:"
head -3 "${MANIFEST}" | sed 's/^/    /'

# ---- Build the per-task wrapper script ----
WRAPPER="${OUTDIR}/_annotate_one_motif.sh"
cat > "${WRAPPER}" << EOF
#!/bin/bash
#SBATCH --output=logs/annotate_%A_%a.out
#SBATCH --error=logs/annotate_%A_%a.err
set -euo pipefail

SCRIPT_DIR="${SCRIPT_DIR}"
MANIFEST="${MANIFEST}"

# Look up this task's motif from manifest (1-based array → manifest is 0-based)
LINE=\$(awk -v idx="\${SLURM_ARRAY_TASK_ID}" '\$1 == idx' "\${MANIFEST}")
if [[ -z "\${LINE}" ]]; then
    echo "[ERROR] No manifest entry for SLURM_ARRAY_TASK_ID=\${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi
MOTIF_ID=\$(echo "\${LINE}" | cut -f2)
MOTIF_BED=\$(echo "\${LINE}" | cut -f3)

echo "Task \${SLURM_ARRAY_TASK_ID}: motif=\${MOTIF_ID}"
echo "BED:    \${MOTIF_BED}"
echo "Date:   \$(date)"

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
echo "Task \${SLURM_ARRAY_TASK_ID} done: \$(date)"
EOF
chmod +x "${WRAPPER}"

ARRAY_SPEC="0-$((N_MOTIFS - 1))%${MAX_CONCURRENT}"

echo ""
echo "Wrapper script:    ${WRAPPER}"
echo ""
echo "Would submit:"
echo "  sbatch \\"
echo "    --array=${ARRAY_SPEC} \\"
echo "    --partition=${PARTITION} \\"
echo "    --time=${TIME} \\"
echo "    --mem=${MEM} \\"
echo "    --cpus-per-task=2 \\"
echo "    --job-name=motif_annot_v2 \\"
echo "    ${WRAPPER}"

if ${APPLY}; then
    echo ""
    JOB_ID=$(sbatch \
        --array="${ARRAY_SPEC}" \
        --partition="${PARTITION}" \
        --time="${TIME}" \
        --mem="${MEM}" \
        --cpus-per-task=2 \
        --job-name=motif_annot_v2 \
        --parsable \
        "${WRAPPER}")
    echo "[OK] Submitted as job ${JOB_ID}"
    echo "Monitor: squeue -j ${JOB_ID}"
    echo "Logs:    logs/annotate_${JOB_ID}_*.{out,err}"
else
    echo ""
    echo "Dry run; re-run with --apply to actually submit."
fi
