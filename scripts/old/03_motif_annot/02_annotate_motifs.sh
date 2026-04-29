#!/bin/bash
# =============================================================================
# 02_annotate_motifs.sh — Submit SLURM array job for motif annotation
#
# Discovers all motif BED files in MOTIF_DIR (from config.sh), creates a file
# list, generates a SLURM array script with parameters substituted, and
# submits it (or dry-runs).
#
# USAGE:
#   bash 02_annotate_motifs.sh                    # dry run
#   bash 02_annotate_motifs.sh --apply            # submit
#   bash 02_annotate_motifs.sh --apply --hybrids="F121-9,BL6xCAST"
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# ---- Parse args (override config defaults) ----
APPLY=false
for arg in "$@"; do
    case "$arg" in
        --apply) APPLY=true ;;
        --hybrids=*) HYBRIDS="${arg#*=}" ;;
        --flank-distances=*) FLANK_DISTANCES="${arg#*=}" ;;
        --dry-run) APPLY=false ;;  # explicit dry-run keyword
    esac
done

mkdir -p "${OUTDIR}/logs"

echo "=============================================="
echo "Motif Annotation — Job Submission"
echo "Date: $(date)"
echo "=============================================="
echo "  Motif dir:        ${MOTIF_DIR}"
echo "  Genome FASTA:     ${GENOME_FA}"
echo "  Preprocessed dir: ${PREPROCESSED_DIR}"
echo "  VCF dir:          ${VCF_DIR}"
echo "  Output dir:       ${OUTDIR}"
echo "  Hybrids:          ${HYBRIDS}"
echo "  Flank distances:  ${FLANK_DISTANCES}"
echo "  Promoter window:  [-${PROMOTER_UPSTREAM}, +${PROMOTER_DOWNSTREAM}]"
echo "  Apply:            ${APPLY}"
echo ""

# ---- Discover motif files ----
FILELIST="${OUTDIR}/motif_file_list.txt"
ls "${MOTIF_DIR}"/*.bed > "${FILELIST}" 2>/dev/null || true
N_FILES=$(wc -l < "${FILELIST}")
echo "Found ${N_FILES} motif BED files."
if [[ ${N_FILES} -eq 0 ]]; then
    echo "ERROR: No .bed files in ${MOTIF_DIR}" >&2
    exit 1
fi
echo "  First 3: $(head -3 ${FILELIST} | xargs -n1 basename | tr '\n' ' ')"
echo "  Last 3:  $(tail -3 ${FILELIST} | xargs -n1 basename | tr '\n' ' ')"
echo ""

# ---- Prerequisite checks ----
echo "Checking prerequisites..."
for f in gencode_vM38_txdb.sqlite tx_metadata.rds tss_granges.rds tss_unique_granges.rds; do
    if [[ ! -f "${PREPROCESSED_DIR}/${f}" ]]; then
        echo "ERROR: Missing preprocessed file: ${PREPROCESSED_DIR}/${f}" >&2
        echo "       Run 01_run_preprocess_resources.sh first." >&2
        exit 1
    fi
done
[[ -f "${GENOME_FA}" ]] || { echo "ERROR: missing ${GENOME_FA}" >&2; exit 1; }
[[ -f "${GENOME_FA}.fai" ]] || { echo "ERROR: missing ${GENOME_FA}.fai (run samtools faidx)" >&2; exit 1; }
[[ -f "${SCRIPT_DIR}/02_annotate_motifs.R" ]] || { echo "ERROR: missing 02_annotate_motifs.R in ${SCRIPT_DIR}" >&2; exit 1; }
echo "  All prerequisites OK."
echo ""

# ---- Generate the array job script ----
JOB_SCRIPT="${OUTDIR}/run_annotate_array.sh"

cat > "${JOB_SCRIPT}" << SLURM_EOF
#!/bin/bash
#SBATCH --job-name=annot_motif
#SBATCH --partition=${DEFAULT_PARTITION}
#SBATCH --time=${DEFAULT_TIME}
#SBATCH --mem=${DEFAULT_MEM}
#SBATCH --cpus-per-task=1
#SBATCH --output=${OUTDIR}/logs/annot_motif_%A_%a.out
#SBATCH --error=${OUTDIR}/logs/annot_motif_%A_%a.err

set -euo pipefail

echo "================================================================"
echo "Motif Annotation — Array Task \${SLURM_ARRAY_TASK_ID}"
echo "Job ID:    \${SLURM_ARRAY_JOB_ID}_\${SLURM_ARRAY_TASK_ID}"
echo "Node:      \$(hostname)"
echo "Date:      \$(date)"
echo "================================================================"

MOTIF_BED=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${FILELIST}")
if [[ -z "\${MOTIF_BED}" ]]; then
    echo "ERROR: No file for array task \${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi
MOTIF_ID=\$(basename "\${MOTIF_BED}" .bed)
echo "Motif: \${MOTIF_ID}"

# Skip if output already exists (idempotent)
OUT_FILE="${OUTDIR}/\${MOTIF_ID}_annotated.tsv.gz"
if [[ -f "\${OUT_FILE}" && -s "\${OUT_FILE}" ]]; then
    echo "Output exists (\${OUT_FILE}); skipping. Delete to re-run."
    exit 0
fi

set +u
source activate ${CONDA_ENV}
set -u

Rscript ${SCRIPT_DIR}/02_annotate_motifs.R \\
    --motif_bed "\${MOTIF_BED}" \\
    --genome_fa ${GENOME_FA} \\
    --preprocessed_dir ${PREPROCESSED_DIR} \\
    --vcf_dir ${VCF_DIR} \\
    --hybrids "${HYBRIDS}" \\
    --flank_distances "${FLANK_DISTANCES}" \\
    --promoter_upstream ${PROMOTER_UPSTREAM} \\
    --promoter_downstream ${PROMOTER_DOWNSTREAM} \\
    --outdir ${OUTDIR}

echo ""
echo "Task \${SLURM_ARRAY_TASK_ID} done: \$(date)"
SLURM_EOF

chmod +x "${JOB_SCRIPT}"
echo "Generated array script: ${JOB_SCRIPT}"
echo ""

# ---- Submit ----
if [[ "${APPLY}" == false ]]; then
    echo "=== DRY RUN ==="
    echo "Would submit:"
    echo "  sbatch --array=1-${N_FILES}%${DEFAULT_MAX_CONCURRENT} ${JOB_SCRIPT}"
    echo ""
    echo "Re-run with --apply to submit."
else
    JOB_ID=$(sbatch --array=1-${N_FILES}%${DEFAULT_MAX_CONCURRENT} \
        "${JOB_SCRIPT}" | awk '{print $NF}')
    echo "=============================================="
    echo "Submitted: ${JOB_ID}"
    echo "Array:     1-${N_FILES} (max ${DEFAULT_MAX_CONCURRENT} concurrent)"
    echo "=============================================="
    echo "Monitor: squeue -u \$(whoami) -j ${JOB_ID}"
    echo "Logs:    ${OUTDIR}/logs/"
fi
