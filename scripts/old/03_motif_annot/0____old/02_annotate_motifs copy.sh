#!/bin/bash
# =============================================================================
# 03_submit_all.sh — Submit SLURM array jobs for motif annotation
#
# This script:
#   1. Discovers all motif archetype BED files
#   2. Creates a file list for the SLURM array
#   3. Submits a SLURM array job that runs 02_annotate_motifs.R per archetype
#
# USAGE:
#   bash 03_submit_all.sh [--dry-run] [--hybrids "F121-9,BL6xCAST"]
#
# ASSUMPTIONS:
#   - 01_preprocess_resources.R has already been run successfully
#   - Conda environment "test_m" exists and has all required packages
#   - All paths below are correct for your setup
# =============================================================================

set -euo pipefail

# ---- Configuration (EDIT THESE) ----
MOTIF_DIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/moods_scan/merged/per_motif"
GENOME_FA="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/genome/mm39.fa"
PREPROCESSED_DIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/annot_motifs/preprocessed"
VCF_DIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources/vcf"
OUTDIR="/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/annot_motifs"
HYBRIDS="F121-9,BL6xCAST"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"  # Directory containing this script

# SLURM settings
PARTITION="short"
TIME="01:00:00"
MEM="16G"
MAX_CONCURRENT=1000  # Max simultaneous array tasks

# ---- Parse optional arguments ----
DRY_RUN=false
for arg in "$@"; do
    case $arg in
        --dry-run)
            DRY_RUN=true
            ;;
        --hybrids=*)
            HYBRIDS="${arg#*=}"
            ;;
    esac
done

# ---- Create directories ----
mkdir -p "${OUTDIR}/logs"

# ---- Discover motif files ----
echo "=============================================="
echo "Motif Annotation Pipeline — Job Submission"
echo "Date: $(date)"
echo "=============================================="
echo ""
echo "Configuration:"
echo "  Motif dir:        ${MOTIF_DIR}"
echo "  Genome FASTA:     ${GENOME_FA}"
echo "  Preprocessed dir: ${PREPROCESSED_DIR}"
echo "  VCF dir:          ${VCF_DIR}"
echo "  Output dir:       ${OUTDIR}"
echo "  Hybrids:          ${HYBRIDS}"
echo "  Script dir:       ${SCRIPT_DIR}"
echo ""

# Create file list for the array
FILELIST="${OUTDIR}/motif_file_list.txt"
ls "${MOTIF_DIR}"/*.bed > "${FILELIST}" 2>/dev/null || true

N_FILES=$(wc -l < "${FILELIST}")
echo "Found ${N_FILES} motif BED files."
echo ""

if [[ ${N_FILES} -eq 0 ]]; then
    echo "ERROR: No .bed files found in ${MOTIF_DIR}"
    exit 1
fi

# Show first and last few files
echo "First 5 files:"
head -5 "${FILELIST}" | while read f; do echo "  $(basename $f)"; done
echo "..."
echo "Last 5 files:"
tail -5 "${FILELIST}" | while read f; do echo "  $(basename $f)"; done
echo ""

# ---- Check prerequisites ----
echo "Checking prerequisites..."

# Check preprocessed files exist
for f in gencode_vM38_txdb.sqlite tx_metadata.rds tss_granges.rds tss_unique_granges.rds; do
    if [[ ! -f "${PREPROCESSED_DIR}/${f}" ]]; then
        echo "ERROR: Missing preprocessed file: ${PREPROCESSED_DIR}/${f}"
        echo "       Run 01_preprocess_resources.R first."
        exit 1
    fi
done
echo "  Preprocessed files: OK"

# Check genome FASTA and index
if [[ ! -f "${GENOME_FA}" ]]; then
    echo "ERROR: Genome FASTA not found: ${GENOME_FA}"
    exit 1
fi
if [[ ! -f "${GENOME_FA}.fai" ]]; then
    echo "ERROR: Genome FASTA index not found: ${GENOME_FA}.fai"
    echo "       Run: samtools faidx ${GENOME_FA}"
    exit 1
fi
echo "  Genome FASTA: OK"

# Check annotation script exists
if [[ ! -f "${SCRIPT_DIR}/02_annotate_motifs.R" ]]; then
    echo "ERROR: Annotation script not found: ${SCRIPT_DIR}/02_annotate_motifs.R"
    exit 1
fi
echo "  Annotation script: OK"
echo ""

# ---- Create SLURM array job script ----
JOB_SCRIPT="${OUTDIR}/run_annotate_array.sh"

cat > "${JOB_SCRIPT}" << 'SLURM_HEREDOC'
#!/bin/bash
#SBATCH --job-name=annot_motif
#SBATCH --partition=PARTITION_PLACEHOLDER
#SBATCH --time=TIME_PLACEHOLDER
#SBATCH --mem=MEM_PLACEHOLDER
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/annot_motif_%A_%a.out
#SBATCH --error=logs/annot_motif_%A_%a.err

set -euo pipefail

echo "================================================================"
echo "Motif Annotation — Array Task ${SLURM_ARRAY_TASK_ID}"
echo "Job ID:    ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Node:      $(hostname)"
echo "Date:      $(date)"
echo "================================================================"
echo ""

# Get the motif file for this array task
FILELIST="FILELIST_PLACEHOLDER"
MOTIF_BED=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${FILELIST}")

if [[ -z "${MOTIF_BED}" ]]; then
    echo "ERROR: No file found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Motif BED file: ${MOTIF_BED}"
echo "Motif ID:       $(basename ${MOTIF_BED} .bed)"
echo ""

# Check if output already exists (for restartability)
MOTIF_ID=$(basename "${MOTIF_BED}" .bed)
OUT_FILE="OUTDIR_PLACEHOLDER/${MOTIF_ID}_annotated.tsv.gz"
if [[ -f "${OUT_FILE}" ]]; then
    echo "Output already exists: ${OUT_FILE}"
    echo "Skipping. Delete the file to re-run."
    exit 0
fi

# Activate conda environment
source activate test_m

# Run the annotation script
Rscript SCRIPT_DIR_PLACEHOLDER/02_annotate_motifs.R \
    --motif_bed "${MOTIF_BED}" \
    --genome_fa GENOME_FA_PLACEHOLDER \
    --preprocessed_dir PREPROCESSED_DIR_PLACEHOLDER \
    --vcf_dir VCF_DIR_PLACEHOLDER \
    --hybrids "HYBRIDS_PLACEHOLDER" \
    --outdir OUTDIR_PLACEHOLDER

echo ""
echo "Array task ${SLURM_ARRAY_TASK_ID} finished: $(date)"
SLURM_HEREDOC

# Replace placeholders with actual values
sed -i "s|PARTITION_PLACEHOLDER|${PARTITION}|g" "${JOB_SCRIPT}"
sed -i "s|TIME_PLACEHOLDER|${TIME}|g" "${JOB_SCRIPT}"
sed -i "s|MEM_PLACEHOLDER|${MEM}|g" "${JOB_SCRIPT}"
sed -i "s|OUTDIR_PLACEHOLDER|${OUTDIR}|g" "${JOB_SCRIPT}"
sed -i "s|FILELIST_PLACEHOLDER|${FILELIST}|g" "${JOB_SCRIPT}"
sed -i "s|SCRIPT_DIR_PLACEHOLDER|${SCRIPT_DIR}|g" "${JOB_SCRIPT}"
sed -i "s|GENOME_FA_PLACEHOLDER|${GENOME_FA}|g" "${JOB_SCRIPT}"
sed -i "s|PREPROCESSED_DIR_PLACEHOLDER|${PREPROCESSED_DIR}|g" "${JOB_SCRIPT}"
sed -i "s|VCF_DIR_PLACEHOLDER|${VCF_DIR}|g" "${JOB_SCRIPT}"
sed -i "s|HYBRIDS_PLACEHOLDER|${HYBRIDS}|g" "${JOB_SCRIPT}"

echo "Generated SLURM array script: ${JOB_SCRIPT}"
echo ""

# ---- Submit or show ----
if [[ "${DRY_RUN}" == true ]]; then
    echo "=== DRY RUN — Would submit: ==="
    echo "  sbatch --array=1-${N_FILES}%${MAX_CONCURRENT} ${JOB_SCRIPT}"
    echo ""
    echo "Array job script contents:"
    echo "---"
    cat "${JOB_SCRIPT}"
    echo "---"
else
    echo "Submitting SLURM array job..."
    JOB_ID=$(sbatch --array=1-${N_FILES}%${MAX_CONCURRENT} "${JOB_SCRIPT}" | awk '{print $NF}')
    echo ""
    echo "=============================================="
    echo "Submitted array job: ${JOB_ID}"
    echo "Array size:          1-${N_FILES}"
    echo "Max concurrent:      ${MAX_CONCURRENT}"
    echo "=============================================="
    echo ""
    echo "Monitor with:"
    echo "  squeue -u \$(whoami) -j ${JOB_ID}"
    echo "  sacct -j ${JOB_ID} --format=JobID,State,Elapsed,MaxRSS"
    echo ""
    echo "Check logs in:"
    echo "  ${OUTDIR}/logs/"
    echo ""
    echo "Output files will appear in:"
    echo "  ${OUTDIR}/"
fi