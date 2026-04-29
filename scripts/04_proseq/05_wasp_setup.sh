#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04_proseq/05_wasp_setup.sh — WASP repo + per-chromosome SNP files
# =============================================================================
#
# WHAT THIS DOES:
#   1. Ensures pysam and numpy are installed in the conda env
#   2. Clones the WASP repo (if not already)
#   3. Generates per-chromosome SNP files via the shared
#      lib/make_wasp_snps.sh (used by both 03_atac/ and 04_proseq/)
#
# CHANGES FROM PREVIOUS VERSION:
#   - The inline VCF-to-WASP-format conversion logic that lived inside the
#     old 06_wasp_setup.sh has been moved to lib/make_wasp_snps.sh and is
#     now shared with 03_atac/. This script becomes a thin wrapper.
#
# SBATCH:
#   sbatch --partition=short --time=01:00:00 --mem=8G --cpus-per-task=4 \
#       -o logs/05_%j.out -e logs/05_%j.err 05_wasp_setup.sh
# =============================================================================

source "config.sh"

step_header "PRO-seq Pipeline Step 05: WASP Setup"

source activate "${CONDA_ENV_NAME}"
create_dirs

# -----------------------------------------------------------------------------
# 1. pysam + numpy
# -----------------------------------------------------------------------------
echo "[INFO] Checking pysam..."
if python -c "import pysam" 2>/dev/null; then
    echo "[INFO]   Already installed: $(python -c 'import pysam; print(pysam.__version__)')"
else
    echo "[INFO]   Installing pysam..."
    conda install -y -n "${CONDA_ENV_NAME}" -c bioconda -c conda-forge pysam
fi
if python -c "import numpy" 2>/dev/null; then
    echo "[INFO] numpy already installed."
else
    echo "[INFO] Installing numpy..."
    conda install -y -n "${CONDA_ENV_NAME}" -c conda-forge numpy
fi

# -----------------------------------------------------------------------------
# 2. WASP repo
# -----------------------------------------------------------------------------
echo ""
echo "[INFO] Checking WASP repo..."
if [[ -d "${WASP_REPO_DIR}" && -f "${WASP_FIND}" ]]; then
    echo "[INFO]   Already cloned: ${WASP_REPO_DIR}"
else
    echo "[INFO]   Cloning WASP from GitHub..."
    mkdir -p "$(dirname "${WASP_REPO_DIR}")"
    git clone https://github.com/bmvdgeijn/WASP.git "${WASP_REPO_DIR}"
fi

for SCRIPT in "${WASP_FIND}" "${WASP_FILTER}"; do
    if [[ ! -f "${SCRIPT}" ]]; then
        echo "[ERROR] Missing WASP script: ${SCRIPT}" >&2
        exit 1
    fi
    echo "  [OK] $(basename "${SCRIPT}")"
done

# -----------------------------------------------------------------------------
# 3. SNP files (delegated to shared lib script)
# -----------------------------------------------------------------------------
echo ""
echo "[INFO] Generating WASP per-chromosome SNP files..."

[[ -f "${ALLELE_VCF}" ]] || \
    { echo "[ERROR] Allele VCF missing: ${ALLELE_VCF}" >&2; exit 1; }

bash "${MAKE_WASP_SNPS_SH}" \
    --vcf "${ALLELE_VCF}" \
    --output_dir "${WASP_SNP_DIR}"

# -----------------------------------------------------------------------------
# 4. Pre-flight check for next step
# -----------------------------------------------------------------------------
echo ""
echo "[PRE-FLIGHT CHECK for 06_wasp_filter.sh]"
ALL_OK=1
[[ -f "${WASP_FIND}"   ]] && echo "  [OK] WASP find_intersecting_snps.py" || { echo "  [MISSING]"; ALL_OK=0; }
[[ -f "${WASP_FILTER}" ]] && echo "  [OK] WASP filter_remapped_reads.py"   || { echo "  [MISSING]"; ALL_OK=0; }

NFILES=$(ls "${WASP_SNP_DIR}"/chr*.snps.txt.gz 2>/dev/null | wc -l)
echo "  [INFO] WASP SNP chromosome files: ${NFILES}"
(( NFILES > 0 )) || ALL_OK=0

[[ -f "${MM39_BT2_IDX}.1.bt2" ]] && echo "  [OK] mm39 bowtie2 index" || { echo "  [MISSING]"; ALL_OK=0; }

echo "  Final BAMs:"
for SAMPLE in "${SAMPLE_ORDER[@]}"; do
    BAM="${BAM_MM39_DIR}/${SAMPLE}_final.bam"
    if [[ -f "${BAM}" ]]; then
        echo "    [OK]    ${SAMPLE}"
    else
        echo "    [MISSING] ${SAMPLE}: ${BAM}"
        ALL_OK=0
    fi
done

(( ALL_OK )) || { echo ""; echo "[ERROR] Pre-flight check failed." >&2; exit 1; }

step_header "Step 05 COMPLETE"
echo "Next: bash 06_wasp_filter.sh"
