#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 05_wasp_setup.sh — WASP repo + per-chromosome SNP files
# =============================================================================
#
# WHAT THIS DOES:
#   1. Ensures pysam is installed in the conda env
#   2. Clones the WASP repo (if not already)
#   3. Generates per-chromosome SNP files in WASP format from the F121-9
#      het VCF, by calling lib/make_wasp_snps.sh
#
# WHY THIS IS A THIN WRAPPER:
#   The actual VCF -> WASP-SNP-files conversion is in lib/make_wasp_snps.sh
#   so it's shared with 04_proseq/06_wasp_setup.sh. This script's only
#   ATAC-specific job is to plug the right paths into the shared utility.
#
# SBATCH RESOURCES:
#   sbatch --partition=short --time=01:00:00 --mem=8G --cpus-per-task=4 \
#       -o logs/05_%j.out -e logs/05_%j.err 05_wasp_setup.sh
# =============================================================================

source "config.sh"

step_header "ATAC Pipeline Step 05: WASP Setup"

source activate "${CONDA_ENV_NAME}"
create_dirs

# -----------------------------------------------------------------------------
# 1. pysam (WASP dependency)
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

# Verify the two scripts we need are present
for SCRIPT in "${WASP_FIND}" "${WASP_FILTER}"; do
    if [[ ! -f "${SCRIPT}" ]]; then
        echo "[ERROR] Missing WASP script: ${SCRIPT}" >&2
        exit 1
    fi
    echo "  [OK] $(basename "${SCRIPT}")"
done

# -----------------------------------------------------------------------------
# 3. SNP files (delegate to shared lib script)
# -----------------------------------------------------------------------------
echo ""
echo "[INFO] Generating WASP per-chromosome SNP files..."

# Verify the input VCF exists
if [[ ! -f "${ALLELE_VCF}" ]]; then
    echo "[ERROR] Allele VCF not found: ${ALLELE_VCF}" >&2
    echo "[ERROR] Was 01_resources/03_convert_vcf_to_ucsc.sh run?" >&2
    exit 1
fi

bash "${MAKE_WASP_SNPS_SH}" \
    --vcf "${ALLELE_VCF}" \
    --output_dir "${WASP_SNP_DIR}"

# -----------------------------------------------------------------------------
# 4. Pre-flight check for the next step
# -----------------------------------------------------------------------------
echo ""
echo "[PRE-FLIGHT for 06_wasp_filter.sh]"
ALL_OK=1
echo "  WASP scripts:"
[[ -f "${WASP_FIND}"   ]] && echo "    [OK] $(basename "${WASP_FIND}")"   || { echo "    [MISSING] $(basename "${WASP_FIND}")";   ALL_OK=0; }
[[ -f "${WASP_FILTER}" ]] && echo "    [OK] $(basename "${WASP_FILTER}")" || { echo "    [MISSING] $(basename "${WASP_FILTER}")"; ALL_OK=0; }

NFILES=$(ls "${WASP_SNP_DIR}"/chr*.snps.txt.gz 2>/dev/null | wc -l)
echo "  WASP SNP files:  ${NFILES} chromosome files"
if (( NFILES == 0 )); then ALL_OK=0; fi

echo "  Bowtie2 index:   ${MM39_BT2_IDX}"
[[ -f "${MM39_BT2_IDX}.1.bt2" ]] && echo "    [OK]" || { echo "    [MISSING]"; ALL_OK=0; }

echo "  Final BAMs:"
for SAMPLE in "${SAMPLE_ORDER[@]}"; do
    BAM="${BAM_FINAL_DIR}/${SAMPLE}_final.bam"
    if [[ -f "${BAM}" ]]; then
        echo "    [OK]    ${SAMPLE}: ${BAM}"
    else
        echo "    [MISSING] ${SAMPLE}: ${BAM}"
        ALL_OK=0
    fi
done

if (( ! ALL_OK )); then
    echo ""
    echo "[ERROR] One or more pre-flight checks failed." >&2
    exit 1
fi

step_header "Step 05 COMPLETE"
echo "Next: bash 06_wasp_filter.sh   (or sbatch as array of 4)"
