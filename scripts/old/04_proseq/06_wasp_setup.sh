#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 06_wasp_setup.sh — Install WASP and prepare SNP input files
# =============================================================================
#
# WHAT THIS DOES:
#   1. Installs pysam into the existing conda environment (WASP dependency)
#   2. Clones the WASP repository from GitHub
#   3. Converts the F121-9 het SNP VCF to WASP's per-chromosome SNP format
#
# WHY WASP:
#   Our BAMs are aligned to mm39, which is effectively a BL6 reference genome.
#   F121-9 cells are 129S1/SvImJ × CAST/EiJ F1 hybrids, so at each het SNP:
#     - The 129S1 allele (close to BL6/reference) maps better → bias
#     - The CAST allele (divergent from reference) may mismap → bias
#   WASP removes reads where mapping depends on which allele is present,
#   giving us allele-agnostic, mapping-bias-free BAMs. This is essential
#   for any allele-specific analysis.
#
# WASP SNP FORMAT:
#   WASP expects a directory with one file per chromosome:
#     chr1.snps.txt.gz    chr2.snps.txt.gz    ...
#   Each file is tab-separated: POSITION(1-based) \t REF_ALLELE \t ALT_ALLELE
#
# SBATCH RESOURCES:
#   sbatch --partition=short --time=1:00:00 --mem=8G --cpus-per-task=4 \
#       -o logs/06_%j.out -e logs/06_%j.err 06_wasp_setup.sh
# =============================================================================

source "config_wasp.sh"

step_header "PRO-seq Pipeline Step 06: WASP Setup"

source activate "${CONDA_ENV_NAME}"
create_dirs
create_wasp_dirs

# =============================================================================
# 1. Install pysam (required by WASP)
# =============================================================================
step_header "STEP 1: Install pysam"

if python -c "import pysam" 2>/dev/null; then
    echo "[INFO] pysam already installed: $(python -c 'import pysam; print(pysam.__version__)')"
else
    echo "[INFO] Installing pysam..."
    conda install -y -c bioconda -c conda-forge pysam
    echo "[INFO] pysam installed: $(python -c 'import pysam; print(pysam.__version__)')"
fi

# Also ensure numpy is available (WASP may need it)
if python -c "import numpy" 2>/dev/null; then
    echo "[INFO] numpy already installed."
else
    echo "[INFO] Installing numpy..."
    conda install -y -c conda-forge numpy
fi

# =============================================================================
# 2. Clone WASP
# =============================================================================
step_header "STEP 2: Clone WASP repository"

if [[ -d "${WASP_DIR}" && -f "${WASP_DIR}/mapping/find_intersecting_snps.py" ]]; then
    echo "[INFO] WASP already cloned: ${WASP_DIR}"
else
    echo "[INFO] Cloning WASP from GitHub..."
    mkdir -p "$(dirname "${WASP_DIR}")"
    git clone https://github.com/bmvdgeijn/WASP.git "${WASP_DIR}"
    echo "[INFO] WASP cloned to: ${WASP_DIR}"
fi

# Verify key scripts exist
for SCRIPT in "mapping/find_intersecting_snps.py" "mapping/filter_remapped_reads.py"; do
    if [[ -f "${WASP_DIR}/${SCRIPT}" ]]; then
        echo "  [OK] ${SCRIPT}"
    else
        echo "  [ERROR] Missing: ${SCRIPT}"
        exit 1
    fi
done

# =============================================================================
# 3. Convert F121-9 het SNP VCF to WASP format
# =============================================================================
step_header "STEP 3: Prepare SNP files for WASP"

# Use UCSC-style VCF (has chr-prefixed chromosome names matching our BAMs)
VCF_INPUT="${VCF_DIR}/F121-9_het_snps.ucsc.vcf.gz"

if [[ ! -f "${VCF_INPUT}" ]]; then
    echo "[WARNING] UCSC VCF not found: ${VCF_INPUT}"
    echo "[INFO] Trying regular VCF..."
    VCF_INPUT="${VCF_DIR}/F121-9_het_snps.vcf.gz"
    if [[ ! -f "${VCF_INPUT}" ]]; then
        echo "[ERROR] No F121-9 het SNP VCF found in ${VCF_DIR}"
        exit 1
    fi
fi
echo "[INFO] Using VCF: ${VCF_INPUT}"

# Count total SNPs
TOTAL_SNPS=$(zcat "${VCF_INPUT}" | grep -v '^#' | wc -l)
echo "[INFO] Total het SNPs in VCF: ${TOTAL_SNPS}"

# Check if SNP files already exist
if [[ -d "${WASP_SNP_DIR}" ]] && ls "${WASP_SNP_DIR}"/chr*.snps.txt.gz 1>/dev/null 2>&1; then
    EXISTING=$(ls "${WASP_SNP_DIR}"/chr*.snps.txt.gz | wc -l)
    echo "[INFO] WASP SNP files already exist (${EXISTING} chromosomes). Skipping conversion."
else
    echo "[INFO] Converting VCF to WASP per-chromosome SNP format..."
    mkdir -p "${WASP_SNP_DIR}"

    # Extract SNPs: VCF POS is already 1-based, which is what WASP expects.
    # We only keep biallelic SNPs (single REF and ALT base) on standard chroms.
    # Format: POS \t REF \t ALT
    #
    # Strategy: extract all SNPs to a temp file, then split by chromosome.
    # This is faster than running zcat + grep for each chromosome.

    TMPFILE="${WASP_SNP_DIR}/_all_snps.tmp"

    zcat "${VCF_INPUT}" \
    | grep -v '^#' \
    | awk -v OFS='\t' '
        # Filter: standard chroms, biallelic SNPs only (single-base REF and ALT)
        length($4) == 1 && length($5) == 1 && $1 ~ /^chr[0-9XY]+$/ {
            print $1, $2, $4, $5
        }
    ' > "${TMPFILE}"

    FILTERED_SNPS=$(wc -l < "${TMPFILE}")
    echo "[INFO] Biallelic SNPs on standard chroms: ${FILTERED_SNPS}"

    # Split by chromosome
    for CHROM in $(cut -f1 "${TMPFILE}" | sort -u); do
        OUTFILE="${WASP_SNP_DIR}/${CHROM}.snps.txt.gz"
        awk -v OFS='\t' -v chrom="${CHROM}" '$1 == chrom {print $2, $3, $4}' "${TMPFILE}" \
        | sort -k1,1n \
        | gzip > "${OUTFILE}"

        N=$(zcat "${OUTFILE}" | wc -l)
        echo "  ${CHROM}: ${N} SNPs"
    done

    rm -f "${TMPFILE}"
    echo "[INFO] SNP files written to: ${WASP_SNP_DIR}"
fi

# =============================================================================
# 4. Sanity checks
# =============================================================================
step_header "Sanity Checks"

echo "WASP SNP directory contents:"
ls -lh "${WASP_SNP_DIR}"/chr*.snps.txt.gz | head -5 | sed 's/^/    /'
echo "    ..."
NFILES=$(ls "${WASP_SNP_DIR}"/chr*.snps.txt.gz | wc -l)
echo "  Total chromosome files: ${NFILES}"

echo ""
echo "Sample SNP entries (chr1, first 5):"
zcat "${WASP_SNP_DIR}/chr1.snps.txt.gz" | head -5 | sed 's/^/    /'

TOTAL_WASP_SNPS=$(zcat "${WASP_SNP_DIR}"/chr*.snps.txt.gz | wc -l)
echo ""
echo "Total SNPs in WASP format: ${TOTAL_WASP_SNPS}"

# Verify key things are in place
echo ""
echo "Pre-flight check for WASP pipeline:"
echo "  [1] WASP scripts:"
echo "      find_intersecting_snps.py: $(ls ${WASP_DIR}/mapping/find_intersecting_snps.py 2>/dev/null && echo OK || echo MISSING)"
echo "      filter_remapped_reads.py:  $(ls ${WASP_DIR}/mapping/filter_remapped_reads.py 2>/dev/null && echo OK || echo MISSING)"
echo "  [2] SNP directory: ${WASP_SNP_DIR} (${NFILES} chrom files)"
echo "  [3] Bowtie2 index: ${MM39_BT2_IDX} ($(ls ${MM39_BT2_IDX}.1.bt2 2>/dev/null && echo OK || echo MISSING))"
echo "  [4] pysam: $(python -c 'import pysam; print(pysam.__version__)' 2>/dev/null || echo MISSING)"
echo "  [5] Final BAMs:"
for SAMPLE in "${SAMPLE_ORDER[@]}"; do
    BAM="${BAM_MM39_DIR}/${SAMPLE}_final.bam"
    echo "      ${SAMPLE}: $(ls ${BAM} 2>/dev/null && echo OK || echo MISSING)"
done

step_header "Step 06 COMPLETE"
echo "WASP is ready. Next: Run 07_wasp_filter.sh"