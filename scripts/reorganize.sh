#!/bin/bash
# ============================================================================
# reorganize.sh — One-time directory reorganization
#
# Moves resource/VCF scripts out of motif_scanning/ into a new 01_resources/
# sibling directory, and renames the existing motif_scanning/ to
# 02_motif_scanning/ for consistency with 03_motif_annot/.
#
# Edit the paths at the top of this file before running. The script uses
# `git mv` if you're inside a git repo (preserves history); otherwise plain
# `mv`. Dry-run by default.
#
# Usage:
#   bash reorganize.sh                # show what would happen
#   bash reorganize.sh --apply        # actually do it
# ============================================================================

set -euo pipefail

# ---- EDIT THESE ----
ROOT="/home/alb1273/pausing_phase_project/scripts"
CURRENT_MOTIF_SCANNING_DIR="${ROOT}/motif_scanning"   # whatever it's named today
NEW_RESOURCES_DIR="${ROOT}/01_resources"
NEW_MOTIF_SCANNING_DIR="${ROOT}/02_motif_scanning"
# 03_motif_annot/ is assumed to already exist as ${ROOT}/03_motif_annot

# ---- Parse args ----
APPLY=false
for arg in "$@"; do
    case "$arg" in
        --apply) APPLY=true ;;
    esac
done

# Decide mv command
if [[ "${APPLY}" == true ]]; then
    if (cd "${CURRENT_MOTIF_SCANNING_DIR}" && git rev-parse --is-inside-work-tree &>/dev/null); then
        MV="git mv"
    else
        MV="mv"
    fi
else
    MV="echo [DRY] mv"
fi

echo "=========================================="
echo "Reorganization plan"
echo "=========================================="
echo "  Root:                       ${ROOT}"
echo "  Current motif_scanning dir: ${CURRENT_MOTIF_SCANNING_DIR}"
echo "  → 01_resources/             ${NEW_RESOURCES_DIR}"
echo "  → 02_motif_scanning/        ${NEW_MOTIF_SCANNING_DIR}"
echo "  → (already exists)          ${ROOT}/03_motif_annot"
echo "  Apply:                      ${APPLY}"
echo "  mv command:                 ${MV}"
echo ""

# ---- Step 1: Create 01_resources/ ----
echo "Step 1: Creating ${NEW_RESOURCES_DIR}/ …"
if [[ "${APPLY}" == true ]]; then
    mkdir -p "${NEW_RESOURCES_DIR}"
fi

# ---- Step 2: Move resource/VCF scripts to 01_resources/ ----
echo ""
echo "Step 2: Moving resource/VCF scripts → 01_resources/"
echo "        (and renaming with numeric prefixes for clarity)"
echo ""

# Old name → new name pairs
declare -a MOVES=(
    "download_resources_mm39.sbatch:01_download_resources.sbatch"
    "extract_f1_het_snps.sbatch:02_extract_f1_het_snps.sbatch"
    "convert_vcf_to_ucsc.sh:03_convert_vcf_to_ucsc.sh"
    "vcf_summary_stats.sh:04_vcf_summary_stats.sh"
)

for pair in "${MOVES[@]}"; do
    old_name="${pair%%:*}"
    new_name="${pair##*:}"
    src="${CURRENT_MOTIF_SCANNING_DIR}/${old_name}"
    dst="${NEW_RESOURCES_DIR}/${new_name}"
    if [[ -f "${src}" ]]; then
        echo "  ${old_name}  →  01_resources/${new_name}"
        ${MV} "${src}" "${dst}"
    else
        echo "  (not found, skipping) ${src}"
    fi
done

# ---- Step 3: Rename motif_scanning/ → 02_motif_scanning/ ----
echo ""
echo "Step 3: Renaming motif_scanning dir → 02_motif_scanning"
if [[ "${CURRENT_MOTIF_SCANNING_DIR}" != "${NEW_MOTIF_SCANNING_DIR}" ]]; then
    if [[ -d "${CURRENT_MOTIF_SCANNING_DIR}" ]]; then
        echo "  ${CURRENT_MOTIF_SCANNING_DIR}  →  ${NEW_MOTIF_SCANNING_DIR}"
        ${MV} "${CURRENT_MOTIF_SCANNING_DIR}" "${NEW_MOTIF_SCANNING_DIR}"
    fi
else
    echo "  (already at target name)"
fi

# ---- Step 4: Suggested deletions ----
echo ""
echo "Step 4: Suggested deletions in 02_motif_scanning/"
echo ""
echo "  03_validate_promoter_enrichment.py    # superseded by 03_qc_report.py"
echo ""
echo "  (Run manually after reviewing.)"

# ---- Step 5: Drop in the new files ----
echo ""
echo "Step 5: Drop in the new files I generated alongside this script:"
echo "  → 01_resources/README.md"
echo "  → 02_motif_scanning/README.md (replaces any existing one)"
echo "  → 02_motif_scanning/02_postprocess.py (replaces existing)"
echo "  → 02_motif_scanning/run_pipeline.sh (replaces existing)"
echo "  → 03_motif_annot/README.md"
echo ""

if [[ "${APPLY}" == false ]]; then
    echo "DRY RUN complete. Re-run with --apply to actually move files."
fi
