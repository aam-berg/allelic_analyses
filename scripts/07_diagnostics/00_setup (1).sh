#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 07_diagnostics/00_setup.sh — Pre-flight input validation
# =============================================================================
#
# Default mode (lenient): reports what's missing as warnings; pipeline can
# still run on whatever subset of motifs is currently annotated.
# Strict mode (set STRICT=1): requires all motifs to be present in 05+06.
#
# Usage:
#   bash 00_setup.sh
#   STRICT=1 bash 00_setup.sh
# =============================================================================

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

step_header "07_diagnostics/00_setup.sh — input validation (STRICT=${STRICT})"

create_dirs

ALL_OK=1
WARN_COUNT=0

check_file() {
    local label="$1" path="$2" required="${3:-required}"
    if [[ -f "${path}" ]]; then
        printf "  [OK]    %-40s %s\n" "${label}:" "${path}"
    else
        if [[ "${required}" == "required" ]]; then
            printf "  [MISS]  %-40s %s\n" "${label}:" "${path}"
            ALL_OK=0
        else
            printf "  [WARN]  %-40s %s (optional)\n" "${label}:" "${path}"
            WARN_COUNT=$((WARN_COUNT + 1))
        fi
    fi
}

count_files() {
    local pattern="$1"
    find "$(dirname "${pattern}")" -maxdepth 1 -name "$(basename "${pattern}")" 2>/dev/null | wc -l
}

# -----------------------------------------------------------------------------
echo ""
echo "--- Reference / metadata ---"
check_file "PWM MEME file"          "${PWM_MEME_FILE}"
check_file "TF metadata"            "${TF_METADATA_FILE}"
check_file "Gene expression"        "${GENE_EXPRESSION_FILE}"
check_file "Genome FASTA"           "${GENOME_FA}"
check_file "F121-9 het VCF"         "${VCF_FILE}"

# -----------------------------------------------------------------------------
echo ""
echo "--- 02_motif_scanning/ raw motif BEDs ---"
N_BEDS=$(count_files "${MOTIF_BED_DIR}/*.bed")
if (( N_BEDS == 0 )); then
    printf "  [MISS]  %-40s %s (no .bed files found)\n" \
        "Motif BED dir:" "${MOTIF_BED_DIR}"
    ALL_OK=0
else
    printf "  [OK]    %-40s %d files\n" "Motif BED dir:" "${N_BEDS}"
fi

# -----------------------------------------------------------------------------
echo ""
echo "--- 05_motif_annot/ annotated tables ---"
N_ANNOT=$(count_files "${MOTIF_ANNOT_DIR}/*_annotated.tsv.gz")
if (( N_ANNOT == 0 )); then
    printf "  [MISS]  %-40s %s\n" "Annotated motifs:" "${MOTIF_ANNOT_DIR}"
    ALL_OK=0
else
    if (( N_ANNOT < N_BEDS )); then
        msg="(${N_ANNOT}/${N_BEDS} motifs — partial)"
        if [[ "${STRICT}" == "1" ]]; then
            printf "  [MISS]  %-40s %d/%d motifs (STRICT mode)\n" \
                "Annotated motifs:" "${N_ANNOT}" "${N_BEDS}"
            ALL_OK=0
        else
            printf "  [PART]  %-40s %d/%d motifs (lenient)\n" \
                "Annotated motifs:" "${N_ANNOT}" "${N_BEDS}"
            WARN_COUNT=$((WARN_COUNT + 1))
        fi
    else
        printf "  [OK]    %-40s %d files\n" "Annotated motifs:" "${N_ANNOT}"
    fi
fi

# -----------------------------------------------------------------------------
echo ""
echo "--- 06_allele_pairs/ pair tables ---"
N_PAIR=$(count_files "${PAIR_TABLE_DIR}/*_pair_table.tsv.gz")
N_INDEX=$(count_files "${PAIR_INDEX_DIR}/*_pair_index.tsv.gz")
if (( N_PAIR == 0 )); then
    printf "  [WARN]  %-40s %s (07 plots needing 06 will be skipped)\n" \
        "Pair tables:" "${PAIR_TABLE_DIR}"
    WARN_COUNT=$((WARN_COUNT + 1))
else
    printf "  [OK]    %-40s %d files\n" "Pair tables:" "${N_PAIR}"
fi
if (( N_INDEX > 0 )); then
    printf "  [OK]    %-40s %d files\n" "Pair index:" "${N_INDEX}"
fi

# -----------------------------------------------------------------------------
echo ""
echo "--- BAMs (for WASP/coverage diagnostics) ---"
for SAMPLE in "${ATAC_SAMPLES[@]}"; do
    check_file "ATAC ${SAMPLE} ref" \
        "${ATAC_ALLELE_BAM_DIR}/${SAMPLE}_ref.bam" optional
    check_file "ATAC ${SAMPLE} alt" \
        "${ATAC_ALLELE_BAM_DIR}/${SAMPLE}_alt.bam" optional
done
for SAMPLE in "${PROSEQ_SAMPLES[@]}"; do
    check_file "PRO-seq ${SAMPLE} ref" \
        "${PROSEQ_ALLELE_BAM_DIR}/${SAMPLE}_ref.bam" optional
    check_file "PRO-seq ${SAMPLE} alt" \
        "${PROSEQ_ALLELE_BAM_DIR}/${SAMPLE}_alt.bam" optional
done
for SAMPLE in "${RNASEQ_SAMPLES[@]}"; do
    check_file "RNA-seq ${SAMPLE} ref" \
        "${RNASEQ_ALLELE_BAM_DIR}/${SAMPLE}_ref.bam" optional
    check_file "RNA-seq ${SAMPLE} alt" \
        "${RNASEQ_ALLELE_BAM_DIR}/${SAMPLE}_alt.bam" optional
done

# -----------------------------------------------------------------------------
echo ""
echo "=========================================="
if (( ALL_OK == 1 )); then
    if (( WARN_COUNT == 0 )); then
        echo "All inputs present. Pipeline ready."
    else
        echo "Required inputs present (${WARN_COUNT} optional/partial)."
        echo "Lenient mode: pipeline will run with available data."
    fi
    echo "=========================================="
    exit 0
else
    echo "[ERROR] Required inputs missing. Cannot proceed."
    echo "=========================================="
    exit 1
fi
