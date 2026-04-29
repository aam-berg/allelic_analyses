#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 06_allele_pairs/00_setup.sh — Pre-flight input validation
# =============================================================================
#
# Verifies that all upstream pipeline outputs needed by 06_allele_pairs/ are
# present and reports anything missing. Does NOT do heavy work — this is a
# fail-fast check before launching SLURM arrays.
#
# Usage:
#   bash 00_setup.sh
# =============================================================================

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

step_header "06_allele_pairs/00_setup.sh — input validation"

create_dirs

ALL_OK=1
WARN_COUNT=0

check_file() {
    local label="$1"
    local path="$2"
    local required="${3:-required}"

    if [[ -f "${path}" ]]; then
        printf "  [OK]    %-40s %s\n" "${label}:" "${path}"
    else
        if [[ "${required}" == "required" ]]; then
            printf "  [MISS]  %-40s %s  (REQUIRED)\n" "${label}:" "${path}"
            ALL_OK=0
        else
            printf "  [WARN]  %-40s %s  (optional)\n" "${label}:" "${path}"
            WARN_COUNT=$((WARN_COUNT + 1))
        fi
    fi
}

check_dir() {
    local label="$1"
    local path="$2"
    local pattern="${3:-*}"
    local required="${4:-required}"

    if [[ -d "${path}" ]]; then
        local n
        n=$(find "${path}" -maxdepth 1 -name "${pattern}" 2>/dev/null | wc -l)
        if (( n > 0 )); then
            printf "  [OK]    %-40s %s  (%d files match %s)\n" \
                "${label}:" "${path}" "${n}" "${pattern}"
        else
            if [[ "${required}" == "required" ]]; then
                printf "  [MISS]  %-40s %s  (no files match %s; REQUIRED)\n" \
                    "${label}:" "${path}" "${pattern}"
                ALL_OK=0
            else
                printf "  [WARN]  %-40s %s  (no files match %s)\n" \
                    "${label}:" "${path}" "${pattern}"
                WARN_COUNT=$((WARN_COUNT + 1))
            fi
        fi
    else
        if [[ "${required}" == "required" ]]; then
            printf "  [MISS]  %-40s %s  (REQUIRED)\n" "${label}:" "${path}"
            ALL_OK=0
        else
            printf "  [WARN]  %-40s %s\n" "${label}:" "${path}"
            WARN_COUNT=$((WARN_COUNT + 1))
        fi
    fi
}

# -----------------------------------------------------------------------------
# Reference resources
# -----------------------------------------------------------------------------
echo ""
echo "--- Reference resources ---"
check_file "mm39 FASTA"             "${GENOME_FA}"
check_file "mm39 FASTA index"       "${GENOME_FA}.fai"
check_file "F121-9 het VCF"         "${VCF_DIR}/F121-9_het_snps.ucsc.vcf.gz"
check_file "F121-9 het VCF index"   "${VCF_DIR}/F121-9_het_snps.ucsc.vcf.gz.tbi"

# -----------------------------------------------------------------------------
# Motif scanning outputs (02_motif_scanning/)
# -----------------------------------------------------------------------------
echo ""
echo "--- 02_motif_scanning/ outputs ---"
check_dir "Motif BED dir"           "${MOTIF_BED_DIR}" "*.bed"
check_file "PWM MEME file"          "${PWM_MEME_FILE}"

# -----------------------------------------------------------------------------
# 05_motif_annot/ outputs
# -----------------------------------------------------------------------------
echo ""
echo "--- 05_motif_annot/ outputs ---"
check_dir "Annotated motifs dir"    "${MOTIF_ANNOT_DIR}" "*_annotated.tsv.gz"
check_file "Splice donors RDS"      "${SPLICE_DONORS_RDS}"
check_file "Splice acceptors RDS"   "${SPLICE_ACCEPTORS_RDS}"
check_file "Genes GR RDS"           "${GENES_GR_RDS}"

# -----------------------------------------------------------------------------
# 03_atac/ outputs
# -----------------------------------------------------------------------------
echo ""
echo "--- 03_atac/ allele-specific BAMs ---"
for SAMPLE in "${ATAC_SAMPLES[@]}"; do
    for ALLELE in ref alt; do
        check_file "ATAC ${SAMPLE} ${ALLELE}" \
            "${ATAC_ALLELE_BAM_DIR}/${SAMPLE}_${ALLELE}.bam"
    done
done

# -----------------------------------------------------------------------------
# 04_proseq/ outputs
# -----------------------------------------------------------------------------
echo ""
echo "--- 04_proseq/ allele-specific BAMs and bigWigs ---"
for SAMPLE in "${PROSEQ_SAMPLES[@]}"; do
    for ALLELE in ref alt; do
        check_file "PRO-seq ${SAMPLE} ${ALLELE} BAM" \
            "${PROSEQ_ALLELE_BAM_DIR}/${SAMPLE}_${ALLELE}.bam"
        # Strand-specific bigwigs (canonical strand-swapped) for profiles.
        # These are produced by 04_proseq/'s 08_make_bigwigs_allele step.
        for STRAND in plus minus; do
            check_file "PRO-seq ${SAMPLE} ${ALLELE} ${STRAND} bw" \
                "${PROSEQ_ALLELE_BIGWIG_DIR}/${SAMPLE}_${ALLELE}_${STRAND}.bw" \
                optional
        done
    done
done

# -----------------------------------------------------------------------------
# 04c_rnaseq/ outputs
# -----------------------------------------------------------------------------
echo ""
echo "--- 04c_rnaseq/ allele-specific BAMs ---"
for SAMPLE in "${RNASEQ_SAMPLES[@]}"; do
    for ALLELE in ref alt; do
        check_file "RNA-seq ${SAMPLE} ${ALLELE} BAM" \
            "${RNASEQ_ALLELE_BAM_DIR}/${SAMPLE}_${ALLELE}.bam"
    done
done

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "=========================================="
if (( ALL_OK == 1 )); then
    if (( WARN_COUNT == 0 )); then
        echo "All checks passed. Pipeline ready to run."
    else
        echo "All required inputs present (${WARN_COUNT} optional missing)."
        echo "Pipeline can run; some outputs will be limited or skipped."
    fi
    echo "=========================================="
    exit 0
else
    echo "[ERROR] Required inputs missing. Pipeline cannot run."
    echo "=========================================="
    exit 1
fi
