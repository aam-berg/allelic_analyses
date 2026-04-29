#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 08_compute_wasp_balance.sh — WASP cascade + per-SNP allelic depth
# =============================================================================
#
# Computes three pieces of data needed by 09_plot_wasp_balance.R:
#
# (A) WASP filter cascade — per assay × replicate, read counts at:
#       total_input    — pre-WASP BAM, all primary mapped reads
#       wasp_passed    — post-WASP BAM, vW:i:1-tagged reads (passed allele filter)
#       allele_split   — sum of {SAMPLE}_ref.bam + {SAMPLE}_alt.bam reads
#     Output: ${WASP_DIR}/wasp_cascade.tsv
#
# (B) Global allele balance — per assay × replicate, ref vs alt total reads.
#     Computed from the allele-specific BAMs.
#     Output: ${WASP_DIR}/allele_balance.tsv
#
# (C) Per-SNP allelic depth — bcftools mpileup at every unique SNP that
#     appears in 06_allele_pairs, across all 10 BAMs. Counts per SNP per
#     sample.
#     Output: ${WASP_DIR}/per_snp_coverage.tsv.gz
#
# USAGE:
#   bash 08_compute_wasp_balance.sh
# =============================================================================

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
step_header "08_compute_wasp_balance.sh — WASP + coverage diagnostics"

mkdir -p "${WASP_DIR}"

# Helpful aliases
SAMTOOLS=$(command -v samtools 2>/dev/null || true)
BCFTOOLS=$(command -v bcftools 2>/dev/null || true)

if [[ -z "${SAMTOOLS}" || -z "${BCFTOOLS}" ]]; then
    echo "[ERROR] samtools and bcftools required. Activate ${CONDA_ENV_TOOLS}."
    exit 1
fi

# -----------------------------------------------------------------------------
# (A) WASP cascade — read counts at each stage
# -----------------------------------------------------------------------------
# We enumerate all BAMs at each stage and read-count via samtools view -c.
# Pre-WASP BAMs may not all be present in scratch (depends on whether 03/04/04c
# kept their intermediates). We're lenient: missing stage = NA.

CASCADE_OUT="${WASP_DIR}/wasp_cascade.tsv"
echo -e "assay\tsample\tstage\tn_reads" > "${CASCADE_OUT}"

count_reads() {
    local bam="$1"
    if [[ -f "${bam}" ]]; then
        "${SAMTOOLS}" view -c -F 256 -F 2048 "${bam}" 2>/dev/null || echo "NA"
    else
        echo "NA"
    fi
}

count_wasp_passed() {
    # Reads with vW:i:1 tag (WASP allele filter passed)
    local bam="$1"
    if [[ -f "${bam}" ]]; then
        "${SAMTOOLS}" view -F 256 -F 2048 -d vW:1 -c "${bam}" 2>/dev/null || echo "NA"
    else
        echo "NA"
    fi
}

count_for_sample() {
    local assay="$1" sample="$2"
    local pre_bam_dir="$3" allele_bam_dir="$4"

    # Pre-WASP (or WASP-tagged but unfiltered) BAM
    local pre_bam=""
    for candidate in "${pre_bam_dir}/${sample}.bam" \
                      "${pre_bam_dir}/${sample}_wasp.bam" \
                      "${pre_bam_dir}/${sample}.sorted.bam"; do
        if [[ -f "${candidate}" ]]; then pre_bam="${candidate}"; break; fi
    done

    local total=$(count_reads "${pre_bam:-/nonexistent}")
    local wasp_pass=$(count_wasp_passed "${pre_bam:-/nonexistent}")

    local ref_bam="${allele_bam_dir}/${sample}_ref.bam"
    local alt_bam="${allele_bam_dir}/${sample}_alt.bam"
    local ref_n=$(count_reads "${ref_bam}")
    local alt_n=$(count_reads "${alt_bam}")
    local split_total
    if [[ "${ref_n}" == "NA" || "${alt_n}" == "NA" ]]; then
        split_total="NA"
    else
        split_total=$((ref_n + alt_n))
    fi

    {
        echo -e "${assay}\t${sample}\ttotal_input\t${total}"
        echo -e "${assay}\t${sample}\twasp_passed\t${wasp_pass}"
        echo -e "${assay}\t${sample}\tallele_split_total\t${split_total}"
        echo -e "${assay}\t${sample}\tref\t${ref_n}"
        echo -e "${assay}\t${sample}\talt\t${alt_n}"
    } >> "${CASCADE_OUT}"
}

echo ""
echo "--- (A) WASP cascade per assay/sample ---"
for s in "${ATAC_SAMPLES[@]}"; do
    echo "  ATAC: ${s}"
    count_for_sample "ATAC" "${s}" "${ATAC_BAM_DIR}" "${ATAC_ALLELE_BAM_DIR}"
done
for s in "${PROSEQ_SAMPLES[@]}"; do
    echo "  PRO-seq: ${s}"
    count_for_sample "PRO-seq" "${s}" "${PROSEQ_BAM_DIR}" "${PROSEQ_ALLELE_BAM_DIR}"
done
for s in "${RNASEQ_SAMPLES[@]}"; do
    echo "  RNA-seq: ${s}"
    count_for_sample "RNA-seq" "${s}" "${RNASEQ_BAM_DIR}" "${RNASEQ_ALLELE_BAM_DIR}"
done

echo "Wrote: ${CASCADE_OUT}"

# -----------------------------------------------------------------------------
# (B) Allele balance (ref vs alt totals) — already in (A) but pivot-friendly
# -----------------------------------------------------------------------------
BALANCE_OUT="${WASP_DIR}/allele_balance.tsv"
echo -e "assay\tsample\tn_ref\tn_alt\tref_frac" > "${BALANCE_OUT}"

awk -F'\t' '
NR == 1 {next}
$3 == "ref" { ref[$1"\t"$2] = $4 }
$3 == "alt" { alt[$1"\t"$2] = $4 }
END {
    for (k in ref) {
        if (ref[k] == "NA" || alt[k] == "NA") {
            print k "\t" ref[k] "\t" alt[k] "\tNA"
        } else {
            tot = ref[k] + alt[k]
            if (tot == 0) {
                print k "\t" ref[k] "\t" alt[k] "\tNA"
            } else {
                printf "%s\t%s\t%s\t%.4f\n", k, ref[k], alt[k], ref[k]/tot
            }
        }
    }
}' "${CASCADE_OUT}" >> "${BALANCE_OUT}"

echo "Wrote: ${BALANCE_OUT}"

# -----------------------------------------------------------------------------
# (C) Per-SNP allelic depth — bcftools mpileup at unique SNPs in 06 pair index
# -----------------------------------------------------------------------------
echo ""
echo "--- (C) Per-SNP allelic depth ---"

UNIQUE_SNP_BED="${WASP_DIR}/unique_snps.bed"
PER_SNP_OUT="${WASP_DIR}/per_snp_coverage.tsv.gz"

# Build unique SNP BED from 06_allele_pairs pair indices
if [[ ! -d "${PAIR_INDEX_DIR}" ]] || \
   [[ -z "$(ls -A "${PAIR_INDEX_DIR}"/*_pair_index.tsv.gz 2>/dev/null)" ]]; then
    echo "  [SKIP] No pair indices found at ${PAIR_INDEX_DIR}; per-SNP coverage skipped."
    echo "  Re-run after 06_allele_pairs/ completes."
    exit 0
fi

echo "  Building unique-SNP BED from pair indices..."
# pair_index headers contain snp_chrom, snp_pos
# Extract chrom+pos columns -> 0-based BED
TMP_SNPS="${WASP_DIR}/_tmp_snps.txt"
> "${TMP_SNPS}"
for f in "${PAIR_INDEX_DIR}"/*_pair_index.tsv.gz; do
    zcat "${f}" | awk -F'\t' '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "snp_chrom") c = i;
                if ($i == "snp_pos")   p = i;
            }
            next
        }
        c && p && $c != "" && $p != "" { print $c "\t" ($p - 1) "\t" $p }
    ' >> "${TMP_SNPS}"
done
sort -u -k1,1 -k2,2n "${TMP_SNPS}" > "${UNIQUE_SNP_BED}"
rm -f "${TMP_SNPS}"

N_SNPS=$(wc -l < "${UNIQUE_SNP_BED}")
echo "  ${N_SNPS} unique SNPs to pile up"

if (( N_SNPS == 0 )); then
    echo "  [SKIP] No SNPs in pair indices."
    exit 0
fi

# bcftools mpileup at every SNP, per BAM. Output: chrom, pos, sample, n_ref, n_alt
PILE_TMP="${WASP_DIR}/_pile_combined.tsv"
> "${PILE_TMP}"
echo -e "chrom\tpos\tassay\tsample\tallele\tn_reads" > "${PILE_TMP}"

run_pileup() {
    local bam="$1" assay="$2" sample="$3" allele="$4"
    if [[ ! -f "${bam}" ]]; then return; fi

    # bcftools mpileup -a AD writes per-position allele depths. We sum across
    # alleles ignoring the call — just need read depth per BAM at each SNP.
    "${BCFTOOLS}" mpileup \
        -f "${GENOME_FA}" \
        -R "${UNIQUE_SNP_BED}" \
        -B -q 1 -Q 13 \
        -a "FORMAT/AD,FORMAT/DP" \
        --no-version --max-depth 1000 \
        "${bam}" 2>/dev/null | \
    awk -v assay="${assay}" -v sample="${sample}" -v allele="${allele}" '
        BEGIN { OFS = "\t" }
        /^#/ { next }
        {
            # Find DP in FORMAT and value in samples
            n = split($9, fmt, ":");
            dp_idx = 0;
            for (i = 1; i <= n; i++) if (fmt[i] == "DP") dp_idx = i;
            if (dp_idx == 0) next;
            split($10, vals, ":");
            dp = vals[dp_idx];
            if (dp == "." || dp == "") dp = 0;
            print $1, $2, assay, sample, allele, dp
        }
    ' >> "${PILE_TMP}"
}

echo ""
echo "  Running mpileup across allele-specific BAMs..."
for s in "${ATAC_SAMPLES[@]}"; do
    run_pileup "${ATAC_ALLELE_BAM_DIR}/${s}_ref.bam" "ATAC" "${s}" "ref"
    run_pileup "${ATAC_ALLELE_BAM_DIR}/${s}_alt.bam" "ATAC" "${s}" "alt"
done
for s in "${PROSEQ_SAMPLES[@]}"; do
    run_pileup "${PROSEQ_ALLELE_BAM_DIR}/${s}_ref.bam" "PRO-seq" "${s}" "ref"
    run_pileup "${PROSEQ_ALLELE_BAM_DIR}/${s}_alt.bam" "PRO-seq" "${s}" "alt"
done
for s in "${RNASEQ_SAMPLES[@]}"; do
    run_pileup "${RNASEQ_ALLELE_BAM_DIR}/${s}_ref.bam" "RNA-seq" "${s}" "ref"
    run_pileup "${RNASEQ_ALLELE_BAM_DIR}/${s}_alt.bam" "RNA-seq" "${s}" "alt"
done

echo "  Compressing per-SNP table..."
gzip -f -c "${PILE_TMP}" > "${PER_SNP_OUT}"
rm -f "${PILE_TMP}"
echo "  Wrote: ${PER_SNP_OUT}"

echo ""
echo "Done."
