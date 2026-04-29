#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04c_rnaseq/04_quantify.sh — featureCounts gene + transcript level + TPM
# =============================================================================
#
# WHAT THIS DOES (per sample):
#   1. Runs featureCounts at the gene level (-t exon -g gene_id) for primary
#      gene quantification.
#   2. Runs featureCounts at the transcript level (-t exon -g transcript_id)
#      for downstream tx-aware analysis (e.g., splice junction usage in
#      06_allele_pairs/).
#   3. Computes per-sample TPM from gene-level counts and saves to a
#      per-sample TPM file. Cross-rep aggregation happens in step 08.
#
# WHY featureCounts (not just STAR's GeneCounts):
#   STAR's --quantMode GeneCounts is a useful sanity check but featureCounts
#   gives more flexibility:
#     - Tunable read assignment behavior (overlap fraction, multi-mappers)
#     - Both gene-level AND transcript-level outputs from a single tool
#     - Standard tool with well-understood behavior across the field
#   We use STAR's output for cross-validation only (in step 08 QC).
#
# WHY -p AND --countReadPairs:
#   -p tells featureCounts the data is paired-end. --countReadPairs counts
#   FRAGMENTS (read pairs) instead of individual reads. For PE RNA-seq this
#   is the right unit of expression measurement.
#
# WHY -s 2:
#   TruSeq Stranded library: R1 maps antisense, so we tell featureCounts
#   the library is reversely-stranded (-s 2).
#
# WHY USE THE FULL UNFILTERED BAM (not WASP-filtered):
#   Gene-level expression should be measured from all valid alignments.
#   WASP filtering would reduce statistical power without benefit because
#   gene counting is not allele-specific. (Allele-specific expression, if
#   needed later, uses the WASP-filtered + allele-split BAMs from step 06.)
#
# PARALLEL EXECUTION:
#   sbatch --array=0-1 --partition=short --time=1:00:00 --mem=8G \
#       --cpus-per-task=8 -o logs/04_%A_%a.out -e logs/04_%A_%a.err \
#       04_quantify.sh
# =============================================================================

source "config.sh"

step_header "RNA-seq Pipeline Step 04: Gene + Transcript Quantification"

source activate "${CONDA_ENV_NAME}"

if ! command -v featureCounts &>/dev/null; then
    echo "[ERROR] featureCounts not found." >&2
    echo "[ERROR] Install: conda install -n ${CONDA_ENV_NAME} -c bioconda subread" >&2
    exit 1
fi
echo "[INFO] featureCounts: $(featureCounts -v 2>&1 | head -1)"

create_dirs
resolve_samples "${1:-}"

# featureCounts wants uncompressed GTF
GTF_USE="${GTF_FILE}"
if [[ "${GTF_FILE}" == *.gz ]]; then
    GTF_DECOMPRESSED="${GTF_FILE%.gz}"
    if [[ ! -f "${GTF_DECOMPRESSED}" ]]; then
        echo "[INFO] Decompressing GTF for featureCounts..."
        gunzip -k "${GTF_FILE}"
    fi
    GTF_USE="${GTF_DECOMPRESSED}"
fi
echo "[INFO] GTF: ${GTF_USE}"

for SAMPLE in "${RUN_SAMPLES[@]}"; do
    step_header "Quantifying ${SAMPLE}"

    BAM="${BAM_STAR_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    [[ -f "${BAM}" ]] || { echo "[ERROR] BAM missing: ${BAM}" >&2; exit 1; }

    GENE_COUNTS="${QUANT_GENE_DIR}/${SAMPLE}_gene_counts.tsv"
    TX_COUNTS="${QUANT_TX_DIR}/${SAMPLE}_tx_counts.tsv"
    GENE_TPM="${QUANT_GENE_DIR}/${SAMPLE}_gene_tpm.tsv"

    # =========================================================================
    # Gene-level quantification
    # =========================================================================
    if [[ -f "${GENE_COUNTS}" ]]; then
        echo "[INFO] Gene counts exist for ${SAMPLE}. Skipping."
    else
        echo "[INFO] Running featureCounts (gene-level)..."
        featureCounts \
            -a "${GTF_USE}" \
            -F GTF \
            -t exon \
            -g gene_id \
            -s "${FEATURECOUNTS_STRAND}" \
            -p --countReadPairs \
            -T "${THREADS}" \
            -o "${GENE_COUNTS}" \
            "${BAM}" \
            2>&1 | tee "${LOG_DIR}/${SAMPLE}_featureCounts_gene.log"
    fi

    # =========================================================================
    # Transcript-level quantification (for downstream tx-aware analyses)
    # =========================================================================
    if [[ -f "${TX_COUNTS}" ]]; then
        echo "[INFO] Tx counts exist for ${SAMPLE}. Skipping."
    else
        echo "[INFO] Running featureCounts (transcript-level)..."
        featureCounts \
            -a "${GTF_USE}" \
            -F GTF \
            -t exon \
            -g transcript_id \
            -s "${FEATURECOUNTS_STRAND}" \
            -p --countReadPairs \
            -T "${THREADS}" \
            -o "${TX_COUNTS}" \
            "${BAM}" \
            2>&1 | tee "${LOG_DIR}/${SAMPLE}_featureCounts_tx.log"
    fi

    # =========================================================================
    # Compute per-sample TPM from gene-level counts
    # =========================================================================
    if [[ -f "${GENE_TPM}" ]]; then
        echo "[INFO] TPM file exists for ${SAMPLE}. Skipping."
    else
        echo "[INFO] Computing TPM..."

        # featureCounts output:
        # Column 1: Geneid; 2: Chr; 3: Start; 4: End; 5: Strand; 6: Length;
        # Column 7: counts (one column per BAM input — we only have one)
        #
        # TPM formula:
        #   rate_i = counts_i / length_i_kb
        #   TPM_i = rate_i / sum(rate) * 1e6

        python3 - "${GENE_COUNTS}" "${GENE_TPM}" << 'PYEOF'
import sys

counts_path, tpm_path = sys.argv[1], sys.argv[2]

# Read featureCounts output: 1 comment line, then header, then data
genes = []
with open(counts_path) as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('Geneid'):
            continue   # header
        parts = line.rstrip('\n').split('\t')
        gene_id = parts[0]
        length = int(parts[5])
        count = int(parts[6])
        genes.append((gene_id, length, count))

# Compute rates: counts per kb of effective length
rates = []
for gene_id, length, count in genes:
    if length > 0:
        rate = count / (length / 1000.0)
    else:
        rate = 0.0
    rates.append(rate)

total_rate = sum(rates)
scale = 1e6 / total_rate if total_rate > 0 else 0

with open(tpm_path, 'w') as out:
    out.write("gene_id\tlength\tcounts\ttpm\n")
    for (gene_id, length, count), rate in zip(genes, rates):
        tpm = rate * scale
        out.write(f"{gene_id}\t{length}\t{count}\t{tpm:.4f}\n")

# Sanity: TPM should sum to 1e6
total_tpm = sum(rate * scale for rate in rates)
print(f"  Genes: {len(genes):,}")
print(f"  Total reads in genes: {sum(c for _, _, c in genes):,}")
print(f"  TPM sum: {total_tpm:.1f} (should be ~1e6)")
covered = sum(1 for _, _, c in genes if c > 0)
print(f"  Genes with counts > 0: {covered:,}")
PYEOF
    fi

    # =========================================================================
    # Per-sample QC stats
    # =========================================================================
    if [[ -f "${GENE_COUNTS}.summary" ]]; then
        ASSIGNED=$(grep "^Assigned" "${GENE_COUNTS}.summary" | awk '{print $2}')
        UNASSIGNED_NOFEATURE=$(grep "^Unassigned_NoFeatures" "${GENE_COUNTS}.summary" | awk '{print $2}')
        UNASSIGNED_AMBIG=$(grep "^Unassigned_Ambiguity" "${GENE_COUNTS}.summary" | awk '{print $2}')
        UNASSIGNED_MULTI=$(grep "^Unassigned_MultiMapping" "${GENE_COUNTS}.summary" | awk '{print $2}')

        TOTAL_FRAG=$((ASSIGNED + UNASSIGNED_NOFEATURE + UNASSIGNED_AMBIG + UNASSIGNED_MULTI))
        if (( TOTAL_FRAG > 0 )); then
            ASSIGN_PCT=$(awk "BEGIN {printf \"%.2f\", ${ASSIGNED}/${TOTAL_FRAG}*100}")
            NOFEAT_PCT=$(awk "BEGIN {printf \"%.2f\", ${UNASSIGNED_NOFEATURE}/${TOTAL_FRAG}*100}")
        else
            ASSIGN_PCT="NA"; NOFEAT_PCT="NA"
        fi

        GENES_NONZERO=$(awk 'NR>2 && $7 > 0' "${GENE_COUNTS}" | wc -l)

        echo ""
        echo "[SANITY CHECK] Quantification for ${SAMPLE}:"
        echo "  Total fragments:         ${TOTAL_FRAG}"
        echo "  Assigned to genes:       ${ASSIGNED} (${ASSIGN_PCT}%)"
        echo "  Unassigned (no feature): ${UNASSIGNED_NOFEATURE} (${NOFEAT_PCT}%)"
        echo "  Unassigned (ambiguous):  ${UNASSIGNED_AMBIG}"
        echo "  Unassigned (multi-map):  ${UNASSIGNED_MULTI}"
        echo "  Genes with counts > 0:   ${GENES_NONZERO}"
        echo "[INTERPRETATION]"
        echo "  Assigned fraction:"
        echo "    >70%:    typical for stranded total RNA"
        echo "    50-70%:  acceptable; rRNA/intronic noise is the usual culprit"
        echo "    <50%:    investigate (wrong strandedness? wrong GTF?)"
        echo "  If the assigned fraction is very low, try -s 0 or -s 1 to see"
        echo "  if the strand setting is wrong."

        cat > "${QC_DIR}/quantify_stats_${SAMPLE}.tsv" << EOF
sample	total_fragments	assigned	assigned_pct	no_feature	no_feature_pct	ambiguous	multi_mapping	genes_nonzero
${SAMPLE}	${TOTAL_FRAG}	${ASSIGNED}	${ASSIGN_PCT}	${UNASSIGNED_NOFEATURE}	${NOFEAT_PCT}	${UNASSIGNED_AMBIG}	${UNASSIGNED_MULTI}	${GENES_NONZERO}
EOF
    fi
done

step_header "Step 04 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output:"
echo "  Gene counts:  ${QUANT_GENE_DIR}/{SAMPLE}_gene_counts.tsv"
echo "  Gene TPM:     ${QUANT_GENE_DIR}/{SAMPLE}_gene_tpm.tsv"
echo "  Tx counts:    ${QUANT_TX_DIR}/{SAMPLE}_tx_counts.tsv"
echo ""
echo "Next: bash 05_make_bigwigs.sh"
