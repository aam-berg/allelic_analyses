#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04_make_bigwigs.sh — Generate strand-specific bigWig files (3' and 5' ends)
# =============================================================================
#
# WHAT THIS DOES:
#   For each replicate, generates 4 bigWig files:
#     {sample}_3prime_plus.bw   : 3' end of nascent RNA, plus strand
#     {sample}_3prime_minus.bw  : 3' end of nascent RNA, minus strand
#     {sample}_5prime_plus.bw   : 5' end of nascent RNA, plus strand
#     {sample}_5prime_minus.bw  : 5' end of nascent RNA, minus strand
#
#   These are single-nucleotide resolution, non-normalized COUNT files.
#
# =============================================================================
# *** CRITICAL: UNDERSTANDING PRO-seq STRAND LOGIC ***
# =============================================================================
#
# The PRO-seq library prep:
#   1. Nuclear run-on: Pol II extends nascent RNA with biotin-NTPs
#   2. RNA extraction + base hydrolysis (fragmentation)
#   3. Biotin enrichment (selects 3'-biotin fragments)
#   4. 3' adapter ligation to RNA 3' end
#   5. Reverse transcription from the 3' adapter → cDNA
#   6. PCR + sequencing of cDNA
#
# Because of reverse transcription:
#   - Sequencing reads = cDNA = REVERSE COMPLEMENT of the RNA
#   - Reads map to the OPPOSITE strand from the original RNA
#
# THE STRAND SWAP:
#   Read on + strand in BAM → RNA was on − strand → MINUS strand signal
#   Read on − strand in BAM → RNA was on + strand → PLUS strand signal
#
# END ASSIGNMENT:
#   5' end of read = 3' end of RNA (Pol II active site — RT starts here)
#   3' end of read = 5' end of RNA (5' end of fragment)
#
# So in bedtools terms:
#   -5 flag = 5' of read = 3' of RNA (Pol II position)
#   -3 flag = 3' of read = 5' of RNA
#
# PUTTING IT TOGETHER for 3' end of RNA (Pol II active site):
#   Reads with FLAG 16 (BAM − strand), bedtools -5  → plus_3prime.bw
#   Reads without FLAG 16 (BAM + strand), bedtools -5 → minus_3prime.bw
#
# And for 5' end of RNA:
#   Reads with FLAG 16, bedtools -3  → plus_5prime.bw
#   Reads without FLAG 16, bedtools -3 → minus_5prime.bw
#
# WHY BOTH 3' AND 5' ENDS?
#   3' end: Most commonly used. Represents active site of RNA Pol II.
#   5' end: Less common but useful for some analyses (TSS mapping, etc.)
#
# PARALLEL EXECUTION:
#   sbatch --array=0-3 --partition=short --time=4:00:00 --mem=16G --cpus-per-task=8 \
#       -o logs/04_%A_%a.out -e logs/04_%A_%a.err 04_make_bigwigs.sh
# =============================================================================

source "config.sh"

step_header "PRO-seq Pipeline Step 04: Generate bigWig Files"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] bedtools version: $(bedtools --version)"
echo "[INFO] bedGraphToBigWig: $(which bedGraphToBigWig)"

create_dirs
resolve_samples "${1:-}"

# Verify chrom sizes file
if [[ ! -f "${MM39_CHROM_SIZES}" ]]; then
    echo "[ERROR] Chromosome sizes file not found: ${MM39_CHROM_SIZES}"
    exit 1
fi

# =============================================================================
# Function: make_bigwig
# =============================================================================
make_bigwig() {
    local BAM="$1"
    local SAM_FLAG="$2"      # -f 16 or -F 16
    local BT_END_FLAG="$3"   # -5 or -3
    local OUTBW="$4"
    local LABEL="$5"

    local TMPBG="${OUTBW%.bw}.bedGraph"

    if [[ -f "${OUTBW}" ]]; then
        echo "  [INFO] ${LABEL}: bigWig already exists. Skipping."
        return 0
    fi

    echo "  [INFO] ${LABEL}: extracting positions..."

    # Pipeline:
    # 1. samtools view: filter reads by strand
    # 2. bedtools genomecov -5/-3 -bg: single-nucleotide coverage
    # 3. Sort bedGraph (required by bedGraphToBigWig)
    # 4. Filter to valid chromosomes in chrom.sizes
    # 5. Convert to bigWig
    samtools view -b ${SAM_FLAG} "${BAM}" \
    | bedtools genomecov -ibam stdin -bg ${BT_END_FLAG} \
    | sort -k1,1 -k2,2n \
    | awk -v OFS='\t' 'NR==FNR {chr[$1]=1; next} ($1 in chr)' "${MM39_CHROM_SIZES}" - \
    > "${TMPBG}"

    local N_POS=$(wc -l < "${TMPBG}")
    local TOTAL_COUNTS=$(awk '{s+=$4} END {print s+0}' "${TMPBG}")
    echo "  [INFO] ${LABEL}: ${N_POS} positions, ${TOTAL_COUNTS} total counts"

    bedGraphToBigWig "${TMPBG}" "${MM39_CHROM_SIZES}" "${OUTBW}"
    rm -f "${TMPBG}"
    echo "  [INFO] ${LABEL}: → ${OUTBW}"
}

# =============================================================================
# Generate bigWigs for each sample
# =============================================================================
for SAMPLE in "${RUN_SAMPLES[@]}"; do
    step_header "Generating bigWigs: ${SAMPLE}"

    FINAL_BAM="${BAM_MM39_DIR}/${SAMPLE}_final.bam"
    if [[ ! -f "${FINAL_BAM}" ]]; then
        echo "[ERROR] Final BAM not found: ${FINAL_BAM}"
        exit 1
    fi

    OUTDIR="${BIGWIG_IND_DIR}"

    # --- 3' end of nascent RNA (Pol II active site) ---
    echo "--- 3' end of nascent RNA (Pol II position) ---"

    # Plus strand: reads on BAM − strand (flag 16), use -5 (5' of read = 3' of RNA)
    make_bigwig "${FINAL_BAM}" "-f 16" "-5" \
        "${OUTDIR}/${SAMPLE}_3prime_plus.bw" "3prime_plus"

    # Minus strand: reads on BAM + strand (NOT flag 16), use -5
    make_bigwig "${FINAL_BAM}" "-F 16" "-5" \
        "${OUTDIR}/${SAMPLE}_3prime_minus.bw" "3prime_minus"

    # --- 5' end of nascent RNA ---
    echo ""
    echo "--- 5' end of nascent RNA ---"

    # Plus strand: reads on BAM − strand (flag 16), use -3 (3' of read = 5' of RNA)
    make_bigwig "${FINAL_BAM}" "-f 16" "-3" \
        "${OUTDIR}/${SAMPLE}_5prime_plus.bw" "5prime_plus"

    # Minus strand: reads on BAM + strand (NOT flag 16), use -3
    make_bigwig "${FINAL_BAM}" "-F 16" "-3" \
        "${OUTDIR}/${SAMPLE}_5prime_minus.bw" "5prime_minus"

    echo ""
    echo "[DONE] ${SAMPLE}: 4 bigWig files"
    ls -lh "${OUTDIR}/${SAMPLE}"_*.bw | sed 's/^/    /'
done

# =============================================================================
# Sanity check: strand balance
# =============================================================================
step_header "SANITY CHECK: Strand balance"

echo "File sizes (plus and minus should be roughly similar):"
for SAMPLE in "${RUN_SAMPLES[@]}"; do
    echo "--- ${SAMPLE} ---"
    for END in 3prime 5prime; do
        PLUS_BW="${BIGWIG_IND_DIR}/${SAMPLE}_${END}_plus.bw"
        MINUS_BW="${BIGWIG_IND_DIR}/${SAMPLE}_${END}_minus.bw"
        if [[ -f "${PLUS_BW}" && -f "${MINUS_BW}" ]]; then
            echo "  ${END}: plus=$(du -h "${PLUS_BW}" | cut -f1), minus=$(du -h "${MINUS_BW}" | cut -f1)"
        fi
    done
done

step_header "Step 04 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output: ${BIGWIG_IND_DIR}/"
echo ""
echo "Next: Run 05_merge_and_qc.sh (after ALL step 04 array tasks finish)"