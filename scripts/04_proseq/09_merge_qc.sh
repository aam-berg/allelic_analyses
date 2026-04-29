#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04_proseq/09_merge_qc.sh — Cross-rep merging + final QC compilation
# =============================================================================
#
# WHAT THIS DOES (single non-array job, runs once at end):
#   1. Merges per-rep regular bigWigs (4 per sample) across replicates by
#      summing counts. Output: merged_{3,5}prime_{plus,minus}.bw (4 files).
#   2. Merges per-rep ref/alt allele bigWigs across replicates. Output:
#      merged_{ref,alt}_3prime_{plus,minus}.bw (4 files).
#   3. Merges per-SNP allele count tables across replicates into one master
#      table at qc/allele_counts_merged.tsv.
#   4. Compiles a single comprehensive pipeline_qc_report.txt that pulls
#      together: read counts, trimming, alignment, WASP, allele-split, and
#      coverage-summary stats from all per-sample TSVs.
#
# WHY ONE COMBINED MERGE STEP:
#   Previously, there were two separate merge/QC scripts (05_merge_and_qc.sh
#   and 09_wasp_qc.sh). They depended on different earlier stages and had
#   overlapping responsibilities. Consolidating into one final step makes
#   the dependency graph cleaner: this step waits for both step 04 (regular
#   bigWigs) AND step 08 (allele bigWigs), then does all merging at once.
#
# SBATCH:
#   sbatch --partition=short --time=2:00:00 --mem=16G --cpus-per-task=4 \
#       -o logs/09_%j.out -e logs/09_%j.err 09_merge_qc.sh
# =============================================================================

source "config.sh"

step_header "PRO-seq Pipeline Step 09: Merge + Final QC"

source activate "${CONDA_ENV_NAME}"
create_dirs

# -----------------------------------------------------------------------------
# Helper: merge a list of bigWigs into a single bigWig (sum of counts)
# -----------------------------------------------------------------------------
merge_bigwigs() {
    local OUTBW="$1"; shift
    local LABEL="$1"; shift
    # remaining args: input bigwig files

    if [[ -f "${OUTBW}" ]]; then
        echo "  [SKIP] ${LABEL}: merged bigWig exists."
        return 0
    fi

    local TMPBG="${OUTBW%.bw}.bedGraph"
    bigWigMerge "$@" "${TMPBG}"
    sort -k1,1 -k2,2n -o "${TMPBG}" "${TMPBG}"
    awk 'NR==FNR {chr[$1]=1; next} ($1 in chr)' "${MM39_CHROM_SIZES}" "${TMPBG}" \
        > "${TMPBG}.f"
    mv "${TMPBG}.f" "${TMPBG}"

    if [[ ! -s "${TMPBG}" ]]; then
        echo "  [WARN] ${LABEL}: empty merge result."
        rm -f "${TMPBG}"
        return 0
    fi

    bedGraphToBigWig "${TMPBG}" "${MM39_CHROM_SIZES}" "${OUTBW}"
    rm -f "${TMPBG}"
    echo "  [OK]   ${LABEL}: -> ${OUTBW}"
}

# =============================================================================
# 1. Verify all per-rep bigWigs exist before merging
# =============================================================================
step_header "Verifying per-rep bigWigs"

ALL_OK=1
for SAMPLE in "${SAMPLE_ORDER[@]}"; do
    for END in 3prime 5prime; do
        for STRAND in plus minus; do
            BW="${BIGWIG_IND_DIR}/${SAMPLE}_${END}_${STRAND}.bw"
            if [[ -f "${BW}" ]]; then
                echo "  [OK]    ${BW}"
            else
                echo "  [MISSING] ${BW}"
                ALL_OK=0
            fi
        done
    done
done
(( ALL_OK )) || \
    { echo "[ERROR] Some per-rep bigWigs missing; run step 04 to completion." >&2; exit 1; }

echo ""
echo "Verifying per-allele bigWigs..."
ALL_OK=1
for SAMPLE in "${SAMPLE_ORDER[@]}"; do
    for ALLELE in ref alt; do
        for STRAND in plus minus; do
            BW="${BIGWIG_ALLELE_DIR}/${SAMPLE}_${ALLELE}_3prime_${STRAND}.bw"
            if [[ -f "${BW}" ]]; then
                echo "  [OK]    ${BW}"
            else
                echo "  [MISSING] ${BW}"
                ALL_OK=0
            fi
        done
    done
done
(( ALL_OK )) || \
    { echo "[ERROR] Some allele bigWigs missing; run step 08 to completion." >&2; exit 1; }

# =============================================================================
# 2. Merge per-rep regular bigWigs (3'/5' × +/-)
# =============================================================================
step_header "Merging per-rep regular bigWigs"

for END in 3prime 5prime; do
    for STRAND in plus minus; do
        BWS=()
        for SAMPLE in "${SAMPLE_ORDER[@]}"; do
            BWS+=("${BIGWIG_IND_DIR}/${SAMPLE}_${END}_${STRAND}.bw")
        done
        OUTBW="${BIGWIG_MERGED_DIR}/merged_${END}_${STRAND}.bw"
        merge_bigwigs "${OUTBW}" "merged_${END}_${STRAND}" "${BWS[@]}"
    done
done

# =============================================================================
# 3. Merge per-rep allele bigWigs (ref/alt × +/-)
# =============================================================================
step_header "Merging per-rep allele bigWigs"

for ALLELE in ref alt; do
    for STRAND in plus minus; do
        BWS=()
        for SAMPLE in "${SAMPLE_ORDER[@]}"; do
            BWS+=("${BIGWIG_ALLELE_DIR}/${SAMPLE}_${ALLELE}_3prime_${STRAND}.bw")
        done
        OUTBW="${BIGWIG_ALLELE_MERGED_DIR}/merged_${ALLELE}_3prime_${STRAND}.bw"
        merge_bigwigs "${OUTBW}" "merged_${ALLELE}_3prime_${STRAND}" "${BWS[@]}"
    done
done

# =============================================================================
# 4. Merge per-SNP allele counts across replicates
# =============================================================================
step_header "Merging per-SNP allele counts across replicates"

MERGED_COUNTS="${QC_DIR}/allele_counts_merged.tsv"

# Pass paths and sample names to a Python heredoc
python3 - "${ALLELE_DIR}" "${MERGED_COUNTS}" "${SAMPLE_ORDER[@]}" << 'PYEOF'
import os, sys
from collections import defaultdict

allele_dir = sys.argv[1]
output_file = sys.argv[2]
samples = sys.argv[3:]

print(f"  Allele dir: {allele_dir}")
print(f"  Samples:    {samples}")

# merged[(chrom, pos)] aggregates per-SNP totals plus per-rep counts
merged = defaultdict(lambda: {
    'ref_allele': '', 'alt_allele': '',
    'ref_total': 0, 'alt_total': 0, 'other_total': 0,
    'per_rep': {}
})

for s in samples:
    f = os.path.join(allele_dir, f"{s}_allele_counts.tsv")
    if not os.path.exists(f):
        print(f"  [WARN] missing: {f}")
        continue
    n = 0
    with open(f) as fh:
        fh.readline()  # header
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 7:
                continue
            chrom, pos = parts[0], int(parts[1])
            ref_a, alt_a = parts[2], parts[3]
            ref_c, alt_c, other_c = int(parts[4]), int(parts[5]), int(parts[6])
            key = (chrom, pos)
            merged[key]['ref_allele'] = ref_a
            merged[key]['alt_allele'] = alt_a
            merged[key]['ref_total'] += ref_c
            merged[key]['alt_total'] += alt_c
            merged[key]['other_total'] += other_c
            merged[key]['per_rep'][s] = (ref_c, alt_c)
            n += 1
    print(f"  {s}: {n:,} SNPs")

print(f"\n  Writing merged counts ({len(merged):,} unique SNPs)...")
with open(output_file, 'w') as out:
    rep_cols = []
    for s in samples:
        rep_cols.extend([f"{s}_ref", f"{s}_alt"])
    out.write("chrom\tpos\tref_allele\talt_allele\t"
              "ref_total\talt_total\tother_total\ttotal\tref_fraction\t"
              + "\t".join(rep_cols) + "\n")

    for (chrom, pos) in sorted(merged.keys(), key=lambda x: (x[0], x[1])):
        d = merged[(chrom, pos)]
        total = d['ref_total'] + d['alt_total'] + d['other_total']
        allelic = d['ref_total'] + d['alt_total']
        ref_frac = f"{d['ref_total']/allelic:.4f}" if allelic > 0 else "NA"
        rep_vals = []
        for s in samples:
            if s in d['per_rep']:
                rep_vals.extend([str(d['per_rep'][s][0]), str(d['per_rep'][s][1])])
            else:
                rep_vals.extend(["0", "0"])
        out.write(f"{chrom}\t{pos}\t{d['ref_allele']}\t{d['alt_allele']}\t"
                  f"{d['ref_total']}\t{d['alt_total']}\t{d['other_total']}\t"
                  f"{total}\t{ref_frac}\t" + "\t".join(rep_vals) + "\n")

# Summary
covered = sum(1 for d in merged.values() if d['ref_total'] + d['alt_total'] > 0)
total_ref = sum(d['ref_total'] for d in merged.values())
total_alt = sum(d['alt_total'] for d in merged.values())
allelic = total_ref + total_alt
print(f"\n  Summary:")
print(f"    SNPs with allelic coverage: {covered:,}")
print(f"    Total ref reads at SNPs:    {total_ref:,}")
print(f"    Total alt reads at SNPs:    {total_alt:,}")
if allelic > 0:
    pct = 100 * total_ref / allelic
    print(f"    Global ref fraction:        {pct:.2f}%  (expect ~50% after WASP)")
for thr in (5, 10, 20):
    n = sum(1 for d in merged.values() if d['ref_total'] + d['alt_total'] >= thr)
    print(f"    SNPs with >={thr} allelic reads: {n:,}")
PYEOF

echo ""
echo "[INFO] Merged allele counts: ${MERGED_COUNTS}"

# =============================================================================
# 5. Compile final pipeline_qc_report.txt
# =============================================================================
step_header "Compiling final QC report"

QC_REPORT="${QC_DIR}/pipeline_qc_report.txt"
{
    echo "============================================================"
    echo "F121-9 PRO-seq Pipeline — QC Report"
    echo "Generated: $(date)"
    echo "Pipeline:  ${SCRIPT_DIR}"
    echo "Output:    ${SCRATCH_DIR}"
    echo "Samples:   ${SAMPLE_ORDER[*]}"
    echo "============================================================"
    echo ""
    echo "Pipeline notes:"
    echo "  - No deduplication (PRO-seq pause sites can have many reads at"
    echo "    the same nucleotide; dedup would collapse real biology)"
    echo "  - rRNA filter: bowtie2 default presets (sensitive)"
    echo "  - Spike-in filter: ${DM6_SPIKEIN_FILTER}"
    echo "  - MAPQ threshold: ${MAPQ_THRESHOLD}"
    echo "  - chrM: removed in step 03"
    echo ""

    echo "--- 1. Raw read counts (step 01) ---"
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        F="${QC_DIR}/read_counts_raw_${SAMPLE}.tsv"
        [[ -f "${F}" ]] && cat "${F}" || echo -e "${SAMPLE}\t?\t?"
    done
    echo ""

    echo "--- 2. Trimming (step 02) ---"
    echo -e "sample\traw_reads\tadapted\ttoo_short\twritten\tpct_adapter\tpct_surviving"
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        F="${QC_DIR}/trimming_stats_${SAMPLE}.tsv"
        [[ -f "${F}" ]] && cat "${F}" || echo -e "${SAMPLE}\t?\t?\t?\t?\t?\t?"
    done
    echo ""

    echo "--- 3. Alignment cascade (step 03) ---"
    H=0
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        F="${QC_DIR}/alignment_stats_${SAMPLE}.tsv"
        if [[ -f "${F}" ]]; then
            (( H == 0 )) && { head -1 "${F}"; H=1; }
            tail -n +2 "${F}"
        fi
    done
    echo ""

    echo "--- 4. WASP mapping bias correction (step 06) ---"
    H=0
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        F="${QC_DIR}/wasp_stats_${SAMPLE}.tsv"
        if [[ -f "${F}" ]]; then
            (( H == 0 )) && { head -1 "${F}"; H=1; }
            tail -n +2 "${F}"
        fi
    done
    echo ""

    echo "--- 5. Allele-specific splitting (step 07) ---"
    H=0
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        F="${QC_DIR}/allele_split_${SAMPLE}.tsv"
        if [[ -f "${F}" ]]; then
            (( H == 0 )) && { head -1 "${F}"; H=1; }
            tail -n +2 "${F}"
        fi
    done
    echo ""

    echo "--- 6. Output bigWigs ---"
    echo "Per-rep (${BIGWIG_IND_DIR}/):"
    ls -lh "${BIGWIG_IND_DIR}"/*.bw 2>/dev/null | awk '{print "  "$NF" ("$5")"}' || echo "  (none)"
    echo ""
    echo "Merged regular (${BIGWIG_MERGED_DIR}/):"
    ls -lh "${BIGWIG_MERGED_DIR}"/*.bw 2>/dev/null | awk '{print "  "$NF" ("$5")"}' || echo "  (none)"
    echo ""
    echo "Per-allele (${BIGWIG_ALLELE_DIR}/):"
    ls -lh "${BIGWIG_ALLELE_DIR}"/*.bw 2>/dev/null | awk '{print "  "$NF" ("$5")"}' || echo "  (none)"
    echo ""
    echo "Merged allele (${BIGWIG_ALLELE_MERGED_DIR}/):"
    ls -lh "${BIGWIG_ALLELE_MERGED_DIR}"/*.bw 2>/dev/null | awk '{print "  "$NF" ("$5")"}' || echo "  (none)"
    echo ""

    echo "--- 7. Per-SNP allele counts (merged) ---"
    echo "  ${MERGED_COUNTS}"
    if [[ -f "${MERGED_COUNTS}" ]]; then
        N=$(($(wc -l < "${MERGED_COUNTS}") - 1))
        echo "  Total SNPs:  ${N}"
    fi
    echo ""

    echo "============================================================"
    echo "End of report. Source TSVs in: ${QC_DIR}/"
    echo "Final BAMs in:                 ${BAM_MM39_DIR}/*_final.bam"
    echo "Allele-split BAMs in:          ${ALLELE_DIR}/"
    echo "============================================================"
} > "${QC_REPORT}"

echo "[INFO] QC report: ${QC_REPORT}"
echo ""
echo "----- Preview (first 80 lines) -----"
head -80 "${QC_REPORT}"
echo "----- ... -----"

step_header "Step 09 COMPLETE — Pipeline finished"
echo ""
echo "Final outputs:"
echo "  Final BAMs:            ${BAM_MM39_DIR}/*_final.bam"
echo "  WASP-filtered BAMs:    ${WASP_FILTERED_DIR}/"
echo "  Allele-split BAMs:     ${ALLELE_DIR}/{SAMPLE}_{ref,alt,nosnp,ambiguous}.bam"
echo "  Per-rep bigWigs:       ${BIGWIG_IND_DIR}/"
echo "  Merged regular:        ${BIGWIG_MERGED_DIR}/merged_{3,5}prime_{plus,minus}.bw"
echo "  Per-allele bigWigs:    ${BIGWIG_ALLELE_DIR}/"
echo "  Merged allele:         ${BIGWIG_ALLELE_MERGED_DIR}/merged_{ref,alt}_3prime_{plus,minus}.bw"
echo "  Per-SNP merged counts: ${MERGED_COUNTS}"
echo "  Final QC report:       ${QC_REPORT}"
echo ""
echo "Downstream consumer: 06_allele_pairs/04_proseq_alleles_at_motifs.R"
