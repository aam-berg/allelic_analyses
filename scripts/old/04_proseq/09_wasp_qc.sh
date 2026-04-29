#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 09_wasp_qc.sh — Compile WASP QC + merge allele count tables across replicates
# =============================================================================
#
# WHAT THIS DOES:
#   1. Consolidates per-sample WASP filtering stats
#   2. Merges per-SNP allele counts across all replicates into one table
#   3. Generates a comprehensive WASP QC report
#
# This script always runs on ALL samples. Do not use with --array.
#
# SBATCH RESOURCES:
#   sbatch --partition=short --time=1:00:00 --mem=8G --cpus-per-task=4 \
#       -o logs/09_%j.out -e logs/09_%j.err 09_wasp_qc.sh
# =============================================================================

source "config_wasp.sh"

step_header "PRO-seq Pipeline Step 09: WASP QC Summary"

source activate "${CONDA_ENV_NAME}"
create_wasp_dirs

# =============================================================================
# 1. Consolidate WASP filtering stats
# =============================================================================
step_header "STEP 1: Consolidate WASP filtering stats"

WASP_COMBINED="${WASP_QC_DIR}/wasp_stats_combined.tsv"

FIRST="${WASP_QC_DIR}/wasp_stats_${SAMPLE_ORDER[0]}.tsv"
if [[ -f "${FIRST}" ]]; then
    head -1 "${FIRST}" > "${WASP_COMBINED}"
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        F="${WASP_QC_DIR}/wasp_stats_${SAMPLE}.tsv"
        if [[ -f "${F}" ]]; then
            tail -1 "${F}" >> "${WASP_COMBINED}"
        else
            echo "[WARNING] Missing: ${F}"
        fi
    done
    echo "[INFO] Combined WASP stats:"
    column -t "${WASP_COMBINED}"
else
    echo "[WARNING] No WASP stats found. Was step 07 completed?"
fi

# =============================================================================
# 2. Merge per-SNP allele counts across replicates
# =============================================================================
step_header "STEP 2: Merge allele counts across replicates"

MERGED_COUNTS="${WASP_QC_DIR}/allele_counts_merged.tsv"

python3 - "${WASP_ALLELE_DIR}" "${MERGED_COUNTS}" ${SAMPLE_ORDER[@]} << 'PYEOF'
import os, sys
from collections import defaultdict

allele_dir = sys.argv[1]
output_file = sys.argv[2]
sample_order = sys.argv[3:]

print(f"  Allele dir: {allele_dir}")
print(f"  Samples: {sample_order}")

merged = defaultdict(lambda: {
    'ref_allele': '', 'alt_allele': '',
    'ref_total': 0, 'alt_total': 0, 'other_total': 0,
    'per_rep': {}
})

for sample in sample_order:
    counts_file = os.path.join(allele_dir, f"{sample}_allele_counts.tsv")
    if not os.path.exists(counts_file):
        print(f"  [WARNING] Missing: {counts_file}")
        continue
    n = 0
    with open(counts_file) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
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
            merged[key]['per_rep'][sample] = (ref_c, alt_c)
            n += 1
    print(f"  {sample}: {n:,} SNPs read")

print(f"\n  Writing merged counts ({len(merged):,} SNPs)...")
with open(output_file, 'w') as f:
    rep_cols = []
    for s in sample_order:
        rep_cols.extend([f"{s}_ref", f"{s}_alt"])
    f.write("chrom\tpos\tref_allele\talt_allele\t"
            "ref_total\talt_total\tother_total\ttotal\tref_fraction\t"
            + "\t".join(rep_cols) + "\n")

    for (chrom, pos) in sorted(merged.keys(), key=lambda x: (x[0], x[1])):
        d = merged[(chrom, pos)]
        total = d['ref_total'] + d['alt_total'] + d['other_total']
        allelic = d['ref_total'] + d['alt_total']
        ref_frac = f"{d['ref_total']/allelic:.4f}" if allelic > 0 else "NA"
        rep_vals = []
        for s in sample_order:
            if s in d['per_rep']:
                rep_vals.extend([str(d['per_rep'][s][0]), str(d['per_rep'][s][1])])
            else:
                rep_vals.extend(["0", "0"])
        f.write(f"{chrom}\t{pos}\t{d['ref_allele']}\t{d['alt_allele']}\t"
                f"{d['ref_total']}\t{d['alt_total']}\t{d['other_total']}\t"
                f"{total}\t{ref_frac}\t" + "\t".join(rep_vals) + "\n")

# Summary
covered = sum(1 for d in merged.values() if d['ref_total'] + d['alt_total'] > 0)
total_ref = sum(d['ref_total'] for d in merged.values())
total_alt = sum(d['alt_total'] for d in merged.values())
total_allelic = total_ref + total_alt
print(f"\n  Summary:")
print(f"    SNPs with allelic coverage: {covered:,}")
print(f"    Total ref reads: {total_ref:,}  |  Total alt reads: {total_alt:,}")
if total_allelic > 0:
    print(f"    Global ref fraction: {100*total_ref/total_allelic:.2f}%")
for threshold in [5, 10, 20]:
    n = sum(1 for d in merged.values() if d['ref_total'] + d['alt_total'] >= threshold)
    print(f"    SNPs with >={threshold} allelic reads: {n:,}")
PYEOF

echo ""
echo "[INFO] Merged allele counts: ${MERGED_COUNTS}"

# =============================================================================
# 3. Final QC report
# =============================================================================
step_header "STEP 3: WASP QC Report"

WASP_REPORT="${WASP_QC_DIR}/wasp_qc_report.txt"

{
    echo "================================================================"
    echo "WASP Allele-Specific Pipeline — QC Report"
    echo "================================================================"
    echo "Generated: $(date)"
    echo ""

    echo "--- WASP FILTERING STATS ---"
    if [[ -f "${WASP_COMBINED}" ]]; then
        column -t "${WASP_COMBINED}"
    fi
    echo ""

    echo "--- ALLELE-SPECIFIC BAMs ---"
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        echo "  ${SAMPLE}:"
        for ALLELE in ref alt nosnp ambiguous; do
            BAM="${WASP_ALLELE_DIR}/${SAMPLE}_${ALLELE}.bam"
            if [[ -f "${BAM}" ]]; then
                COUNT=$(samtools view -c "${BAM}" 2>/dev/null || echo "?")
                echo "    ${ALLELE}: ${COUNT} reads"
            fi
        done
    done
    echo ""

    echo "--- ALLELE-SPECIFIC bigWigs ---"
    ls -lh "${WASP_ALLELE_BW_DIR}"/*.bw 2>/dev/null | awk '{print "  "$NF" ("$5")"}' \
        || echo "  (none found)"
    echo ""

    echo "--- DOWNSTREAM ANALYSIS GUIDE ---"
    echo ""
    echo "For your TFBS motif-strength analysis, the key workflow is:"
    echo ""
    echo "1. Load annotated motif file (e.g. AC0007_annotated.tsv.gz)"
    echo ""
    echo "2. Filter to TFBSs of interest:"
    echo "   is_intragenic=TRUE, is_promoter=FALSE,"
    echo "   atac_overlap=TRUE, snp_F121-9_overlap=TRUE"
    echo ""
    echo "3. For each TFBS with het SNP(s), look up allele counts:"
    echo "   Merged counts: ${MERGED_COUNTS}"
    echo "   Per-replicate: ${WASP_ALLELE_DIR}/*_allele_counts.tsv"
    echo ""
    echo "4. Determine which allele gives stronger motif:"
    echo "   Use the snp_F121-9_details column (chr:pos:REF>ALT) to identify"
    echo "   the SNP, then score the motif with each allele. The 'score' column"
    echo "   in the annotated file gives the reference-allele MOODS score."
    echo ""
    echo "5. Compare PRO-seq signal at stronger vs weaker motif alleles,"
    echo "   pooling across TFBSs within each motif archetype."
    echo ""

} > "${WASP_REPORT}"

cat "${WASP_REPORT}"

step_header "WASP PIPELINE COMPLETE"
echo ""
echo "FINAL OUTPUT SUMMARY"
echo ""
echo "WASP-filtered BAMs:            ${WASP_FILTERED_DIR}/*_wasp.bam"
echo "Allele-specific BAMs:          ${WASP_ALLELE_DIR}/*_{ref,alt}.bam"
echo "Allele-specific bigWigs:       ${WASP_ALLELE_BW_DIR}/"
echo "Per-SNP allele counts (merged): ${MERGED_COUNTS}"
echo "Per-SNP allele counts (per-rep): ${WASP_ALLELE_DIR}/*_allele_counts.tsv"
echo "QC report:                     ${WASP_REPORT}"