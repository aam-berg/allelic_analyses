#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04c_rnaseq/08_merge_qc.sh — Cross-rep merging + QC + expression summary
# =============================================================================
#
# WHAT THIS DOES (single non-array job, runs once at end):
#   1. Merges per-rep strand-specific bigWigs across replicates into
#      merged_{plus,minus}.bw (sum of counts).
#   2. Merges per-rep allele bigWigs across replicates into
#      merged_{ref,alt}_{plus,minus}.bw.
#   3. Merges per-SNP allele count tables across replicates into
#      qc/allele_counts_merged.tsv.
#   4. **Builds gene_expression_summary.tsv** — the canonical expression
#      table consumed by 05_motif_annot/02f_gene_expression.R. One row per
#      gene_id, with per-rep counts/TPM and across-rep summary statistics.
#   5. Compiles a comprehensive pipeline_qc_report.txt.
#
# THE EXPRESSION TABLE (the most important output of this step):
#   For each gene_id present in any per-sample featureCounts output:
#     gene_id          GENCODE gene ID (e.g., ENSMUSG00000051951.7)
#     gene_name        Gene symbol (joined from GTF)
#     gene_strand      Strand from GTF ('+' or '-')
#     length           Sum of exonic length used by featureCounts
#     {SAMPLE}_counts  Per-rep raw count column (one per sample)
#     {SAMPLE}_tpm     Per-rep TPM column (one per sample)
#     counts_sum       Sum of counts across replicates
#     tpm_mean         Mean of per-rep TPM
#     tpm_sd           SD of per-rep TPM
#     tpm_log2_mean    log2(tpm_mean + 1)  for downstream "expressed" filtering
#
#   Strand-aware joining in 05_motif_annot/ uses gene_strand: a (motif, +)
#   gets joined with genes on '+', and so on for each motif orientation.
#
# WHY MERGE PER-SNP ALLELE COUNTS HERE:
#   Same as in 03_atac/ and 04_proseq/. Cross-rep aggregation makes a single
#   table that 06_allele_pairs/ can reference without re-reading per-rep
#   files.
#
# SBATCH:
#   sbatch --partition=short --time=2:00:00 --mem=16G --cpus-per-task=4 \
#       -o logs/08_%j.out -e logs/08_%j.err 08_merge_qc.sh
# =============================================================================

source "config.sh"

step_header "RNA-seq Pipeline Step 08: Merge + QC + Expression Summary"

source activate "${CONDA_ENV_NAME}"
create_dirs

# -----------------------------------------------------------------------------
# Helper: merge a list of bigWigs into one (sum of counts)
# -----------------------------------------------------------------------------
merge_bigwigs() {
    local OUTBW="$1"; shift
    local LABEL="$1"; shift

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
# 1. Verify per-rep bigWigs exist
# =============================================================================
step_header "Verifying per-rep bigWigs"
ALL_OK=1
for SAMPLE in "${SAMPLE_ORDER[@]}"; do
    for STRAND in plus minus; do
        BW="${BIGWIG_IND_DIR}/${SAMPLE}_${STRAND}.bw"
        if [[ -f "${BW}" ]]; then
            echo "  [OK]    ${BW}"
        else
            echo "  [MISSING] ${BW}"
            ALL_OK=0
        fi
    done
done
(( ALL_OK )) || \
    { echo "[ERROR] Some per-rep bigWigs missing; run step 05 to completion." >&2; exit 1; }

echo ""
echo "Verifying per-allele bigWigs..."
ALL_OK=1
for SAMPLE in "${SAMPLE_ORDER[@]}"; do
    for ALLELE in ref alt; do
        for STRAND in plus minus; do
            BW="${BIGWIG_ALLELE_DIR}/${SAMPLE}_${ALLELE}_${STRAND}.bw"
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
    { echo "[ERROR] Some allele bigWigs missing; run step 07 to completion." >&2; exit 1; }

# =============================================================================
# 2. Merge per-rep regular bigWigs (plus, minus)
# =============================================================================
step_header "Merging per-rep regular bigWigs"
for STRAND in plus minus; do
    BWS=()
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        BWS+=("${BIGWIG_IND_DIR}/${SAMPLE}_${STRAND}.bw")
    done
    OUTBW="${BIGWIG_MERGED_DIR}/merged_${STRAND}.bw"
    merge_bigwigs "${OUTBW}" "merged_${STRAND}" "${BWS[@]}"
done

# =============================================================================
# 3. Merge per-rep allele bigWigs (ref/alt × plus/minus)
# =============================================================================
step_header "Merging per-rep allele bigWigs"
for ALLELE in ref alt; do
    for STRAND in plus minus; do
        BWS=()
        for SAMPLE in "${SAMPLE_ORDER[@]}"; do
            BWS+=("${BIGWIG_ALLELE_DIR}/${SAMPLE}_${ALLELE}_${STRAND}.bw")
        done
        OUTBW="${BIGWIG_ALLELE_MERGED_DIR}/merged_${ALLELE}_${STRAND}.bw"
        merge_bigwigs "${OUTBW}" "merged_${ALLELE}_${STRAND}" "${BWS[@]}"
    done
done

# =============================================================================
# 4. Build gene_expression_summary.tsv (the table 05_motif_annot will read)
# =============================================================================
step_header "Building gene expression summary"

EXPR_OUT="${EXPRESSION_DIR}/gene_expression_summary.tsv"

# We need an uncompressed GTF to extract gene_name and gene_strand.
GTF_USE="${GTF_FILE}"
if [[ "${GTF_FILE}" == *.gz ]]; then
    GTF_USE="${GTF_FILE%.gz}"
    [[ -f "${GTF_USE}" ]] || gunzip -k "${GTF_FILE}"
fi

# Pass per-sample count files + GTF to a Python heredoc; the script reads
# featureCounts gene-level outputs and assembles a single gene-by-sample
# table with summary statistics.
python3 - "${QUANT_GENE_DIR}" "${GTF_USE}" "${EXPR_OUT}" "${SAMPLE_ORDER[@]}" << 'PYEOF'
import sys
import os
import math
from collections import defaultdict

quant_dir   = sys.argv[1]
gtf_path    = sys.argv[2]
out_path    = sys.argv[3]
samples     = sys.argv[4:]

print(f"  Quant dir: {quant_dir}")
print(f"  GTF:       {gtf_path}")
print(f"  Samples:   {samples}")
print(f"  Output:    {out_path}")

# ---- 1. Parse GTF for gene_id -> (gene_name, gene_strand) ----
# featureCounts uses gene_id from the GTF, so we match on that.
print("\n  Parsing GTF for gene_name and gene_strand...")
gene_info = {}
with open(gtf_path) as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9 or parts[2] != 'gene':
            continue
        strand = parts[6]
        attrs = parts[8]
        # Extract gene_id and gene_name
        gene_id = None
        gene_name = None
        for kv in attrs.split(';'):
            kv = kv.strip()
            if kv.startswith('gene_id "'):
                gene_id = kv[len('gene_id "'):-1]
            elif kv.startswith('gene_name "'):
                gene_name = kv[len('gene_name "'):-1]
        if gene_id is not None:
            gene_info[gene_id] = (gene_name or '', strand)
print(f"    Parsed {len(gene_info):,} gene entries from GTF")

# ---- 2. Read per-sample featureCounts gene counts ----
# featureCounts output (gene-level):
#   line starts with '#' = comment
#   header line: Geneid  Chr  Start  End  Strand  Length  <bam_path>
#   then one row per gene
# We index by gene_id.

# Per-sample dict: gene_id -> count
sample_counts = {s: {} for s in samples}
gene_lengths = {}   # gene_id -> length (use first sample's value)

for s in samples:
    counts_path = os.path.join(quant_dir, f"{s}_gene_counts.tsv")
    if not os.path.exists(counts_path):
        print(f"  [WARN] missing: {counts_path}")
        continue
    n = 0
    with open(counts_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('Geneid'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 7:
                continue
            gene_id = parts[0]
            length  = int(parts[5])
            count   = int(parts[6])
            sample_counts[s][gene_id] = count
            if gene_id not in gene_lengths:
                gene_lengths[gene_id] = length
            n += 1
    print(f"  {s}: {n:,} genes, total counts {sum(sample_counts[s].values()):,}")

# ---- 3. Compute per-sample TPM ----
# TPM_i = (count_i / length_i_kb) / sum(count_j / length_j_kb) * 1e6
sample_tpm = {s: {} for s in samples}
for s in samples:
    rates = {}
    for g, c in sample_counts[s].items():
        L = gene_lengths.get(g, 0)
        rates[g] = (c / (L / 1000.0)) if L > 0 else 0.0
    total_rate = sum(rates.values())
    if total_rate > 0:
        scale = 1e6 / total_rate
        for g, r in rates.items():
            sample_tpm[s][g] = r * scale
    else:
        for g in rates:
            sample_tpm[s][g] = 0.0

# ---- 4. Assemble output table ----
# Universe of genes = union across samples
all_genes = set()
for s in samples:
    all_genes.update(sample_counts[s].keys())
all_genes = sorted(all_genes)
print(f"\n  Total unique genes across samples: {len(all_genes):,}")

def mean(xs): return sum(xs)/len(xs) if xs else 0.0
def sd(xs):
    if len(xs) < 2: return 0.0
    m = mean(xs)
    return math.sqrt(sum((x-m)**2 for x in xs) / (len(xs)-1))

with open(out_path, 'w') as out:
    # Header
    cols = ['gene_id', 'gene_name', 'gene_strand', 'length']
    for s in samples:
        cols.append(f"{s}_counts")
    for s in samples:
        cols.append(f"{s}_tpm")
    cols += ['counts_sum', 'tpm_mean', 'tpm_sd', 'tpm_log2_mean']
    out.write('\t'.join(cols) + '\n')

    n_with_info  = 0
    n_expressed  = 0   # tpm_mean >= 1
    for g in all_genes:
        gname, gstrand = gene_info.get(g, ('', '.'))
        if gname or gstrand != '.':
            n_with_info += 1
        L = gene_lengths.get(g, 0)
        per_rep_counts = [sample_counts[s].get(g, 0)   for s in samples]
        per_rep_tpm    = [sample_tpm[s].get(g, 0.0)    for s in samples]
        c_sum   = sum(per_rep_counts)
        t_mean  = mean(per_rep_tpm)
        t_sd    = sd(per_rep_tpm)
        t_log2  = math.log2(t_mean + 1)
        if t_mean >= 1:
            n_expressed += 1

        row = [g, gname, gstrand, str(L)]
        row += [str(c) for c in per_rep_counts]
        row += [f"{t:.4f}" for t in per_rep_tpm]
        row += [str(c_sum), f"{t_mean:.4f}", f"{t_sd:.4f}", f"{t_log2:.4f}"]
        out.write('\t'.join(row) + '\n')

# Sanity print
print(f"\n  Wrote: {out_path}")
print(f"  Genes total:                 {len(all_genes):,}")
print(f"  Genes with GTF metadata:     {n_with_info:,}")
print(f"  Genes 'expressed' (TPM>=1):  {n_expressed:,}")
print(f"  Per-sample TPM should sum to ~1e6:")
for s in samples:
    print(f"    {s}:  {sum(sample_tpm[s].values()):,.0f}")
PYEOF

echo ""
echo "[INFO] gene_expression_summary.tsv ready at: ${EXPR_OUT}"

# =============================================================================
# 5. Merge per-SNP allele counts across replicates
# =============================================================================
step_header "Merging per-SNP allele counts across replicates"

MERGED_COUNTS="${QC_DIR}/allele_counts_merged.tsv"

python3 - "${ALLELE_DIR}" "${MERGED_COUNTS}" "${SAMPLE_ORDER[@]}" << 'PYEOF'
import os, sys
from collections import defaultdict

allele_dir = sys.argv[1]
output_file = sys.argv[2]
samples = sys.argv[3:]

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
        fh.readline()
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

covered = sum(1 for d in merged.values() if d['ref_total'] + d['alt_total'] > 0)
total_ref = sum(d['ref_total'] for d in merged.values())
total_alt = sum(d['alt_total'] for d in merged.values())
allelic = total_ref + total_alt
print(f"\n  Summary:")
print(f"    SNPs with allelic coverage: {covered:,}")
print(f"    Total ref reads at SNPs:    {total_ref:,}")
print(f"    Total alt reads at SNPs:    {total_alt:,}")
if allelic > 0:
    print(f"    Global ref fraction:        {100*total_ref/allelic:.2f}%")
for thr in (5, 10, 20):
    n = sum(1 for d in merged.values() if d['ref_total'] + d['alt_total'] >= thr)
    print(f"    SNPs with >={thr} allelic reads: {n:,}")
PYEOF

# =============================================================================
# 6. Compile final QC report
# =============================================================================
step_header "Compiling final QC report"

QC_REPORT="${QC_DIR}/pipeline_qc_report.txt"
{
    echo "============================================================"
    echo "F121-9 RNA-seq Pipeline — QC Report"
    echo "Generated: $(date)"
    echo "Pipeline:  ${SCRIPT_DIR}"
    echo "Output:    ${SCRATCH_DIR}"
    echo "Samples:   ${SAMPLE_ORDER[*]}"
    echo "============================================================"
    echo ""
    echo "Pipeline notes:"
    echo "  - Aligner: STAR 2-pass with WASP via vW SAMtag"
    echo "  - Library: TruSeq Stranded Total RNA (reverse-stranded)"
    echo "  - Quantification: featureCounts -s 2 -p --countReadPairs"
    echo "  - WASP filter: keep vW:i:1 + reads with no vW tag"
    echo ""

    echo "--- 1. Raw read counts (step 01) ---"
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        F="${QC_DIR}/read_counts_raw_${SAMPLE}.tsv"
        [[ -f "${F}" ]] && cat "${F}" || echo -e "${SAMPLE}\t?\t?"
    done
    echo ""

    echo "--- 2. Trimming (step 02) ---"
    echo -e "sample\traw_pairs\ttoo_short\twritten\tpct_short\tpct_surviving"
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        F="${QC_DIR}/trimming_stats_${SAMPLE}.tsv"
        [[ -f "${F}" ]] && cat "${F}" || echo -e "${SAMPLE}\t?\t?\t?\t?\t?"
    done
    echo ""

    echo "--- 3. STAR alignment (step 03) ---"
    H=0
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        F="${QC_DIR}/star_align_${SAMPLE}.tsv"
        if [[ -f "${F}" ]]; then
            (( H == 0 )) && { head -1 "${F}"; H=1; }
            tail -n +2 "${F}"
        fi
    done
    echo ""

    echo "--- 4. Quantification (step 04) ---"
    H=0
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        F="${QC_DIR}/quantify_stats_${SAMPLE}.tsv"
        if [[ -f "${F}" ]]; then
            (( H == 0 )) && { head -1 "${F}"; H=1; }
            tail -n +2 "${F}"
        fi
    done
    echo ""

    echo "--- 5. WASP tag filter (step 06) ---"
    H=0
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        F="${QC_DIR}/wasp_filter_${SAMPLE}.tsv"
        if [[ -f "${F}" ]]; then
            (( H == 0 )) && { head -1 "${F}"; H=1; }
            tail -n +2 "${F}"
        fi
    done
    echo ""

    echo "--- 6. Allele splitting (step 06) ---"
    H=0
    for SAMPLE in "${SAMPLE_ORDER[@]}"; do
        F="${QC_DIR}/allele_split_${SAMPLE}.tsv"
        if [[ -f "${F}" ]]; then
            (( H == 0 )) && { head -1 "${F}"; H=1; }
            tail -n +2 "${F}"
        fi
    done
    echo ""

    echo "--- 7. Gene expression summary table ---"
    if [[ -f "${EXPR_OUT}" ]]; then
        N=$(($(wc -l < "${EXPR_OUT}") - 1))
        EXPRESSED=$(awk -F'\t' 'NR>1 && $(NF-2) >= 1' "${EXPR_OUT}" | wc -l)
        echo "  File:                 ${EXPR_OUT}"
        echo "  Total genes:          ${N}"
        echo "  Expressed (TPM>=1):   ${EXPRESSED}"
        echo "  Top 5 by TPM mean:"
        # tpm_mean is the third-from-last column
        awk -F'\t' 'NR>1 {print $1"\t"$2"\t"$3"\t"$(NF-2)}' "${EXPR_OUT}" \
            | sort -k4,4nr | head -5 \
            | awk -F'\t' '{printf "    %s (%s, %s)  TPM=%s\n", $1, $2, $3, $4}'
    fi
    echo ""

    echo "--- 8. Per-SNP allele counts (merged) ---"
    if [[ -f "${MERGED_COUNTS}" ]]; then
        N=$(($(wc -l < "${MERGED_COUNTS}") - 1))
        echo "  File:        ${MERGED_COUNTS}"
        echo "  Total SNPs:  ${N}"
    fi
    echo ""

    echo "--- 9. Output bigWigs ---"
    echo "Per-rep stranded (${BIGWIG_IND_DIR}/):"
    ls -lh "${BIGWIG_IND_DIR}"/*.bw 2>/dev/null | awk '{print "  "$NF" ("$5")"}' || echo "  (none)"
    echo ""
    echo "Merged stranded (${BIGWIG_MERGED_DIR}/):"
    ls -lh "${BIGWIG_MERGED_DIR}"/*.bw 2>/dev/null | awk '{print "  "$NF" ("$5")"}' || echo "  (none)"
    echo ""
    echo "Per-allele stranded (${BIGWIG_ALLELE_DIR}/):"
    ls -lh "${BIGWIG_ALLELE_DIR}"/*.bw 2>/dev/null | awk '{print "  "$NF" ("$5")"}' || echo "  (none)"
    echo ""
    echo "Merged allele stranded (${BIGWIG_ALLELE_MERGED_DIR}/):"
    ls -lh "${BIGWIG_ALLELE_MERGED_DIR}"/*.bw 2>/dev/null | awk '{print "  "$NF" ("$5")"}' || echo "  (none)"
    echo ""

    echo "============================================================"
    echo "End of report."
    echo "Source TSVs:                  ${QC_DIR}/"
    echo "STAR BAMs (unfiltered):       ${BAM_STAR_DIR}/"
    echo "WASP-filtered BAMs:           ${BAM_FINAL_DIR}/"
    echo "Allele-split BAMs:            ${ALLELE_DIR}/"
    echo "Splice junctions:             ${BAM_STAR_DIR}/{SAMPLE}_SJ.out.tab"
    echo "============================================================"
} > "${QC_REPORT}"

echo "[INFO] QC report: ${QC_REPORT}"
echo ""
echo "----- Preview (first 80 lines) -----"
head -80 "${QC_REPORT}"
echo "----- ... -----"

step_header "Step 08 COMPLETE — Pipeline finished"
echo ""
echo "Final outputs:"
echo "  STAR BAMs:                    ${BAM_STAR_DIR}/"
echo "  WASP-filtered BAMs:           ${BAM_FINAL_DIR}/"
echo "  Allele-split BAMs:            ${ALLELE_DIR}/{SAMPLE}_{ref,alt,nosnp,ambiguous}.bam"
echo "  Splice junctions:             ${BAM_STAR_DIR}/{SAMPLE}_SJ.out.tab"
echo "  Per-rep stranded bigwigs:     ${BIGWIG_IND_DIR}/"
echo "  Merged stranded bigwigs:      ${BIGWIG_MERGED_DIR}/"
echo "  Per-allele stranded bigwigs:  ${BIGWIG_ALLELE_DIR}/"
echo "  Merged allele stranded:       ${BIGWIG_ALLELE_MERGED_DIR}/"
echo "  Per-SNP merged counts:        ${MERGED_COUNTS}"
echo "  >>> Gene expression summary:  ${EXPR_OUT}"
echo "  Final QC report:              ${QC_REPORT}"
echo ""
echo "Downstream consumers:"
echo "  05_motif_annot/02f_gene_expression.R ← gene_expression_summary.tsv"
echo "  06_allele_pairs/05_splice_jx_alleles.R ← allele_specific BAMs + SJ.out.tab"
