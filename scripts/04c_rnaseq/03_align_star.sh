#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04c_rnaseq/03_align_star.sh — STAR alignment with WASP + 2-pass + GeneCounts
# =============================================================================
#
# WHAT THIS DOES (per sample):
#   1. Aligns trimmed PE reads to mm39 with STAR
#   2. Uses --twopassMode Basic for novel splice junction discovery
#   3. Adds vW SAM tags via --waspOutputMode SAMtag (for allele filtering)
#   4. Produces stranded gene counts via --quantMode GeneCounts
#   5. Outputs sorted+indexed BAM, SJ.out.tab, ReadsPerGene.out.tab
#
# WHY STAR'S BUILT-IN WASP:
#   Unlike the bowtie2-based pipelines (03_atac/, 04_proseq/) which require
#   an explicit find-intersecting-snps -> remap -> filter cycle, STAR has
#   native WASP support. It tags each read with a vW:i:N integer based on
#   whether the alignment changes when the allele is swapped. We can filter
#   later by simple SAM-tag selection in step 06.
#
# vW TAG VALUES (STAR docs):
#   vW:i:1  passed WASP test (alignment unchanged with allele swap)
#   vW:i:2  multi-mapping when allele swapped (mapping bias)
#   vW:i:3  different number of mismatches with allele swap
#   vW:i:4  ungapped alignment with different start position
#   vW:i:5  different soft-clipping on right side
#   vW:i:6  read overlaps a SNP, but VCF genotype is missing
#   vW:i:7  read overlaps a SNP, ref/alt indistinguishable
#   no tag  read does not overlap any SNP (kept as-is)
#
# WHY 2-PASS:
#   STAR's first pass discovers novel junctions; second pass aligns reads
#   using both the discovered junctions AND the GTF-supplied junctions from
#   the index. Better sensitivity for novel splicing events. The original
#   GSE200699 paper used --outFilterType BySJout for similar reasons.
#
# WHY --quantMode GeneCounts:
#   STAR can produce a stranded gene-count matrix as a byproduct of
#   alignment. The output ReadsPerGene.out.tab has 4 columns:
#     gene_id, unstranded_count, forward_stranded_count, reverse_stranded_count
#   For TruSeq Stranded reverse libraries we read column 4. This is
#   essentially free (no extra alignment passes) and lets us cross-check
#   featureCounts output in step 04.
#
# WHY ALLELE_VCF_STAR (not the original parental VCF):
#   STAR's --varVCFfile uses the first sample column and looks for 0/1
#   genotypes. The source parental VCF has two homozygous parent samples
#   (every genotype is 1/1 or 0/0) so STAR finds no het SNPs and aborts.
#   Step 00 synthesizes ALLELE_VCF_STAR — a single-sample VCF with one
#   record per F1-het site, all genotyped 0/1 — which is what STAR needs.
#
# RESOURCES:
#   STAR with mm39 + sjdb + 2-pass + WASP needs ~32 GB peak RAM.
#   Wall time: 4-8 hr per sample for ~20-40M PE reads.
#
# PARALLEL EXECUTION:
#   sbatch --array=0-1 --partition=short --time=12:00:00 --mem=64G \
#       --cpus-per-task=8 -o logs/03_%A_%a.out -e logs/03_%A_%a.err \
#       03_align_star.sh
# =============================================================================

source "config.sh"

step_header "RNA-seq Pipeline Step 03: STAR Alignment + WASP + GeneCounts"

source activate "${CONDA_ENV_NAME}"
echo "[INFO] STAR:     $(STAR --version 2>&1 | head -1)"
echo "[INFO] samtools: $(samtools --version | head -1)"

create_dirs
resolve_samples "${1:-}"

# Pre-flight checks
[[ -f "${STAR_INDEX_DIR}/SAindex" ]] || \
    { echo "[ERROR] STAR index missing. Run 00_setup_references.sh first." >&2; exit 1; }
[[ -f "${ALLELE_VCF_STAR}" ]] || \
    { echo "[ERROR] STAR VCF missing: ${ALLELE_VCF_STAR}" >&2;
      echo "[ERROR] Run 00_setup_references.sh first." >&2; exit 1; }

for SAMPLE in "${RUN_SAMPLES[@]}"; do
    step_header "Aligning ${SAMPLE}"

    R1="${TRIMMED_DIR}/${SAMPLE}_R1_trimmed.fastq.gz"
    R2="${TRIMMED_DIR}/${SAMPLE}_R2_trimmed.fastq.gz"

    [[ -f "${R1}" && -f "${R2}" ]] || \
        { echo "[ERROR] Trimmed FASTQs missing for ${SAMPLE}" >&2; exit 1; }

    SAMPLE_STAR_PREFIX="${BAM_STAR_DIR}/${SAMPLE}_"
    BAM_OUT="${SAMPLE_STAR_PREFIX}Aligned.sortedByCoord.out.bam"
    SJ_OUT="${SAMPLE_STAR_PREFIX}SJ.out.tab"
    GENE_COUNTS="${SAMPLE_STAR_PREFIX}ReadsPerGene.out.tab"
    LOG_FINAL="${SAMPLE_STAR_PREFIX}Log.final.out"

    if [[ -f "${BAM_OUT}" && -f "${BAM_OUT}.bai" && -f "${SJ_OUT}" && -f "${GENE_COUNTS}" ]]; then
        echo "[INFO] STAR outputs exist for ${SAMPLE}. Skipping alignment."
    else
        echo "[INFO] Running STAR (2-pass + WASP) for ${SAMPLE}..."

        # Per-sample tmp dir (STAR is picky about pre-existing dirs;
        # remove and let STAR create it)
        STAR_TMP="${SCRATCH_DIR}/tmp/star_${SAMPLE}"
        rm -rf "${STAR_TMP}"

        STAR \
            --runMode alignReads \
            --runThreadN "${THREADS}" \
            --genomeDir "${STAR_INDEX_DIR}" \
            --readFilesIn "${R1}" "${R2}" \
            --readFilesCommand zcat \
            --outFileNamePrefix "${SAMPLE_STAR_PREFIX}" \
            --outTmpDir "${STAR_TMP}" \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMstrandField intronMotif \
            --outSAMattributes NH HI AS nM MD vW \
            --outFilterMismatchNmax "${STAR_OUTFILTER_MISMATCH}" \
            --outMultimapperOrder Random \
            --outFilterMultimapNmax "${STAR_OUTFILTER_MULTIMAP}" \
            --outFilterType "${STAR_OUTFILTER_TYPE}" \
            --alignSJoverhangMin "${STAR_ALIGN_SJ_OVERHANG}" \
            --outFilterIntronMotifs "${STAR_INTRON_MOTIFS}" \
            --twopassMode Basic \
            --varVCFfile "${ALLELE_VCF_STAR}" \
            --waspOutputMode SAMtag \
            --quantMode GeneCounts \
            --outBAMsortingThreadN 4

        # Index the sorted BAM
        samtools index "${BAM_OUT}"

        # Move STAR's auxiliary log files to the LOG_DIR for tidiness
        for ext in Log.out Log.progress.out Log.final.out; do
            F="${SAMPLE_STAR_PREFIX}${ext}"
            [[ -f "${F}" ]] && cp "${F}" "${LOG_DIR}/${SAMPLE}_STAR_${ext}"
        done
    fi

    # =========================================================================
    # Sanity checks
    # =========================================================================
    if [[ ! -f "${LOG_FINAL}" ]]; then
        echo "[WARN] STAR Log.final.out missing; skipping detailed stats parse."
    else
        # Parse STAR final log
        INPUT_READS=$(grep "Number of input reads" "${LOG_FINAL}" | awk -F'\t' '{print $2}' | tr -d ' ')
        UNIQ_MAP=$(grep "Uniquely mapped reads number" "${LOG_FINAL}" | awk -F'\t' '{print $2}' | tr -d ' ')
        UNIQ_PCT=$(grep "Uniquely mapped reads %" "${LOG_FINAL}" | awk -F'\t' '{print $2}' | tr -d ' %')
        MULTI_MAP=$(grep "Number of reads mapped to multiple loci" "${LOG_FINAL}" | awk -F'\t' '{print $2}' | tr -d ' ')
        MULTI_PCT=$(grep "% of reads mapped to multiple loci" "${LOG_FINAL}" | awk -F'\t' '{print $2}' | tr -d ' %')
        SPLICE=$(grep "Number of splices: Total" "${LOG_FINAL}" | awk -F'\t' '{print $2}' | tr -d ' ')
        ANNOT_SPLICE=$(grep "Number of splices: Annotated" "${LOG_FINAL}" | awk -F'\t' '{print $2}' | tr -d ' ')

        echo ""
        echo "[SANITY CHECK] STAR alignment for ${SAMPLE}:"
        echo "  Input read pairs:                  ${INPUT_READS}"
        echo "  Uniquely mapped:                   ${UNIQ_MAP} (${UNIQ_PCT}%)"
        echo "  Multi-mapped:                      ${MULTI_MAP} (${MULTI_PCT}%)"
        echo "  Total splices:                     ${SPLICE}"
        echo "  Annotated splices:                 ${ANNOT_SPLICE}"
        echo "[INTERPRETATION]"
        echo "  Unique mapping >75%: good"
        echo "  Unique mapping <60%: investigate"
    fi

    # =========================================================================
    # WASP tag distribution (full BAM)
    # =========================================================================
    # ⚠ pipefail safety: the previous version was
    #     samtools view BAM | head -1000000 | awk ... | sort
    # When `head` closes after 1M lines, `samtools view` is still writing
    # and gets SIGPIPE (exit 141). Under `set -o pipefail` that kills the
    # script silently — visible symptom is that the WASP distribution
    # prints but no later output (strand check, QC TSV) appears.
    # Fix: drop `head` entirely. awk processes every alignment, samtools
    # view runs to completion, no SIGPIPE possible. ~60s extra per sample
    # for ~50M-read BAMs; negligible vs. STAR runtime.
    echo ""
    echo "[WASP TAG DISTRIBUTION] (full BAM ${BAM_OUT}):"
    samtools view "${BAM_OUT}" \
    | awk '
        {
            wasp_tag = "no_vW_tag"
            for (i=12; i<=NF; i++) {
                if ($i ~ /^vW:i:/) { wasp_tag = $i; break }
            }
            counts[wasp_tag]++
            total++
        }
        END {
            for (k in counts) {
                printf "  %-20s %12d  (%5.2f%%)\n", k, counts[k], 100*counts[k]/total
            }
        }
    ' | sort

    echo ""
    echo "[INFO] vW:i:1 = WASP-pass (kept for allele analysis)"
    echo "       no_vW_tag = read does not overlap any SNP (kept)"
    echo "       vW:i:2-7 = WASP-fail variants (filtered out in step 06)"

    # =========================================================================
    # Strand check on STAR's gene counts
    # =========================================================================
    if [[ -f "${GENE_COUNTS}" ]]; then
        echo ""
        echo "[SANITY CHECK] STAR ReadsPerGene strand inference:"
        # Skip first 4 rows (N_unmapped, N_multimapping, N_noFeature, N_ambiguous)
        TOTAL_UNSTR=$(awk 'NR>4 {s+=$2} END {print s+0}' "${GENE_COUNTS}")
        TOTAL_FWD=$(awk 'NR>4 {s+=$3} END {print s+0}' "${GENE_COUNTS}")
        TOTAL_REV=$(awk 'NR>4 {s+=$4} END {print s+0}' "${GENE_COUNTS}")
        echo "  Reads in genes (unstranded col):  ${TOTAL_UNSTR}"
        echo "  Reads in genes (forward col):     ${TOTAL_FWD}"
        echo "  Reads in genes (reverse col):     ${TOTAL_REV}"
        echo "  For TruSeq Stranded reverse: column 4 (reverse) should be"
        echo "  much LARGER than column 3 (forward). If similar, library may"
        echo "  not be reverse-stranded; if reversed, library is forward-stranded."
    fi

    # Per-sample QC TSV
    cat > "${QC_DIR}/star_align_${SAMPLE}.tsv" << EOF
sample	input_reads	unique_mapped	unique_pct	multi_mapped	multi_pct	splices_total	splices_annotated	reads_in_genes_unstranded	reads_in_genes_forward	reads_in_genes_reverse
${SAMPLE}	${INPUT_READS:-NA}	${UNIQ_MAP:-NA}	${UNIQ_PCT:-NA}	${MULTI_MAP:-NA}	${MULTI_PCT:-NA}	${SPLICE:-NA}	${ANNOT_SPLICE:-NA}	${TOTAL_UNSTR:-NA}	${TOTAL_FWD:-NA}	${TOTAL_REV:-NA}
EOF

done

step_header "Step 03 COMPLETE for: ${RUN_SAMPLES[*]}"
echo "Output:"
echo "  STAR BAMs:           ${BAM_STAR_DIR}/{SAMPLE}_Aligned.sortedByCoord.out.bam"
echo "  Splice junctions:    ${BAM_STAR_DIR}/{SAMPLE}_SJ.out.tab"
echo "  Gene counts (STAR):  ${BAM_STAR_DIR}/{SAMPLE}_ReadsPerGene.out.tab"
echo "  QC stats:            ${QC_DIR}/star_align_*.tsv"
echo ""
echo "Next: bash 04_quantify.sh   (and bash 05_make_bigwigs.sh in parallel)"
