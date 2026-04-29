#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 04c_rnaseq/00_setup_references.sh — STAR index for mm39 + vM38 GTF
# =============================================================================
#
# WHAT THIS DOES:
#   1. Verifies mm39 FASTA is present (REUSE from earlier pipelines)
#   2. Builds STAR index with sjdb annotations from vM38 GTF
#   3. Verifies VCF for STAR's WASP is present
#   4. Synthesizes a STAR-compatible single-sample F1 VCF
#   5. Builds WASP per-chromosome SNP files for step 06
#
# WHY STAR (NOT BOWTIE2):
#   RNA-seq reads span splice junctions. STAR is splice-aware; bowtie2 isn't.
#
# DESIGN — REUSE FIRST:
#   mm39 FASTA already exists from 03_atac/ or 04_proseq/. We reuse it.
#   The STAR index is RNA-seq-specific (different from bowtie2 indices) and
#   gets built fresh — but only once, since it's also reused on subsequent
#   pipeline runs. ~1-2 hr build, ~30 GB output, ~30 GB peak RAM.
#
# WHY --sjdbOverhang 99:
#   STAR docs recommend (max read length - 1). Reads are trimmed to ~100 nt
#   in the original GSE200699 protocol; we'll match that.
#
# SBATCH RESOURCES:
#   sbatch --partition=short --time=4:00:00 --mem=64G --cpus-per-task=8 \
#       -o logs/00_%j.out -e logs/00_%j.err 00_setup_references.sh
# =============================================================================

source "config.sh"

step_header "RNA-seq Pipeline Step 00: Reference Setup"

source activate "${CONDA_ENV_NAME}"

# Verify STAR is installed
if ! command -v STAR &>/dev/null; then
    echo "[ERROR] STAR not found in conda env ${CONDA_ENV_NAME}." >&2
    echo "[ERROR] Install with: conda install -n ${CONDA_ENV_NAME} -c bioconda star" >&2
    exit 1
fi

if ! command -v featureCounts &>/dev/null; then
    echo "[WARN] featureCounts not found. Step 04 will fail until installed." >&2
    echo "[WARN] Install with: conda install -n ${CONDA_ENV_NAME} -c bioconda subread" >&2
fi

echo "[INFO] STAR:           $(STAR --version 2>&1 | head -1)"
echo "[INFO] featureCounts:  $(featureCounts -v 2>&1 | head -1 || echo 'NOT INSTALLED')"
echo "[INFO] samtools:       $(samtools --version | head -1)"

create_dirs

# -----------------------------------------------------------------------------
# 1. mm39 FASTA must be present (built by 03_atac/ or 04_proseq/)
# -----------------------------------------------------------------------------
if [[ ! -f "${MM39_FA}" ]]; then
    echo "[ERROR] mm39 FASTA missing: ${MM39_FA}" >&2
    echo "[ERROR] Run 03_atac/00_setup_references.sh or 04_proseq/00_setup_references.sh first," >&2
    echo "[ERROR] which downloads and indexes the mm39 reference." >&2
    exit 1
fi
echo "[OK]  mm39 FASTA: ${MM39_FA} ($(du -h "${MM39_FA}" | cut -f1))"

# -----------------------------------------------------------------------------
# 2. vM38 GTF must be present (built by 01_resources/)
# -----------------------------------------------------------------------------
if [[ ! -f "${GTF_FILE}" ]]; then
    echo "[ERROR] vM38 GTF missing: ${GTF_FILE}" >&2
    echo "[ERROR] This should have been provisioned by 01_resources/." >&2
    exit 1
fi
echo "[OK]  vM38 GTF: ${GTF_FILE} ($(du -h "${GTF_FILE}" | cut -f1))"

# STAR can't read .gz GTF directly; decompress if needed
GTF_USE="${GTF_FILE}"
if [[ "${GTF_FILE}" == *.gz ]]; then
    GTF_DECOMPRESSED="${GTF_FILE%.gz}"
    if [[ ! -f "${GTF_DECOMPRESSED}" ]]; then
        echo "[INFO] Decompressing GTF for STAR (keeping .gz)..."
        gunzip -k "${GTF_FILE}"
    fi
    GTF_USE="${GTF_DECOMPRESSED}"
    echo "[INFO] Using decompressed GTF: ${GTF_USE}"
fi

# -----------------------------------------------------------------------------
# 3. STAR index
# -----------------------------------------------------------------------------
if [[ -f "${STAR_INDEX_DIR}/SAindex" ]]; then
    echo ""
    echo "[INFO] STAR index already present at ${STAR_INDEX_DIR}. Reusing."
    echo "       To rebuild, remove the directory first."
    # ls then awk-truncate to avoid a `| head` pipe (SIGPIPE-with-pipefail trap)
    ls -lh "${STAR_INDEX_DIR}/" | awk 'NR<=10 {print "    "$0}'
else
    echo ""
    echo "[INFO] Building STAR index (1-2 hours; needs ~30-50 GB RAM)..."
    mkdir -p "${STAR_INDEX_DIR}"

    # --genomeSAindexNbases 14 is the default for mm39 (calculated as
    # min(14, log2(genome_size)/2 - 1) ≈ 14 for ~2.7 Gb mouse).
    STAR \
        --runMode genomeGenerate \
        --runThreadN "${THREADS}" \
        --genomeDir "${STAR_INDEX_DIR}" \
        --genomeFastaFiles "${MM39_FA}" \
        --sjdbGTFfile "${GTF_USE}" \
        --sjdbOverhang "${SJDB_OVERHANG}" \
        --outFileNamePrefix "${STAR_INDEX_DIR}/_buildlog_"

    echo "[OK]  STAR index built at ${STAR_INDEX_DIR}"
fi

# -----------------------------------------------------------------------------
# 4. VCF for STAR's WASP — verify present, decompress if needed
# -----------------------------------------------------------------------------
echo ""
if [[ ! -f "${ALLELE_VCF}" ]]; then
    echo "[ERROR] Allele VCF missing: ${ALLELE_VCF}" >&2
    exit 1
fi
echo "[OK]  Source allele VCF: ${ALLELE_VCF}"

# We always need the uncompressed VCF available — both for the synthetic
# F1 build below and as the canonical decompressed copy other tools may use.
if [[ ! -f "${ALLELE_VCF_UNZIPPED}" ]]; then
    echo "[INFO] Decompressing source VCF (keeping .gz)..."
    gunzip -kc "${ALLELE_VCF}" > "${ALLELE_VCF_UNZIPPED}"
    echo "[OK]  Decompressed: ${ALLELE_VCF_UNZIPPED}"
else
    echo "[OK]  Decompressed VCF already present: ${ALLELE_VCF_UNZIPPED}"
fi

# -----------------------------------------------------------------------------
# 5. Build STAR-compatible single-sample F1 het VCF
# -----------------------------------------------------------------------------
# The source VCF has TWO PARENTAL samples (129S1_SvImJ + CAST_EiJ). Each
# parent is scored as homozygous at every site (1/1 or 0/0). STAR's
# --varVCFfile uses the first sample column and looks for 0/1 genotypes;
# with the parental layout it finds zero het sites and aborts:
#     "could not find any SNPs in VCF file"
#
# We synthesize a single-sample F1 VCF here: every site where the parents
# are OPPOSITE homozygotes (so the F1 is heterozygous) becomes a 0/1 record
# for a notional "F121_9" sample. As a side benefit we drop the giant CSQ
# INFO blocks, shrinking the file from ~20 GB to a few hundred MB.
#
# This file is STAR-only. lib/make_wasp_snps.sh and split_by_allele.py
# continue to use the original parental ALLELE_VCF — they key off
# (chrom, pos, ref, alt) and don't care about how genotypes are encoded.
echo ""
if [[ ! -f "${ALLELE_VCF_STAR}" ]]; then
    echo "[INFO] Building synthetic F1 single-sample VCF for STAR..."
    echo "       Source: ${ALLELE_VCF_UNZIPPED}"
    echo "       Output: ${ALLELE_VCF_STAR}"
    echo "       (this can take a few minutes for a large parental VCF)"

    awk 'BEGIN{OFS="\t"; kept=0; nonbi=0; concord=0; missing=0}
        /^##/ { print; next }
        /^#CHROM/ {
            # Replace the two parent sample columns with one F1 column
            printf "%s", $1
            for (i=2; i<=9; i++) printf "\t%s", $i
            print "\tF121_9"
            next
        }
        {
            # Biallelic SNPs only: REF and ALT both single chars, no comma
            if (length($4)!=1 || length($5)!=1 || $5 ~ /,/) { nonbi++; next }
            split($10, p1, ":"); gt1=p1[1]; gsub(/\|/, "/", gt1)
            split($11, p2, ":"); gt2=p2[1]; gsub(/\|/, "/", gt2)
            # Skip missing genotypes
            if (gt1 ~ /\./ || gt2 ~ /\./) { missing++; next }
            # Keep only sites where parents are opposite homozygotes
            if ((gt1=="0/0" && gt2=="1/1") || (gt1=="1/1" && gt2=="0/0")) {
                # Strip INFO (saves enormous space; STAR ignores INFO anyway)
                printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\tGT\t0/1\n", \
                    $1, $2, $3, $4, $5, $6, $7
                kept++
            } else if (gt1==gt2) { concord++ }
        }
        END {
            printf "  Kept (F1 het sites):         %d\n", kept     > "/dev/stderr"
            printf "  Skipped non-biallelic/indel: %d\n", nonbi    > "/dev/stderr"
            printf "  Skipped parent-concordant:   %d\n", concord  > "/dev/stderr"
            printf "  Skipped missing genotypes:   %d\n", missing  > "/dev/stderr"
        }' "${ALLELE_VCF_UNZIPPED}" > "${ALLELE_VCF_STAR}.tmp"
    mv "${ALLELE_VCF_STAR}.tmp" "${ALLELE_VCF_STAR}"

    NREC=$(grep -vc "^#" "${ALLELE_VCF_STAR}" || true)
    echo "[OK]  STAR VCF: ${ALLELE_VCF_STAR}"
    echo "      Records: ${NREC}    Size: $(du -h "${ALLELE_VCF_STAR}" | cut -f1)"
    if (( NREC == 0 )); then
        echo "[ERROR] Synthetic VCF is empty — check parent column order/genotypes." >&2
        echo "[ERROR] Inspect with:" >&2
        echo "[ERROR]   grep '^#CHROM' ${ALLELE_VCF_UNZIPPED}" >&2
        echo "[ERROR]   grep -v '^#' ${ALLELE_VCF_UNZIPPED} | head -3" >&2
        exit 1
    fi
else
    NREC=$(grep -vc "^#" "${ALLELE_VCF_STAR}" || true)
    echo "[OK]  STAR VCF already present: ${ALLELE_VCF_STAR} (${NREC} records)"
fi

# Chromosome-naming sanity check (run unconditionally so it's visible on re-runs).
#
# ⚠ pipefail safety: the previous version of this check used
#   `grep ... | cut | head -100000 | sort -u | head -5 | tr ...`
# Once `head -5` closes its stdin, upstream `sort` (and possibly `grep`/`cut`)
# write to a broken pipe, get SIGPIPE, exit 141. With `set -o pipefail` that
# kills the script silently right here. The fix below uses `awk` with `exit`
# (clean stdout close, no upstream broken pipe) and a `head FILE` form (head
# reads the file directly, no upstream pipe at all).
VCF_CHR_SAMPLE=$(awk '!/^#/ && !($1 in seen) {seen[$1]=1; print $1; if (++c>=5) exit}' \
                     "${ALLELE_VCF_STAR}" | tr '\n' ' ')
GENOME_CHR_SAMPLE=$(head -5 "${MM39_CHROM_SIZES}" | cut -f1 | tr '\n' ' ')
echo "[INFO] First chroms in STAR VCF:    ${VCF_CHR_SAMPLE}"
echo "[INFO] First chroms in mm39:        ${GENOME_CHR_SAMPLE}"
if ! echo "${VCF_CHR_SAMPLE}" | grep -q "chr"; then
    echo "[WARN] STAR VCF chroms don't have 'chr' prefix but mm39 does."
    echo "[WARN] STAR will silently match nothing. Fix the source VCF naming."
fi

# -----------------------------------------------------------------------------
# 6. WASP SNP files (for lib/split_by_allele.py in step 06)
# -----------------------------------------------------------------------------
echo ""
echo "[INFO] Building WASP per-chromosome SNP files (needed by step 06)..."
if ls "${WASP_SNP_DIR}"/chr*.snps.txt.gz 1>/dev/null 2>&1; then
    NCHR=$(ls "${WASP_SNP_DIR}"/chr*.snps.txt.gz | wc -l)
    echo "[OK]  WASP SNP files exist (${NCHR} chromosome files)"
else
    bash "${MAKE_WASP_SNPS_SH}" \
        --vcf "${ALLELE_VCF}" \
        --output_dir "${WASP_SNP_DIR}"
fi

# -----------------------------------------------------------------------------
# 7. Final pre-flight summary
# -----------------------------------------------------------------------------
echo ""
echo "[PRE-FLIGHT CHECK]"
ALL_OK=1
for F in "${MM39_FA}" "${GTF_USE}" "${ALLELE_VCF_STAR}" \
         "${STAR_INDEX_DIR}/SAindex" "${MM39_CHROM_SIZES}"; do
    if [[ -f "${F}" ]]; then
        echo "  [OK]    ${F}"
    else
        echo "  [MISSING] ${F}"
        ALL_OK=0
    fi
done

NCHR=$(ls "${WASP_SNP_DIR}"/chr*.snps.txt.gz 2>/dev/null | wc -l)
if (( NCHR > 0 )); then
    echo "  [OK]    ${WASP_SNP_DIR}/chr*.snps.txt.gz (${NCHR} files)"
else
    echo "  [MISSING] WASP SNP files"
    ALL_OK=0
fi

(( ALL_OK )) || { echo ""; echo "[ERROR] Pre-flight check failed." >&2; exit 1; }

step_header "Step 00 COMPLETE"
echo "Next: bash 01_download_fastq.sh   (or sbatch as array of ${#SAMPLE_ORDER[@]})"
