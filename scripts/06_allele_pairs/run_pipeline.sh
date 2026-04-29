#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 06_allele_pairs/run_pipeline.sh — Orchestrator
# =============================================================================
#
# Stages (each runs as a SLURM array of size N_motifs):
#   01 build_pair_index        per-motif: build (motif, SNP) pair index
#   02 score_pwm               per-motif: PWM ref/alt for body-SNPs
#   03 count_atac              per-motif: ATAC ref/alt counts at motif body
#   04 count_proseq_window     per-motif: PRO-seq pause + gene body counts
#   04b extract_proseq_profile per-motif: ±500 bp per-base PRO-seq profile
#   04c compute_pause_indices  per-motif: pause-index columns
#   05 extract_rnaseq_junc     per-motif: junction use/skip counts at 8 cells
#   05b compute_psi            per-motif: PSI per allele per cell
#   06 assemble_pair_table     per-motif: master per-pair table
#   07 pair_qc                 single:    cross-motif QC report
#
# DEPENDENCY GRAPH:
#   01 → 02
#   01 → 03 ─┐
#   01 → 04 → 04c
#   01 → 04b
#   01 → 05 → 05b
#   {02, 03, 04c, 05b} → 06
#   06 → 07
#
# USAGE:
#   bash run_pipeline.sh                           # dry run (default)
#   bash run_pipeline.sh --apply                   # actually submit
#   bash run_pipeline.sh --apply --start-from 04   # skip 01-03 (already done)
# =============================================================================

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

SUBMIT=false
START_FROM=01
while [[ $# -gt 0 ]]; do
    case "$1" in
        --apply|--submit)   SUBMIT=true; shift ;;
        --dry-run)          SUBMIT=false; shift ;;
        --start-from)       START_FROM="$2"; shift 2 ;;
        -h|--help)
            grep -E '^# ' "$0" | sed 's/^# //' | head -40
            exit 0 ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

create_dirs

# -----------------------------------------------------------------------------
# Discover motifs
# -----------------------------------------------------------------------------
mapfile -t BED_FILES < <(find "${MOTIF_BED_DIR}" -maxdepth 1 -name "*.bed" -type f | sort)
N_MOTIFS=${#BED_FILES[@]}

if (( N_MOTIFS == 0 )); then
    echo "[ERROR] No .bed files in ${MOTIF_BED_DIR}" >&2
    exit 1
fi

# Build manifest (motif_id → BED path), task indices 0..N-1
MANIFEST="${OUTDIR}/_motif_manifest.txt"
> "${MANIFEST}"
for i in "${!BED_FILES[@]}"; do
    bed="${BED_FILES[$i]}"
    motif_id=$(basename "${bed}" .bed)
    motif_id_short="${motif_id%%:*}"   # strip trailing colon-suffix if any
    echo -e "${i}\t${motif_id_short}\t${bed}" >> "${MANIFEST}"
done

ARRAY_SPEC="0-$((N_MOTIFS - 1))%${DEFAULT_MAX_CONCURRENT}"

echo "============================================================"
echo "06_allele_pairs/ pipeline"
echo "============================================================"
echo "  Mode:           $( ${SUBMIT} && echo APPLY || echo DRY-RUN )"
echo "  Start from:     ${START_FROM}"
echo "  Motifs:         ${N_MOTIFS}"
echo "  Manifest:       ${MANIFEST}"
echo "  Array spec:     ${ARRAY_SPEC}"
echo "  Output:         ${OUTDIR}"
echo ""

# -----------------------------------------------------------------------------
# Helper: write per-stage wrapper script that reads manifest and runs cmd.
# Args: $1 = stage_id, $2 = stage_name, $3 = command_template (uses
# {MOTIF_ID} and {MOTIF_BED} placeholders).
# -----------------------------------------------------------------------------
make_stage_wrapper() {
    local stage_id="$1"
    local stage_name="$2"
    local cmd_tpl="$3"

    local wrapper="${OUTDIR}/_stage_${stage_id}.sh"
    cat > "${wrapper}" << EOF
#!/usr/bin/env bash
#SBATCH --output=${OUTDIR}/logs/${stage_id}_%A_%a.out
#SBATCH --error=${OUTDIR}/logs/${stage_id}_%A_%a.err
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

LINE=\$(awk -v idx="\${SLURM_ARRAY_TASK_ID}" '\$1 == idx' "${MANIFEST}")
[[ -n "\${LINE}" ]] || { echo "[ERROR] No manifest entry for task \${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
MOTIF_ID=\$(echo "\${LINE}" | cut -f2)
MOTIF_BED=\$(echo "\${LINE}" | cut -f3)

echo "Stage ${stage_id} (${stage_name}) — motif \${MOTIF_ID}"
echo "Date: \$(date)"

${cmd_tpl}

echo "Stage ${stage_id} done: \$(date)"
EOF
    chmod +x "${wrapper}"
    echo "${wrapper}"
}

# -----------------------------------------------------------------------------
# Submit helper (with dependency)
# -----------------------------------------------------------------------------
SUBMITTED_JOBID="NA"
submit_array() {
    local stage_label="$1"
    local sbatch_args="$2"
    local script_path="$3"
    local deps="${4:-}"

    echo ""
    echo "--- ${stage_label} ---"
    echo "  sbatch ${sbatch_args} ${deps} ${script_path}"

    if ${SUBMIT}; then
        out=$(sbatch ${sbatch_args} ${deps} --parsable "${script_path}")
        SUBMITTED_JOBID="${out%%;*}"
        echo "  -> Submitted as ${SUBMITTED_JOBID}"
    else
        SUBMITTED_JOBID="DRY_$(date +%N)_${stage_label}"
        echo "  (dry run)"
    fi
}

should_run() {
    local s="$1"
    [[ "${s}" > "${START_FROM}" || "${s}" == "${START_FROM}" ]]
}

base_sbatch="--array=${ARRAY_SPEC} --partition=${DEFAULT_PARTITION} \
--time=${DEFAULT_TIME} --mem=${DEFAULT_MEM} \
--cpus-per-task=${DEFAULT_THREADS}"

# -----------------------------------------------------------------------------
# Stage 01 — build pair index (R script)
# -----------------------------------------------------------------------------
JOB_01=""
if should_run 01; then
    SCRIPT=$(make_stage_wrapper "01" "build_pair_index" "$(cat <<EOC
set +u
source activate "${CONDA_ENV_R}"
set -u

OUT="${PAIR_INDEX_DIR}/\${MOTIF_ID}_pair_index.tsv.gz"
if [[ -f "\${OUT}" ]]; then
    echo "[INFO] Pair index exists; skipping."
    exit 0
fi

Rscript "$(dirname "${BASH_SOURCE[0]}")/01_build_pair_index.R" \\
    --motif_bed "\${MOTIF_BED}" \\
    --vcf "${VCF_DIR}/${HYBRID}_het_snps.ucsc.vcf.gz" \\
    --max_flank_bp ${MAX_PAIR_FLANK_BP} \\
    --motif_id "\${MOTIF_ID}" \\
    --outdir "${PAIR_INDEX_DIR}"
EOC
)")
    submit_array "01_build_pair_index" "${base_sbatch} --job-name=ap01" "${SCRIPT}"
    JOB_01="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# Stage 02 — PWM scoring
# -----------------------------------------------------------------------------
JOB_02=""
if should_run 02; then
    DEPS=""
    [[ -n "${JOB_01}" ]] && DEPS="--dependency=afterok:${JOB_01}"
    SCRIPT=$(make_stage_wrapper "02" "score_pwm" "$(cat <<EOC
set +u
source activate "${CONDA_ENV_R}"
set -u

PI="${PAIR_INDEX_DIR}/\${MOTIF_ID}_pair_index.tsv.gz"
[[ -f "\${PI}" ]] || { echo "[ERROR] Missing pair index"; exit 1; }

OUT="${PWM_DIR}/\${MOTIF_ID}_pwm_scores.tsv.gz"
if [[ -f "\${OUT}" ]]; then
    echo "[INFO] PWM output exists; skipping."
    exit 0
fi

Rscript "$(dirname "${BASH_SOURCE[0]}")/02_score_pwm_per_allele.R" \\
    --pair_index "\${PI}" \\
    --pwm_meme "${PWM_MEME_FILE}" \\
    --genome_fa "${GENOME_FA}" \\
    --motif_id "\${MOTIF_ID}" \\
    --pseudocount ${PWM_PSEUDOCOUNT} \\
    --outdir "${PWM_DIR}"
EOC
)")
    submit_array "02_score_pwm" "${base_sbatch} --job-name=ap02" "${SCRIPT}" "${DEPS}"
    JOB_02="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# Stage 03 — ATAC counts
# -----------------------------------------------------------------------------
JOB_03=""
if should_run 03; then
    DEPS=""
    [[ -n "${JOB_01}" ]] && DEPS="--dependency=afterok:${JOB_01}"
    SCRIPT=$(make_stage_wrapper "03" "count_atac" "$(cat <<EOC
PI="${PAIR_INDEX_DIR}/\${MOTIF_ID}_pair_index.tsv.gz"
[[ -f "\${PI}" ]] || { echo "[ERROR] Missing pair index"; exit 1; }

bash "$(dirname "${BASH_SOURCE[0]}")/03_count_atac_alleles.sh" \\
    --pair_index "\${PI}" \\
    --motif_id "\${MOTIF_ID}" \\
    --outdir "${ATAC_COUNTS_DIR}"
EOC
)")
    submit_array "03_count_atac" "${base_sbatch} --job-name=ap03" "${SCRIPT}" "${DEPS}"
    JOB_03="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# Stage 04 — PRO-seq window counts
# -----------------------------------------------------------------------------
JOB_04=""
if should_run 04; then
    DEPS=""
    [[ -n "${JOB_01}" ]] && DEPS="--dependency=afterok:${JOB_01}"
    SCRIPT=$(make_stage_wrapper "04" "count_proseq" "$(cat <<EOC
PI="${PAIR_INDEX_DIR}/\${MOTIF_ID}_pair_index.tsv.gz"
[[ -f "\${PI}" ]] || { echo "[ERROR] Missing pair index"; exit 1; }

bash "$(dirname "${BASH_SOURCE[0]}")/04_count_proseq_alleles.sh" \\
    --pair_index "\${PI}" \\
    --motif_id "\${MOTIF_ID}" \\
    --outdir "${PROSEQ_COUNTS_DIR}"
EOC
)")
    submit_array "04_count_proseq" "${base_sbatch} --job-name=ap04" "${SCRIPT}" "${DEPS}"
    JOB_04="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# Stage 04b — PRO-seq profile extraction
# -----------------------------------------------------------------------------
JOB_04B=""
if should_run 04b; then
    DEPS=""
    [[ -n "${JOB_01}" ]] && DEPS="--dependency=afterok:${JOB_01}"
    SCRIPT=$(make_stage_wrapper "04b" "extract_proseq_profile" "$(cat <<EOC
PI="${PAIR_INDEX_DIR}/\${MOTIF_ID}_pair_index.tsv.gz"
[[ -f "\${PI}" ]] || { echo "[ERROR] Missing pair index"; exit 1; }

bash "$(dirname "${BASH_SOURCE[0]}")/04b_extract_proseq_profile.sh" \\
    --pair_index "\${PI}" \\
    --motif_id "\${MOTIF_ID}" \\
    --outdir "${PROSEQ_PROFILE_DIR}"
EOC
)")
    submit_array "04b_extract_proseq_profile" "${base_sbatch} --job-name=ap04b --time=08:00:00" "${SCRIPT}" "${DEPS}"
    JOB_04B="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# Stage 04c — pause indices (depends on 04)
# -----------------------------------------------------------------------------
JOB_04C=""
if should_run 04c; then
    DEPS=""
    [[ -n "${JOB_04}" ]] && DEPS="--dependency=afterok:${JOB_04}"
    SCRIPT=$(make_stage_wrapper "04c" "compute_pause_indices" "$(cat <<EOC
set +u
source activate "${CONDA_ENV_R}"
set -u

WC="${PROSEQ_COUNTS_DIR}/\${MOTIF_ID}_proseq_window_counts.tsv.gz"
[[ -f "\${WC}" ]] || { echo "[ERROR] Missing PRO-seq window counts"; exit 1; }

Rscript "$(dirname "${BASH_SOURCE[0]}")/04c_compute_pause_indices.R" \\
    --window_counts "\${WC}" \\
    --pause_window_half_bp ${PAUSE_WINDOW_HALF_BP} \\
    --gene_body_flank_bp ${GENE_BODY_FLANK_BP} \\
    --motif_id "\${MOTIF_ID}" \\
    --outdir "${PAUSE_DIR}"
EOC
)")
    submit_array "04c_compute_pause_indices" \
        "--array=${ARRAY_SPEC} --partition=${DEFAULT_PARTITION} \
--time=01:00:00 --mem=8G --cpus-per-task=2 --job-name=ap04c" \
        "${SCRIPT}" "${DEPS}"
    JOB_04C="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# Stage 05 — RNA-seq junction counts
# -----------------------------------------------------------------------------
JOB_05=""
if should_run 05; then
    DEPS=""
    [[ -n "${JOB_01}" ]] && DEPS="--dependency=afterok:${JOB_01}"
    SCRIPT=$(make_stage_wrapper "05" "extract_rnaseq_junctions" "$(cat <<EOC
PI="${PAIR_INDEX_DIR}/\${MOTIF_ID}_pair_index.tsv.gz"
[[ -f "\${PI}" ]] || { echo "[ERROR] Missing pair index"; exit 1; }

ANNOT="${MOTIF_ANNOT_DIR}/\${MOTIF_ID}_annotated.tsv.gz"

bash "$(dirname "${BASH_SOURCE[0]}")/05_extract_rnaseq_junctions.sh" \\
    --pair_index "\${PI}" \\
    --motif_id "\${MOTIF_ID}" \\
    --motif_annot_tsv "\${ANNOT}" \\
    --outdir "${RNASEQ_JUNC_DIR}"
EOC
)")
    submit_array "05_extract_rnaseq_junctions" "${base_sbatch} --job-name=ap05 --time=06:00:00" "${SCRIPT}" "${DEPS}"
    JOB_05="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# Stage 05b — PSI per allele
# -----------------------------------------------------------------------------
JOB_05B=""
if should_run 05b; then
    DEPS=""
    [[ -n "${JOB_05}" ]] && DEPS="--dependency=afterok:${JOB_05}"
    SCRIPT=$(make_stage_wrapper "05b" "compute_psi" "$(cat <<EOC
set +u
source activate "${CONDA_ENV_R}"
set -u

JC="${RNASEQ_JUNC_DIR}/\${MOTIF_ID}_rnaseq_junctions.tsv.gz"
[[ -f "\${JC}" ]] || { echo "[ERROR] Missing junction counts"; exit 1; }

Rscript "$(dirname "${BASH_SOURCE[0]}")/05b_compute_psi_per_allele.R" \\
    --junction_counts "\${JC}" \\
    --psi_min_total_reads ${PSI_MIN_TOTAL_READS} \\
    --motif_id "\${MOTIF_ID}" \\
    --outdir "${PSI_DIR}"
EOC
)")
    submit_array "05b_compute_psi" \
        "--array=${ARRAY_SPEC} --partition=${DEFAULT_PARTITION} \
--time=01:00:00 --mem=8G --cpus-per-task=2 --job-name=ap05b" \
        "${SCRIPT}" "${DEPS}"
    JOB_05B="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# Stage 06 — assemble pair table
# -----------------------------------------------------------------------------
JOB_06=""
if should_run 06; then
    DEPS_LIST=()
    [[ -n "${JOB_02}" ]] && DEPS_LIST+=("${JOB_02}")
    [[ -n "${JOB_03}" ]] && DEPS_LIST+=("${JOB_03}")
    [[ -n "${JOB_04C}" ]] && DEPS_LIST+=("${JOB_04C}")
    [[ -n "${JOB_05B}" ]] && DEPS_LIST+=("${JOB_05B}")
    DEPS=""
    if (( ${#DEPS_LIST[@]} > 0 )); then
        DEPS="--dependency=afterok:$(IFS=:; echo "${DEPS_LIST[*]}")"
    fi

    SCRIPT=$(make_stage_wrapper "06" "assemble_pair_table" "$(cat <<EOC
set +u
source activate "${CONDA_ENV_R}"
set -u

Rscript "$(dirname "${BASH_SOURCE[0]}")/06_assemble_pair_table.R" \\
    --pair_index "${PAIR_INDEX_DIR}/\${MOTIF_ID}_pair_index.tsv.gz" \\
    --pwm_scores "${PWM_DIR}/\${MOTIF_ID}_pwm_scores.tsv.gz" \\
    --atac_counts "${ATAC_COUNTS_DIR}/\${MOTIF_ID}_atac_counts.tsv.gz" \\
    --proseq_window "${PROSEQ_COUNTS_DIR}/\${MOTIF_ID}_proseq_window_counts.tsv.gz" \\
    --pause_indices "${PAUSE_DIR}/\${MOTIF_ID}_pause_indices.tsv.gz" \\
    --psi_wide "${PSI_DIR}/\${MOTIF_ID}_psi_wide.tsv.gz" \\
    --motif_id "\${MOTIF_ID}" \\
    --outdir "${PAIR_TABLE_DIR}"
EOC
)")
    submit_array "06_assemble_pair_table" \
        "--array=${ARRAY_SPEC} --partition=${DEFAULT_PARTITION} \
--time=01:00:00 --mem=8G --cpus-per-task=2 --job-name=ap06" \
        "${SCRIPT}" "${DEPS}"
    JOB_06="${SUBMITTED_JOBID}"
fi

# -----------------------------------------------------------------------------
# Stage 07 — single QC job
# -----------------------------------------------------------------------------
JOB_07=""
if should_run 07; then
    DEPS=""
    [[ -n "${JOB_06}" ]] && DEPS="--dependency=afterok:${JOB_06}"

    QC_WRAPPER="${OUTDIR}/_stage_07.sh"
    cat > "${QC_WRAPPER}" << EOF
#!/usr/bin/env bash
#SBATCH --output=${OUTDIR}/logs/07_%j.out
#SBATCH --error=${OUTDIR}/logs/07_%j.err
set -euo pipefail
set +u
source activate "${CONDA_ENV_R}"
set -u

Rscript "$(dirname "${BASH_SOURCE[0]}")/07_pair_qc.R" \\
    --pair_table_dir "${PAIR_TABLE_DIR}" \\
    --outdir "${OUTDIR}/qc"
EOF
    chmod +x "${QC_WRAPPER}"

    echo ""
    echo "--- 07_pair_qc ---"
    echo "  sbatch --partition=${ASSEMBLY_PARTITION} --time=${ASSEMBLY_TIME} \\"
    echo "         --mem=${ASSEMBLY_MEM} --cpus-per-task=2 --job-name=ap07 \\"
    echo "         ${DEPS} ${QC_WRAPPER}"
    if ${SUBMIT}; then
        out=$(sbatch --partition=${ASSEMBLY_PARTITION} --time=${ASSEMBLY_TIME} \
            --mem=${ASSEMBLY_MEM} --cpus-per-task=2 --job-name=ap07 \
            ${DEPS} --parsable "${QC_WRAPPER}")
        JOB_07="${out%%;*}"
        echo "  -> Submitted as ${JOB_07}"
    else
        JOB_07="DRY_07"
        echo "  (dry run)"
    fi
fi

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "Submission summary"
echo "============================================================"
printf "  %-24s %s\n" "01 build_pair_index:"   "${JOB_01:-(skipped)}"
printf "  %-24s %s\n" "02 score_pwm:"          "${JOB_02:-(skipped)}"
printf "  %-24s %s\n" "03 count_atac:"         "${JOB_03:-(skipped)}"
printf "  %-24s %s\n" "04 count_proseq:"       "${JOB_04:-(skipped)}"
printf "  %-24s %s\n" "04b extract_profile:"   "${JOB_04B:-(skipped)}"
printf "  %-24s %s\n" "04c pause_indices:"     "${JOB_04C:-(skipped)}"
printf "  %-24s %s\n" "05 extract_junctions:"  "${JOB_05:-(skipped)}"
printf "  %-24s %s\n" "05b compute_psi:"       "${JOB_05B:-(skipped)}"
printf "  %-24s %s\n" "06 assemble_pair:"      "${JOB_06:-(skipped)}"
printf "  %-24s %s\n" "07 pair_qc:"            "${JOB_07:-(skipped)}"

echo ""
if ${SUBMIT}; then
    echo "Monitor: squeue -u \$(whoami)"
    echo "Logs:    ${OUTDIR}/logs/"
else
    echo "Dry run complete. Re-run with --apply to actually submit."
fi
