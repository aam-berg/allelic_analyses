#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# 07_diagnostics/run_pipeline.sh — SLURM orchestrator
# =============================================================================
#
# Order of stages:
#   00      validate inputs
#   01      per-motif summary (SLURM array, fanned out)
#   02      aggregate per-motif summaries -> archetype_summary.tsv
#   03-06   plotting from aggregated tables (single jobs)
#   07      SNP-distance ATAC diagnostic (single)
#   08      WASP cascade + per-SNP coverage compute (single, heavy)
#   09      WASP plots (single, depends on 08)
#   10      powered-archetype plots (single, depends on 06 pair tables)
#
# Usage:
#   bash run_pipeline.sh                  # dry run
#   bash run_pipeline.sh --apply          # actually submit
#   bash run_pipeline.sh --apply --start-from 03   # skip earlier stages
#   STRICT=1 bash run_pipeline.sh --apply  # require all motifs present
#
# Stage selection:
#   --skip 08,09     Comma-separated list of stages to skip.
#   --start-from N   Start from stage N (skip everything earlier).
# =============================================================================

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

APPLY=0
START_FROM=00
SKIP=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --apply) APPLY=1; shift ;;
        --start-from) START_FROM="$2"; shift 2 ;;
        --skip) SKIP="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

stage_skipped() {
    local stage="$1"
    [[ ",${SKIP}," == *",${stage},"* ]]
}

stage_too_early() {
    local stage="$1"
    [[ "${stage}" < "${START_FROM}" ]]
}

submit_or_dry() {
    local label="$1"; shift
    if (( APPLY == 1 )); then
        echo "  >>> ${label}: $*"
        eval "$@"
    else
        echo "  [DRY] ${label}: $*"
    fi
}

create_dirs

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LIB="${SCRIPT_DIR}/lib_diagnostics.R"

# -----------------------------------------------------------------------------
# Stage 00: Setup / validate inputs
# -----------------------------------------------------------------------------
if ! stage_too_early "00" && ! stage_skipped "00"; then
    step_header "Stage 00 — input validation"
    if (( APPLY == 1 )); then
        bash "${SCRIPT_DIR}/00_setup.sh"
    else
        echo "  [DRY] bash 00_setup.sh"
    fi
fi

# -----------------------------------------------------------------------------
# Stage 01: per-motif summary (SLURM array)
# -----------------------------------------------------------------------------
if ! stage_too_early "01" && ! stage_skipped "01"; then
    step_header "Stage 01 — per-motif summary"

    # Build motif list from 05_motif_annot outputs
    MOTIF_LIST="${LOGS_DIR}/motif_list.txt"
    if (( APPLY == 1 )); then
        find "${MOTIF_ANNOT_DIR}" -maxdepth 1 -name "*_annotated.tsv.gz" \
            -printf "%f\n" | sed 's/_annotated\.tsv\.gz$//' | sort > "${MOTIF_LIST}"
        N_MOTIFS=$(wc -l < "${MOTIF_LIST}")
        echo "  Built motif list (${N_MOTIFS} entries): ${MOTIF_LIST}"
    else
        N_MOTIFS=10
        echo "  [DRY] would build motif list at ${MOTIF_LIST}"
    fi

    # SLURM array script
    ARRAY_SH="${LOGS_DIR}/01_array_task.sh"
    cat > "${ARRAY_SH}" << EOF
#!/bin/bash
#SBATCH --partition=${DEFAULT_PARTITION}
#SBATCH --time=${DEFAULT_TIME}
#SBATCH --mem=${DEFAULT_MEM}
#SBATCH --cpus-per-task=${DEFAULT_THREADS}
#SBATCH --output=${LOGS_DIR}/01_%A_%a.out
#SBATCH --error=${LOGS_DIR}/01_%A_%a.err

set -euo pipefail
source "${SCRIPT_DIR}/config.sh"
source activate ${CONDA_ENV_R}

MOTIF_ID=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${MOTIF_LIST}")
ANNOT_TSV="${MOTIF_ANNOT_DIR}/\${MOTIF_ID}_annotated.tsv.gz"

# Get motif width by reading PWM file (do once via R inside the script)
MOTIF_WIDTH=\$(Rscript -e "
source('${LIB}')
w <- get_archetype_pwm_widths('${PWM_MEME_FILE}')
mw <- w[motif_id == '\${MOTIF_ID}', motif_width]
cat(if (length(mw) == 0) NA else mw)
")

Rscript "${SCRIPT_DIR}/01_build_summary.R" \\
    --annot_tsv "\${ANNOT_TSV}" \\
    --motif_id "\${MOTIF_ID}" \\
    --motif_width "\${MOTIF_WIDTH}" \\
    --apply_expression_filter TRUE \\
    --expression_tpm_min ${EXPRESSION_TPM_MIN} \\
    --standard_chroms "${STANDARD_CHROMS}" \\
    --gene_types_include "${GENE_TYPES_INCLUDE}" \\
    --hybrid "${HYBRID}" \\
    --outdir "${SUMMARY_DIR}" \\
    --lib_path "${LIB}"
EOF
    chmod +x "${ARRAY_SH}"

    if (( APPLY == 1 )); then
        JID01=$(sbatch --parsable \
            --array=1-${N_MOTIFS}%${DEFAULT_MAX_CONCURRENT} \
            --job-name=07d_01 \
            "${ARRAY_SH}")
        echo "  Submitted 01: jobid=${JID01}"
        DEP_01="--dependency=afterok:${JID01}"
    else
        echo "  [DRY] sbatch --array=1-N --parsable ${ARRAY_SH}"
        DEP_01=""
    fi
else
    DEP_01=""
fi

# -----------------------------------------------------------------------------
# Helper: submit a single-job R script
# -----------------------------------------------------------------------------
submit_r_step() {
    local stage="$1" name="$2" script="$3" mem="$4" time_lim="$5" deps="$6"
    shift 6
    local args="$*"

    local job_sh="${LOGS_DIR}/${stage}_${name}.sh"
    cat > "${job_sh}" << EOF
#!/bin/bash
#SBATCH --partition=${DEFAULT_PARTITION}
#SBATCH --time=${time_lim}
#SBATCH --mem=${mem}
#SBATCH --cpus-per-task=${DEFAULT_THREADS}
#SBATCH --output=${LOGS_DIR}/${stage}_${name}_%j.out
#SBATCH --error=${LOGS_DIR}/${stage}_${name}_%j.err

set -euo pipefail
source "${SCRIPT_DIR}/config.sh"
source activate ${CONDA_ENV_R}

Rscript "${script}" ${args} --lib_path "${LIB}"
EOF
    chmod +x "${job_sh}"

    if (( APPLY == 1 )); then
        local jid
        jid=$(sbatch --parsable ${deps} --job-name="07d_${stage}" "${job_sh}")
        echo "  Submitted ${stage} (${name}): jobid=${jid}"
        echo "${jid}"
    else
        echo "  [DRY] sbatch ${deps} ${job_sh}"
        echo ""
    fi
}

submit_sh_step() {
    local stage="$1" name="$2" script="$3" mem="$4" time_lim="$5" deps="$6"

    local job_sh="${LOGS_DIR}/${stage}_${name}.sh"
    cat > "${job_sh}" << EOF
#!/bin/bash
#SBATCH --partition=${HEAVY_PARTITION}
#SBATCH --time=${time_lim}
#SBATCH --mem=${mem}
#SBATCH --cpus-per-task=4
#SBATCH --output=${LOGS_DIR}/${stage}_${name}_%j.out
#SBATCH --error=${LOGS_DIR}/${stage}_${name}_%j.err

set -euo pipefail
source activate ${CONDA_ENV_TOOLS}
bash "${script}"
EOF
    chmod +x "${job_sh}"

    if (( APPLY == 1 )); then
        local jid
        jid=$(sbatch --parsable ${deps} --job-name="07d_${stage}" "${job_sh}")
        echo "  Submitted ${stage} (${name}): jobid=${jid}"
        echo "${jid}"
    else
        echo "  [DRY] sbatch ${deps} ${job_sh}"
        echo ""
    fi
}

# -----------------------------------------------------------------------------
# Stage 02: aggregate
# -----------------------------------------------------------------------------
if ! stage_too_early "02" && ! stage_skipped "02"; then
    step_header "Stage 02 — aggregate summaries"
    JID02=$(submit_r_step "02" "aggregate" \
        "${SCRIPT_DIR}/02_aggregate_summary.R" \
        "${DEFAULT_MEM}" "${DEFAULT_TIME}" \
        "${DEP_01}" \
        --summary_dir "${SUMMARY_DIR}" \
        --pwm_meme "${PWM_MEME_FILE}" \
        --metadata "${TF_METADATA_FILE}" \
        --expression "${GENE_EXPRESSION_FILE}" \
        --outdir "${TABLES_DIR}")
    DEP_02=$([[ -n "${JID02}" ]] && echo "--dependency=afterok:${JID02}" || echo "")
else
    DEP_02=""
fi

# -----------------------------------------------------------------------------
# Stages 03-06: aggregated-table plots
# -----------------------------------------------------------------------------
if ! stage_too_early "03" && ! stage_skipped "03"; then
    step_header "Stage 03 — filter cascade plots"
    submit_r_step "03" "cascade" \
        "${SCRIPT_DIR}/03_plot_filter_cascade.R" \
        "${DEFAULT_MEM}" "${DEFAULT_TIME}" \
        "${DEP_02}" \
        --archetype_summary "${TABLES_DIR}/archetype_summary.tsv" \
        --plots_dir "${PLOTS_DIR}" \
        --min_motif_width "${MIN_MOTIF_WIDTH}" >/dev/null
fi

if ! stage_too_early "04" && ! stage_skipped "04"; then
    step_header "Stage 04 — qualifying histograms"
    submit_r_step "04" "qualifying" \
        "${SCRIPT_DIR}/04_plot_qualifying_histogram.R" \
        "${DEFAULT_MEM}" "${DEFAULT_TIME}" \
        "${DEP_02}" \
        --qualifying_table "${TABLES_DIR}/qualifying_per_archetype.tsv" \
        --plots_dir "${PLOTS_DIR}" \
        --min_motif_width "${MIN_MOTIF_WIDTH}" >/dev/null
fi

if ! stage_too_early "05" && ! stage_skipped "05"; then
    step_header "Stage 05 — top archetypes panel"
    submit_r_step "05" "top_archetypes" \
        "${SCRIPT_DIR}/05_plot_top_archetypes.R" \
        "${DEFAULT_MEM}" "${DEFAULT_TIME}" \
        "${DEP_02}" \
        --qualifying_table "${TABLES_DIR}/qualifying_per_archetype.tsv" \
        --pwm_meme "${PWM_MEME_FILE}" \
        --plots_dir "${PLOTS_DIR}" \
        --top_n_list "${TOP_N_LIST}" \
        --min_motif_width "${MIN_MOTIF_WIDTH}" >/dev/null
fi

if ! stage_too_early "06" && ! stage_skipped "06"; then
    step_header "Stage 06 — PWM score distributions"
    submit_r_step "06" "pwm_dists" \
        "${SCRIPT_DIR}/06_plot_pwm_score_dists.R" \
        "${DEFAULT_MEM}" "${DEFAULT_TIME}" \
        "${DEP_02}" \
        --summary_dir "${SUMMARY_DIR}" \
        --qualifying_table "${TABLES_DIR}/qualifying_per_archetype.tsv" \
        --plots_dir "${PLOTS_DIR}" \
        --top_n 25 \
        --min_motif_width "${MIN_MOTIF_WIDTH}" >/dev/null
fi

# -----------------------------------------------------------------------------
# Stage 07: SNP-distance diagnostic
# -----------------------------------------------------------------------------
if ! stage_too_early "07" && ! stage_skipped "07"; then
    step_header "Stage 07 — SNP-distance ATAC diagnostic"
    submit_r_step "07" "snp_dist" \
        "${SCRIPT_DIR}/07_plot_snp_distance_diag.R" \
        "${HEAVY_MEM}" "${DEFAULT_TIME}" \
        "" \
        --pair_table_dir "${PAIR_TABLE_DIR}" \
        --plots_dir "${PLOTS_DIR}" \
        --tables_dir "${TABLES_DIR}" \
        --bin_edges "${ATAC_DIST_BIN_EDGES}" \
        --n_permutations "${ATAC_DIAGNOSTIC_N_PERMUTATIONS}" \
        --min_motif_width "${MIN_MOTIF_WIDTH}" \
        --pwm_meme "${PWM_MEME_FILE}" >/dev/null
fi

# -----------------------------------------------------------------------------
# Stage 08: WASP cascade + per-SNP coverage compute (heavy)
# -----------------------------------------------------------------------------
if ! stage_too_early "08" && ! stage_skipped "08"; then
    step_header "Stage 08 — WASP balance + coverage compute"
    JID08=$(submit_sh_step "08" "wasp" \
        "${SCRIPT_DIR}/08_compute_wasp_balance.sh" \
        "${HEAVY_MEM}" "${HEAVY_TIME}" "")
    DEP_08=$([[ -n "${JID08}" ]] && echo "--dependency=afterok:${JID08}" || echo "")
else
    DEP_08=""
fi

# -----------------------------------------------------------------------------
# Stage 09: WASP plots
# -----------------------------------------------------------------------------
if ! stage_too_early "09" && ! stage_skipped "09"; then
    step_header "Stage 09 — WASP/balance plots"
    submit_r_step "09" "wasp_plot" \
        "${SCRIPT_DIR}/09_plot_wasp_balance.R" \
        "${DEFAULT_MEM}" "${DEFAULT_TIME}" \
        "${DEP_08}" \
        --wasp_dir "${WASP_DIR}" \
        --plots_dir "${PLOTS_DIR}" >/dev/null
fi

# -----------------------------------------------------------------------------
# Stage 10: powered archetypes
# -----------------------------------------------------------------------------
if ! stage_too_early "10" && ! stage_skipped "10"; then
    step_header "Stage 10 — powered archetypes"
    submit_r_step "10" "powered" \
        "${SCRIPT_DIR}/10_plot_powered_archetypes.R" \
        "${HEAVY_MEM}" "${DEFAULT_TIME}" \
        "${DEP_02}" \
        --pair_table_dir "${PAIR_TABLE_DIR}" \
        --qualifying_table "${TABLES_DIR}/qualifying_per_archetype.tsv" \
        --plots_dir "${PLOTS_DIR}" \
        --tables_dir "${TABLES_DIR}" \
        --power_atac_min "${POWER_ATAC_MIN_PER_ALLELE}" \
        --power_proseq_pause_min "${POWER_PROSEQ_PAUSE_MIN_PER_ALLELE}" \
        --power_rnaseq_min "${POWER_RNASEQ_AT_JUNCTION_MIN_PER_ALLELE}" \
        --min_motif_width "${MIN_MOTIF_WIDTH}" \
        --pwm_meme "${PWM_MEME_FILE}" >/dev/null
fi

echo ""
echo "============================================================"
if (( APPLY == 1 )); then
    echo "All stages submitted. Monitor with: squeue -u \$USER -n '07d_*'"
else
    echo "Dry run complete. Run with --apply to actually submit."
fi
echo "============================================================"
