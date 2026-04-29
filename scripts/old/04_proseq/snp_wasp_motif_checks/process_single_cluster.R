#!/usr/bin/env Rscript
# ==============================================================================
# process_single_cluster.R — Process one motif archetype for SNP coverage
# ==============================================================================
# Usage:
#   Rscript process_single_cluster.R <cluster_id>
#   Rscript process_single_cluster.R AC0002
#
# Or via SLURM array (see run_motif_snp_coverage.sh):
#   Uses SLURM_ARRAY_TASK_ID to pick cluster from metadata.
# ==============================================================================

library(data.table)

# ---- PATHS (edit these) ----
ANNOT_DIR          <- "/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/annot_motifs"
ALLELE_COUNTS_FILE <- "/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/proseq/wasp/qc/allele_counts_merged.tsv"
METADATA_FILE      <- "/home/alb1273/pausing_phase_project/resources/metadata.tsv"
OUTPUT_DIR         <- "/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/motif_snp_coverage"

# ---- SPEED FLAG ----
# Set TRUE to ONLY process motifs inside ATAC peaks (much faster).
# Set FALSE to process all intragenic motifs (reports both with_atac and no_atac).
ATAC_ONLY <- TRUE

# Suffix output dir so ATAC-only results don't overwrite full runs
if (ATAC_ONLY) OUTPUT_DIR <- paste0(OUTPUT_DIR, "_atac_only")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Resolve cluster ID ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
    cluster_id <- args[1]
} else {
    # SLURM array mode: pick cluster by task ID
    task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA))
    if (is.na(task_id)) stop("No cluster_id argument and no SLURM_ARRAY_TASK_ID set.")
    meta <- fread(METADATA_FILE, select = "cluster")
    clusters <- sort(unique(meta$cluster))
    # SLURM arrays are 0-indexed in our submission script
    idx <- task_id + 1L
    if (idx < 1 || idx > length(clusters)) {
        stop(sprintf("SLURM_ARRAY_TASK_ID=%d out of range (0-%d)", task_id, length(clusters) - 1))
    }
    cluster_id <- clusters[idx]
}

cat(sprintf("[INFO] Processing cluster: %s\n", cluster_id))

outfile <- file.path(OUTPUT_DIR, paste0(cluster_id, "_snp_coverage.tsv"))
if (file.exists(outfile)) {
    cat(sprintf("[INFO] Output already exists: %s. Skipping.\n", outfile))
    quit(status = 0)
}

# ---- Load allele counts and build keyed lookup ----
cat("[INFO] Loading allele counts...\n")
ac <- fread(ALLELE_COUNTS_FILE)
ac[, total_allelic := ref_total + alt_total]

# Identify per-rep ref/alt columns (e.g. WT_PROseq_rep1_ref, WT_PROseq_rep1_alt)
rep_cols <- grep("_(ref|alt)$", colnames(ac), value = TRUE)
# Exclude the merged totals
rep_cols <- rep_cols[!rep_cols %in% c("ref_total", "alt_total")]
cat(sprintf("[INFO] Per-rep allele columns: %s\n", paste(rep_cols, collapse = ", ")))

# Compute min(ref, alt, across all reps) per SNP
# This is the most conservative metric: the weakest allele in the weakest rep
if (length(rep_cols) > 0) {
    ac[, min_allele_any_rep := do.call(pmin, .SD), .SDcols = rep_cols]
} else {
    # Fallback if no per-rep columns: min of merged ref/alt
    ac[, min_allele_any_rep := pmin(ref_total, alt_total)]
}

ac_lookup <- ac[, .(chrom, pos, total_allelic, min_allele_any_rep)]
setkey(ac_lookup, chrom, pos)

cat(sprintf("[INFO] Allele counts: %s SNPs\n", format(nrow(ac_lookup), big.mark = ",")))
cat(sprintf("[INFO] Depth metrics: 'merged_total' (ref+alt merged) and 'min_allele_rep' (min across all rep×allele)\n"))

# ---- Helper: parse SNP detail string → max depth for BOTH metrics ----
# Returns a named list: list(merged_total = X, min_allele_rep = Y)
# "max" across SNPs in the string, separately for each metric.
get_max_depths_both <- function(snp_str) {
    if (is.na(snp_str) || snp_str == "") return(list(merged_total = 0L, min_allele_rep = 0L))
    parts <- strsplit(snp_str, "\\|")[[1]]
    max_mt <- 0L
    max_mar <- 0L
    for (p in parts) {
        fields <- strsplit(p, ":")[[1]]
        if (length(fields) >= 2) {
            chr <- fields[1]
            pos <- as.integer(fields[2])
            row <- ac_lookup[.(chr, pos), nomatch = 0L]
            if (nrow(row) > 0) {
                if (row$total_allelic[1] > max_mt)       max_mt  <- row$total_allelic[1]
                if (row$min_allele_any_rep[1] > max_mar)  max_mar <- row$min_allele_any_rep[1]
            }
        }
    }
    return(list(merged_total = max_mt, min_allele_rep = max_mar))
}

# Vectorized: returns a data.table with two columns
get_max_depths_vec <- function(snp_col) {
    out <- lapply(snp_col, get_max_depths_both)
    data.table(
        merged_total   = vapply(out, `[[`, integer(1), "merged_total"),
        min_allele_rep = vapply(out, `[[`, integer(1), "min_allele_rep")
    )
}

# ---- Distance columns and depth thresholds ----
DIST_COLS <- list(
    "0bp"   = list(overlap = "snp_F121-9_overlap",           details = "snp_F121-9_details"),
    "10bp"  = list(overlap = "snp_F121-9_flank10bp_overlap",  details = "snp_F121-9_flank10bp_details"),
    "25bp"  = list(overlap = "snp_F121-9_flank25bp_overlap",  details = "snp_F121-9_flank25bp_details"),
    "50bp"  = list(overlap = "snp_F121-9_flank50bp_overlap",  details = "snp_F121-9_flank50bp_details"),
    "100bp" = list(overlap = "snp_F121-9_flank100bp_overlap", details = "snp_F121-9_flank100bp_details"),
    "200bp" = list(overlap = "snp_F121-9_flank200bp_overlap", details = "snp_F121-9_flank200bp_details"),
    "400bp" = list(overlap = "snp_F121-9_flank400bp_overlap", details = "snp_F121-9_flank400bp_details")
)

MIN_DEPTHS <- c(1, 5, 10, 20, 50)

# ---- Load annotated motif file ----
fpath <- file.path(ANNOT_DIR, paste0(cluster_id, "_annotated.tsv.gz"))
if (!file.exists(fpath)) {
    fpath <- file.path(ANNOT_DIR, paste0(cluster_id, "_annotated.tsv"))
}
if (!file.exists(fpath)) {
    cat(sprintf("[WARN] Annotated file not found for %s. Writing empty output.\n", cluster_id))
    empty <- data.table(
        cluster = character(), total_motifs = integer(),
        n_gene_filter = integer(), atac_filter = character(),
        n_base = integer(), distance = character(),
        min_depth = integer(), n_with_snp = integer()
    )
    fwrite(empty, outfile, sep = "\t")
    quit(status = 0)
}

cat(sprintf("[INFO] Reading %s\n", fpath))
dt <- fread(fpath, showProgress = FALSE)
total_motifs <- nrow(dt)
cat(sprintf("[INFO] Total motifs: %d\n", total_motifs))

# ---- Base filters ----
dt[, passes_gene_filter := (
    is_intragenic == TRUE &
    is_promoter == FALSE &
    !is.na(intragenic_gene_types) &
    grepl("protein_coding|lncRNA", intragenic_gene_types)
)]

atac_col1 <- "atac_4DNFIAEQI3RP_overlap"
atac_col2 <- "atac_4DNFIZNPOOZN_overlap"
if (atac_col1 %in% names(dt) & atac_col2 %in% names(dt)) {
    dt[, passes_atac := (get(atac_col1) == TRUE | get(atac_col2) == TRUE)]
} else {
    dt[, passes_atac := FALSE]
}

n_gene_filter      <- sum(dt$passes_gene_filter)
n_gene_atac_filter <- sum(dt$passes_gene_filter & dt$passes_atac)
cat(sprintf("[INFO] Pass gene filter: %d | + ATAC: %d\n", n_gene_filter, n_gene_atac_filter))

# ---- Determine which ATAC levels to report ----
atac_levels <- if (ATAC_ONLY) "with_atac" else c("no_atac", "with_atac")

DEPTH_METRICS <- c("merged_total", "min_allele_rep")

# ---- Apply early filter ----
if (ATAC_ONLY) {
    dt_filt <- dt[passes_gene_filter == TRUE & passes_atac == TRUE]
    cat(sprintf("[INFO] ATAC_ONLY mode: %d motifs after early ATAC filter\n", nrow(dt_filt)))
} else {
    dt_filt <- dt[passes_gene_filter == TRUE]
}

# ---- Process ----
results_list <- list()

if (nrow(dt_filt) == 0) {
    for (dist_name in names(DIST_COLS)) {
        for (min_d in MIN_DEPTHS) {
            for (atac_req in atac_levels) {
                for (dm in DEPTH_METRICS) {
                    results_list[[length(results_list) + 1]] <- data.table(
                        cluster = cluster_id, total_motifs = total_motifs,
                        n_gene_filter = n_gene_filter, atac_filter = atac_req,
                        depth_metric = dm,
                        n_base = 0L, distance = dist_name,
                        min_depth = min_d, n_with_snp = 0L
                    )
                }
            }
        }
    }
} else {
    # Compute max depth for each distance window (both metrics)
    for (dist_name in names(DIST_COLS)) {
        det_col <- DIST_COLS[[dist_name]]$details
        ovl_col <- DIST_COLS[[dist_name]]$overlap
        mt_col  <- paste0("mt_", dist_name)   # merged_total
        mar_col <- paste0("mar_", dist_name)   # min_allele_rep

        dt_filt[, (mt_col)  := 0L]
        dt_filt[, (mar_col) := 0L]

        if (det_col %in% names(dt_filt)) {
            has_overlap <- which(dt_filt[[ovl_col]] == TRUE)
            if (length(has_overlap) > 0) {
                both <- get_max_depths_vec(dt_filt[[det_col]][has_overlap])
                dt_filt[has_overlap, (mt_col)  := both$merged_total]
                dt_filt[has_overlap, (mar_col) := both$min_allele_rep]
            }
        }
    }

    # Count qualifying motifs for each combination
    for (dist_name in names(DIST_COLS)) {
        mt_col  <- paste0("mt_", dist_name)
        mar_col <- paste0("mar_", dist_name)

        for (min_d in MIN_DEPTHS) {
            # Precompute pass vectors for each metric
            pass_mt  <- dt_filt[[mt_col]]  >= min_d
            pass_mar <- dt_filt[[mar_col]] >= min_d

            for (dm in DEPTH_METRICS) {
                has_snp <- if (dm == "merged_total") pass_mt else pass_mar

                if (ATAC_ONLY) {
                    results_list[[length(results_list) + 1]] <- data.table(
                        cluster = cluster_id, total_motifs = total_motifs,
                        n_gene_filter = n_gene_filter, atac_filter = "with_atac",
                        depth_metric = dm,
                        n_base = nrow(dt_filt), distance = dist_name,
                        min_depth = min_d, n_with_snp = sum(has_snp)
                    )
                } else {
                    results_list[[length(results_list) + 1]] <- data.table(
                        cluster = cluster_id, total_motifs = total_motifs,
                        n_gene_filter = n_gene_filter, atac_filter = "no_atac",
                        depth_metric = dm,
                        n_base = nrow(dt_filt), distance = dist_name,
                        min_depth = min_d, n_with_snp = sum(has_snp)
                    )
                    results_list[[length(results_list) + 1]] <- data.table(
                        cluster = cluster_id, total_motifs = total_motifs,
                        n_gene_filter = n_gene_filter, atac_filter = "with_atac",
                        depth_metric = dm,
                        n_base = n_gene_atac_filter, distance = dist_name,
                        min_depth = min_d, n_with_snp = sum(has_snp & dt_filt$passes_atac)
                    )
                }
            }
        }
    }
}

# ---- Save ----
out <- rbindlist(results_list)
fwrite(out, outfile, sep = "\t")
cat(sprintf("[DONE] %s → %s (%d rows)\n", cluster_id, outfile, nrow(out)))
