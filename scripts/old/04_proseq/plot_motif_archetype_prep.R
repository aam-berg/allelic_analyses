#!/usr/bin/env Rscript
# =============================================================================
# plot_motif_archetype_prep.R — Preparatory analysis for TFBS × PRO-seq study
# =============================================================================
#
# Generates summary plots characterizing the annotated motif archetypes,
# focusing on the subset relevant for allele-specific PRO-seq analysis:
#   - Non-promoter, intragenic (protein-coding/lncRNA), chromatin-accessible
#     TFBSs with het SNPs at various proximity thresholds
#
# KEY FEATURES:
#   - Parses MEME-format PWM file to extract motif width and information content
#   - Filters archetypes by minimum width (default: 8 bp) and optionally by IC
#   - Cumulative SNP proximity analysis: overlap, ≤10bp, ≤25bp, ≤50bp, ≤100bp
#   - Multiple plot versions at each SNP proximity threshold
#
# USAGE:
#   Rscript plot_motif_archetype_prep.R \
#       --annot_dir /n/scratch/users/a/alb1273/pausing_phase_project_intermediates/annot_motifs \
#       --metadata /home/alb1273/pausing_phase_project/resources/metadata.tsv \
#       --meme_file /home/alb1273/pausing_phase_project/resources/consensus_pwms.meme \
#       --min_width 8 \
#       --min_ic 0 \
#       --outdir /n/scratch/users/a/alb1273/pausing_phase_project_intermediates/annot_motifs/overview_plots
#
# DEPENDENCIES:
#   R packages: data.table, ggplot2, patchwork, scales, viridis
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(ggplot2)
    library(patchwork)
    library(scales)
    library(viridis)
})

# =============================================================================
# 0. Parse arguments
# =============================================================================

option_list <- list(
    make_option("--annot_dir", type = "character",
                default = "/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/annot_motifs",
                help = "Directory with annotated motif files (*_annotated.tsv.gz)"),
    make_option("--metadata", type = "character",
                default = "/home/alb1273/pausing_phase_project/resources/metadata.tsv",
                help = "Motif metadata file (motif_id, cluster, tf_name, family_name)"),
    make_option("--meme_file", type = "character",
                default = "/home/alb1273/pausing_phase_project/resources/consensus_pwms.meme",
                help = "MEME-format PWM file for archetype widths and IC"),
    make_option("--min_width", type = "integer", default = 8,
                help = "Minimum motif width (bp) to include archetype [default: %default]"),
    make_option("--min_ic", type = "double", default = 0,
                help = "Minimum total information content (bits) [default: %default]"),
    make_option("--outdir", type = "character",
                default = "/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/annot_motifs/prep_plots",
                help = "Output directory"),
    make_option("--atac_cols", type = "character",
                default = "atac_4DNFIAEQI3RP_overlap,atac_4DNFIZNPOOZN_overlap",
                help = "Comma-separated ATAC overlap column names"),
    make_option("--snp_col", type = "character",
                default = "snp_F121-9_overlap",
                help = "SNP overlap column name"),
    make_option("--snp_detail_col", type = "character",
                default = "snp_F121-9_details",
                help = "SNP detail column name"),
    make_option("--snp_count_col", type = "character",
                default = "snp_F121-9_count",
                help = "SNP count column name"),
    make_option("--snp_prefix", type = "character",
                default = "snp_F121-9",
                help = "SNP column prefix for flanking columns")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
atac_cols <- trimws(strsplit(opt$atac_cols, ",")[[1]])

cat("================================================================\n")
cat("Motif Archetype Preparatory Analysis\n")
cat("================================================================\n\n")
cat("Parameters:\n")
cat("  annot_dir:  ", opt$annot_dir, "\n")
cat("  metadata:   ", opt$metadata, "\n")
cat("  meme_file:  ", opt$meme_file, "\n")
cat("  min_width:  ", opt$min_width, " bp\n")
cat("  min_ic:     ", opt$min_ic, " bits\n")
cat("  outdir:     ", opt$outdir, "\n\n")

# =============================================================================
# === CUSTOMIZE === Global theme and layout settings
# =============================================================================

BASE_SIZE       <- 13
TITLE_SIZE      <- 16
SUBTITLE_SIZE   <- 13
AXIS_TEXT_SIZE  <- 13
AXIS_TITLE_SIZE <- 16
LEGEND_SIZE     <- 12
BAR_LABEL_SIZE  <- 3.8
FUNNEL_LABEL_SIZE <- 4.2

FIG_WIDTH  <- 14
FIG_HEIGHT <- 11
TOP_N_ARCHETYPES <- 40

# Cumulative SNP proximity thresholds (must be subset of flank distances)
SNP_CUMUL_DISTS <- c(100, 50, 25, 10, 0)  # 0 = direct overlap
SNP_CUMUL_LABELS <- c(
    "0"   = "Overlaps motif",
    "10"  = "\u226410 bp",
    "25"  = "\u226425 bp",
    "50"  = "\u226450 bp",
    "100" = "\u2264100 bp"
)
SNP_CUMUL_COLORS <- c(
    "Overlaps motif" = "#d62728",
    "\u226410 bp"    = "#e6550d",
    "\u226425 bp"    = "#fd8d3c",
    "\u226450 bp"    = "#fdae6b",
    "\u2264100 bp"   = "#fdd0a2"
)

# Colors for filter categories in funnel plots
FILTER_COLORS <- c(
    "All TFBSs"                       = "#bdbdbd",
    "Intragenic"                      = "#a1d99b",
    "Protein-coding / lncRNA"         = "#74c476",
    "Non-promoter"                    = "#238b45",
    "Accessible (ATAC)"               = "#fd8d3c",
    "Het SNP overlaps motif"          = "#d62728"
)

THEME <- theme_bw(base_size = BASE_SIZE) +
    theme(
        plot.title    = element_text(size = TITLE_SIZE, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = SUBTITLE_SIZE, color = "grey40"),
        axis.text     = element_text(size = AXIS_TEXT_SIZE),
        axis.title    = element_text(size = AXIS_TITLE_SIZE),
        legend.text   = element_text(size = LEGEND_SIZE),
        legend.title  = element_text(size = LEGEND_SIZE, face = "bold"),
        strip.text    = element_text(size = BASE_SIZE),
        panel.grid.minor = element_blank()
    )

# =============================================================================
# 1. Load metadata
# =============================================================================

cat("[1] Loading metadata...\n")
meta <- fread(opt$metadata)
cat("  ", nrow(meta), " motif entries\n")

arch_family <- unique(meta[, .(cluster, family_name)])
arch_family_map <- arch_family[, .(
    families = paste(unique(family_name), collapse = "|")
), by = cluster]

arch_tf_map <- meta[, .(
    tf_names = paste(unique(tf_name), collapse = "|")
), by = cluster]

# =============================================================================
# 1b. Parse MEME file
# =============================================================================

cat("[1b] Parsing MEME file for motif widths and information content...\n")

if (!file.exists(opt$meme_file)) {
    stop("MEME file not found: ", opt$meme_file)
}

parse_meme <- function(meme_path) {
    lines <- readLines(meme_path)
    results <- list()
    i <- 1
    while (i <= length(lines)) {
        line <- trimws(lines[i])
        if (startsWith(line, "MOTIF ")) {
            tokens <- strsplit(sub("^MOTIF\\s+", "", line), ":")[[1]]
            arch_id <- tokens[1]
            i <- i + 1
            while (i <= length(lines) && !grepl("^letter-probability matrix", lines[i])) {
                i <- i + 1
            }
            if (i > length(lines)) break
            header <- lines[i]
            width <- as.integer(sub(".*w=\\s*(\\d+).*", "\\1", header))
            pwm <- matrix(0, nrow = width, ncol = 4)
            row_idx <- 0
            i <- i + 1
            while (row_idx < width && i <= length(lines)) {
                row_line <- trimws(lines[i])
                if (nchar(row_line) > 0 && !startsWith(row_line, "MOTIF") &&
                    !startsWith(row_line, "letter")) {
                    vals <- as.numeric(strsplit(row_line, "\\s+")[[1]])
                    if (length(vals) == 4) {
                        row_idx <- row_idx + 1
                        pwm[row_idx, ] <- vals
                    }
                }
                i <- i + 1
            }
            ic_per_pos <- apply(pwm, 1, function(p) {
                p <- pmax(p, 1e-10)
                2 + sum(p * log2(p))
            })
            total_ic <- sum(ic_per_pos)
            results[[arch_id]] <- data.table(
                archetype       = arch_id,
                motif_width     = as.integer(width),
                total_ic        = round(total_ic, 3),
                mean_ic_per_pos = round(total_ic / width, 3)
            )
        } else {
            i <- i + 1
        }
    }
    rbindlist(results)
}

pwm_all <- parse_meme(opt$meme_file)
cat("  Parsed ", nrow(pwm_all), " archetypes from MEME file\n")
cat("  Width range: ", min(pwm_all$motif_width), "–",
    max(pwm_all$motif_width), " bp\n")

pwm_filtered <- copy(pwm_all)
n_before <- nrow(pwm_filtered)
if (opt$min_width > 0) {
    pwm_filtered <- pwm_filtered[motif_width >= opt$min_width]
    cat("  After min_width >= ", opt$min_width, ": ",
        nrow(pwm_filtered), " / ", n_before, " archetypes\n")
}
if (opt$min_ic > 0) {
    n_before_ic <- nrow(pwm_filtered)
    pwm_filtered <- pwm_filtered[total_ic >= opt$min_ic]
    cat("  After min_ic >= ", opt$min_ic, ": ",
        nrow(pwm_filtered), " / ", n_before_ic, " archetypes\n")
}
VALID_ARCHETYPES <- pwm_filtered$archetype
cat("  Archetypes passing PWM filters: ", length(VALID_ARCHETYPES), "\n\n")

# =============================================================================
# 2. Load and summarize all annotated motif files
# =============================================================================

cat("[2] Loading annotated motif files...\n")

annot_files <- list.files(opt$annot_dir, pattern = "_annotated\\.tsv\\.gz$",
                           full.names = TRUE)
cat("  Found ", length(annot_files), " annotated archetype files\n")
if (length(annot_files) == 0) stop("No annotated files found in ", opt$annot_dir)

annot_files <- annot_files[
    sub("_annotated\\.tsv\\.gz$", "", basename(annot_files)) %in% VALID_ARCHETYPES
]
cat("  After PWM filter: ", length(annot_files), " files\n")

# Build flanking column names needed for cumulative SNP analysis
flank_overlap_cols <- paste0(opt$snp_prefix, "_flank",
                              c(10, 25, 50, 100), "bp_overlap")
flank_count_cols <- paste0(opt$snp_prefix, "_flank",
                            c(10, 25, 50, 100), "bp_count")

# Columns to load from each file
load_cols <- c("chrom", "start", "end", "motif_id", "score", "strand",
               "is_intragenic", "is_promoter", "intragenic_gene_types",
               atac_cols, opt$snp_col, opt$snp_count_col,
               flank_overlap_cols, flank_count_cols)

archetype_summary <- list()
score_data <- list()
snp_count_data <- list()
all_tfbs_scores <- list()

# Also collect cumulative SNP count data for each threshold
snp_count_cumul_data <- list()

for (f in annot_files) {
    arch_id <- sub("_annotated\\.tsv\\.gz$", "", basename(f))
    if ((which(annot_files == f) %% 50) == 1) {
        cat("  Processing ", arch_id, " (",
            which(annot_files == f), "/", length(annot_files), ")\n")
    }

    dt <- fread(f, select = load_cols)

    n_total <- nrow(dt)
    if (n_total == 0) {
        archetype_summary[[arch_id]] <- data.table(
            archetype = arch_id, n_total = 0, n_intragenic = 0,
            n_intragenic_pcl = 0, n_intragenic_nonpromoter = 0,
            n_accessible = 0,
            n_snp_overlap = 0, n_snp_10bp = 0, n_snp_25bp = 0,
            n_snp_50bp = 0, n_snp_100bp = 0,
            n_noatac_snp_overlap = 0, n_noatac_snp_10bp = 0,
            n_noatac_snp_25bp = 0, n_noatac_snp_50bp = 0,
            n_noatac_snp_100bp = 0
        )
        next
    }

    # --- Progressive filters ---
    dt_intra <- dt[is_intragenic == TRUE]
    n_intragenic <- nrow(dt_intra)

    # NEW: protein-coding / lncRNA filter
    # A motif passes if intragenic_gene_types contains "protein_coding" or "lncRNA"
    dt_pcl <- dt_intra[
        grepl("protein_coding", intragenic_gene_types, fixed = TRUE) |
        grepl("lncRNA", intragenic_gene_types, fixed = TRUE)
    ]
    n_intragenic_pcl <- nrow(dt_pcl)

    dt_nonprom <- dt_pcl[is_promoter == FALSE]
    n_nonprom <- nrow(dt_nonprom)

    # Accessible
    if (length(atac_cols) > 0 && all(atac_cols %in% names(dt_nonprom))) {
        dt_nonprom[, accessible := Reduce(`|`, .SD), .SDcols = atac_cols]
    } else {
        dt_nonprom[, accessible := FALSE]
    }
    dt_acc <- dt_nonprom[accessible == TRUE]
    n_accessible <- nrow(dt_acc)

    # --- Cumulative SNP proximity on accessible subset ---
    snp_col_safe <- opt$snp_col
    flank10_col <- paste0(opt$snp_prefix, "_flank10bp_overlap")
    flank25_col <- paste0(opt$snp_prefix, "_flank25bp_overlap")
    flank50_col <- paste0(opt$snp_prefix, "_flank50bp_overlap")
    flank100_col <- paste0(opt$snp_prefix, "_flank100bp_overlap")

    # Compute cumulative flags on accessible subset
    if (nrow(dt_acc) > 0) {
        dt_acc[, snp_within_0   := get(snp_col_safe) == TRUE]
        dt_acc[, snp_within_10  := snp_within_0 | (get(flank10_col) == TRUE)]
        dt_acc[, snp_within_25  := snp_within_10 | (get(flank25_col) == TRUE)]
        dt_acc[, snp_within_50  := snp_within_25 | (get(flank50_col) == TRUE)]
        dt_acc[, snp_within_100 := snp_within_50 | (get(flank100_col) == TRUE)]
    }

    n_snp_overlap <- if (nrow(dt_acc) > 0) sum(dt_acc$snp_within_0) else 0L
    n_snp_10bp    <- if (nrow(dt_acc) > 0) sum(dt_acc$snp_within_10) else 0L
    n_snp_25bp    <- if (nrow(dt_acc) > 0) sum(dt_acc$snp_within_25) else 0L
    n_snp_50bp    <- if (nrow(dt_acc) > 0) sum(dt_acc$snp_within_50) else 0L
    n_snp_100bp   <- if (nrow(dt_acc) > 0) sum(dt_acc$snp_within_100) else 0L

    # --- Same cumulative SNP on non-promoter WITHOUT ATAC ---
    if (nrow(dt_nonprom) > 0) {
        dt_nonprom[, snp_within_0   := get(snp_col_safe) == TRUE]
        dt_nonprom[, snp_within_10  := snp_within_0 | (get(flank10_col) == TRUE)]
        dt_nonprom[, snp_within_25  := snp_within_10 | (get(flank25_col) == TRUE)]
        dt_nonprom[, snp_within_50  := snp_within_25 | (get(flank50_col) == TRUE)]
        dt_nonprom[, snp_within_100 := snp_within_50 | (get(flank100_col) == TRUE)]
    }

    n_noatac_snp_overlap <- if (nrow(dt_nonprom) > 0) sum(dt_nonprom$snp_within_0) else 0L
    n_noatac_snp_10bp    <- if (nrow(dt_nonprom) > 0) sum(dt_nonprom$snp_within_10) else 0L
    n_noatac_snp_25bp    <- if (nrow(dt_nonprom) > 0) sum(dt_nonprom$snp_within_25) else 0L
    n_noatac_snp_50bp    <- if (nrow(dt_nonprom) > 0) sum(dt_nonprom$snp_within_50) else 0L
    n_noatac_snp_100bp   <- if (nrow(dt_nonprom) > 0) sum(dt_nonprom$snp_within_100) else 0L

    archetype_summary[[arch_id]] <- data.table(
        archetype = arch_id,
        n_total = n_total,
        n_intragenic = n_intragenic,
        n_intragenic_pcl = n_intragenic_pcl,
        n_intragenic_nonpromoter = n_nonprom,
        n_accessible = n_accessible,
        n_snp_overlap = n_snp_overlap,
        n_snp_10bp = n_snp_10bp,
        n_snp_25bp = n_snp_25bp,
        n_snp_50bp = n_snp_50bp,
        n_snp_100bp = n_snp_100bp,
        n_noatac_snp_overlap = n_noatac_snp_overlap,
        n_noatac_snp_10bp = n_noatac_snp_10bp,
        n_noatac_snp_25bp = n_noatac_snp_25bp,
        n_noatac_snp_50bp = n_noatac_snp_50bp,
        n_noatac_snp_100bp = n_noatac_snp_100bp
    )

    # Score data for directly overlapping SNP TFBSs
    if (n_snp_overlap > 0) {
        score_data[[arch_id]] <- data.table(
            archetype = arch_id, score = dt_acc[snp_within_0 == TRUE]$score
        )
        snp_count_data[[arch_id]] <- data.table(
            archetype = arch_id,
            snp_count = dt_acc[snp_within_0 == TRUE][[opt$snp_count_col]]
        )
    }

    # Cumulative SNP count data at each threshold (for SNP-per-TFBS plots)
    for (dist_label in c("0", "10", "25", "50", "100")) {
        col_flag <- paste0("snp_within_", dist_label)
        if (nrow(dt_acc) > 0 && sum(dt_acc[[col_flag]]) > 0) {
            ready_dt <- dt_acc[get(col_flag) == TRUE]
            # For cumulative SNP count, sum all relevant overlap columns
            # E.g., within 25bp = direct_count + flank10bp_count + flank25bp_count
            if (dist_label == "0") {
                total_snps <- ready_dt[[opt$snp_count_col]]
            } else {
                # Sum direct + all flanking bins up to this distance
                count_cols <- opt$snp_count_col
                dists_to_include <- c(10, 25, 50, 100)
                dists_to_include <- dists_to_include[dists_to_include <= as.integer(dist_label)]
                for (d in dists_to_include) {
                    count_cols <- c(count_cols,
                                    paste0(opt$snp_prefix, "_flank", d, "bp_count"))
                }
                # Need to read those columns - they should already be loaded
                # Actually they might not all be in dt_acc. Let's be safe:
                available_count_cols <- intersect(count_cols, names(ready_dt))
                if (length(available_count_cols) > 0) {
                    total_snps <- rowSums(ready_dt[, ..available_count_cols])
                } else {
                    total_snps <- ready_dt[[opt$snp_count_col]]
                }
            }
            key <- paste0(arch_id, "_", dist_label)
            snp_count_cumul_data[[key]] <- data.table(
                archetype = arch_id,
                snp_dist = dist_label,
                snp_count = total_snps
            )
        }
    }

    all_tfbs_scores[[arch_id]] <- data.table(
        archetype = arch_id, score = dt$score
    )
}

cat("  Done processing ", length(annot_files), " archetypes\n\n")

# Combine
arch_sum <- rbindlist(archetype_summary)
scores_ready <- rbindlist(score_data)
scores_all   <- rbindlist(all_tfbs_scores)
snp_counts   <- rbindlist(snp_count_data)
snp_counts_cumul <- rbindlist(snp_count_cumul_data)

# Merge with metadata and PWM info
arch_sum <- merge(arch_sum, arch_family_map, by.x = "archetype", by.y = "cluster",
                   all.x = TRUE)
arch_sum <- merge(arch_sum, arch_tf_map, by.x = "archetype", by.y = "cluster",
                   all.x = TRUE)
arch_sum <- merge(arch_sum, pwm_filtered, by = "archetype", all.x = TRUE)

# Backward compat columns
arch_sum[, n_analysis_ready := n_snp_overlap]
arch_sum[, n_no_atac_constraint := n_noatac_snp_overlap]

arch_sum <- arch_sum[order(-n_snp_overlap)]

# Save summary table
summary_out <- file.path(opt$outdir, "archetype_filter_summary.tsv")
fwrite(arch_sum, summary_out, sep = "\t")
cat("[INFO] Archetype summary table: ", summary_out, "\n")
cat("  Total archetypes:        ", nrow(arch_sum), "\n")
cat("  With SNP overlap:        ", sum(arch_sum$n_snp_overlap > 0), "\n")
cat("  With SNP \u226410bp:        ", sum(arch_sum$n_snp_10bp > 0), "\n")
cat("  With SNP \u226425bp:        ", sum(arch_sum$n_snp_25bp > 0), "\n")
cat("  With SNP \u226450bp:        ", sum(arch_sum$n_snp_50bp > 0), "\n")
cat("  With SNP \u2264100bp:       ", sum(arch_sum$n_snp_100bp > 0), "\n")
cat("  Total overlap TFBSs:     ", comma(sum(arch_sum$n_snp_overlap)), "\n")
cat("  Total \u2264100bp TFBSs:    ", comma(sum(arch_sum$n_snp_100bp)), "\n\n")

# =============================================================================
# HELPER: Create a display label from archetype + tf_names
# =============================================================================

make_display <- function(dt) {
    paste0(dt$archetype, " (", dt$tf_names, ")")
}

# =============================================================================
# 3. PLOT: Global filtering funnel (with protein-coding/lncRNA step)
# =============================================================================

cat("[3] Building global filtering funnel...\n")

total_all     <- sum(arch_sum$n_total)
total_intra   <- sum(arch_sum$n_intragenic)
total_pcl     <- sum(arch_sum$n_intragenic_pcl)
total_nonprom <- sum(arch_sum$n_intragenic_nonpromoter)
total_acc     <- sum(arch_sum$n_accessible)
total_snp_ov  <- sum(arch_sum$n_snp_overlap)

global_funnel <- data.table(
    step = c("All TFBSs", "Intragenic",
             "Protein-coding /\nlncRNA",
             "Non-promoter",
             "Accessible\n(ATAC)",
             "Het SNP\noverlaps motif"),
    count = c(total_all, total_intra, total_pcl, total_nonprom,
              total_acc, total_snp_ov)
)
global_funnel[, step := factor(step, levels = step)]
global_funnel[, pct := count / total_all * 100]
global_funnel[, pct_label := paste0(comma(count), "\n(",
                                     sprintf("%.2f%%", pct), ")")]

p_global_funnel <- ggplot(global_funnel, aes(x = step, y = count, fill = step)) +
    geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
    geom_text(aes(label = pct_label), vjust = -0.3, size = FUNNEL_LABEL_SIZE,
              fontface = "bold") +
    scale_fill_manual(values = c("#bdbdbd", "#a1d99b", "#74c476",
                                  "#238b45", "#fd8d3c", "#d62728")) +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0, 0.15))) +
    labs(title = "Global Filtering Funnel (All Archetypes Combined)",
         subtitle = paste0(comma(total_all), " total TFBSs across ",
                           nrow(arch_sum), " archetypes (width \u2265 ",
                           opt$min_width, " bp) \u2192 ",
                           comma(total_snp_ov), " with het SNP overlapping motif"),
         x = NULL, y = "Number of TFBSs") +
    THEME +
    theme(axis.text.x = element_text(size = AXIS_TEXT_SIZE, lineheight = 0.9))

# =============================================================================
# 3b. PLOT: Zoomed funnel — cumulative SNP proximity after accessible
# =============================================================================

cat("[3b] Building zoomed SNP proximity funnel...\n")

total_snp_100 <- sum(arch_sum$n_snp_100bp)
total_snp_50  <- sum(arch_sum$n_snp_50bp)
total_snp_25  <- sum(arch_sum$n_snp_25bp)
total_snp_10  <- sum(arch_sum$n_snp_10bp)

zoomed_funnel <- data.table(
    step = c("Accessible\n(ATAC)",
             "Het SNP\n\u2264100 bp",
             "Het SNP\n\u226450 bp",
             "Het SNP\n\u226425 bp",
             "Het SNP\n\u226410 bp",
             "Het SNP\noverlaps motif"),
    count = c(total_acc, total_snp_100, total_snp_50,
              total_snp_25, total_snp_10, total_snp_ov)
)
zoomed_funnel[, step := factor(step, levels = step)]
zoomed_funnel[, pct_of_acc := count / total_acc * 100]
zoomed_funnel[, pct_label := paste0(comma(count), "\n(",
                                     sprintf("%.2f%%", pct_of_acc),
                                     ")")]

p_zoomed_funnel <- ggplot(zoomed_funnel, aes(x = step, y = count, fill = step)) +
    geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
    geom_text(aes(label = pct_label), vjust = -0.3, size = FUNNEL_LABEL_SIZE - 0.3,
              fontface = "bold") +
    scale_fill_manual(values = c("#fd8d3c", "#fdd0a2", "#fdae6b",
                                  "#fd8d3c", "#e6550d", "#d62728")) +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0, 0.18))) +
    labs(title = "SNP Proximity Breakdown (Accessible TFBSs)",
         subtitle = paste0("Starting from ", comma(total_acc),
                           " accessible TFBSs \u2014 cumulative het SNP proximity"),
         x = NULL, y = "Number of TFBSs") +
    THEME +
    theme(axis.text.x = element_text(size = AXIS_TEXT_SIZE, lineheight = 0.9))

# =============================================================================
# 4. Per-archetype funnel plots (standard + detailed)
# =============================================================================

cat("[4] Building per-archetype funnel plots...\n")

# --- Standard funnel (counts) ---
top_arch <- head(arch_sum[n_snp_overlap > 0], TOP_N_ARCHETYPES)

funnel_long <- melt(top_arch,
    id.vars = c("archetype", "tf_names"),
    measure.vars = c("n_total", "n_intragenic", "n_intragenic_pcl",
                      "n_intragenic_nonpromoter", "n_accessible", "n_snp_overlap"),
    variable.name = "filter_level", value.name = "count"
)
funnel_long[, filter_label := fcase(
    filter_level == "n_total",                    "All TFBSs",
    filter_level == "n_intragenic",               "Intragenic",
    filter_level == "n_intragenic_pcl",           "Protein-coding / lncRNA",
    filter_level == "n_intragenic_nonpromoter",   "Non-promoter",
    filter_level == "n_accessible",               "Accessible (ATAC)",
    filter_level == "n_snp_overlap",              "Het SNP overlaps motif"
)]
funnel_long[, filter_label := factor(filter_label,
    levels = names(FILTER_COLORS))]

funnel_long[, display := make_display(.SD), .SDcols = c("archetype", "tf_names")]
arch_order <- top_arch[order(n_snp_overlap)]$archetype
display_order <- paste0(arch_order, " (",
                         top_arch[match(arch_order, archetype)]$tf_names, ")")
funnel_long[, display := factor(display, levels = display_order)]

p_funnel_counts <- ggplot(funnel_long,
                           aes(x = display, y = count, fill = filter_label)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8) +
    scale_fill_manual(values = FILTER_COLORS, name = "Filter") +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0, 0.05))) +
    coord_flip() +
    labs(title = "Filtering Funnel: TFBS Counts per Archetype",
         subtitle = paste0("Top ", nrow(top_arch),
                           " archetypes by het-SNP-overlapping count (width \u2265 ",
                           opt$min_width, " bp)"),
         x = NULL, y = "Number of TFBSs") +
    THEME +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 7))

# --- Standard funnel (proportions) ---
funnel_long[, pct := count / count[filter_label == "All TFBSs"] * 100,
             by = archetype]

p_funnel_pct <- ggplot(funnel_long[filter_label != "All TFBSs"],
             aes(x = display, y = pct, fill = filter_label)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8) +
    scale_fill_manual(values = FILTER_COLORS[-1], name = "Filter") +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.05))) +
    coord_flip() +
    labs(title = "Filtering Funnel: Proportion of All TFBSs",
         subtitle = "Fraction of each archetype's TFBSs passing each filter",
         x = NULL, y = "Percentage of all TFBSs") +
    THEME +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 7))

# --- Detailed funnel (counts) with cumulative SNP proximity ---
DETAIL_FILTER_COLORS <- c(
    "All TFBSs"                       = "#bdbdbd",
    "Intragenic"                      = "#a1d99b",
    "Protein-coding / lncRNA"         = "#74c476",
    "Non-promoter"                    = "#238b45",
    "Accessible (ATAC)"               = "#fd8d3c",
    "Het SNP \u2264100 bp"            = "#fdd0a2",
    "Het SNP \u226450 bp"             = "#fdae6b",
    "Het SNP \u226425 bp"             = "#fd8d3c",
    "Het SNP \u226410 bp"             = "#e6550d",
    "Het SNP overlaps motif"          = "#d62728"
)

detail_long <- melt(top_arch,
    id.vars = c("archetype", "tf_names"),
    measure.vars = c("n_total", "n_intragenic", "n_intragenic_pcl",
                      "n_intragenic_nonpromoter", "n_accessible",
                      "n_snp_100bp", "n_snp_50bp", "n_snp_25bp",
                      "n_snp_10bp", "n_snp_overlap"),
    variable.name = "filter_level", value.name = "count"
)
detail_long[, filter_label := fcase(
    filter_level == "n_total",                    "All TFBSs",
    filter_level == "n_intragenic",               "Intragenic",
    filter_level == "n_intragenic_pcl",           "Protein-coding / lncRNA",
    filter_level == "n_intragenic_nonpromoter",   "Non-promoter",
    filter_level == "n_accessible",               "Accessible (ATAC)",
    filter_level == "n_snp_100bp",                "Het SNP \u2264100 bp",
    filter_level == "n_snp_50bp",                 "Het SNP \u226450 bp",
    filter_level == "n_snp_25bp",                 "Het SNP \u226425 bp",
    filter_level == "n_snp_10bp",                 "Het SNP \u226410 bp",
    filter_level == "n_snp_overlap",              "Het SNP overlaps motif"
)]
detail_long[, filter_label := factor(filter_label,
    levels = names(DETAIL_FILTER_COLORS))]

detail_long[, display := paste0(archetype, " (", tf_names, ")")]
detail_long[, display := factor(display, levels = display_order)]

p_funnel_detail_counts <- ggplot(detail_long,
                                  aes(x = display, y = count, fill = filter_label)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.85) +
    scale_fill_manual(values = DETAIL_FILTER_COLORS, name = "Filter") +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0, 0.05))) +
    coord_flip() +
    labs(title = "Detailed Filtering Funnel: TFBS Counts (with SNP proximity)",
         subtitle = paste0("Top ", nrow(top_arch),
                           " archetypes — cumulative het SNP proximity thresholds"),
         x = NULL, y = "Number of TFBSs") +
    THEME +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 7),
          legend.key.size = unit(0.35, "cm"))

# --- Detailed funnel (proportions) ---
detail_long[, pct := count / count[filter_label == "All TFBSs"] * 100,
             by = archetype]

p_funnel_detail_pct <- ggplot(detail_long[filter_label != "All TFBSs"],
             aes(x = display, y = pct, fill = filter_label)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.85) +
    scale_fill_manual(values = DETAIL_FILTER_COLORS[-1], name = "Filter") +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.05))) +
    coord_flip() +
    labs(title = "Detailed Filtering Funnel: Proportions (with SNP proximity)",
         subtitle = "Fraction of each archetype's TFBSs passing each filter",
         x = NULL, y = "Percentage of all TFBSs") +
    THEME +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 7),
          legend.key.size = unit(0.35, "cm"))

# =============================================================================
# 5. ATAC constraint comparison (standard + detailed)
# =============================================================================

cat("[5] Building ATAC constraint comparison plots...\n")

# --- Standard version (overlap only) ---
atac_compare <- arch_sum[n_noatac_snp_overlap > 0][order(-n_noatac_snp_overlap)]
atac_compare <- head(atac_compare, TOP_N_ARCHETYPES)

atac_long <- melt(atac_compare,
    id.vars = c("archetype", "tf_names"),
    measure.vars = c("n_noatac_snp_overlap", "n_snp_overlap"),
    variable.name = "filter", value.name = "count"
)
atac_long[, filter_label := fifelse(
    filter == "n_snp_overlap",
    "With ATAC + het SNP overlaps motif",
    "Without ATAC + het SNP overlaps motif"
)]
atac_long[, display := paste0(archetype, " (", tf_names, ")")]
atac_order <- atac_compare[order(n_noatac_snp_overlap)]$archetype
atac_display_order <- paste0(atac_order, " (",
                              atac_compare[match(atac_order, archetype)]$tf_names, ")")
atac_long[, display := factor(display, levels = atac_display_order)]

p_atac <- ggplot(atac_long, aes(x = display, y = count, fill = filter_label)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    scale_fill_manual(values = c("With ATAC + het SNP overlaps motif" = "#d62728",
                                  "Without ATAC + het SNP overlaps motif" = "#ff9896"),
                       name = NULL) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
    coord_flip() +
    labs(title = "Effect of Chromatin Accessibility Constraint",
         subtitle = "TFBSs with het SNP overlapping motif: with vs without ATAC requirement",
         x = NULL, y = "Number of TFBSs") +
    THEME +
    theme(legend.position = "bottom", axis.text.y = element_text(size = 7))

# --- Detailed version (all SNP proximity thresholds) ---
atac_detail_measures <- c(
    "n_noatac_snp_100bp", "n_snp_100bp",
    "n_noatac_snp_50bp",  "n_snp_50bp",
    "n_noatac_snp_25bp",  "n_snp_25bp",
    "n_noatac_snp_10bp",  "n_snp_10bp",
    "n_noatac_snp_overlap", "n_snp_overlap"
)

atac_detail_long <- melt(atac_compare,
    id.vars = c("archetype", "tf_names"),
    measure.vars = atac_detail_measures,
    variable.name = "filter", value.name = "count"
)
atac_detail_long[, atac_status := fifelse(
    grepl("noatac", filter), "Without ATAC", "With ATAC"
)]
atac_detail_long[, snp_thresh := fcase(
    grepl("100bp", filter), "\u2264100 bp",
    grepl("50bp", filter),  "\u226450 bp",
    grepl("25bp", filter),  "\u226425 bp",
    grepl("10bp", filter),  "\u226410 bp",
    grepl("overlap", filter), "Overlaps motif"
)]
atac_detail_long[, label := paste0(atac_status, " + ", snp_thresh)]
atac_detail_long[, display := paste0(archetype, " (", tf_names, ")")]
atac_detail_long[, display := factor(display, levels = atac_display_order)]

p_atac_detail <- ggplot(atac_detail_long,
                          aes(x = display, y = count, fill = label)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.85) +
    scale_fill_manual(
        values = c(
            "With ATAC + \u2264100 bp" = "#fdd0a2",
            "Without ATAC + \u2264100 bp" = "#c6dbef",
            "With ATAC + \u226450 bp" = "#fdae6b",
            "Without ATAC + \u226450 bp" = "#9ecae1",
            "With ATAC + \u226425 bp" = "#fd8d3c",
            "Without ATAC + \u226425 bp" = "#6baed6",
            "With ATAC + \u226410 bp" = "#e6550d",
            "Without ATAC + \u226410 bp" = "#3182bd",
            "With ATAC + Overlaps motif" = "#d62728",
            "Without ATAC + Overlaps motif" = "#08519c"
        ), name = NULL) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
    coord_flip() +
    labs(title = "ATAC Constraint × SNP Proximity (Detailed)",
         subtitle = "With vs without ATAC at each cumulative SNP distance threshold",
         x = NULL, y = "Number of TFBSs") +
    THEME +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 7),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 8))

# =============================================================================
# 6. Top archetypes ranked — one version per SNP proximity threshold
# =============================================================================

cat("[6] Building top archetype ranking plots...\n")

snp_rank_versions <- list(
    list(col = "n_snp_overlap", label = "Het SNP overlaps motif", suffix = "overlap"),
    list(col = "n_snp_10bp",    label = "Het SNP \u226410 bp from motif", suffix = "10bp"),
    list(col = "n_snp_25bp",    label = "Het SNP \u226425 bp from motif", suffix = "25bp"),
    list(col = "n_snp_50bp",    label = "Het SNP \u226450 bp from motif", suffix = "50bp"),
    list(col = "n_snp_100bp",   label = "Het SNP \u2264100 bp from motif", suffix = "100bp")
)

p_top_ranked <- list()
p_heatmaps  <- list()

for (ver in snp_rank_versions) {
    col_name <- ver$col
    label    <- ver$label
    suffix   <- ver$suffix

    top_r <- head(arch_sum[get(col_name) > 0][order(-get(col_name))],
                  TOP_N_ARCHETYPES)
    if (nrow(top_r) == 0) next

    top_r[, display := paste0(archetype, " (", tf_names, ")")]
    top_r[, display := factor(display,
        levels = top_r[order(get(col_name))]$display)]

    p_top_ranked[[suffix]] <- ggplot(
        top_r, aes(x = display, y = get(col_name))) +
        geom_bar(stat = "identity", width = 0.7, fill = "#4a90d9") +
        scale_y_continuous(labels = comma,
                           expand = expansion(mult = c(0, 0.05))) +
        coord_flip() +
        labs(title = paste0("Top Archetypes: ", label),
             subtitle = paste0("Intragenic, protein-coding/lncRNA, non-promoter, ",
                               "accessible (width \u2265 ", opt$min_width, " bp)"),
             x = NULL, y = "Number of qualifying TFBSs") +
        THEME +
        theme(axis.text.y = element_text(size = 7))

    # --- Pass-rate heatmap for this threshold ---
    top_h <- head(arch_sum[get(col_name) > 0][order(-get(col_name))], 30)
    if (nrow(top_h) == 0) next

    top_h[, `:=`(
        pct_intragenic = n_intragenic / n_total * 100,
        pct_pcl = n_intragenic_pcl / n_total * 100,
        pct_nonpromoter = n_intragenic_nonpromoter / n_total * 100,
        pct_accessible = n_accessible / n_total * 100,
        pct_ready = get(col_name) / n_total * 100
    )]

    heat_l <- melt(top_h,
        id.vars = c("archetype", "tf_names"),
        measure.vars = c("pct_intragenic", "pct_pcl", "pct_nonpromoter",
                          "pct_accessible", "pct_ready"),
        variable.name = "filter", value.name = "pct"
    )
    heat_l[, filter_label := fcase(
        filter == "pct_intragenic",  "Intragenic",
        filter == "pct_pcl",         "PC/lncRNA",
        filter == "pct_nonpromoter", "Non-promoter",
        filter == "pct_accessible",  "Accessible",
        filter == "pct_ready",       label
    )]
    heat_l[, filter_label := factor(filter_label,
        levels = c("Intragenic", "PC/lncRNA", "Non-promoter",
                    "Accessible", label))]

    heat_l[, display := paste0(archetype, " (", tf_names, ")")]
    arch_h_order <- top_h[order(get(col_name))]$archetype
    heat_display_order <- paste0(arch_h_order, " (",
        top_h[match(arch_h_order, archetype)]$tf_names, ")")
    heat_l[, display := factor(display, levels = heat_display_order)]

    p_heatmaps[[suffix]] <- ggplot(heat_l,
                                    aes(x = filter_label, y = display, fill = pct)) +
        geom_tile(color = "white", linewidth = 0.5) +
        geom_text(aes(label = sprintf("%.1f", pct)), size = 2.8) +
        scale_fill_viridis_c(option = "magma", direction = -1,
                              name = "% of all TFBSs", limits = c(0, NA)) +
        labs(title = paste0("Filter Pass Rates: ", label),
             subtitle = "Percentage of all TFBSs passing each filter level",
             x = NULL, y = NULL) +
        THEME +
        theme(axis.text.y = element_text(size = 7),
              axis.text.x = element_text(angle = 30, hjust = 1, size = 9))
}

# =============================================================================
# 7. MOODS score distributions (uniform color, not rainbow)
# =============================================================================

cat("[7] Building MOODS score distribution plots...\n")

MIN_FOR_SCORE_DIST <- 20
score_arch_counts <- scores_ready[, .N, by = archetype][order(-N)]
top_score_arch <- score_arch_counts[N >= MIN_FOR_SCORE_DIST]
top_score_arch <- head(top_score_arch, 30)

scores_plot <- scores_ready[archetype %in% top_score_arch$archetype]
scores_plot <- merge(scores_plot, arch_tf_map, by.x = "archetype",
                      by.y = "cluster", all.x = TRUE)
scores_plot[, display := paste0(archetype, " (", tf_names, ")")]

med_scores <- scores_plot[, .(med = median(score)), by = display]
scores_plot[, display := factor(display,
    levels = med_scores[order(med)]$display)]

p_moods <- ggplot(scores_plot, aes(x = display, y = score)) +
    geom_boxplot(fill = "#6baed6", color = "#2171b5",
                 outlier.size = 0.5, outlier.alpha = 0.3) +
    coord_flip() +
    labs(title = "MOODS Score Distribution (Het-SNP-Overlapping TFBSs)",
         subtitle = paste0("Archetypes with \u2265", MIN_FOR_SCORE_DIST,
                           " intragenic, non-promoter, accessible, SNP-overlapping TFBSs"),
         x = NULL, y = "MOODS score") +
    THEME +
    theme(axis.text.y = element_text(size = 7))

# Score comparison (all vs analysis-ready)
rep_archetypes <- head(top_score_arch$archetype, 12)
scores_compare <- rbind(
    scores_all[archetype %in% rep_archetypes][, subset := "All TFBSs"],
    scores_ready[archetype %in% rep_archetypes][, subset := "SNP-overlapping"]
)
scores_compare <- merge(scores_compare, arch_tf_map, by.x = "archetype",
                         by.y = "cluster", all.x = TRUE)
scores_compare[, display := paste0(archetype, "\n(", tf_names, ")")]

p_score_compare <- ggplot(scores_compare,
                           aes(x = score, fill = subset, color = subset)) +
    geom_density(alpha = 0.4, linewidth = 0.5) +
    facet_wrap(~ display, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = c("All TFBSs" = "#bdbdbd",
                                  "SNP-overlapping" = "#d62728"), name = NULL) +
    scale_color_manual(values = c("All TFBSs" = "#636363",
                                   "SNP-overlapping" = "#d62728"), name = NULL) +
    labs(title = "MOODS Score: All TFBSs vs SNP-Overlapping Subset",
         subtitle = "Checking whether filtering biases the score distribution",
         x = "MOODS score", y = "Density") +
    THEME +
    theme(legend.position = "bottom", strip.text = element_text(size = 9))

# =============================================================================
# 8. SNP count distributions — multiple versions per proximity threshold
# =============================================================================

cat("[8] Building SNP count distribution plots...\n")

p_snp_count <- list()

snp_dist_labels <- c(
    "0"   = "Het SNP overlaps motif",
    "10"  = "Het SNP \u226410 bp from motif",
    "25"  = "Het SNP \u226425 bp from motif",
    "50"  = "Het SNP \u226450 bp from motif",
    "100" = "Het SNP \u2264100 bp from motif"
)

for (dist_val in c("0", "10", "25", "50", "100")) {
    dist_data <- snp_counts_cumul[snp_dist == dist_val]
    if (nrow(dist_data) == 0) next

    n_total_tfbs <- nrow(dist_data)
    pct_one <- round(100 * mean(dist_data$snp_count == 1), 1)

    # Bar counts for labeling
    count_tab <- dist_data[, .N, by = snp_count][order(snp_count)]
    count_tab[, pct := round(N / sum(N) * 100, 1)]

    p_snp_count[[dist_val]] <- ggplot(dist_data, aes(x = factor(snp_count))) +
        geom_bar(fill = "#d62728", width = 0.7) +
        geom_text(
            data = count_tab,
            aes(x = factor(snp_count), y = N,
                label = paste0(comma(N), "\n(", pct, "%)")),
            vjust = -0.2, size = BAR_LABEL_SIZE, fontface = "bold"
        ) +
        scale_y_continuous(labels = comma,
                           expand = expansion(mult = c(0, 0.15))) +
        labs(title = paste0("Het SNPs per Qualifying TFBS: ",
                            snp_dist_labels[dist_val]),
             subtitle = paste0(comma(n_total_tfbs), " TFBSs; ",
                               pct_one, "% with exactly 1 SNP"),
             x = "Number of het SNPs per TFBS",
             y = "Number of TFBSs") +
        THEME
}

# =============================================================================
# 9. Width and IC distribution plots (increased font sizes)
# =============================================================================

cat("[9] Building motif width and IC distribution plots...\n")

pwm_all[, passed_filter := archetype %in% VALID_ARCHETYPES]

p_width_hist <- ggplot(pwm_all, aes(x = motif_width, fill = passed_filter)) +
    geom_histogram(binwidth = 1, color = "white", linewidth = 0.3) +
    geom_vline(xintercept = opt$min_width - 0.5, linetype = "dashed",
               color = "red", linewidth = 0.6) +
    annotate("text", x = opt$min_width + 0.3, y = Inf,
             label = paste0("min_width = ", opt$min_width),
             hjust = 0, vjust = 1.5, size = 4.5, color = "red",
             fontface = "bold") +
    scale_fill_manual(values = c("TRUE" = "#4a90d9", "FALSE" = "#d9d9d9"),
                       name = NULL,
                       labels = c("TRUE" = "Included", "FALSE" = "Excluded")) +
    scale_x_continuous(breaks = seq(0, 30, by = 2)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(title = "Motif Archetype Width Distribution",
         subtitle = paste0("All ", nrow(pwm_all), " archetypes; ",
                           sum(pwm_all$passed_filter), " pass width \u2265 ",
                           opt$min_width, " bp filter"),
         x = "Motif width (bp)", y = "Number of archetypes") +
    THEME +
    theme(legend.position = c(0.85, 0.85),
          legend.text = element_text(size = LEGEND_SIZE))

p_ic_hist <- ggplot(pwm_filtered, aes(x = total_ic)) +
    geom_histogram(bins = 40, fill = "#e6550d", color = "white", linewidth = 0.3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(title = "Total Information Content Distribution",
         subtitle = paste0("Archetypes passing width filter (n = ",
                           nrow(pwm_filtered), "); range: ",
                           round(min(pwm_filtered$total_ic), 1), " \u2013 ",
                           round(max(pwm_filtered$total_ic), 1), " bits"),
         x = "Total information content (bits)", y = "Number of archetypes") +
    THEME

p_ic_per_pos <- ggplot(pwm_filtered, aes(x = mean_ic_per_pos)) +
    geom_histogram(bins = 40, fill = "#756bb1", color = "white", linewidth = 0.3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(title = "Mean IC per Position",
         subtitle = paste0("Median = ",
                           round(median(pwm_filtered$mean_ic_per_pos), 2),
                           " bits/pos (0 = uniform, 2 = perfect)"),
         x = "Mean IC per position (bits)", y = "Number of archetypes") +
    THEME

# Width vs IC scatter
pwm_scatter <- merge(pwm_filtered, arch_sum[, .(archetype, n_snp_overlap)],
                      by = "archetype", all.x = TRUE)
pwm_scatter[is.na(n_snp_overlap), n_snp_overlap := 0]
pwm_scatter[, has_ready := n_snp_overlap > 0]

p_width_vs_ic <- ggplot(pwm_scatter,
                         aes(x = motif_width, y = total_ic,
                             size = pmax(n_snp_overlap, 1),
                             color = has_ready)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("FALSE" = "#bdbdbd", "TRUE" = "#d62728"),
                        name = "Has qualifying\nTFBSs",
                        labels = c("No", "Yes")) +
    scale_size_continuous(name = "Qualifying\nTFBSs", range = c(1, 8),
                           trans = "sqrt",
                           breaks = c(1, 10, 50, 100, 500)) +
    scale_x_continuous(breaks = seq(0, 30, by = 2)) +
    labs(title = "Width vs Information Content",
         subtitle = paste0("Each point = one archetype (n = ",
                           nrow(pwm_scatter), ")"),
         x = "Motif width (bp)", y = "Total IC (bits)") +
    THEME +
    theme(legend.position = "right")

# =============================================================================
# 10. Histogram of archetypes by number of qualifying TFBSs
# =============================================================================

cat("[10] Building archetype yield histograms...\n")

p_arch_yield <- list()

for (ver in snp_rank_versions) {
    col_name <- ver$col
    label    <- ver$label
    suffix   <- ver$suffix

    yield_data <- arch_sum[get(col_name) > 0]
    yield_data <- data.table(count = yield_data[[col_name]])

    if (nrow(yield_data) == 0) next

    # Create sensible bins
    max_val <- max(yield_data$count)
    if (max_val <= 50) {
        bw <- 5
    } else if (max_val <= 200) {
        bw <- 10
    } else if (max_val <= 500) {
        bw <- 25
    } else {
        bw <- 50
    }

    n_arch <- nrow(yield_data)
    med_val <- median(yield_data$count)

    p_arch_yield[[suffix]] <- ggplot(yield_data, aes(x = count)) +
        geom_histogram(binwidth = bw, fill = "#4a90d9", color = "white",
                       linewidth = 0.3) +
        geom_vline(xintercept = med_val, linetype = "dashed",
                   color = "#d62728", linewidth = 0.6) +
        annotate("text", x = med_val, y = Inf,
                 label = paste0("median = ", round(med_val)),
                 hjust = -0.1, vjust = 1.5, size = 4, color = "#d62728",
                 fontface = "bold") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
        labs(title = paste0("Archetype Yield Distribution: ", label),
             subtitle = paste0(n_arch, " archetypes with \u22651 qualifying TFBS; ",
                               "total = ", comma(sum(yield_data$count)),
                               " TFBSs"),
             x = "Number of qualifying TFBSs per archetype",
             y = "Number of archetypes") +
        THEME
}

# =============================================================================
# 11. Assemble and save
# =============================================================================

cat("[11] Saving PDF...\n")

pdf_path <- file.path(opt$outdir, "motif_archetype_prep.pdf")

PREP_FIG_W <- 15
PREP_FIG_H <- 12

pdf(pdf_path, width = PREP_FIG_W, height = PREP_FIG_H)

# Page 1: Global funnel
print(p_global_funnel)

# Page 2: Zoomed SNP proximity funnel
print(p_zoomed_funnel)

# Page 3: Per-archetype funnel (counts)
print(p_funnel_counts)

# Page 4: Per-archetype funnel (proportions)
print(p_funnel_pct)

# Page 5: Detailed per-archetype funnel (counts)
print(p_funnel_detail_counts)

# Page 6: Detailed per-archetype funnel (proportions)
print(p_funnel_detail_pct)

# Page 7: ATAC constraint comparison
print(p_atac)

# Page 8: ATAC constraint detailed
print(p_atac_detail)

# Page 9+: Top ranked + heatmap for each SNP threshold
for (suffix in c("overlap", "10bp", "25bp", "50bp", "100bp")) {
    if (!is.null(p_top_ranked[[suffix]])) {
        print(p_top_ranked[[suffix]])
    }
    if (!is.null(p_heatmaps[[suffix]])) {
        print(p_heatmaps[[suffix]])
    }
}

# MOODS scores
print(p_moods)

# Score comparison
print(p_score_compare)

# SNP count distributions for each threshold
for (dist_val in c("0", "10", "25", "50", "100")) {
    if (!is.null(p_snp_count[[dist_val]])) {
        print(p_snp_count[[dist_val]])
    }
}

# Width and IC distributions
print((p_width_hist + p_ic_hist) / (p_ic_per_pos + p_width_vs_ic))

# Archetype yield histograms
for (suffix in c("overlap", "10bp", "25bp", "50bp", "100bp")) {
    if (!is.null(p_arch_yield[[suffix]])) {
        print(p_arch_yield[[suffix]])
    }
}

dev.off()

cat("\n[DONE] PDF saved: ", pdf_path, "\n")

# =============================================================================
# 12. Save individual PNGs
# =============================================================================

cat("[INFO] Saving individual PNGs...\n")
png_dir <- file.path(opt$outdir, "png")
dir.create(png_dir, showWarnings = FALSE)

PNG_DPI <- 300

# Core plots
core_plots <- list(
    "01_global_funnel"          = p_global_funnel,
    "02_zoomed_snp_funnel"      = p_zoomed_funnel,
    "03_funnel_counts"          = p_funnel_counts,
    "04_funnel_proportions"     = p_funnel_pct,
    "05_funnel_detail_counts"   = p_funnel_detail_counts,
    "06_funnel_detail_pct"      = p_funnel_detail_pct,
    "07_atac_comparison"        = p_atac,
    "08_atac_detail"            = p_atac_detail,
    "20_moods_scores"           = p_moods,
    "21_score_comparison"       = p_score_compare,
    "30_width_distribution"     = p_width_hist,
    "31_total_ic_dist"          = p_ic_hist,
    "32_ic_per_position"        = p_ic_per_pos,
    "33_width_vs_ic"            = p_width_vs_ic
)

for (nm in names(core_plots)) {
    w <- if (grepl("funnel|ranked|heatmap|comparison|detail", nm)) 14 else 11
    h <- if (grepl("funnel|ranked|heatmap|detail", nm)) 11 else 8
    ggsave(file.path(png_dir, paste0(nm, ".png")),
           plot = core_plots[[nm]], width = w, height = h, dpi = PNG_DPI)
}

# Versioned plots (per SNP threshold)
for (suffix in c("overlap", "10bp", "25bp", "50bp", "100bp")) {
    if (!is.null(p_top_ranked[[suffix]])) {
        ggsave(file.path(png_dir, paste0("10_top_ranked_", suffix, ".png")),
               plot = p_top_ranked[[suffix]], width = 14, height = 11, dpi = PNG_DPI)
    }
    if (!is.null(p_heatmaps[[suffix]])) {
        ggsave(file.path(png_dir, paste0("11_passrate_heatmap_", suffix, ".png")),
               plot = p_heatmaps[[suffix]], width = 14, height = 11, dpi = PNG_DPI)
    }
    if (!is.null(p_arch_yield[[suffix]])) {
        ggsave(file.path(png_dir, paste0("40_archetype_yield_", suffix, ".png")),
               plot = p_arch_yield[[suffix]], width = 11, height = 8, dpi = PNG_DPI)
    }
}

for (dist_val in c("0", "10", "25", "50", "100")) {
    if (!is.null(p_snp_count[[dist_val]])) {
        ggsave(file.path(png_dir, paste0("22_snp_count_", dist_val, "bp.png")),
               plot = p_snp_count[[dist_val]], width = 11, height = 8, dpi = PNG_DPI)
    }
}

cat("[DONE] PNGs saved: ", png_dir, "\n\n")

# =============================================================================
# Console summary
# =============================================================================

cat("================================================================\n")
cat("Quick Summary\n")
cat("================================================================\n")
cat(sprintf("  MEME archetypes parsed:       %s\n", comma(nrow(pwm_all))))
cat(sprintf("  Excluded by width < %d bp:    %s\n",
            opt$min_width, comma(sum(!pwm_all$passed_filter))))
cat(sprintf("  Archetypes after PWM filter:  %s\n", comma(nrow(arch_sum))))
cat(sprintf("  Width range (included):       %d \u2013 %d bp\n",
            min(arch_sum$motif_width, na.rm = TRUE),
            max(arch_sum$motif_width, na.rm = TRUE)))
cat(sprintf("  Total TFBSs:                  %s\n", comma(total_all)))
cat(sprintf("  Intragenic:                   %s (%.1f%%)\n",
            comma(total_intra), 100 * total_intra / total_all))
cat(sprintf("  Protein-coding/lncRNA:        %s (%.1f%%)\n",
            comma(total_pcl), 100 * total_pcl / total_all))
cat(sprintf("  Non-promoter:                 %s (%.1f%%)\n",
            comma(total_nonprom), 100 * total_nonprom / total_all))
cat(sprintf("  Accessible:                   %s (%.1f%%)\n",
            comma(total_acc), 100 * total_acc / total_all))
cat(sprintf("  Het SNP overlaps motif:       %s (%.2f%%)\n",
            comma(total_snp_ov), 100 * total_snp_ov / total_all))
cat(sprintf("  Het SNP \u226410 bp:              %s (%.2f%%)\n",
            comma(sum(arch_sum$n_snp_10bp)),
            100 * sum(arch_sum$n_snp_10bp) / total_all))
cat(sprintf("  Het SNP \u226425 bp:              %s (%.2f%%)\n",
            comma(sum(arch_sum$n_snp_25bp)),
            100 * sum(arch_sum$n_snp_25bp) / total_all))
cat(sprintf("  Het SNP \u226450 bp:              %s (%.2f%%)\n",
            comma(sum(arch_sum$n_snp_50bp)),
            100 * sum(arch_sum$n_snp_50bp) / total_all))
cat(sprintf("  Het SNP \u2264100 bp:             %s (%.2f%%)\n",
            comma(sum(arch_sum$n_snp_100bp)),
            100 * sum(arch_sum$n_snp_100bp) / total_all))
cat(sprintf("\n  Archetypes w/ \u22651 overlap:    %s\n",
            sum(arch_sum$n_snp_overlap > 0)))
cat(sprintf("  Archetypes w/ \u226520 overlap:   %s\n",
            sum(arch_sum$n_snp_overlap >= 20)))
cat(sprintf("  Archetypes w/ \u226550 overlap:   %s\n",
            sum(arch_sum$n_snp_overlap >= 50)))
cat(sprintf("  Archetypes w/ \u2265100 overlap:  %s\n",
            sum(arch_sum$n_snp_overlap >= 100)))
cat("================================================================\n")