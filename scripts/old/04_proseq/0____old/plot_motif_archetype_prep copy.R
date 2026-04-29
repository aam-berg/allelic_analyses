#!/usr/bin/env Rscript
# =============================================================================
# plot_motif_archetype_prep.R — Preparatory analysis for TFBS × PRO-seq study
# =============================================================================
#
# Generates summary plots characterizing the annotated motif archetypes,
# focusing on the subset relevant for allele-specific PRO-seq analysis:
#   - Non-promoter, intragenic, chromatin-accessible TFBSs with het SNPs
#
# NEW IN THIS VERSION:
#   - Parses MEME-format PWM file to extract motif width and information content
#   - Filters archetypes by minimum width (default: 8 bp) and optionally by IC
#   - Adds plots for width/IC distributions and width-vs-IC scatter
#
# PLOTS GENERATED:
#   Page 1: Global filtering funnel + TF family breakdown
#   Page 2: Per-archetype filtering funnel (absolute counts)
#   Page 3: Per-archetype filtering funnel (proportions)
#   Page 4: ATAC constraint comparison
#   Page 5: Top archetypes ranked + pass-rate heatmap
#   Page 6: MOODS score boxplots
#   Page 7: Score comparison (all vs analysis-ready)
#   Page 8: SNP count distribution
#   Page 9: Motif width/IC distributions + width vs IC scatter   ← NEW
#
# USAGE:
#   Rscript plot_motif_archetype_prep.R \
#       --annot_dir /path/to/annot_motifs \
#       --metadata /path/to/metadata.tsv \
#       --meme_file /path/to/consensus_pwms.meme \
#       --min_width 8 \
#       --min_ic 0 \
#       --outdir /path/to/output
#
# DEPENDENCIES:
#   R packages: data.table, ggplot2, patchwork, scales, viridis
#
# CUSTOMIZATION:
#   Search for "=== CUSTOMIZE ===" for adjustable settings.
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
                help = "MEME-format PWM file for archetype widths and IC [default: %default]"),
    make_option("--min_width", type = "integer", default = 8,
                help = "Minimum motif width (bp) to include archetype [default: %default]. Set to 0 to disable."),
    make_option("--min_ic", type = "double", default = 0,
                help = "Minimum total information content (bits) to include archetype [default: %default]. Set to 0 to disable."),
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
                help = "SNP count column name")
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
# === CUSTOMIZE === Global theme settings
# =============================================================================

BASE_SIZE       <- 11       # Base font size for body text
TITLE_SIZE      <- 13       # Plot title size
SUBTITLE_SIZE   <- 10       # Subtitle / annotation size
AXIS_TEXT_SIZE  <- 9        # Axis tick labels
LEGEND_SIZE     <- 9        # Legend text

FIG_WIDTH  <- 13            # PDF page width (inches)
FIG_HEIGHT <- 10            # PDF page height (inches)

# Colors for filter categories in the funnel plots
FILTER_COLORS <- c(
    "All TFBSs"                     = "#bdbdbd",
    "Intragenic"                    = "#74c476",
    "Intragenic, non-promoter"      = "#238b45",
    "...+ accessible"               = "#fd8d3c",
    "...+ het SNP (analysis-ready)" = "#d62728"
)

# === CUSTOMIZE === Top N archetypes to highlight in ranked bar charts
TOP_N_ARCHETYPES <- 40

THEME <- theme_bw(base_size = BASE_SIZE) +
    theme(
        plot.title    = element_text(size = TITLE_SIZE, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = SUBTITLE_SIZE, color = "grey40"),
        axis.text     = element_text(size = AXIS_TEXT_SIZE),
        axis.title    = element_text(size = BASE_SIZE),
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

# Build archetype → TF family mapping
arch_family <- unique(meta[, .(cluster, family_name)])
arch_family_map <- arch_family[, .(
    families = paste(unique(family_name), collapse = "|")
), by = cluster]

# Archetype → TF name(s)
arch_tf_map <- meta[, .(
    tf_names = paste(unique(tf_name), collapse = "|")
), by = cluster]

# =============================================================================
# 1b. Parse MEME file → extract motif width and information content
# =============================================================================
#
# The MEME file contains position weight matrices (PWMs) for each archetype.
# We parse it to get:
#   - motif_width:     number of positions in the PWM
#   - total_ic:        total information content in bits
#   - mean_ic_per_pos: IC per position (total_ic / width)
#
# Information Content (IC) per position:
#   IC_j = 2 + sum_b [ p(b,j) * log2(p(b,j)) ]
#   where b = {A,C,G,T} and background is uniform (0.25 each).
#   Range: 0 bits (uniform) to 2 bits (perfectly conserved).
#   Total IC = sum across all positions.
#   A 10-bp motif with perfect conservation = 20 bits total IC.
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

        # Look for MOTIF lines: "MOTIF AC0001:DLX/LHX:Homeodomain ..."
        if (startsWith(line, "MOTIF ")) {
            # Extract archetype ID = first colon-delimited token after "MOTIF "
            # e.g. "MOTIF AC0001:DLX/LHX:Homeodomain ..." → "AC0001"
            tokens <- strsplit(sub("^MOTIF\\s+", "", line), ":")[[1]]
            arch_id <- tokens[1]

            # Find the "letter-probability matrix" line
            i <- i + 1
            while (i <= length(lines) && !grepl("^letter-probability matrix", lines[i])) {
                i <- i + 1
            }
            if (i > length(lines)) break

            # Parse header: "letter-probability matrix: alength= 4 w= 6 ..."
            header <- lines[i]
            width <- as.integer(sub(".*w=\\s*(\\d+).*", "\\1", header))

            # Read PWM rows (next 'width' non-empty numeric lines)
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

            # Compute IC per position: IC_j = 2 + sum(p * log2(p))
            ic_per_pos <- apply(pwm, 1, function(p) {
                p <- pmax(p, 1e-10)  # avoid log(0)
                2 + sum(p * log2(p))
            })
            total_ic <- sum(ic_per_pos)
            mean_ic <- total_ic / width

            results[[arch_id]] <- data.table(
                archetype       = arch_id,
                motif_width     = as.integer(width),
                total_ic        = round(total_ic, 3),
                mean_ic_per_pos = round(mean_ic, 3)
            )
        } else {
            i <- i + 1
        }
    }
    rbindlist(results)
}

# pwm_all: ALL archetypes from the MEME file (before any filtering)
# We keep this for the "before vs after" width distribution plots.
pwm_all <- parse_meme(opt$meme_file)
cat("  Parsed ", nrow(pwm_all), " archetypes from MEME file\n")
cat("  Width range:    ", min(pwm_all$motif_width), " – ",
    max(pwm_all$motif_width), " bp\n")
cat("  Total IC range: ", round(min(pwm_all$total_ic), 1), " – ",
    round(max(pwm_all$total_ic), 1), " bits\n")
cat("  Mean IC/pos:    ", round(min(pwm_all$mean_ic_per_pos), 2), " – ",
    round(max(pwm_all$mean_ic_per_pos), 2), " bits/pos\n")

# --- Apply width and IC filters ---
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
    cat("  After min_ic >= ", opt$min_ic, ":    ",
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

if (length(annot_files) == 0) {
    stop("No annotated files found in ", opt$annot_dir)
}

# Filter to only archetypes passing the PWM width/IC filters
annot_files <- annot_files[
    sub("_annotated\\.tsv\\.gz$", "", basename(annot_files)) %in% VALID_ARCHETYPES
]
cat("  After PWM filter (width >= ", opt$min_width,
    if (opt$min_ic > 0) paste0(", IC >= ", opt$min_ic) else "",
    "): ", length(annot_files), " files\n")

# For each archetype, compute counts at each filter level.
archetype_summary <- list()
score_data <- list()       # MOODS scores for analysis-ready TFBSs
snp_count_data <- list()   # SNP counts for analysis-ready TFBSs
all_tfbs_scores <- list()  # all MOODS scores (for comparison)

for (f in annot_files) {
    arch_id <- sub("_annotated\\.tsv\\.gz$", "", basename(f))
    if ((which(annot_files == f) %% 50) == 1) {
        cat("  Processing ", arch_id, " (",
            which(annot_files == f), "/", length(annot_files), ")\n")
    }

    dt <- fread(f, select = c("chrom", "start", "end", "motif_id", "score",
                               "strand", "is_intragenic", "is_promoter",
                               atac_cols, opt$snp_col, opt$snp_count_col))

    n_total <- nrow(dt)
    if (n_total == 0) {
        archetype_summary[[arch_id]] <- data.table(
            archetype = arch_id, n_total = 0, n_intragenic = 0,
            n_intragenic_nonpromoter = 0, n_accessible = 0,
            n_analysis_ready = 0, n_no_atac_constraint = 0
        )
        next
    }

    # --- Progressive filters ---
    dt_intra <- dt[is_intragenic == TRUE]
    n_intragenic <- nrow(dt_intra)

    dt_nonprom <- dt_intra[is_promoter == FALSE]
    n_nonprom <- nrow(dt_nonprom)

    if (length(atac_cols) > 0 && all(atac_cols %in% names(dt_nonprom))) {
        dt_nonprom[, accessible := Reduce(`|`, .SD), .SDcols = atac_cols]
    } else {
        dt_nonprom[, accessible := FALSE]
    }
    dt_acc <- dt_nonprom[accessible == TRUE]
    n_accessible <- nrow(dt_acc)

    snp_col_safe <- opt$snp_col
    dt_ready <- dt_acc[get(snp_col_safe) == TRUE]
    n_ready <- nrow(dt_ready)

    dt_noatac <- dt_nonprom[get(snp_col_safe) == TRUE]
    n_noatac <- nrow(dt_noatac)

    archetype_summary[[arch_id]] <- data.table(
        archetype = arch_id,
        n_total = n_total,
        n_intragenic = n_intragenic,
        n_intragenic_nonpromoter = n_nonprom,
        n_accessible = n_accessible,
        n_analysis_ready = n_ready,
        n_no_atac_constraint = n_noatac
    )

    if (n_ready > 0) {
        score_data[[arch_id]] <- data.table(
            archetype = arch_id, score = dt_ready$score
        )
        snp_count_data[[arch_id]] <- data.table(
            archetype = arch_id, snp_count = dt_ready[[opt$snp_count_col]]
        )
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

# Merge with metadata and PWM info
arch_sum <- merge(arch_sum, arch_family_map, by.x = "archetype", by.y = "cluster",
                   all.x = TRUE)
arch_sum <- merge(arch_sum, arch_tf_map, by.x = "archetype", by.y = "cluster",
                   all.x = TRUE)
arch_sum <- merge(arch_sum, pwm_filtered, by = "archetype", all.x = TRUE)

arch_sum <- arch_sum[order(-n_analysis_ready)]

# Save summary table
summary_out <- file.path(opt$outdir, "archetype_filter_summary.tsv")
fwrite(arch_sum, summary_out, sep = "\t")
cat("[INFO] Archetype summary table: ", summary_out, "\n")
cat("  Total archetypes:      ", nrow(arch_sum), "\n")
cat("  Width range:           ", min(arch_sum$motif_width, na.rm = TRUE), " – ",
    max(arch_sum$motif_width, na.rm = TRUE), " bp\n")
cat("  Width filter applied:  >= ", opt$min_width, " bp\n")
if (opt$min_ic > 0) {
    cat("  IC filter applied:     >= ", opt$min_ic, " bits\n")
} else {
    cat("  IC filter:             none (min_ic = 0)\n")
}
cat("  With analysis-ready:   ",
    sum(arch_sum$n_analysis_ready > 0), "\n")
cat("  Total analysis-ready:  ",
    comma(sum(arch_sum$n_analysis_ready)), " TFBSs\n")
cat("  Without ATAC constraint: ",
    comma(sum(arch_sum$n_no_atac_constraint)), " TFBSs\n\n")

# =============================================================================
# 3. PLOT 1 & 2: Filtering funnel — top archetypes
# =============================================================================

cat("[3] Building filtering funnel plots...\n")

top_arch <- head(arch_sum[n_analysis_ready > 0], TOP_N_ARCHETYPES)

funnel_long <- melt(top_arch,
    id.vars = c("archetype", "tf_names"),
    measure.vars = c("n_total", "n_intragenic", "n_intragenic_nonpromoter",
                      "n_accessible", "n_analysis_ready"),
    variable.name = "filter_level", value.name = "count"
)

funnel_long[, filter_label := fcase(
    filter_level == "n_total",                    "All TFBSs",
    filter_level == "n_intragenic",               "Intragenic",
    filter_level == "n_intragenic_nonpromoter",   "Intragenic, non-promoter",
    filter_level == "n_accessible",               "...+ accessible",
    filter_level == "n_analysis_ready",           "...+ het SNP (analysis-ready)"
)]
funnel_long[, filter_label := factor(filter_label,
    levels = c("All TFBSs", "Intragenic", "Intragenic, non-promoter",
               "...+ accessible", "...+ het SNP (analysis-ready)"))]

funnel_long[, display := paste0(archetype, " (", tf_names, ")")]
arch_order <- top_arch[order(n_analysis_ready)]$archetype
display_order <- paste0(arch_order, " (",
                         top_arch[match(arch_order, archetype)]$tf_names, ")")
funnel_long[, display := factor(display, levels = display_order)]

p1 <- ggplot(funnel_long, aes(x = display, y = count, fill = filter_label)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8) +
    scale_fill_manual(values = FILTER_COLORS, name = "Filter") +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0, 0.05))) +
    coord_flip() +
    labs(title = "A. Filtering Funnel: TFBS Counts per Archetype",
         subtitle = paste0("Top ", nrow(top_arch),
                           " archetypes by analysis-ready count",
                           " (width \u2265 ", opt$min_width, " bp)"),
         x = NULL, y = "Number of TFBSs") +
    THEME +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 7))

# Proportion version
funnel_long[, pct := count / count[filter_label == "All TFBSs"] * 100,
             by = archetype]

p2 <- ggplot(funnel_long[filter_label != "All TFBSs"],
             aes(x = display, y = pct, fill = filter_label)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8) +
    scale_fill_manual(values = FILTER_COLORS[-1], name = "Filter") +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.05))) +
    coord_flip() +
    labs(title = "B. Filtering Funnel: Proportion of All TFBSs",
         subtitle = "What fraction of each archetype's TFBSs pass each filter",
         x = NULL, y = "Percentage of all TFBSs") +
    THEME +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 7))

# =============================================================================
# 4. PLOT 3: With vs without ATAC constraint
# =============================================================================

cat("[4] Building ATAC constraint comparison...\n")

atac_compare <- arch_sum[n_no_atac_constraint > 0][order(-n_no_atac_constraint)]
atac_compare <- head(atac_compare, TOP_N_ARCHETYPES)

atac_long <- melt(atac_compare,
    id.vars = c("archetype", "tf_names"),
    measure.vars = c("n_no_atac_constraint", "n_analysis_ready"),
    variable.name = "filter", value.name = "count"
)
atac_long[, filter_label := fifelse(
    filter == "n_analysis_ready",
    "With ATAC (accessible + SNP)",
    "Without ATAC (intragenic + non-promoter + SNP)"
)]
atac_long[, display := paste0(archetype, " (", tf_names, ")")]
atac_order <- atac_compare[order(n_no_atac_constraint)]$archetype
atac_display_order <- paste0(atac_order, " (",
                              atac_compare[match(atac_order, archetype)]$tf_names, ")")
atac_long[, display := factor(display, levels = atac_display_order)]

ATAC_COMPARE_COLORS <- c(
    "With ATAC (accessible + SNP)" = "#d62728",
    "Without ATAC (intragenic + non-promoter + SNP)" = "#ff9896"
)

p3 <- ggplot(atac_long, aes(x = display, y = count, fill = filter_label)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    scale_fill_manual(values = ATAC_COMPARE_COLORS, name = NULL) +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0, 0.05))) +
    coord_flip() +
    labs(title = "C. Effect of Chromatin Accessibility Constraint",
         subtitle = "How many TFBSs with SNPs are lost by requiring ATAC overlap",
         x = NULL, y = "Number of TFBSs") +
    THEME +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 7))

# =============================================================================
# 5. PLOT 4: MOODS score distributions
# =============================================================================

cat("[5] Building MOODS score distribution plots...\n")

# === CUSTOMIZE === Minimum TFBSs for score distribution plot
MIN_FOR_SCORE_DIST <- 20

score_arch_counts <- scores_ready[, .N, by = archetype][order(-N)]
top_score_arch <- score_arch_counts[N >= MIN_FOR_SCORE_DIST]
cat("  Archetypes with >=", MIN_FOR_SCORE_DIST,
    " analysis-ready TFBSs: ", nrow(top_score_arch), "\n")

top_score_arch <- head(top_score_arch, 30)

scores_plot <- scores_ready[archetype %in% top_score_arch$archetype]
scores_plot <- merge(scores_plot, arch_tf_map, by.x = "archetype",
                      by.y = "cluster", all.x = TRUE)
scores_plot[, display := paste0(archetype, " (", tf_names, ")")]

med_scores <- scores_plot[, .(med = median(score)), by = display]
scores_plot[, display := factor(display,
    levels = med_scores[order(med)]$display)]

p4 <- ggplot(scores_plot, aes(x = display, y = score, fill = display)) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3, show.legend = FALSE) +
    scale_fill_viridis_d(option = "turbo", guide = "none") +
    coord_flip() +
    labs(title = "D. MOODS Score Distribution (Analysis-Ready TFBSs)",
         subtitle = paste0("Archetypes with \u2265", MIN_FOR_SCORE_DIST,
                           " intragenic, non-promoter, accessible TFBSs with het SNPs"),
         x = NULL, y = "MOODS score") +
    THEME +
    theme(axis.text.y = element_text(size = 7))

# =============================================================================
# 6. PLOT 5: Comparison of MOODS scores — all TFBSs vs analysis-ready
# =============================================================================

cat("[6] Building score comparison plot...\n")

rep_archetypes <- head(top_score_arch$archetype, 12)

scores_compare <- rbind(
    scores_all[archetype %in% rep_archetypes][, subset := "All TFBSs"],
    scores_ready[archetype %in% rep_archetypes][, subset := "Analysis-ready"]
)
scores_compare <- merge(scores_compare, arch_tf_map, by.x = "archetype",
                         by.y = "cluster", all.x = TRUE)
scores_compare[, display := paste0(archetype, "\n(", tf_names, ")")]

p5 <- ggplot(scores_compare, aes(x = score, fill = subset, color = subset)) +
    geom_density(alpha = 0.4, linewidth = 0.5) +
    facet_wrap(~ display, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = c("All TFBSs" = "#bdbdbd",
                                  "Analysis-ready" = "#d62728"),
                       name = NULL) +
    scale_color_manual(values = c("All TFBSs" = "#636363",
                                   "Analysis-ready" = "#d62728"),
                        name = NULL) +
    labs(title = "E. MOODS Score: All TFBSs vs Analysis-Ready Subset",
         subtitle = "Checking whether filtering biases the score distribution (ideally similar shapes)",
         x = "MOODS score", y = "Density") +
    THEME +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 8))

# =============================================================================
# 7. PLOT 6: SNPs per TFBS
# =============================================================================

cat("[7] Building SNP count distribution...\n")

p6 <- ggplot(snp_counts, aes(x = factor(snp_count))) +
    geom_bar(fill = "#d62728", width = 0.7) +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0, 0.05))) +
    labs(title = "F. Het SNPs per Analysis-Ready TFBS",
         subtitle = paste0("Distribution across all ", comma(nrow(snp_counts)),
                           " analysis-ready TFBSs (",
                           round(100 * mean(snp_counts$snp_count == 1), 1),
                           "% with exactly 1 SNP)"),
         x = "Number of F121-9 het SNPs per TFBS",
         y = "Number of TFBSs") +
    THEME

# =============================================================================
# 8. PLOT 7: Top archetypes ranked by analysis-ready count
# =============================================================================

cat("[8] Building top archetype ranking...\n")

top_ranked <- head(arch_sum[n_analysis_ready > 0], TOP_N_ARCHETYPES)
top_ranked <- merge(top_ranked, arch_tf_map, by.x = "archetype",
                     by.y = "cluster", all.x = TRUE, suffixes = c("", ".dup"))
top_ranked[, display := paste0(archetype, " (", tf_names, ")")]
top_ranked[, display := factor(display,
    levels = top_ranked[order(n_analysis_ready)]$display)]

p7 <- ggplot(top_ranked, aes(x = display, y = n_analysis_ready,
                               fill = families)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_viridis_d(option = "plasma", name = "TF family") +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0, 0.05))) +
    coord_flip() +
    labs(title = "G. Top Archetypes by Analysis-Ready TFBS Count",
         subtitle = paste0("Intragenic, non-promoter, accessible, with het SNP",
                           " (width \u2265 ", opt$min_width, " bp)"),
         x = NULL, y = "Number of analysis-ready TFBSs") +
    THEME +
    theme(axis.text.y = element_text(size = 7),
          legend.position = "right")

# =============================================================================
# 8b. PLOT 7b: TF family breakdown
# =============================================================================

cat("[8b] Building TF family breakdown...\n")

family_counts <- arch_sum[n_analysis_ready > 0 & !is.na(families),
                           .(total_ready = sum(n_analysis_ready),
                             n_archetypes = .N),
                           by = families]
family_counts <- family_counts[order(-total_ready)]
family_top <- head(family_counts, 15)
family_top[, display := factor(families, levels = rev(families))]

p8 <- ggplot(family_top, aes(x = display, y = total_ready, fill = display)) +
    geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
    geom_text(aes(label = paste0("(", n_archetypes, " arch.)")),
              hjust = -0.1, size = 3, color = "grey30") +
    scale_fill_viridis_d(option = "mako", guide = "none") +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0, 0.15))) +
    coord_flip() +
    labs(title = "H. Analysis-Ready TFBSs by TF Family",
         subtitle = "Summed across all archetypes in each family",
         x = NULL, y = "Total analysis-ready TFBSs") +
    THEME

# =============================================================================
# 8c. PLOT 7c: Filtering pass-rate heatmap
# =============================================================================

cat("[8c] Building pass-rate heatmap...\n")

top_heat <- head(arch_sum[n_analysis_ready > 0], 30)
top_heat[, `:=`(
    pct_intragenic = n_intragenic / n_total * 100,
    pct_nonpromoter = n_intragenic_nonpromoter / n_total * 100,
    pct_accessible = n_accessible / n_total * 100,
    pct_ready = n_analysis_ready / n_total * 100
)]

heat_long <- melt(top_heat,
    id.vars = c("archetype", "tf_names"),
    measure.vars = c("pct_intragenic", "pct_nonpromoter",
                      "pct_accessible", "pct_ready"),
    variable.name = "filter", value.name = "pct"
)
heat_long[, filter_label := fcase(
    filter == "pct_intragenic",  "Intragenic",
    filter == "pct_nonpromoter", "Non-promoter",
    filter == "pct_accessible",  "Accessible",
    filter == "pct_ready",       "Analysis-ready"
)]
heat_long[, filter_label := factor(filter_label,
    levels = c("Intragenic", "Non-promoter", "Accessible", "Analysis-ready"))]

heat_long[, display := paste0(archetype, " (", tf_names, ")")]
arch_heat_order <- top_heat[order(n_analysis_ready)]$archetype
heat_display_order <- paste0(arch_heat_order, " (",
    top_heat[match(arch_heat_order, archetype)]$tf_names, ")")
heat_long[, display := factor(display, levels = heat_display_order)]

p_heatmap <- ggplot(heat_long, aes(x = filter_label, y = display, fill = pct)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f", pct)), size = 2.5) +
    scale_fill_viridis_c(option = "magma", direction = -1,
                          name = "% of all TFBSs",
                          limits = c(0, NA)) +
    labs(title = "I. Filter Pass Rates by Archetype",
         subtitle = "Percentage of all TFBSs passing each filter level",
         x = NULL, y = NULL) +
    THEME +
    theme(axis.text.y = element_text(size = 7),
          axis.text.x = element_text(angle = 30, hjust = 1))

# =============================================================================
# 9. PLOT 8: Motif width and IC distributions (NEW)
# =============================================================================

cat("[9] Building motif width and IC distribution plots...\n")

# --- 9a. Width distribution: all archetypes in MEME file, highlighting the
#     filter cutoff and which archetypes are included vs excluded.
pwm_all[, passed_filter := archetype %in% VALID_ARCHETYPES]

p_width_hist <- ggplot(pwm_all, aes(x = motif_width, fill = passed_filter)) +
    geom_histogram(binwidth = 1, color = "white", linewidth = 0.3) +
    geom_vline(xintercept = opt$min_width - 0.5, linetype = "dashed",
               color = "red", linewidth = 0.6) +
    annotate("text", x = opt$min_width + 0.3, y = Inf,
             label = paste0("min_width = ", opt$min_width),
             hjust = 0, vjust = 1.5, size = 3.5, color = "red") +
    scale_fill_manual(values = c("TRUE" = "#4a90d9", "FALSE" = "#d9d9d9"),
                       name = NULL,
                       labels = c("TRUE" = "Included", "FALSE" = "Excluded")) +
    scale_x_continuous(breaks = seq(0, 30, by = 2)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(title = "J. Motif Archetype Width Distribution",
         subtitle = paste0("All ", nrow(pwm_all), " archetypes; ",
                           sum(pwm_all$passed_filter), " pass width \u2265 ",
                           opt$min_width, " bp filter (",
                           sum(!pwm_all$passed_filter), " excluded)"),
         x = "Motif width (bp)", y = "Number of archetypes") +
    THEME +
    theme(legend.position = c(0.85, 0.85))

# --- 9b. Total IC distribution (for archetypes passing width filter)
p_ic_hist <- ggplot(pwm_filtered, aes(x = total_ic)) +
    geom_histogram(bins = 40, fill = "#e6550d", color = "white",
                   linewidth = 0.3) +
    {if (opt$min_ic > 0) geom_vline(xintercept = opt$min_ic, linetype = "dashed",
                                     color = "red", linewidth = 0.6)} +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(title = "K. Total Information Content Distribution",
         subtitle = paste0("Archetypes passing width filter (n = ",
                           nrow(pwm_filtered), "); range: ",
                           round(min(pwm_filtered$total_ic), 1), " \u2013 ",
                           round(max(pwm_filtered$total_ic), 1), " bits",
                           if (opt$min_ic > 0) paste0("; IC filter \u2265 ",
                                                       opt$min_ic) else ""),
         x = "Total information content (bits)", y = "Number of archetypes") +
    THEME

# --- 9c. Mean IC per position
p_ic_per_pos <- ggplot(pwm_filtered, aes(x = mean_ic_per_pos)) +
    geom_histogram(bins = 40, fill = "#756bb1", color = "white",
                   linewidth = 0.3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(title = "L. Mean IC per Position",
         subtitle = paste0("Average specificity per position (0 = uniform, 2 = perfect); ",
                           "median = ",
                           round(median(pwm_filtered$mean_ic_per_pos), 2),
                           " bits/pos"),
         x = "Mean IC per position (bits)", y = "Number of archetypes") +
    THEME

# --- 9d. Width vs IC scatter, with analysis-ready TFBSs highlighted
pwm_scatter <- merge(pwm_filtered, arch_sum[, .(archetype, n_analysis_ready)],
                      by = "archetype", all.x = TRUE)
pwm_scatter[is.na(n_analysis_ready), n_analysis_ready := 0]
pwm_scatter[, has_ready := n_analysis_ready > 0]

p_width_vs_ic <- ggplot(pwm_scatter,
                         aes(x = motif_width, y = total_ic,
                             size = pmax(n_analysis_ready, 1),
                             color = has_ready)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("FALSE" = "#bdbdbd", "TRUE" = "#d62728"),
                        name = "Has analysis-\nready TFBSs",
                        labels = c("No", "Yes")) +
    scale_size_continuous(name = "Analysis-ready\nTFBSs", range = c(1, 8),
                           trans = "sqrt",
                           breaks = c(1, 10, 50, 100, 500)) +
    scale_x_continuous(breaks = seq(0, 30, by = 2)) +
    labs(title = "M. Width vs Information Content",
         subtitle = paste0("Each point = one archetype (n = ",
                           nrow(pwm_scatter),
                           "); red = has analysis-ready TFBSs"),
         x = "Motif width (bp)", y = "Total IC (bits)") +
    THEME +
    theme(legend.position = "right")

# --- 9e. Width distribution: all archetypes vs analysis-ready archetypes
width_compare <- rbind(
    pwm_filtered[, .(archetype, motif_width, group = "All archetypes (passing filter)")],
    pwm_filtered[archetype %in% arch_sum[n_analysis_ready > 0]$archetype,
                 .(archetype, motif_width, group = "With analysis-ready TFBSs")]
)

p_width_compare <- ggplot(width_compare, aes(x = motif_width, fill = group)) +
    geom_histogram(binwidth = 1, position = "dodge", color = "white",
                   linewidth = 0.2) +
    scale_fill_manual(values = c("All archetypes (passing filter)" = "#bdbdbd",
                                  "With analysis-ready TFBSs" = "#d62728"),
                       name = NULL) +
    scale_x_continuous(breaks = seq(0, 30, by = 2)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(title = "N. Width: All Archetypes vs Those With Analysis-Ready TFBSs",
         subtitle = "Does motif width influence whether we find analysis-ready TFBSs?",
         x = "Motif width (bp)", y = "Number of archetypes") +
    THEME +
    theme(legend.position = "bottom")

# =============================================================================
# 10. PLOT 9: Global filtering funnel
# =============================================================================

cat("[10] Building global filtering funnel...\n")

total_all <- sum(arch_sum$n_total)
total_intra <- sum(arch_sum$n_intragenic)
total_nonprom <- sum(arch_sum$n_intragenic_nonpromoter)
total_acc <- sum(arch_sum$n_accessible)
total_ready <- sum(arch_sum$n_analysis_ready)
total_noatac <- sum(arch_sum$n_no_atac_constraint)

global_funnel <- data.table(
    step = c("All TFBSs", "Intragenic", "Intragenic +\nnon-promoter",
             "...+ accessible", "...+ het SNP\n(analysis-ready)"),
    count = c(total_all, total_intra, total_nonprom, total_acc, total_ready)
)
global_funnel[, step := factor(step, levels = step)]
global_funnel[, pct := count / total_all * 100]
global_funnel[, pct_label := paste0(comma(count), "\n(",
                                     sprintf("%.1f%%", pct), ")")]

p9 <- ggplot(global_funnel, aes(x = step, y = count, fill = step)) +
    geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
    geom_text(aes(label = pct_label), vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = unname(FILTER_COLORS)) +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0, 0.12))) +
    labs(title = "O. Global Filtering Funnel (All Archetypes Combined)",
         subtitle = paste0(comma(total_all), " total TFBSs across ",
                           nrow(arch_sum), " archetypes (width \u2265 ",
                           opt$min_width, " bp) \u2192 ",
                           comma(total_ready), " analysis-ready"),
         x = NULL, y = "Number of TFBSs") +
    THEME +
    theme(axis.text.x = element_text(size = 9))

# =============================================================================
# 11. Assemble and save
# =============================================================================

cat("[11] Saving PDF...\n")

pdf_path <- file.path(opt$outdir, "motif_archetype_prep.pdf")

# === CUSTOMIZE === PDF dimensions
PREP_FIG_W <- 14
PREP_FIG_H <- 12

pdf(pdf_path, width = PREP_FIG_W, height = PREP_FIG_H)

# Page 1: Global funnel + family breakdown
print(p9 / p8 + plot_layout(heights = c(1, 1.2)))

# Page 2: Per-archetype funnel (counts)
print(p1)

# Page 3: Per-archetype funnel (proportions)
print(p2)

# Page 4: ATAC constraint comparison
print(p3)

# Page 5: Top archetypes ranked + heatmap
print(p7 / p_heatmap + plot_layout(heights = c(1, 1.2)))

# Page 6: MOODS score distributions
print(p4)

# Page 7: Score comparison (all vs analysis-ready)
print(p5)

# Page 8: SNP count distribution
print(p6 + plot_spacer() + plot_layout(widths = c(2, 1)))

# Page 9 (NEW): Width and IC distributions
print((p_width_hist + p_ic_hist) / (p_ic_per_pos + p_width_compare))

# Page 10 (NEW): Width vs IC scatter + width comparison
print(p_width_vs_ic)

dev.off()

cat("\n[DONE] PDF saved: ", pdf_path, "\n")

# Also save individual PNGs
cat("[INFO] Saving individual PNGs...\n")
png_dir <- file.path(opt$outdir, "png")
dir.create(png_dir, showWarnings = FALSE)

# === CUSTOMIZE === PNG resolution
PNG_DPI <- 300

plots <- list(
    "01_global_funnel"      = p9,
    "02_family_breakdown"   = p8,
    "03_funnel_counts"      = p1,
    "04_funnel_proportions" = p2,
    "05_atac_comparison"    = p3,
    "06_top_ranked"         = p7,
    "07_passrate_heatmap"   = p_heatmap,
    "08_moods_scores"       = p4,
    "09_score_comparison"   = p5,
    "10_snp_count"          = p6,
    "11_width_distribution" = p_width_hist,
    "12_total_ic_dist"      = p_ic_hist,
    "13_ic_per_position"    = p_ic_per_pos,
    "14_width_vs_ic"        = p_width_vs_ic,
    "15_width_comparison"   = p_width_compare
)

for (nm in names(plots)) {
    w <- if (grepl("funnel|ranked|heatmap|comparison", nm)) 12 else 10
    h <- if (grepl("funnel|ranked|heatmap", nm)) 10 else 7
    ggsave(file.path(png_dir, paste0(nm, ".png")),
           plot = plots[[nm]], width = w, height = h, dpi = PNG_DPI)
}

cat("[DONE] PNGs saved: ", png_dir, "\n\n")

# Console summary
cat("================================================================\n")
cat("Quick Summary\n")
cat("================================================================\n")
cat(sprintf("  MEME archetypes parsed:       %s\n", comma(nrow(pwm_all))))
cat(sprintf("  Excluded by width < %d bp:    %s\n",
            opt$min_width, comma(sum(!pwm_all$passed_filter))))
if (opt$min_ic > 0) {
    cat(sprintf("  Excluded by IC < %.1f bits:   (applied)\n", opt$min_ic))
}
cat(sprintf("  Archetypes after PWM filter:  %s\n", comma(nrow(arch_sum))))
cat(sprintf("  Width range (included):       %d – %d bp\n",
            min(arch_sum$motif_width, na.rm = TRUE),
            max(arch_sum$motif_width, na.rm = TRUE)))
cat(sprintf("  Total TFBSs:                  %s\n", comma(total_all)))
cat(sprintf("  Intragenic:                   %s (%.1f%%)\n",
            comma(total_intra), 100*total_intra/total_all))
cat(sprintf("  Intragenic + non-promoter:    %s (%.1f%%)\n",
            comma(total_nonprom), 100*total_nonprom/total_all))
cat(sprintf("  ...+ accessible:              %s (%.1f%%)\n",
            comma(total_acc), 100*total_acc/total_all))
cat(sprintf("  ...+ het SNP (analysis-ready): %s (%.1f%%)\n",
            comma(total_ready), 100*total_ready/total_all))
cat(sprintf("  Without ATAC constraint:       %s (%.1f%%)\n",
            comma(total_noatac), 100*total_noatac/total_all))
cat(sprintf("\n  Archetypes with >=1 ready:    %s\n",
            sum(arch_sum$n_analysis_ready > 0)))
cat(sprintf("  Archetypes with >=20 ready:   %s\n",
            sum(arch_sum$n_analysis_ready >= 20)))
cat(sprintf("  Archetypes with >=50 ready:   %s\n",
            sum(arch_sum$n_analysis_ready >= 50)))
cat(sprintf("  Archetypes with >=100 ready:  %s\n",
            sum(arch_sum$n_analysis_ready >= 100)))
cat("================================================================\n")