#!/usr/bin/env Rscript
# ==============================================================================
# merge_and_plot.R — Merge per-cluster TSVs + generate summary plots
# ==============================================================================
# Run after all SLURM array jobs complete:
#   Rscript merge_and_plot.R
# ==============================================================================

library(data.table)
library(ggplot2)
library(scales)
library(gridExtra)

# ---- Paths (edit these) ----
OUTPUT_DIR    <- "/home/alb1273/pausing_phase_project/results/motif_snp_coverage"
METADATA_FILE <- "/home/alb1273/pausing_phase_project/resources/metadata.tsv"
PLOT_DIR      <- "plots_proseq_wasp"
MERGED_FILE   <- file.path(OUTPUT_DIR, "ALL_clusters_snp_coverage.tsv")

dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

PNG_DPI <- 300
THEME <- theme_bw(base_size = 12) +
    theme(
        plot.title    = element_text(size = 14, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = 11, color = "grey40"),
        axis.text     = element_text(size = 10),
        axis.title    = element_text(size = 12),
        legend.text   = element_text(size = 10),
        legend.title  = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_blank()
    )

# ==============================================================================
# 1. MERGE
# ==============================================================================

per_cluster_files <- list.files(OUTPUT_DIR, pattern = "_snp_coverage\\.tsv$",
                                 full.names = TRUE)
# Exclude the merged file itself if it exists
per_cluster_files <- per_cluster_files[!grepl("ALL_clusters", per_cluster_files)]

cat(sprintf("[INFO] Found %d per-cluster files in %s\n",
            length(per_cluster_files), OUTPUT_DIR))

# Check for missing clusters
meta <- fread(METADATA_FILE)
expected <- sort(unique(meta$cluster))
found <- gsub("_snp_coverage\\.tsv$", "", basename(per_cluster_files))
missing <- setdiff(expected, found)

if (length(missing) > 0) {
    cat(sprintf("[WARN] %d clusters missing output:\n", length(missing)))
    cat(paste("  ", head(missing, 20), collapse = "\n"), "\n")
    if (length(missing) > 20) cat(sprintf("  ... and %d more\n", length(missing) - 20))
} else {
    cat("[INFO] All clusters accounted for.\n")
}

# Read and merge
cat("[INFO] Merging...\n")
res <- rbindlist(lapply(per_cluster_files, fread), fill = TRUE)
fwrite(res, MERGED_FILE, sep = "\t")
cat(sprintf("[SAVED] %s (%s rows, %d clusters)\n",
            MERGED_FILE, format(nrow(res), big.mark = ","),
            length(unique(res$cluster))))

# ==============================================================================
# 2. ENRICH WITH TF NAMES
# ==============================================================================

cluster_info <- unique(meta[, .(cluster, tf_name, family_name)])
cluster_info <- cluster_info[!duplicated(cluster)]
res <- merge(res, cluster_info, by = "cluster", all.x = TRUE)

# Factor ordering for distances
dist_order <- c("0bp", "10bp", "25bp", "50bp", "100bp", "200bp", "400bp")
res[, distance := factor(distance, levels = dist_order)]

MIN_DEPTHS <- c(1, 5, 10, 20, 50)

# ==============================================================================
# 3. SUMMARY PRINTOUTS
# ==============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("CLUSTERS WITH >= Y MOTIF-SNP PAIRS\n")
cat("  (intragenic, non-promoter, protein-coding/lncRNA)\n")
cat(strrep("=", 70), "\n")

motif_thresholds <- c(1, 5, 10, 20, 50, 100, 200, 500)

for (dist_name in c("0bp", "25bp", "100bp")) {
    for (atac_val in c("no_atac", "with_atac")) {
        atac_label <- ifelse(atac_val == "with_atac", " + ATAC", "")
        cat(sprintf("\n  Distance = %s, %s%s\n", dist_name,
                    ifelse(atac_val == "no_atac", "no ATAC filter", "ATAC required"),
                    ""))
        cat(sprintf("  %15s", "min_depth →"))
        cat(sprintf("%8d", MIN_DEPTHS))
        cat("\n")
        cat(sprintf("  %15s", ""))
        cat(paste(rep("--------", length(MIN_DEPTHS)), collapse = ""))
        cat("\n")

        sub <- res[distance == dist_name & atac_filter == atac_val]
        for (mt in motif_thresholds) {
            cat(sprintf("  ≥%4d motifs  ", mt))
            for (md in MIN_DEPTHS) {
                n_cl <- sum(sub[min_depth == md, n_with_snp] >= mt)
                cat(sprintf("%8d", n_cl))
            }
            cat("\n")
        }
    }
}

# ---- Top 20 at key settings ----
cat("\n")
cat(strrep("=", 70), "\n")
cat("TOP 20 CLUSTERS: SNP within motif (0bp), depth >= 10, no ATAC\n")
cat(strrep("=", 70), "\n\n")

top20 <- res[distance == "0bp" & atac_filter == "no_atac" & min_depth == 10
             ][order(-n_with_snp)][1:20]
for (i in seq_len(nrow(top20))) {
    cat(sprintf("  %3d. %-8s %-20s %5d motifs with SNP (of %5d intragenic)\n",
                i, top20$cluster[i],
                ifelse(is.na(top20$tf_name[i]), "", paste0("(", top20$tf_name[i], ")")),
                top20$n_with_snp[i], top20$n_gene_filter[i]))
}

# ==============================================================================
# 4. PLOTS
# ==============================================================================

# ========================================================================
# PLOT 1: Heatmap — # clusters with >= Y motifs, by depth × distance
# ========================================================================

for (mt in c(20, 50, 100)) {
    heat_rows <- list()
    for (atac_val in c("no_atac", "with_atac")) {
        for (dist_name in dist_order) {
            for (md in MIN_DEPTHS) {
                sub <- res[distance == dist_name & atac_filter == atac_val & min_depth == md]
                n_cl <- sum(sub$n_with_snp >= mt)
                heat_rows[[length(heat_rows) + 1]] <- data.table(
                    atac     = atac_val,
                    distance = dist_name,
                    min_depth = md,
                    n_clusters = n_cl
                )
            }
        }
    }
    heat_df <- rbindlist(heat_rows)
    heat_df[, distance := factor(distance, levels = dist_order)]
    heat_df[, min_depth_label := factor(paste0("≥", min_depth),
                                         levels = paste0("≥", MIN_DEPTHS))]
    heat_df[, atac_label := fifelse(atac == "with_atac",
                                     "ATAC-accessible only", "All (no ATAC filter)")]

    p_heat <- ggplot(heat_df, aes(x = min_depth_label, y = distance, fill = n_clusters)) +
        geom_tile(color = "white", linewidth = 0.5) +
        geom_text(aes(label = n_clusters), size = 3.5) +
        facet_wrap(~ atac_label) +
        scale_fill_gradient(low = "white", high = "#2166ac", name = "# archetypes") +
        labs(title = sprintf("Archetypes with ≥ %d Qualifying TFBS-SNP Pairs", mt),
             subtitle = "Intragenic, non-promoter, protein-coding/lncRNA motifs",
             x = "Min allelic read depth at SNP",
             y = "Max distance: motif to SNP") +
        THEME +
        theme(panel.grid = element_blank())

    fname <- sprintf("17a_cluster_snp_heatmap_ge%d.png", mt)
    ggsave(file.path(PLOT_DIR, fname), plot = p_heat,
           width = 10, height = 5, dpi = PNG_DPI, bg = "white")
    cat(sprintf("[SAVED] %s\n", fname))
}

# ========================================================================
# PLOT 2: Distribution of motif-SNP counts across clusters (histograms)
# ========================================================================

combo_df <- res[atac_filter == "no_atac" &
                distance %in% c("0bp", "25bp", "100bp") &
                min_depth %in% c(1, 10, 50)]
combo_df[, facet_label := sprintf("≤ %s, depth ≥ %d",
                                   as.character(distance), min_depth)]
combo_df[, facet_label := factor(facet_label,
    levels = unique(combo_df[order(distance, min_depth), facet_label]))]

p_hist <- ggplot(combo_df[n_with_snp > 0], aes(x = n_with_snp)) +
    geom_histogram(bins = 50, fill = "#2166ac", color = "white", linewidth = 0.2) +
    facet_wrap(~ facet_label, scales = "free", ncol = 3) +
    scale_x_log10(labels = comma) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = "Distribution of Qualifying TFBS-SNP Pairs per Archetype",
         subtitle = "Intragenic, non-promoter, protein-coding/lncRNA (no ATAC filter). Log-scale x-axis.",
         x = "# motifs with SNP at given depth (per archetype)",
         y = "# archetypes") +
    THEME

ggsave(file.path(PLOT_DIR, "17b_motif_snp_count_distributions.png"),
       plot = p_hist, width = 12, height = 8, dpi = PNG_DPI, bg = "white")
cat("[SAVED] 17b_motif_snp_count_distributions.png\n")

# ========================================================================
# PLOT 3: Cumulative curves — # clusters with >= Y motifs
# ========================================================================

cum_rows <- list()
y_range <- 1:500
for (dist_name in c("0bp", "25bp", "100bp")) {
    for (md in c(5, 10, 20)) {
        sub <- res[distance == dist_name & atac_filter == "no_atac" & min_depth == md]
        for (y in y_range) {
            cum_rows[[length(cum_rows) + 1]] <- data.table(
                distance  = dist_name,
                min_depth = md,
                min_motifs = y,
                n_clusters = sum(sub$n_with_snp >= y)
            )
        }
    }
}
cum_df <- rbindlist(cum_rows)
cum_df[, combo_label := sprintf("%s, depth ≥ %d", distance, min_depth)]

p_cum <- ggplot(cum_df, aes(x = min_motifs, y = n_clusters, color = combo_label)) +
    geom_line(linewidth = 0.8) +
    scale_x_log10(labels = comma) +
    scale_y_continuous(labels = comma) +
    labs(title = "Cumulative: Archetypes with ≥ Y Qualifying TFBS-SNP Pairs",
         subtitle = "Each line: a distance × depth combination (no ATAC filter)",
         x = "Minimum # qualifying motif-SNP pairs (Y)",
         y = "# archetypes with ≥ Y pairs",
         color = "Setting") +
    THEME +
    theme(legend.position = "right")

ggsave(file.path(PLOT_DIR, "17c_cluster_cumulative_curves.png"),
       plot = p_cum, width = 9, height = 5.5, dpi = PNG_DPI, bg = "white")
cat("[SAVED] 17c_cluster_cumulative_curves.png\n")

# ========================================================================
# PLOT 4: ATAC vs no-ATAC comparison
# ========================================================================

for (md in c(1, 10)) {
    comp_df <- dcast(res[distance == "0bp" & min_depth == md],
                      cluster + tf_name ~ atac_filter,
                      value.var = "n_with_snp")

    p_atac <- ggplot(comp_df[no_atac > 0], aes(x = no_atac, y = with_atac)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
        geom_point(alpha = 0.4, size = 1.5, color = "#2166ac") +
        scale_x_log10(labels = comma) +
        scale_y_log10(labels = comma) +
        labs(title = "Effect of ATAC Filter on Qualifying Motif-SNP Pairs",
             subtitle = sprintf("SNP within motif (0bp), depth ≥ %d. Each dot = one archetype.", md),
             x = "# qualifying pairs (no ATAC filter)",
             y = "# qualifying pairs (ATAC-accessible only)") +
        THEME

    fname <- sprintf("17d_atac_filter_comparison_depth%d.png", md)
    ggsave(file.path(PLOT_DIR, fname), plot = p_atac,
           width = 6, height = 5.5, dpi = PNG_DPI, bg = "white")
    cat(sprintf("[SAVED] %s\n", fname))
}

# ========================================================================
# PLOT 5: Top archetypes ranked bar chart
# ========================================================================

for (dist_name in c("0bp", "25bp")) {
    top_n <- 30
    top_df <- res[distance == dist_name & atac_filter == "no_atac" & min_depth == 10
                  ][order(-n_with_snp)][1:top_n]
    top_df[, label := paste0(cluster,
                              ifelse(is.na(tf_name), "", paste0(" (", tf_name, ")")))]
    top_df[, label := factor(label, levels = rev(label))]

    p_top <- ggplot(top_df, aes(x = label, y = n_with_snp)) +
        geom_bar(stat = "identity", fill = "#2166ac", width = 0.7) +
        geom_text(aes(label = format(n_with_snp, big.mark = ",")),
                  hjust = -0.1, size = 3) +
        coord_flip() +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        labs(title = sprintf("Top %d Archetypes: TFBS-SNP Pairs (%s)", top_n, dist_name),
             subtitle = sprintf("Depth ≥ 10, intragenic/non-promoter (no ATAC filter)"),
             x = NULL, y = "# qualifying motif-SNP pairs") +
        THEME

    fname <- sprintf("17e_top_archetypes_%s.png", dist_name)
    ggsave(file.path(PLOT_DIR, fname), plot = p_top,
           width = 9, height = 7, dpi = PNG_DPI, bg = "white")
    cat(sprintf("[SAVED] %s\n", fname))
}

# ========================================================================
# PLOT 6: Distance effect — how much does extending the window help?
# ========================================================================

# For each cluster at depth >= 10, show n_with_snp at each distance
dist_curve <- res[atac_filter == "no_atac" & min_depth == 10]
dist_agg <- dist_curve[, .(median_n = median(n_with_snp),
                            q25 = quantile(n_with_snp, 0.25),
                            q75 = quantile(n_with_snp, 0.75)),
                        by = distance]

p_dist <- ggplot(dist_agg, aes(x = distance, y = median_n, group = 1)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#2166ac", alpha = 0.2) +
    geom_line(linewidth = 1, color = "#2166ac") +
    geom_point(size = 3, color = "#2166ac") +
    labs(title = "Effect of Distance Window on TFBS-SNP Pair Count",
         subtitle = "Median (line) ± IQR (ribbon) across all archetypes. Depth ≥ 10, no ATAC.",
         x = "Max distance: motif to SNP",
         y = "# qualifying pairs per archetype (median)") +
    THEME

ggsave(file.path(PLOT_DIR, "17f_distance_effect.png"),
       plot = p_dist, width = 8, height = 5, dpi = PNG_DPI, bg = "white")
cat("[SAVED] 17f_distance_effect.png\n")

cat("\n[DONE] All plots saved to:", PLOT_DIR, "\n")