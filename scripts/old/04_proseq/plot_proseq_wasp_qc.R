#!/usr/bin/env Rscript
# =============================================================================
# plot_proseq_wasp_qc.R — QC visualization for PRO-seq + WASP pipeline
# =============================================================================
#
# Generates a multi-page PDF with publication-quality QC plots covering:
#   1. Read filtering cascade (trimming → rRNA → spike-in → alignment → MAPQ)
#   2. Adapter content and survival rates
#   3. rRNA and spike-in contamination
#   4. WASP mapping bias removal
#   5. Allele balance (ref vs alt) per sample and globally
#   6. Per-SNP allele coverage distribution
#   7. Per-SNP reference fraction distribution
#
# USAGE:
#   Rscript plot_proseq_wasp_qc.R \
#       --qc_dir /path/to/proseq/qc \
#       --wasp_qc_dir /path/to/proseq/wasp/qc \
#       --allele_dir /path/to/proseq/wasp/allele_specific \
#       --outdir /path/to/output
#
# CUSTOMIZATION:
#   Search for "=== CUSTOMIZE ===" comments below to find adjustable settings
#   for font sizes, colors, figure dimensions, etc.
#
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(ggplot2)
    library(patchwork)  # for combining plots
    library(scales)     # for number formatting
})

# =============================================================================
# 0. Parse arguments
# =============================================================================

option_list <- list(
    make_option("--qc_dir", type = "character",
                default = "/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/proseq/qc",
                help = "PRO-seq QC directory [default: %default]"),
    make_option("--wasp_qc_dir", type = "character",
                default = "/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/proseq/wasp/qc",
                help = "WASP QC directory [default: %default]"),
    make_option("--allele_dir", type = "character",
                default = "/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/proseq/wasp/allele_specific",
                help = "Allele-specific output directory [default: %default]"),
    make_option("--outdir", type = "character",
                default = "/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/proseq/qc/plots",
                help = "Output directory for plots [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

cat("================================================================\n")
cat("PRO-seq + WASP QC Visualization\n")
cat("================================================================\n\n")

# =============================================================================
# === CUSTOMIZE === Global theme settings
# =============================================================================
# Adjust these to change the look of ALL plots:

BASE_SIZE     <- 12      # Base font size for axis text, labels, etc.
TITLE_SIZE    <- 14      # Plot title font size
SUBTITLE_SIZE <- 11      # Subtitle / annotation size
AXIS_TEXT_SIZE <- 10     # Axis tick labels
LEGEND_SIZE   <- 10      # Legend text

# Figure dimensions (inches) for the output PDF
FIG_WIDTH  <- 12
FIG_HEIGHT <- 9

# Color palettes
# --- For samples (4 replicates) ---
SAMPLE_COLORS <- c(
    "WT_PROseq_rep1" = "#1b9e77",
    "WT_PROseq_rep2" = "#d95f02",
    "WT_PROseq_rep3" = "#7570b3",
    "WT_PROseq_rep4" = "#e7298a"
)

# --- For alleles ---
ALLELE_COLORS <- c("ref (129S1)" = "#2166ac", "alt (CAST)" = "#b2182b")

# --- For filtering cascade categories ---
CASCADE_COLORS <- c(
    "Final reads"      = "#2166ac",
    "chrM removed"     = "#92c5de",
    "Multi-mappers"    = "#d1e5f0",
    "Unmapped (mm39)"  = "#f4a582",
    "Spike-in (dm6)"   = "#d6604d",
    "rRNA"             = "#b2182b",
    "Too short"        = "#878787",
    "Low quality"      = "#bababa"
)

# Global ggplot theme
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

# Short sample labels for x-axis
SHORT_LABELS <- c(
    "WT_PROseq_rep1" = "Rep 1",
    "WT_PROseq_rep2" = "Rep 2",
    "WT_PROseq_rep3" = "Rep 3",
    "WT_PROseq_rep4" = "Rep 4"
)

# =============================================================================
# 1. Load data
# =============================================================================

cat("[1] Loading QC data...\n")

# --- Trimming stats ---
trim_file <- file.path(opt$qc_dir, "trimming_stats.tsv")
if (!file.exists(trim_file)) stop("Trimming stats not found: ", trim_file)
trim <- fread(trim_file)
cat("  Trimming stats: ", nrow(trim), " samples\n")

# --- Alignment stats ---
align_file <- file.path(opt$qc_dir, "alignment_stats.tsv")
if (!file.exists(align_file)) stop("Alignment stats not found: ", align_file)
align <- fread(align_file)
cat("  Alignment stats: ", nrow(align), " samples\n")

# --- WASP stats ---
wasp_file <- file.path(opt$wasp_qc_dir, "wasp_stats_combined.tsv")
if (!file.exists(wasp_file)) stop("WASP stats not found: ", wasp_file)
wasp <- fread(wasp_file)
cat("  WASP stats: ", nrow(wasp), " samples\n")

# --- Merged allele counts ---
allele_file <- file.path(opt$wasp_qc_dir, "allele_counts_merged.tsv")
if (!file.exists(allele_file)) stop("Allele counts not found: ", allele_file)
allele_counts <- fread(allele_file)
cat("  Allele counts: ", nrow(allele_counts), " SNPs\n\n")

# =============================================================================
# 2. FIGURE 1: Read filtering cascade
# =============================================================================

cat("[2] Building filtering cascade plot...\n")

# Compute the number of reads lost at each step
cascade <- data.table(
    sample = align$sample,
    raw_reads = trim$raw_reads,
    too_short = trim$reads_too_short,
    # "low quality" = raw - too_short - surviving (quality-trimmed out)
    low_quality = trim$raw_reads - trim$reads_too_short - trim$reads_after_trim,
    rRNA = align$rRNA_reads,
    spike_in = align$dm6_reads,
    unmapped_mm39 = (align$mm39_input - align$mm39_aligned),
    multi_mappers = align$mm39_multi,
    chrM = align$chrM_removed,
    final = align$final_reads
)

# Reshape for stacked bar
cascade_long <- melt(cascade, id.vars = "sample",
                      measure.vars = c("final", "chrM", "multi_mappers",
                                       "unmapped_mm39", "spike_in", "rRNA",
                                       "too_short", "low_quality"),
                      variable.name = "category", value.name = "reads")

# Pretty labels
cascade_long[, category_label := fcase(
    category == "final",         "Final reads",
    category == "chrM",          "chrM removed",
    category == "multi_mappers", "Multi-mappers",
    category == "unmapped_mm39", "Unmapped (mm39)",
    category == "spike_in",      "Spike-in (dm6)",
    category == "rRNA",          "rRNA",
    category == "too_short",     "Too short",
    category == "low_quality",   "Low quality"
)]

# Order categories so "Final reads" is at the bottom of the stack
cascade_long[, category_label := factor(category_label,
    levels = rev(c("Final reads", "chrM removed", "Multi-mappers",
                   "Unmapped (mm39)", "Spike-in (dm6)", "rRNA",
                   "Too short", "Low quality")))]

cascade_long[, sample_short := SHORT_LABELS[sample]]

p1 <- ggplot(cascade_long, aes(x = sample_short, y = reads / 1e6,
                                fill = category_label)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = CASCADE_COLORS, name = "Category") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = "A. Read Filtering Cascade",
         subtitle = "Breakdown of reads lost at each processing step",
         x = NULL, y = "Reads (millions)") +
    THEME +
    theme(legend.position = "right")

# Also make a percentage version
cascade_long[, pct := reads / sum(reads) * 100, by = sample]

p1b <- ggplot(cascade_long, aes(x = sample_short, y = pct,
                                 fill = category_label)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = CASCADE_COLORS, name = "Category") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01)),
                       labels = function(x) paste0(x, "%")) +
    labs(title = "B. Read Filtering Cascade (%)",
         subtitle = "Proportion of reads at each step",
         x = NULL, y = "Percentage of raw reads") +
    THEME +
    theme(legend.position = "right")

# =============================================================================
# 3. FIGURE 2: Adapter content & survival rates
# =============================================================================

cat("[3] Building trimming summary plots...\n")

# Adapter content
trim[, sample_short := SHORT_LABELS[sample]]

p2a <- ggplot(trim, aes(x = sample_short, y = pct_adapter,
                         fill = sample)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_hline(yintercept = c(30, 60), linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    annotate("text", x = 0.5, y = 62, label = "Expected range",
             hjust = 0, size = 3, color = "grey40") +
    scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
    scale_y_continuous(limits = c(0, 75),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = "C. Adapter Content",
         subtitle = "30-60% is normal for PRO-seq (short nascent RNA fragments)",
         x = NULL, y = "Reads with adapter (%)") +
    THEME

# Survival rate
p2b <- ggplot(trim, aes(x = sample_short, y = pct_surviving,
                         fill = sample)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
    scale_y_continuous(limits = c(0, 100),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = "D. Post-Trim Survival Rate",
         subtitle = "Reads passing quality + length filters",
         x = NULL, y = "Surviving reads (%)") +
    THEME

# =============================================================================
# 4. FIGURE 3: rRNA and spike-in contamination
# =============================================================================

cat("[4] Building contamination plots...\n")

align[, sample_short := SHORT_LABELS[sample]]

p3a <- ggplot(align, aes(x = sample_short, y = rRNA_pct, fill = sample)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_hline(yintercept = c(5, 20), linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    annotate("text", x = 0.5, y = 21, label = "Typical range (5-20%)",
             hjust = 0, size = 3, color = "grey40") +
    scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
    scale_y_continuous(limits = c(0, 25),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = "E. rRNA Contamination",
         subtitle = "Pol I transcription products (filtered before genome alignment)",
         x = NULL, y = "rRNA reads (%)") +
    THEME

p3b <- ggplot(align, aes(x = sample_short, y = dm6_pct, fill = sample)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_hline(yintercept = c(1, 10), linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    annotate("text", x = 0.5, y = 11, label = "Typical range (1-10%)",
             hjust = 0, size = 3, color = "grey40") +
    scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
    scale_y_continuous(limits = c(0, 18),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = "F. Drosophila Spike-In Rate",
         subtitle = "dm6-mapping reads (used for normalization)",
         x = NULL, y = "Spike-in reads (%)") +
    THEME

# =============================================================================
# 5. FIGURE 4: Alignment rate + strand balance
# =============================================================================

cat("[5] Building alignment summary plots...\n")

p4a <- ggplot(align, aes(x = sample_short, y = mm39_pct, fill = sample)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_hline(yintercept = 70, linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    annotate("text", x = 0.5, y = 71.5, label = "Good threshold (>70%)",
             hjust = 0, size = 3, color = "grey40") +
    scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
    scale_y_continuous(limits = c(0, 100),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = "G. mm39 Alignment Rate",
         subtitle = "Percentage of reads mapping to mouse genome",
         x = NULL, y = "Alignment rate (%)") +
    THEME

p4b <- ggplot(align, aes(x = sample_short, y = strand_ratio, fill = sample)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_hline(yintercept = 1.0, linetype = "solid", color = "black",
               linewidth = 0.3) +
    geom_hline(yintercept = c(0.8, 1.2), linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
    scale_y_continuous(limits = c(0.7, 1.3),
                       expand = expansion(mult = c(0.02, 0.02))) +
    labs(title = "H. Strand Balance",
         subtitle = "Plus/minus ratio (expect ~1.0; 0.8-1.2 acceptable)",
         x = NULL, y = "+/\u2212 strand ratio") +
    THEME

# =============================================================================
# 6. FIGURE 5: Final read counts
# =============================================================================

cat("[6] Building final read count plot...\n")

p5 <- ggplot(align, aes(x = sample_short, y = final_reads / 1e6,
                         fill = sample)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_hline(yintercept = c(5, 10, 20), linetype = "dashed", color = "grey50",
               linewidth = 0.3) +
    annotate("text", x = 4.5, y = c(5, 10, 20),
             label = c("Adequate", "Good", "Excellent"),
             hjust = 1, vjust = -0.5, size = 3, color = "grey40") +
    scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = "I. Final Read Counts (Pre-WASP)",
         subtitle = "After all filtering steps",
         x = NULL, y = "Reads (millions)") +
    THEME

# =============================================================================
# 7. FIGURE 6: WASP filtering
# =============================================================================

cat("[7] Building WASP filtering plots...\n")

wasp[, sample_short := SHORT_LABELS[sample]]

# WASP removal rate (of SNP-overlapping reads)
p6a <- ggplot(wasp, aes(x = sample_short, y = snp_removed_pct, fill = sample)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_hline(yintercept = c(5, 15), linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    annotate("text", x = 0.5, y = 16, label = "Typical range (5-15%)",
             hjust = 0, size = 3, color = "grey40") +
    scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
    scale_y_continuous(limits = c(0, 20),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = "J. WASP Mapping Bias Removal",
         subtitle = "% of SNP-overlapping reads removed for mapping bias",
         x = NULL, y = "Reads removed (%)") +
    THEME

# Total read loss from WASP
p6b <- ggplot(wasp, aes(x = sample_short, y = total_removed_pct, fill = sample)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
    scale_y_continuous(limits = c(0, 5),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = "K. Total WASP Read Loss",
         subtitle = "% of all reads lost (most reads don't overlap SNPs)",
         x = NULL, y = "Total reads removed (%)") +
    THEME

# =============================================================================
# 8. FIGURE 7: Allele balance per sample
# =============================================================================

cat("[8] Building allele balance plots...\n")

# Read per-sample allele counts from the allele-specific BAMs
# (from the WASP QC report data)
allele_bam_counts <- data.table(
    sample = wasp$sample,
    sample_short = SHORT_LABELS[wasp$sample]
)

# Read from allele count files to get totals
for (i in seq_len(nrow(allele_bam_counts))) {
    s <- allele_bam_counts$sample[i]
    cf <- file.path(opt$allele_dir, paste0(s, "_allele_counts.tsv"))
    if (file.exists(cf)) {
        ac <- fread(cf)
        allele_bam_counts[i, ref_reads := sum(ac$ref_count)]
        allele_bam_counts[i, alt_reads := sum(ac$alt_count)]
    }
}

allele_bam_counts[, total := ref_reads + alt_reads]
allele_bam_counts[, ref_pct := ref_reads / total * 100]

allele_long <- melt(allele_bam_counts, id.vars = c("sample", "sample_short"),
                     measure.vars = c("ref_reads", "alt_reads"),
                     variable.name = "allele", value.name = "reads")
allele_long[, allele_label := fifelse(allele == "ref_reads",
                                       "ref (129S1)", "alt (CAST)")]

p7a <- ggplot(allele_long, aes(x = sample_short, y = reads / 1e6,
                                fill = allele_label)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    scale_fill_manual(values = ALLELE_COLORS, name = "Allele") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = "L. Allele-Specific Read Counts",
         subtitle = "Total reads assigned to each allele at het SNPs",
         x = NULL, y = "Reads (millions)") +
    THEME

# Ref fraction per sample
p7b <- ggplot(allele_bam_counts, aes(x = sample_short, y = ref_pct,
                                      fill = sample)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_hline(yintercept = 50, linetype = "solid", color = "black",
               linewidth = 0.4) +
    geom_hline(yintercept = c(45, 55), linetype = "dashed", color = "grey50",
               linewidth = 0.3) +
    scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
    scale_y_continuous(limits = c(40, 60),
                       expand = expansion(mult = c(0.02, 0.02))) +
    labs(title = "M. Global Reference Fraction",
         subtitle = "Expect ~50% after WASP; slight ref bias (~53%) is common in F1s",
         x = NULL, y = "Reference allele (%)") +
    THEME

# =============================================================================
# 9. FIGURE 8: Per-SNP coverage distribution
# =============================================================================

cat("[9] Building per-SNP coverage plots...\n")

# Filter to SNPs with at least 1 allelic read
allele_covered <- allele_counts[ref_total + alt_total > 0]
allele_covered[, total_allelic := ref_total + alt_total]

# Coverage distribution (log scale)
p8a <- ggplot(allele_covered, aes(x = total_allelic)) +
    geom_histogram(bins = 80, fill = "#2166ac", color = "white",
                   linewidth = 0.2) +
    scale_x_log10(labels = comma) +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0, 0.05))) +
    labs(title = "N. Per-SNP Allelic Coverage (Merged Across Replicates)",
         subtitle = paste0(
             comma(nrow(allele_covered)), " SNPs with \u22651 allelic read; ",
             "median = ", median(allele_covered$total_allelic), " reads"),
         x = "Total allelic reads per SNP (log scale)",
         y = "Number of SNPs") +
    THEME

# Ref fraction distribution (for well-covered SNPs)
allele_well_covered <- allele_covered[total_allelic >= 10]

p8b <- ggplot(allele_well_covered,
              aes(x = as.numeric(ref_fraction))) +
    geom_histogram(bins = 50, fill = "#2166ac", color = "white",
                   linewidth = 0.2) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "red",
               linewidth = 0.5) +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0, 0.05))) +
    labs(title = "O. Reference Allele Fraction Distribution",
         subtitle = paste0(
             "SNPs with \u226510 total reads (n = ", comma(nrow(allele_well_covered)),
             "); red line = 0.5 (no bias)"),
         x = "Reference allele fraction",
         y = "Number of SNPs") +
    THEME

# =============================================================================
# 10. FIGURE 9: Summary statistics table as a plot
# =============================================================================

cat("[10] Building summary table...\n")

# Create a text-based summary
summary_text <- paste0(
    "PRO-seq + WASP Pipeline Summary (F121-9 WT mESCs, mm39)\n",
    "─────────────────────────────────────────────────────────\n",
    sprintf("  Raw reads:         %s – %s per replicate\n",
            comma(min(trim$raw_reads)), comma(max(trim$raw_reads))),
    sprintf("  Adapter content:   %.1f – %.1f%% (normal for PRO-seq)\n",
            min(trim$pct_adapter), max(trim$pct_adapter)),
    sprintf("  rRNA content:      %.1f – %.1f%%\n",
            min(align$rRNA_pct), max(align$rRNA_pct)),
    sprintf("  Spike-in (dm6):    %.1f – %.1f%%\n",
            min(align$dm6_pct), max(align$dm6_pct)),
    sprintf("  mm39 alignment:    %.1f – %.1f%%\n",
            min(align$mm39_pct), max(align$mm39_pct)),
    sprintf("  Final reads:       %s – %s\n",
            comma(min(align$final_reads)), comma(max(align$final_reads))),
    sprintf("  Strand ratio:      %.3f – %.3f\n",
            min(align$strand_ratio), max(align$strand_ratio)),
    "─────────────────────────────────────────────────────────\n",
    sprintf("  WASP bias removal: %.1f – %.1f%% of SNP reads\n",
            min(wasp$snp_removed_pct), max(wasp$snp_removed_pct)),
    sprintf("  WASP total loss:   %.2f – %.2f%% of all reads\n",
            min(wasp$total_removed_pct), max(wasp$total_removed_pct)),
    sprintf("  Global ref frac:   %.1f%% (expect ~50%%)\n",
            mean(allele_bam_counts$ref_pct)),
    sprintf("  SNPs with ≥10x:    %s / %s covered\n",
            comma(nrow(allele_well_covered)), comma(nrow(allele_covered)))
)

p_summary <- ggplot() +
    annotate("text", x = 0, y = 0, label = summary_text,
             hjust = 0, vjust = 0.5, size = 3.8, family = "mono") +
    theme_void() +
    labs(title = "P. Pipeline Summary")

# =============================================================================
# 11. Assemble and save PDF
# =============================================================================

cat("[11] Saving PDF...\n")

pdf_path <- file.path(opt$outdir, "proseq_wasp_qc.pdf")

pdf(pdf_path, width = FIG_WIDTH, height = FIG_HEIGHT)

# Page 1: Filtering cascade
print(p1 + p1b + plot_layout(guides = "collect") &
          theme(legend.position = "bottom"))

# Page 2: Trimming + contamination
print((p2a + p2b) / (p3a + p3b))

# Page 3: Alignment + strand + read counts
print((p4a + p4b) / (p5 + plot_spacer()))

# Page 4: WASP filtering
print((p6a + p6b) / (p7a + p7b))

# Page 5: Per-SNP distributions
print((p8a + p8b) / p_summary)

dev.off()

cat("\n[DONE] PDF saved: ", pdf_path, "\n")

# Also save individual PNGs for easy access
cat("[INFO] Saving individual PNGs...\n")
png_dir <- file.path(opt$outdir, "png")
dir.create(png_dir, showWarnings = FALSE)

# === CUSTOMIZE === PNG resolution (dpi) and dimensions
PNG_DPI    <- 300
PNG_WIDTH  <- 10
PNG_HEIGHT <- 6

plots <- list(
    "01_filtering_cascade" = p1,
    "02_filtering_cascade_pct" = p1b,
    "03_adapter_content" = p2a,
    "04_survival_rate" = p2b,
    "05_rRNA_contamination" = p3a,
    "06_spikein_rate" = p3b,
    "07_alignment_rate" = p4a,
    "08_strand_balance" = p4b,
    "09_final_read_counts" = p5,
    "10_wasp_bias_removal" = p6a,
    "11_wasp_total_loss" = p6b,
    "12_allele_read_counts" = p7a,
    "13_global_ref_fraction" = p7b,
    "14_snp_coverage_dist" = p8a,
    "15_ref_fraction_dist" = p8b
)

for (nm in names(plots)) {
    ggsave(file.path(png_dir, paste0(nm, ".png")),
           plot = plots[[nm]], width = PNG_WIDTH, height = PNG_HEIGHT,
           dpi = PNG_DPI)
}

cat("[DONE] Individual PNGs saved: ", png_dir, "\n")
cat("\nAll done!\n")