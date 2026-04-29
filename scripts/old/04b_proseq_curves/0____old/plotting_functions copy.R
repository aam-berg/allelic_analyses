#!/usr/bin/env Rscript
# =============================================================================
# plotting_functions.R — Plotting functions for PRO-seq × TFBS metaprofiles (v3)
#
# LAYOUT (for each figure):
#
#   +----------------------------------------------+
#   | Co-directional panel         RNAPII 5'→3' →  |
#   |   "Mean Normalized PRO-seq"                   |
#   +----------------------------------------------+
#   | PWM logo: original motif (e.g. AAACC)         |
#   +----------------------------------------------+
#   | Anti-directional panel       RNAPII 5'→3' →  |
#   |   "Mean Normalized PRO-seq"                   |
#   +----------------------------------------------+
#   | PWM logo: reverse complement (e.g. GGTTT)     |
#   +----------------------------------------------+
#
# Both panels show RNAPII moving left → right.
# Both panels use ALL motif instances (same n).
#
# profiles_dt$encounter = "Co-directional", "Anti-directional"
# names(counts)         = "co_directional", "anti_directional"
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggseqlogo)
  library(grid)
  library(gridExtra)
  library(patchwork)
  library(data.table)
})


# =============================================================================
# FONT SIZE CONSTANTS
# =============================================================================

TITLE_FONT_SIZE        <- 22
SUBTITLE_FONT_SIZE     <- 16
SUBPLOT_TITLE_SIZE     <- 16
AXIS_LABEL_SIZE        <- 14
AXIS_TICK_SIZE         <- 12
RNAPII_LABEL_SIZE      <- 13
SNP_STATS_SIZE         <- 10
INSTANCE_COUNT_SIZE    <- 12
LEGEND_TEXT_SIZE        <- 13
LEGEND_TITLE_SIZE      <- 14

ARROW_HEAD_LENGTH_CM   <- 0.5
ARROW_LINE_WIDTH_MM    <- 0.8
ARROW_LINE_LENGTH_FRAC <- 0.25


# =============================================================================
# ENCOUNTER TABLE — maps encounter names to display properties + PWM transform
# =============================================================================

ENC_TABLE <- list(
  list(enc_label   = "Co-directional",
       counts_key  = "co_directional",
       title       = "Co-directional",
       pwm_fn_name = "identity"),       # original PWM

  list(enc_label   = "Anti-directional",
       counts_key  = "anti_directional",
       title       = "Anti-directional",
       pwm_fn_name = "revcomp_pwm")     # reverse complement PWM
)


# =============================================================================
# HELPER: Apply a PWM transformation by name
# =============================================================================

apply_pwm_fn <- function(pwm_matrix, fn_name) {
  if (fn_name == "identity") return(pwm_matrix)
  fn <- match.fun(fn_name)
  fn(pwm_matrix)
}


# =============================================================================
# HELPER: Convert PWM to ggseqlogo-compatible matrix
# =============================================================================

pwm_to_seqlogo_matrix <- function(pwm_matrix) {
  if (is.null(colnames(pwm_matrix))) {
    colnames(pwm_matrix) <- c("A", "C", "G", "T")
  }
  m <- t(pwm_matrix)
  rownames(m) <- c("A", "C", "G", "T")
  m
}


# =============================================================================
# HELPER: Create sequence logo ggplot from a PWM matrix
# =============================================================================

make_logo_plot <- function(pwm_matrix) {
  mat <- pwm_to_seqlogo_matrix(pwm_matrix)
  p <- ggseqlogo(mat, method = "bits") +
    theme_void() +
    theme(
      axis.text   = element_blank(),
      axis.title  = element_blank(),
      plot.margin = margin(0, 2, 0, 2)
    )
  return(p)
}


# =============================================================================
# HELPER: Format cumulative SNP stats as multi-line string
# =============================================================================

format_snp_stats_text <- function(snp_stats) {
  if (is.null(snp_stats) || nrow(snp_stats) == 0) return("")
  lines <- character(0)
  for (i in seq_len(nrow(snp_stats))) {
    if (!is.na(snp_stats$count[i])) {
      lines <- c(lines, sprintf(
        "SNP %s: %d/%d (%.1f%%)",
        snp_stats$label[i],
        snp_stats$count[i],
        snp_stats$total[i],
        100 * snp_stats$count[i] / snp_stats$total[i]
      ))
    }
  }
  paste(lines, collapse = "\n")
}


# =============================================================================
# make_single_panel() — one metaprofile subplot
# =============================================================================

make_single_panel <- function(
    curves_list,
    counts_list,
    panel_title,
    shared_ylim = NULL,
    colors = NULL,
    show_legend = FALSE,
    snp_stats = NULL,
    motif_half_width = 6
) {

  plot_dt <- rbindlist(lapply(names(curves_list), function(nm) {
    d <- copy(curves_list[[nm]])
    d[, curve_label := nm]
    d
  }))

  if (nrow(plot_dt) == 0) {
    return(ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0(panel_title, "\nNo data"),
               size = 5, color = "grey50") +
      theme_void())
  }

  single_curve <- (length(curves_list) == 1)

  if (is.null(colors)) {
    if (single_curve) {
      colors <- setNames("black", names(curves_list)[1])
    } else {
      pal <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
      colors <- setNames(pal[seq_along(curves_list)], names(curves_list))
    }
  }

  x_lo <- min(plot_dt$position, na.rm = TRUE)
  x_hi <- max(plot_dt$position, na.rm = TRUE)
  x_span <- x_hi - x_lo

  if (!is.null(shared_ylim)) {
    y_lo <- shared_ylim[1]; y_hi <- shared_ylim[2]
  } else {
    y_lo <- min(0, min(plot_dt$mean_signal, na.rm = TRUE))
    y_hi <- max(plot_dt$mean_signal, na.rm = TRUE) * 1.08
  }
  y_span <- y_hi - y_lo

  p <- ggplot(plot_dt, aes(x = position, y = mean_signal,
                            color = curve_label, group = curve_label)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = colors)

  # Motif boundary lines
  p <- p +
    geom_vline(xintercept = -motif_half_width, linetype = "dashed",
               color = "red", linewidth = 0.6) +
    geom_vline(xintercept =  motif_half_width, linetype = "dashed",
               color = "red", linewidth = 0.6)

  p <- p +
    labs(x = "Position relative to motif center (bp)",
         y = "Mean Normalized PRO-seq",
         title = panel_title) +
    theme_bw(base_size = AXIS_TICK_SIZE) +
    theme(
      plot.title   = element_text(size = SUBPLOT_TITLE_SIZE, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = AXIS_LABEL_SIZE),
      axis.title.y = element_text(size = AXIS_LABEL_SIZE),
      axis.text.x  = element_text(size = AXIS_TICK_SIZE),
      axis.text.y  = element_text(size = AXIS_TICK_SIZE),
      legend.text  = element_text(size = LEGEND_TEXT_SIZE),
      legend.title = element_text(size = LEGEND_TITLE_SIZE),
      plot.margin  = margin(8, 8, 4, 8)
    )

  if (single_curve || !show_legend) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + theme(legend.position = "bottom") +
      guides(color = guide_legend(title = NULL))
  }

  if (!is.null(shared_ylim)) {
    p <- p + coord_cartesian(ylim = shared_ylim)
  }

  # RNAPII arrow (always pointing right)
  y_arrow <- y_hi - y_span * 0.06
  xa <- x_lo + x_span * 0.05
  xb <- xa  + x_span * ARROW_LINE_LENGTH_FRAC

  p <- p +
    annotate("segment",
      x = xa, xend = xb, y = y_arrow, yend = y_arrow,
      arrow = arrow(length = unit(ARROW_HEAD_LENGTH_CM, "cm"), type = "closed"),
      linewidth = ARROW_LINE_WIDTH_MM, color = "black") +
    annotate("text",
      x = (xa + xb) / 2, y = y_arrow + y_span * 0.04,
      label = "RNAPII 5'\u21923'",
      size = RNAPII_LABEL_SIZE / .pt,
      fontface = "bold")

  # Instance count
  if (single_curve) {
    count_text <- sprintf("n = %s", format(counts_list[[1]], big.mark = ","))
  } else {
    count_text <- paste(sapply(names(counts_list), function(nm) {
      sprintf("%s: n = %s", nm, format(counts_list[[nm]], big.mark = ","))
    }), collapse = "\n")
  }
  p <- p + annotate("text",
    x = x_hi - x_span * 0.02,
    y = y_lo + y_span * 0.08,
    label = count_text, hjust = 1, vjust = 0,
    size = INSTANCE_COUNT_SIZE / .pt, color = "grey30")

  # SNP stats
  if (!is.null(snp_stats) && nrow(snp_stats) > 0) {
    snp_text <- format_snp_stats_text(snp_stats)
    if (nchar(snp_text) > 0) {
      p <- p + annotate("text",
        x = x_lo + x_span * 0.02,
        y = y_arrow - y_span * 0.06,
        label = snp_text, hjust = 0, vjust = 1,
        size = SNP_STATS_SIZE / .pt,
        color = "grey40", lineheight = 0.9)
    }
  }

  return(p)
}


# =============================================================================
# make_encounter_plot() — 2-panel plot with PWM logos
# =============================================================================

make_encounter_plot <- function(
    profiles_dt,
    counts,
    pwm_matrix,
    motif_id,
    bw_label,
    filter_label,
    motif_half_width,
    tf_subtitle = "",
    snp_stats = NULL,
    shared_ylim = NULL
) {

  panels <- list()
  logos  <- list()

  for (enc in ENC_TABLE) {
    key <- enc$enc_label

    prof <- profiles_dt[encounter == key]

    if (nrow(prof) == 0) {
      panels[[key]] <- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste0(enc$title, "\nInsufficient data"),
                 size = 5, color = "grey50") +
        theme_void()
    } else {
      n_inst <- counts[enc$counts_key]
      if (is.na(n_inst)) n_inst <- 0

      panels[[key]] <- make_single_panel(
        curves_list      = setNames(list(prof), filter_label),
        counts_list      = setNames(list(n_inst), filter_label),
        panel_title      = enc$title,
        shared_ylim      = shared_ylim,
        colors           = setNames("black", filter_label),
        show_legend      = FALSE,
        snp_stats        = snp_stats,
        motif_half_width = motif_half_width
      )
    }

    oriented_pwm <- apply_pwm_fn(pwm_matrix, enc$pwm_fn_name)
    logos[[key]] <- make_logo_plot(oriented_pwm)
  }

  # Assemble: panel / logo / panel / logo
  k <- sapply(ENC_TABLE, function(e) e$enc_label)

  combined <- (
    panels[[k[1]]] /
    logos[[k[1]]]  /
    panels[[k[2]]] /
    logos[[k[2]]]
  ) +
    plot_layout(heights = c(4, 1, 4, 1))

  combined <- combined +
    plot_annotation(
      title    = sprintf("%s \u2014 %s \u2014 %s", motif_id, filter_label, bw_label),
      subtitle = tf_subtitle,
      theme    = theme(
        plot.title    = element_text(size = TITLE_FONT_SIZE, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = SUBTITLE_FONT_SIZE, hjust = 0.5, color = "grey30")
      )
    )

  return(combined)
}


# =============================================================================
# make_overlay_encounter() — Multi-curve overlay, 2-panel with logos
# =============================================================================

make_overlay_encounter <- function(
    profiles_list,
    counts_list,
    pwm_matrix,
    motif_id,
    bw_label,
    preset_names,
    motif_half_width,
    tf_subtitle = "",
    snp_stats_list = NULL,
    shared_ylim = NULL,
    colors = NULL
) {

  if (is.null(colors)) {
    pal <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00")
    colors <- setNames(pal[seq_along(preset_names)], preset_names)
  }

  panels <- list()
  logos  <- list()

  for (enc in ENC_TABLE) {
    key <- enc$enc_label

    curves <- list()
    cnts   <- list()
    for (pn in preset_names) {
      prof <- profiles_list[[pn]][encounter == key]
      if (nrow(prof) > 0) {
        curves[[pn]] <- prof
        n <- counts_list[[pn]][enc$counts_key]
        cnts[[pn]] <- if (!is.na(n)) n else nrow(prof)
      }
    }

    if (length(curves) == 0) {
      panels[[key]] <- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste0(enc$title, "\nNo data"),
                 size = 5, color = "grey50") +
        theme_void()
    } else {
      snp_st <- if (!is.null(snp_stats_list)) snp_stats_list[[preset_names[1]]] else NULL
      panels[[key]] <- make_single_panel(
        curves_list      = curves,
        counts_list      = cnts,
        panel_title      = enc$title,
        shared_ylim      = shared_ylim,
        colors           = colors[names(curves)],
        show_legend      = TRUE,
        snp_stats        = snp_st,
        motif_half_width = motif_half_width
      )
    }

    oriented_pwm <- apply_pwm_fn(pwm_matrix, enc$pwm_fn_name)
    logos[[key]] <- make_logo_plot(oriented_pwm)
  }

  k <- sapply(ENC_TABLE, function(e) e$enc_label)

  combined <- (
    panels[[k[1]]] /
    logos[[k[1]]]  /
    panels[[k[2]]] /
    logos[[k[2]]]
  ) +
    plot_layout(heights = c(4, 1, 4, 1))

  overlay_label <- paste(preset_names, collapse = " vs ")
  combined <- combined +
    plot_annotation(
      title    = sprintf("%s \u2014 %s \u2014 %s", motif_id, overlay_label, bw_label),
      subtitle = tf_subtitle,
      theme    = theme(
        plot.title    = element_text(size = TITLE_FONT_SIZE, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = SUBTITLE_FONT_SIZE, hjust = 0.5, color = "grey30")
      )
    )

  return(combined)
}
