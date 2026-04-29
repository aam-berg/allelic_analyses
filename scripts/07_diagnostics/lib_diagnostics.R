#!/usr/bin/env Rscript
# =============================================================================
# lib_diagnostics.R — Shared utilities for the 07_diagnostics pipeline
# =============================================================================
#
# Sourced by every R script in 07_diagnostics/. Provides:
#   - parse_meme_pwms(): minimal MEME parser
#   - pwm_consensus(), revcomp_pwm()
#   - get_archetype_pwm_widths(): map motif_id -> PWM width
#   - load_tf_metadata(): TF mappings
#   - load_gene_expression(): gene-level TPM
#   - compute_max_tpm_per_archetype(): mean-of-reps then max-of-mapped-TFs
#   - DIAG_THEME: consistent ggplot theme
#   - safe_fread(): silent fread with NULL on missing
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(scales)
})

# -----------------------------------------------------------------------------
# MEME PWM parsing
# -----------------------------------------------------------------------------
parse_meme_pwms <- function(meme_file) {
    lines <- readLines(meme_file)
    pwms <- list()
    motif_starts <- grep("^MOTIF\\s+", lines)
    for (i in seq_along(motif_starts)) {
        ms <- motif_starts[i]
        block_end <- if (i < length(motif_starts)) motif_starts[i + 1] - 1 else length(lines)
        block <- lines[ms:block_end]

        full_name <- strsplit(trimws(sub("^MOTIF\\s+", "", lines[ms])), "\\s+")[[1]][1]
        motif_id  <- sub(":.*$", "", full_name)

        lpm <- grep("^letter-probability matrix", block)
        if (length(lpm) == 0) next
        w <- as.integer(sub(".*w=\\s*([0-9]+).*", "\\1", block[lpm[1]]))
        if (is.na(w)) next

        rows <- block[(lpm[1] + 1):(lpm[1] + w)]
        mat <- do.call(rbind, lapply(rows, function(x)
            as.numeric(strsplit(trimws(x), "\\s+")[[1]])))
        if (ncol(mat) != 4 || nrow(mat) != w) next
        colnames(mat) <- c("A", "C", "G", "T")
        pwms[[motif_id]] <- list(name_full = full_name, width = w, ppm = mat)
    }
    pwms
}

pwm_consensus <- function(ppm) {
    bases <- c("A", "C", "G", "T")
    paste0(bases[apply(ppm, 1, which.max)], collapse = "")
}

revcomp_pwm <- function(ppm) {
    rc <- ppm[nrow(ppm):1, c("T", "G", "C", "A"), drop = FALSE]
    colnames(rc) <- c("A", "C", "G", "T")
    rc
}

get_archetype_pwm_widths <- function(meme_file) {
    pwms <- parse_meme_pwms(meme_file)
    data.table(
        motif_id = names(pwms),
        motif_width = vapply(pwms, function(x) x$width, integer(1))
    )
}

# -----------------------------------------------------------------------------
# TF metadata
# -----------------------------------------------------------------------------
load_tf_metadata <- function(metadata_file) {
    if (!file.exists(metadata_file)) {
        warning("TF metadata file not found: ", metadata_file)
        return(data.table())
    }
    md <- fread(metadata_file)
    # Standard columns: motif_id, cluster, source_id, tf_name, family_name, motif_type, PMID
    # 'cluster' is the archetype ID (e.g. AC0001) — this is what we join on.
    md
}

# -----------------------------------------------------------------------------
# Gene expression
# -----------------------------------------------------------------------------
load_gene_expression <- function(expr_file) {
    if (!file.exists(expr_file)) {
        warning("Gene expression file not found: ", expr_file)
        return(data.table())
    }
    fread(expr_file)
}

#' For each archetype, compute max TPM across mapped TFs.
#' TPM = mean of replicates (already in tpm_mean column); max across TFs in cluster.
#'
#' @return data.table(archetype, n_mapped_tfs, max_tpm_across_mapped, max_tf_name)
compute_max_tpm_per_archetype <- function(metadata_dt, expression_dt) {
    if (nrow(metadata_dt) == 0 || nrow(expression_dt) == 0) {
        return(data.table())
    }
    # metadata 'cluster' = archetype ID; 'tf_name' = the TF (matches gene_name)
    md <- metadata_dt[, .(motif_id = cluster, tf_name)]
    md <- unique(md[!is.na(tf_name) & tf_name != ""])

    # Join on gene_name. Multiple genes may match a TF name; take the max.
    expr <- expression_dt[, .(gene_name, tpm_mean)]
    # If duplicate gene_names exist, take max
    expr <- expr[, .(tpm_mean = max(tpm_mean, na.rm = TRUE)), by = gene_name]

    joined <- expr[md, on = c(gene_name = "tf_name"), allow.cartesian = TRUE]
    setnames(joined, "gene_name", "tf_name")

    out <- joined[, .(
        n_mapped_tfs = .N,
        n_with_tpm = sum(!is.na(tpm_mean)),
        max_tpm_across_mapped = if (all(is.na(tpm_mean))) NA_real_ else max(tpm_mean, na.rm = TRUE),
        max_tf_name = if (all(is.na(tpm_mean))) NA_character_
                      else tf_name[which.max(replace(tpm_mean, is.na(tpm_mean), -Inf))]
    ), by = motif_id]
    out
}

# -----------------------------------------------------------------------------
# Plot theme — compact and readable
# -----------------------------------------------------------------------------
DIAG_THEME <- theme_bw(base_size = 11) +
    theme(
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        panel.grid.minor = element_blank()
    )

# Save plot at consistent dimensions (compact + readable)
save_diag_plot <- function(p, path, width = 8, height = 5, dpi = 200) {
    ggsave(filename = path, plot = p, width = width, height = height, dpi = dpi)
    # Also save PNG version alongside
    if (grepl("\\.pdf$", path)) {
        png_path <- sub("\\.pdf$", ".png", path)
        ggsave(filename = png_path, plot = p, width = width, height = height, dpi = dpi)
    }
}

# -----------------------------------------------------------------------------
# Resilient fread (returns empty DT instead of erroring on missing files)
# -----------------------------------------------------------------------------
safe_fread <- function(path, ...) {
    if (!file.exists(path)) return(data.table())
    tryCatch(fread(path, ...), error = function(e) {
        warning("safe_fread failed on ", path, ": ", conditionMessage(e))
        data.table()
    })
}

# -----------------------------------------------------------------------------
# Convenience: count unique bp covered by a set of (chrom, start_1b, end_1b) rows
# (assumes start/end are 1-based inclusive). Returns total bp after reduce.
# -----------------------------------------------------------------------------
count_unique_bp <- function(chrom, start_1b, end_1b) {
    if (length(chrom) == 0) return(0L)
    suppressPackageStartupMessages(library(GenomicRanges))
    gr <- GRanges(seqnames = chrom,
                  ranges = IRanges(start = start_1b, end = end_1b))
    sum(width(reduce(gr)))
}

# -----------------------------------------------------------------------------
# Filter cascade specification (consumed by 01_build_summary.R)
# -----------------------------------------------------------------------------
# Each step: list(label, fn)
# fn takes a data.table and returns the row indices passing the filter.
# Steps are applied cumulatively.
build_cascade_steps <- function(standard_chroms, gene_types_include,
                                  expression_tpm_min,
                                  apply_expression_filter = TRUE) {
    steps <- list(
        list(
            label = "all",
            fn = function(dt) seq_len(nrow(dt))
        ),
        list(
            label = "standard_chroms",
            fn = function(dt) which(dt$chrom %in% standard_chroms)
        ),
        list(
            label = "intragenic_sense",
            fn = function(dt) which(isTRUE_safe(dt$intragenic_sense))
        ),
        list(
            label = "pc_or_lncrna",
            fn = function(dt) {
                gt <- dt$intragenic_gene_types_sense
                pat <- paste(gene_types_include, collapse = "|")
                which(!is.na(gt) & grepl(pat, gt))
            }
        )
    )
    if (apply_expression_filter) {
        steps <- c(steps, list(
            list(
                label = "expressed",
                fn = function(dt) {
                    tpm <- dt$expression_tpm_max_sense
                    which(!is.na(tpm) & tpm >= expression_tpm_min)
                }
            )
        ))
    }
    steps <- c(steps, list(
        list(
            label = "not_promoter",
            fn = function(dt) which(!isTRUE_safe(dt$is_promoter_sense))
        ),
        list(
            label = "atac_accessible",
            fn = function(dt) which(isTRUE_safe(dt$atac_consensus_overlap))
        )
    ))
    steps
}

# Vectorized "is TRUE" that handles NA, character, logical
isTRUE_safe <- function(x) {
    if (is.logical(x)) {
        ifelse(is.na(x), FALSE, x)
    } else {
        # Handle character "TRUE"/"FALSE" or NA
        x_chr <- as.character(x)
        ifelse(is.na(x_chr), FALSE, x_chr %in% c("TRUE", "true", "T", "1"))
    }
}
