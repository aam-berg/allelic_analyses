#!/usr/bin/env Rscript
# =============================================================================
# functions.R — Helper functions for PRO-seq × TFBS pausing analysis (v3)
#
# v3: Two-encounter model
#   - Co-directional:  RNAPII reads motif consensus 5'→3'
#   - Anti-directional: RNAPII reads reverse complement 5'→3'
#
# Each encounter type combines signal from BOTH + and − strand motif
# instances, using the appropriate bigwig and flip for each.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(patchwork)
})

# =============================================================================
# 1. MEME PWM PARSING
# =============================================================================

parse_meme_pwms <- function(meme_file) {
  lines <- readLines(meme_file)
  pwms <- list()
  i <- 1
  while (i <= length(lines)) {
    if (grepl("^MOTIF ", lines[i])) {
      motif_full <- sub("^MOTIF\\s+", "", lines[i])
      motif_id <- sub(":.*", "", strsplit(motif_full, "\\s+")[[1]][1])

      i <- i + 1
      while (i <= length(lines) && !grepl("^letter-probability matrix", lines[i])) {
        i <- i + 1
      }
      if (i > length(lines)) break

      header <- lines[i]
      w <- as.integer(sub(".*w=\\s*(\\d+).*", "\\1", header))

      mat <- matrix(0, nrow = w, ncol = 4)
      colnames(mat) <- c("A", "C", "G", "T")
      for (j in 1:w) {
        i <- i + 1
        vals <- as.numeric(strsplit(trimws(lines[i]), "\\s+")[[1]])
        mat[j, ] <- vals
      }
      pwms[[motif_id]] <- mat
    }
    i <- i + 1
  }
  return(pwms)
}

pwm_consensus <- function(pwm_matrix) {
  bases <- c("A", "C", "G", "T")
  paste0(bases[apply(pwm_matrix, 1, which.max)], collapse = "")
}

revcomp_pwm <- function(pwm_matrix) {
  # Reverse complement: reverse rows, swap A<->T and C<->G columns
  rc <- pwm_matrix[nrow(pwm_matrix):1, c("T", "G", "C", "A"), drop = FALSE]
  colnames(rc) <- c("A", "C", "G", "T")
  return(rc)
}


# =============================================================================
# 2. TFBS FILTERING
# =============================================================================

filter_tfbs <- function(dt, filters) {
  d <- copy(dt)

  if (!is.null(filters$chroms_exclude)) {
    d <- d[!chrom %in% filters$chroms_exclude]
  }
  if (!is.null(filters$is_intragenic)) {
    d <- d[is_intragenic == filters$is_intragenic]
  }
  if (!is.null(filters$is_promoter)) {
    d <- d[is_promoter == filters$is_promoter]
  }
  if (isTRUE(filters$atac_accessible)) {
    d <- d[atac_4DNFIAEQI3RP_overlap == TRUE | atac_4DNFIZNPOOZN_overlap == TRUE]
  }
  if (!is.null(filters$gene_type_include)) {
    if (!is.null(d$intragenic_gene_types)) {
      pattern <- paste(filters$gene_type_include, collapse = "|")
      d <- d[!is.na(intragenic_gene_types) & grepl(pattern, intragenic_gene_types)]
    }
  }
  if (!is.null(filters$gene_type_exclude)) {
    if (!is.null(d$intragenic_gene_types)) {
      pattern <- paste(filters$gene_type_exclude, collapse = "|")
      d <- d[is.na(intragenic_gene_types) | !grepl(pattern, intragenic_gene_types)]
    }
  }
  if (!is.null(filters$genic_region_include)) {
    pattern <- paste(filters$genic_region_include, collapse = "|")
    d <- d[!is.na(genic_regions) & grepl(pattern, genic_regions)]
  }
  if (!is.null(filters$genic_region_exclude)) {
    pattern <- paste(filters$genic_region_exclude, collapse = "|")
    d <- d[is.na(genic_regions) | !grepl(pattern, genic_regions)]
  }
  if (isTRUE(filters$no_snp)) {
    d <- d[`snp_F121-9_overlap` == FALSE & `snp_BL6xCAST_overlap` == FALSE]
  }
  if (!is.null(filters$min_score)) {
    d <- d[score >= filters$min_score]
  }

  return(d)
}


# =============================================================================
# 3. SIGNAL EXTRACTION FROM BIGWIG
# =============================================================================

extract_signal_matrix <- function(bw_path, regions_gr, window_width) {
  bw <- BigWigFile(bw_path)
  si <- seqinfo(bw)

  n <- length(regions_gr)
  mat <- matrix(0, nrow = n, ncol = window_width)

  chroms <- as.character(seqnames(regions_gr))
  unique_chroms <- unique(chroms)

  for (chr in unique_chroms) {
    idx <- which(chroms == chr)
    if (length(idx) == 0) next

    if (!chr %in% seqnames(si)) {
      mat[idx, ] <- NA
      next
    }

    chr_len <- seqlengths(si)[chr]
    sub_gr <- regions_gr[idx]

    s <- pmax(start(sub_gr), 1L)
    e <- pmin(end(sub_gr), chr_len)

    for (j in seq_along(idx)) {
      if (s[j] > e[j]) {
        mat[idx[j], ] <- NA
        next
      }

      query_gr <- GRanges(chr, IRanges(s[j], e[j]))
      sig <- tryCatch(
        import(bw, selection = query_gr, as = "NumericList"),
        error = function(e) NULL
      )

      if (is.null(sig) || length(sig) == 0) {
        mat[idx[j], ] <- NA
        next
      }

      vals <- unlist(sig)
      offset_start <- s[j] - start(sub_gr)[j] + 1L
      offset_end   <- offset_start + length(vals) - 1L

      if (offset_start >= 1 && offset_end <= window_width && length(vals) > 0) {
        mat[idx[j], offset_start:offset_end] <- vals
      }
    }
  }

  return(mat)
}

extract_signal_matrix_v2 <- function(bw_path, regions_dt, flank_bp) {
  window_width <- 2L * flank_bp + 1L

  gr <- GRanges(
    seqnames = regions_dt$chrom,
    ranges   = IRanges(
      start = regions_dt$center - flank_bp,
      end   = regions_dt$center + flank_bp
    )
  )

  return(extract_signal_matrix(bw_path, gr, window_width))
}


# =============================================================================
# 4. NORMALIZATION
# =============================================================================

normalize_matrix <- function(mat, method = "minmax") {
  if (method == "none") return(mat)

  for (i in seq_len(nrow(mat))) {
    row <- mat[i, ]
    if (all(is.na(row))) next

    if (method == "minmax") {
      rmin <- min(row, na.rm = TRUE)
      rmax <- max(row, na.rm = TRUE)
      if (rmax > rmin) {
        mat[i, ] <- (row - rmin) / (rmax - rmin)
      } else {
        mat[i, ] <- 0
      }
    } else if (method == "zscore") {
      mu <- mean(row, na.rm = TRUE)
      sd_ <- sd(row, na.rm = TRUE)
      if (sd_ > 0) {
        mat[i, ] <- (row - mu) / sd_
      } else {
        mat[i, ] <- 0
      }
    }
  }
  return(mat)
}


# =============================================================================
# 5. ENCOUNTER DEFINITIONS
# =============================================================================
#
# There are exactly 2 physically distinct RNAPII-motif encounters:
#
# 1. CO-DIRECTIONAL: RNAPII reads the motif consensus sequence (5'→3')
#    - From + strand motif sites: + RNAPII (L→R), use plus bigwig, no flip
#    - From − strand motif sites: − RNAPII (R→L), use minus bigwig, flip
#
# 2. ANTI-DIRECTIONAL: RNAPII reads the reverse complement (5'→3')
#    - From + strand motif sites: − RNAPII (R→L), use minus bigwig, flip
#    - From − strand motif sites: + RNAPII (L→R), use plus bigwig, no flip
#
# Flip rule: flip whenever the bigwig is "minus" (RNAPII moves R→L in genome
# coords; we flip to show RNAPII moving L→R in the plot).
#
# In both panels, RNAPII is displayed moving LEFT → RIGHT.
# Left of motif center = RNAPII upstream (approaching).
# Right of motif center = RNAPII downstream (departing).

ENCOUNTER_DEFS <- list(
  co_directional = list(
    label       = "Co-directional",
    description = "RNAPII reads motif consensus (5'\u21923')",
    contributions = list(
      list(motif_strand = "+", bw_strand = "plus",  flip = FALSE),
      list(motif_strand = "-", bw_strand = "minus", flip = TRUE)
    ),
    pwm_fn_name = "identity"     # show original PWM (e.g. AAACC)
  ),
  anti_directional = list(
    label       = "Anti-directional",
    description = "RNAPII reads reverse complement (5'\u21923')",
    contributions = list(
      list(motif_strand = "+", bw_strand = "minus", flip = TRUE),
      list(motif_strand = "-", bw_strand = "plus",  flip = FALSE)
    ),
    pwm_fn_name = "revcomp_pwm"  # show reverse complement (e.g. GGTTT)
  )
)


# =============================================================================
# 6. MAIN EXTRACTION + AGGREGATION FUNCTION
# =============================================================================

compute_metaprofiles <- function(filtered_dt, bw_plus, bw_minus, flank_bp,
                                 norm_method = "minmax", min_instances = 20) {
  # Compute mean normalized PRO-seq metaprofiles for each encounter type.
  #
  # Each encounter combines signal from BOTH + and − strand motif instances,
  # using the appropriate bigwig file and flip operation.
  #
  # Returns a list with:
  #   $profiles : data.table with columns: encounter, position, mean_signal, se_signal, n_instances
  #   $matrices : list of combined raw signal matrices per encounter
  #   $counts   : named vector of total instance counts per encounter

  profiles_list <- list()
  matrices_list <- list()
  counts <- c()

  for (enc_name in names(ENCOUNTER_DEFS)) {
    enc <- ENCOUNTER_DEFS[[enc_name]]

    # Collect signal matrices from both strand contributions
    mat_parts <- list()

    for (contrib in enc$contributions) {
      sub_dt <- filtered_dt[strand == contrib$motif_strand]
      n_sub <- nrow(sub_dt)

      if (n_sub == 0) next

      # Compute center of each motif
      sub_dt[, center := as.integer(floor((start + end) / 2))]

      # Select bigwig
      bw_path <- if (contrib$bw_strand == "plus") bw_plus else bw_minus

      message(sprintf("    [%s] %s-strand motifs (n=%d) × %s bigwig%s",
                      enc$label, contrib$motif_strand, n_sub,
                      contrib$bw_strand,
                      ifelse(contrib$flip, " (flip)", "")))

      # Extract signal
      mat <- extract_signal_matrix_v2(bw_path, sub_dt, flank_bp)

      # Absolute value (minus bigwig values are typically negative)
      mat <- abs(mat)

      # Flip if RNAPII moves R→L (to standardize to L→R display)
      if (contrib$flip) {
        mat <- mat[, ncol(mat):1, drop = FALSE]
      }

      mat_parts[[length(mat_parts) + 1]] <- mat
    }

    # Combine matrices from both strand contributions
    if (length(mat_parts) == 0) {
      counts[enc_name] <- 0L
      matrices_list[[enc_name]] <- NULL
      next
    }

    combined_mat <- do.call(rbind, mat_parts)
    n_total <- nrow(combined_mat)
    counts[enc_name] <- n_total

    message(sprintf("    [%s] Combined: %d total instances", enc$label, n_total))

    if (n_total < min_instances) {
      message(sprintf("    [%s] Only %d instances (< %d minimum), skipping.",
                      enc$label, n_total, min_instances))
      matrices_list[[enc_name]] <- combined_mat
      next
    }

    # Store raw combined matrix
    matrices_list[[enc_name]] <- combined_mat

    # Normalize
    mat_norm <- normalize_matrix(combined_mat, method = norm_method)

    # Remove rows that are all NA
    valid_rows <- rowSums(is.na(mat_norm)) < ncol(mat_norm)
    mat_valid <- mat_norm[valid_rows, , drop = FALSE]
    n_valid <- nrow(mat_valid)

    if (n_valid < min_instances) {
      message(sprintf("    [%s] Only %d valid instances after NA removal, skipping.",
                      enc$label, n_valid))
      next
    }

    mean_sig <- colMeans(mat_valid, na.rm = TRUE)
    se_sig   <- apply(mat_valid, 2, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

    positions <- seq(-flank_bp, flank_bp, by = 1L)

    profiles_list[[enc_name]] <- data.table(
      encounter   = enc$label,
      position    = positions,
      mean_signal = mean_sig,
      se_signal   = se_sig,
      n_instances = n_valid
    )
  }

  profiles <- rbindlist(profiles_list, fill = TRUE)

  return(list(
    profiles = profiles,
    matrices = matrices_list,
    counts   = counts
  ))
}


# =============================================================================
# 7. PWM RE-DERIVATION FROM FILTERED INSTANCES
# =============================================================================

rederive_pwm <- function(filtered_dt, motif_width) {
  seqs <- filtered_dt$motif_sequence
  seqs <- seqs[!is.na(seqs) & nchar(seqs) == motif_width]

  if (length(seqs) == 0) return(NULL)

  mat <- matrix(0, nrow = motif_width, ncol = 4)
  colnames(mat) <- c("A", "C", "G", "T")

  for (i in seq_len(motif_width)) {
    bases <- substr(seqs, i, i)
    tab <- table(factor(toupper(bases), levels = c("A", "C", "G", "T")))
    mat[i, ] <- tab / sum(tab)
  }

  return(mat)
}