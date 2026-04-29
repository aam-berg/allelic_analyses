#!/usr/bin/env Rscript
# =============================================================================
# functions.R — Helper functions for PRO-seq × TFBS pausing analysis
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
  # Returns a named list: list(AC0001 = matrix, AC0002 = matrix, ...)
  # Each matrix has columns A, C, G, T and rows = positions
  lines <- readLines(meme_file)
  pwms <- list()
  i <- 1
  while (i <= length(lines)) {
    if (grepl("^MOTIF ", lines[i])) {
      # Parse motif ID from "MOTIF AC0001:DLX/LHX:Homeodomain ..."
      motif_full <- sub("^MOTIF\\s+", "", lines[i])
      motif_id <- sub(":.*", "", strsplit(motif_full, "\\s+")[[1]][1])

      # Find the letter-probability matrix line
      i <- i + 1
      while (i <= length(lines) && !grepl("^letter-probability matrix", lines[i])) {
        i <- i + 1
      }
      if (i > length(lines)) break

      # Parse w= from header
      header <- lines[i]
      w <- as.integer(sub(".*w=\\s*(\\d+).*", "\\1", header))

      # Read w rows of probabilities
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
  # Get consensus sequence from a PWM matrix
  bases <- c("A", "C", "G", "T")
  paste0(bases[apply(pwm_matrix, 1, which.max)], collapse = "")
}

revcomp_pwm <- function(pwm_matrix) {
  # Reverse complement a PWM: reverse rows, swap A<->T and C<->G columns
  rc <- pwm_matrix[nrow(pwm_matrix):1, c("T", "G", "C", "A"), drop = FALSE]
  colnames(rc) <- c("A", "C", "G", "T")
  return(rc)
}

reverse_pwm <- function(pwm_matrix) {
  # Reverse a PWM (reverse row order, keep columns)
  pwm_matrix[nrow(pwm_matrix):1, , drop = FALSE]
}

complement_pwm <- function(pwm_matrix) {
  # Complement a PWM (swap A<->T and C<->G columns, keep row order)
  comp <- pwm_matrix[, c("T", "G", "C", "A"), drop = FALSE]
  colnames(comp) <- c("A", "C", "G", "T")
  return(comp)
}


# =============================================================================
# 2. TFBS FILTERING
# =============================================================================

filter_tfbs <- function(dt, filters) {
  # Apply a filter preset (named list) to an annotated TFBS data.table.
  # Returns filtered data.table.

  d <- copy(dt)

  # Chromosome exclusion
  if (!is.null(filters$chroms_exclude)) {
    d <- d[!chrom %in% filters$chroms_exclude]
  }

  # Intragenic
  if (!is.null(filters$is_intragenic)) {
    d <- d[is_intragenic == filters$is_intragenic]
  }

  # Promoter
  if (!is.null(filters$is_promoter)) {
    d <- d[is_promoter == filters$is_promoter]
  }

  # ATAC accessibility (either dataset)
  if (isTRUE(filters$atac_accessible)) {
    d <- d[atac_4DNFIAEQI3RP_overlap == TRUE | atac_4DNFIZNPOOZN_overlap == TRUE]
  }

  # Gene type include (grep patterns, OR logic)
  if (!is.null(filters$gene_type_include)) {
    # For intragenic: check intragenic_gene_types
    # For promoter: we don't have gene_type in promoter columns directly,
    # so use intragenic_gene_types if is_intragenic, else skip
    if (!is.null(d$intragenic_gene_types)) {
      pattern <- paste(filters$gene_type_include, collapse = "|")
      d <- d[!is.na(intragenic_gene_types) & grepl(pattern, intragenic_gene_types)]
    }
  }

  # Gene type exclude
  if (!is.null(filters$gene_type_exclude)) {
    if (!is.null(d$intragenic_gene_types)) {
      pattern <- paste(filters$gene_type_exclude, collapse = "|")
      d <- d[is.na(intragenic_gene_types) | !grepl(pattern, intragenic_gene_types)]
    }
  }

  # Genic region include
  if (!is.null(filters$genic_region_include)) {
    pattern <- paste(filters$genic_region_include, collapse = "|")
    d <- d[!is.na(genic_regions) & grepl(pattern, genic_regions)]
  }

  # Genic region exclude
  if (!is.null(filters$genic_region_exclude)) {
    pattern <- paste(filters$genic_region_exclude, collapse = "|")
    d <- d[is.na(genic_regions) | !grepl(pattern, genic_regions)]
  }

  # SNP exclusion
  if (isTRUE(filters$no_snp)) {
    d <- d[`snp_F121-9_overlap` == FALSE & `snp_BL6xCAST_overlap` == FALSE]
  }

  # Minimum motif score
  if (!is.null(filters$min_score)) {
    d <- d[score >= filters$min_score]
  }

  return(d)
}

# Deduplicate TFBS to unique genomic intervals (avoid counting same position twice
# when it appears on both strands). Keep both strand entries for orientation analysis.
deduplicate_positions <- function(dt) {
  # Each position appears twice (+ and - strand). We keep both but mark them
  # as the same "site" for later grouping.
  dt[, site_id := paste(chrom, start, end, sep = "_")]
  return(dt)
}


# =============================================================================
# 3. SIGNAL EXTRACTION FROM BIGWIG
# =============================================================================

extract_signal_matrix <- function(bw_path, regions_gr, window_width) {
  # Extract per-base signal from a bigWig file for a set of genomic regions.
  #
  # Args:
  #   bw_path     : path to bigWig file
  #   regions_gr  : GRanges object with regions (all same width = window_width)
  #   window_width: expected width of each region in bp
  #
  # Returns:
  #   numeric matrix (n_regions × window_width), NAs for out-of-bounds

  bw <- BigWigFile(bw_path)
  si <- seqinfo(bw)

  n <- length(regions_gr)
  mat <- matrix(0, nrow = n, ncol = window_width)

  # Process by chromosome for efficiency
  chroms <- as.character(seqnames(regions_gr))
  unique_chroms <- unique(chroms)

  for (chr in unique_chroms) {
    idx <- which(chroms == chr)
    if (length(idx) == 0) next

    # Check chromosome exists in bigwig
    if (!chr %in% seqnames(si)) {
      mat[idx, ] <- NA
      next
    }

    chr_len <- seqlengths(si)[chr]
    sub_gr <- regions_gr[idx]

    # Clip to chromosome bounds
    s <- pmax(start(sub_gr), 1L)
    e <- pmin(end(sub_gr), chr_len)

    # Import signal for each region
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

      # Convert NumericList to vector
      vals <- unlist(sig)

      # Place into the correct positions of the row
      # Account for clipping at boundaries
      offset_start <- s[j] - start(sub_gr)[j] + 1L
      offset_end   <- offset_start + length(vals) - 1L

      if (offset_start >= 1 && offset_end <= window_width && length(vals) > 0) {
        mat[idx[j], offset_start:offset_end] <- vals
      }
    }
  }

  return(mat)
}

# More efficient batch extraction using coverage
extract_signal_matrix_v2 <- function(bw_path, regions_dt, flank_bp) {
  # Batch extraction using rtracklayer import + coverage approach.
  #
  # Args:
  #   bw_path    : path to bigWig file
  #   regions_dt : data.table with chrom, center columns
  #   flank_bp   : bp to extend on each side of center
  #
  # Returns:
  #   numeric matrix (n_regions × (2*flank_bp + 1))

  window_width <- 2L * flank_bp + 1L

  # Build query GRanges
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
  # Normalize each row (instance) of a signal matrix.
  #
  # Methods:
  #   "minmax"  : scale each row to [0, 1]
  #   "zscore"  : z-score per row
  #   "none"    : no normalization

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
# 5. ORIENTATION LOGIC
# =============================================================================
#
# For a motif with PWM consensus reading e.g. AAACC (5'→3' on + strand):
#
# There are 4 encounter orientations defined by (motif_strand × PRO-seq_strand):
#
# ┌─────────────────────────────────────────────────────────────────────────────┐
# │ Motif     PRO-seq    What RNAPII         Orientation   RNAPII dir   Flip   │
# │ strand    strand     "reads" on          label         in plot      signal │
# │                      + strand (in                                          │
# │                      RNAPII direction)                                     │
# ├─────────────────────────────────────────────────────────────────────────────┤
# │   +       + (L→R)    AAACC (L→R)         Default       L→R          No     │
# │   +       - (R→L)    CCAAA (R→L→flip)    Reverse       L→R          Yes    │
# │   -       - (R→L)    GGTTT (R→L)         RevComp       L→R          Yes    │
# │   -       + (L→R)    TTTGG (L→R→flip)    Complement    L→R          No     │
# └─────────────────────────────────────────────────────────────────────────────┘
#
# Explanation:
# - Motif on + strand: genomic + strand reads AAACC left→right.
# - Motif on - strand: genomic + strand reads GGTTT left→right (the revcomp of PWM).
#
# For the PLOT, we want:
# - Upper left:  Default    — motif shown as AAACC, RNAPII arrow →
# - Upper right: Reverse    — motif shown as CCAAA, RNAPII arrow →
# - Lower left:  RevComp    — motif shown as GGTTT, RNAPII arrow ←
# - Lower right: Complement — motif shown as TTTGG, RNAPII arrow ←
#
# "Flip" means we reverse the signal vector so that RNAPII's direction of travel
# appears consistently in the plot:
# - Upper row: RNAPII always shown moving L→R (5'→3')
# - Lower row: RNAPII always shown moving R→L (5'→3')
#
# For "Reverse" (motif+, proseq-): RNAPII is actually R→L, so we flip to show L→R.
#   After flip, the motif (AAACC) appears reversed (CCAAA). ✓
# For "Complement" (motif-, proseq+): RNAPII is actually L→R, so we flip to show R→L.
#   After flip, GGTTT on + strand appears as TTTGG. ✓

ORIENTATION_DEFS <- list(
  default = list(
    label        = "Default",
    motif_strand = "+",
    proseq_bw    = "plus",    # which bigwig strand to use
    flip_signal  = FALSE,
    rnapii_arrow = "right",   # arrow direction in plot
    quadrant     = "UL",
    description  = "RNAPII co-directional with motif (same strand, same direction)"
  ),
  reverse = list(
    label        = "Reverse",
    motif_strand = "+",
    proseq_bw    = "minus",
    flip_signal  = TRUE,
    rnapii_arrow = "right",
    quadrant     = "UR",
    description  = "RNAPII anti-directional on same strand (approaches motif from 3' end)"
  ),
  revcomp = list(
    label        = "RevComp",
    motif_strand = "-",
    proseq_bw    = "minus",
    flip_signal  = FALSE,
    rnapii_arrow = "left",
    quadrant     = "LL",
    description  = "RNAPII on opposite strand, approaching from motif 3' side"
  ),
  complement = list(
    label        = "Complement",
    motif_strand = "-",
    proseq_bw    = "plus",
    flip_signal  = TRUE,
    rnapii_arrow = "left",
    quadrant     = "LR",
    description  = "RNAPII on opposite strand, approaching from motif 5' side"
  )
)


# =============================================================================
# 6. MAIN EXTRACTION + AGGREGATION FUNCTION
# =============================================================================

compute_metaprofiles <- function(filtered_dt, bw_plus, bw_minus, flank_bp,
                                 norm_method = "minmax", min_instances = 20) {
  # For a filtered TFBS data.table and a pair of bigWig files, compute mean

  # normalized PRO-seq metaprofiles for each of the 4 orientations.
  #
  # Returns a list with:
  #   $profiles : data.table with columns: orientation, position, mean_signal, se_signal, n
  #   $matrices : list of raw signal matrices per orientation (for re-plotting)
  #   $counts   : named vector of instance counts per orientation

  results <- list()
  profiles_list <- list()
  matrices_list <- list()
  counts <- c()

  for (ori_name in names(ORIENTATION_DEFS)) {
    ori <- ORIENTATION_DEFS[[ori_name]]

    # Subset to the correct motif strand
    sub_dt <- filtered_dt[strand == ori$motif_strand]

    n_instances <- nrow(sub_dt)
    counts[ori_name] <- n_instances

    if (n_instances < min_instances) {
      message(sprintf("  [%s] Only %d instances (< %d minimum), skipping.",
                      ori$label, n_instances, min_instances))
      matrices_list[[ori_name]] <- NULL
      next
    }

    # Compute center of each motif
    sub_dt[, center := as.integer(floor((start + end) / 2))]

    # Select bigwig
    bw_path <- if (ori$proseq_bw == "plus") bw_plus else bw_minus

    message(sprintf("  [%s] Extracting signal for %d instances from %s ...",
                    ori$label, n_instances, basename(bw_path)))

    # Extract signal matrix
    mat <- extract_signal_matrix_v2(bw_path, sub_dt, flank_bp)

    # For minus-strand bigwigs, values are typically negative; take absolute value
    mat <- abs(mat)

    # Flip if needed
    if (ori$flip_signal) {
      mat <- mat[, ncol(mat):1, drop = FALSE]
    }

    # Store raw matrix
    matrices_list[[ori_name]] <- mat

    # Normalize
    mat_norm <- normalize_matrix(mat, method = norm_method)

    # Compute mean and SE across instances (column-wise)
    # Remove rows that are all NA
    valid_rows <- rowSums(is.na(mat_norm)) < ncol(mat_norm)
    mat_valid <- mat_norm[valid_rows, , drop = FALSE]
    n_valid <- nrow(mat_valid)

    if (n_valid < min_instances) {
      message(sprintf("  [%s] Only %d valid instances after NA removal, skipping.",
                      ori$label, n_valid))
      next
    }

    mean_sig <- colMeans(mat_valid, na.rm = TRUE)
    se_sig   <- apply(mat_valid, 2, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

    window_width <- 2L * flank_bp + 1L
    positions <- seq(-flank_bp, flank_bp, by = 1L)

    profiles_list[[ori_name]] <- data.table(
      orientation = ori$label,
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
  # Re-derive PWM from actual motif_sequence of filtered instances.
  # Returns a PWM matrix (positions × 4 bases).

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

