#!/usr/bin/env Rscript
# =============================================================================
# 02_score_pwm_per_allele.R — PWM ref/alt scoring for body-SNPs
# =============================================================================
#
# For each (motif, SNP) pair where SNP falls in the motif body:
#   1. Extract motif-spanning genomic sequence in transcription direction
#   2. Build the REF version (genome as-is) and ALT version (SNP swapped)
#   3. Score both against the motif's PWM
#   4. Output: pwm_score_ref, pwm_score_alt, delta_pwm_score (=alt-ref),
#              motif_sequence_ref, motif_sequence_alt
#
# For SNPs in flank (not in motif body), all PWM columns = NA. We add a
# `pwm_scoring_applicable` boolean for downstream filtering.
#
# PWM SCORING:
#   Standard log-odds against uniform background, with pseudocount.
#   For motif with width W and per-base probabilities P[A,C,G,T] at each
#   position i, score for sequence S of length W:
#     score(S) = sum_i log2( (P[S_i, i] + ps) / (0.25 + ps) )
#   Higher = better match.
#
# WHY MAX-OVER-WINDOW INSTEAD OF FIXED-POSITION:
#   The motif BED has the motif scanned by MOODS at a specific position +
#   strand. The motif sequence is exactly (motif_start, motif_end) of width
#   matching the PWM. So the score at this specific position IS the score.
#   We don't slide; we just use the fixed window.
#
# STRAND HANDLING:
#   getSeq() respects strand; for "-" motifs, it returns reverse-complement.
#   The PWM is in motif-orientation by convention (MEME); so we score the
#   getSeq() output directly without further flipping.
#
# USAGE:
#   Rscript 02_score_pwm_per_allele.R \
#     --pair_index /path/to/AC0001_pair_index.tsv.gz \
#     --pwm_meme /path/to/consensus_pwms.meme \
#     --genome_fa /path/to/mm39.fa \
#     --motif_id AC0001 \
#     --pseudocount 0.01 \
#     --outdir /path/to/pwm_scores/
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(Biostrings)
    library(Rsamtools)
    library(data.table)
})

option_list <- list(
    make_option("--pair_index",  type = "character"),
    make_option("--pwm_meme",    type = "character"),
    make_option("--genome_fa",   type = "character"),
    make_option("--motif_id",    type = "character"),
    make_option("--pseudocount", type = "double", default = 0.01),
    make_option("--outdir",      type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

cat("============================================================\n")
cat("02_score_pwm_per_allele.R — motif:", opt$motif_id, "\n")
cat("============================================================\n\n")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Load pair index
# =============================================================================
pair_dt <- fread(opt$pair_index)
if (nrow(pair_dt) == 0) {
    cat("[INFO] Empty pair index; writing empty PWM scores file.\n")
    fwrite(data.table(), file.path(opt$outdir,
        paste0(opt$motif_id, "_pwm_scores.tsv.gz")), sep = "\t",
        compress = "gzip")
    quit(save = "no", status = 0)
}
cat("Pairs read:", nrow(pair_dt), "\n")

# =============================================================================
# Parse MEME file → list of named PWMs (one per motif)
# =============================================================================
parse_meme <- function(path) {
    lines <- readLines(path)
    # Locate MOTIF blocks
    motif_starts <- grep("^MOTIF\\s+", lines)
    if (length(motif_starts) == 0) stop("No MOTIF lines in ", path)

    pwms <- list()
    for (i in seq_along(motif_starts)) {
        ms_line <- lines[motif_starts[i]]
        # Format: "MOTIF AC0001:DLX/LHX:Homeodomain  AC0001:..."
        # The first token after "MOTIF" is the motif name in MEME.
        name_full <- strsplit(trimws(sub("^MOTIF\\s+", "", ms_line)),
                               "\\s+")[[1]][1]
        # Use prefix before ":" as motif_id (e.g., "AC0001")
        motif_id_short <- sub(":.*$", "", name_full)

        # Find the letter-probability matrix line within this MOTIF block
        block_end <- if (i < length(motif_starts)) motif_starts[i + 1] - 1
                     else length(lines)
        block <- lines[motif_starts[i]:block_end]

        lpm_idx <- grep("^letter-probability matrix", block)
        if (length(lpm_idx) == 0) {
            warning("MOTIF ", name_full, " has no letter-probability matrix; skipping.")
            next
        }
        lpm_line <- block[lpm_idx[1]]
        # Parse "w= 6" from the LPM header
        w_match <- regmatches(lpm_line, regexpr("w=\\s*[0-9]+", lpm_line))
        if (length(w_match) == 0) {
            warning("Cannot parse w= for motif ", name_full, "; skipping.")
            next
        }
        W <- as.integer(sub("w=\\s*", "", w_match))

        # Read W rows starting at lpm_idx[1] + 1
        mat_rows <- block[(lpm_idx[1] + 1):(lpm_idx[1] + W)]
        # Each row has 4 numbers (A C G T)
        mat <- do.call(rbind, lapply(mat_rows, function(x) {
            as.numeric(strsplit(trimws(x), "\\s+")[[1]])
        }))
        if (ncol(mat) != 4 || nrow(mat) != W) {
            warning("Bad matrix dims for motif ", name_full, "; skipping.")
            next
        }
        colnames(mat) <- c("A", "C", "G", "T")
        pwms[[motif_id_short]] <- list(name_full = name_full, w = W, ppm = mat)
    }
    pwms
}

cat("Parsing MEME PWMs from:", opt$pwm_meme, "\n")
pwms <- parse_meme(opt$pwm_meme)
cat("  Parsed", length(pwms), "motifs from MEME file\n")

if (!opt$motif_id %in% names(pwms)) {
    cat("[ERROR] motif_id '", opt$motif_id, "' not in MEME file.\n", sep = "")
    cat("        Available:", paste(head(names(pwms), 10), collapse = ", "),
        "...\n")
    quit(save = "no", status = 1)
}

ppm <- pwms[[opt$motif_id]]$ppm
W <- pwms[[opt$motif_id]]$w
cat("  Selected motif:", opt$motif_id, "(width =", W, ")\n\n")

# Build log-odds PWM (background uniform)
ps <- opt$pseudocount
log_odds <- log2((ppm + ps) / (0.25 + ps))

# =============================================================================
# Score each pair
# =============================================================================
cat("Scoring pairs...\n")

# Verify motif_width matches PWM width for body-SNP rows
body_idx <- which(pair_dt$snp_in_motif_body)
if (length(body_idx) > 0) {
    width_mismatch <- sum(pair_dt[body_idx]$motif_width != W)
    if (width_mismatch > 0) {
        cat("[WARN]", width_mismatch, "body pairs have motif_width !=",
            W, "(PWM width). PWM scores will be NA for those.\n")
    }
}

# Initialize output columns
pair_dt[, pwm_scoring_applicable := snp_in_motif_body & (motif_width == W)]
pair_dt[, motif_sequence_ref := NA_character_]
pair_dt[, motif_sequence_alt := NA_character_]
pair_dt[, pwm_score_ref      := NA_real_]
pair_dt[, pwm_score_alt      := NA_real_]
pair_dt[, delta_pwm_score    := NA_real_]

idx <- which(pair_dt$pwm_scoring_applicable)
if (length(idx) == 0) {
    cat("  No pairs eligible for PWM scoring (no body-SNPs of correct width).\n")
} else {
    cat("  Eligible body pairs:", length(idx), "\n")

    # Open FASTA
    fa <- FaFile(opt$genome_fa)

    # Build motif-extent GRanges for body pairs (in transcription direction)
    motif_gr <- GRanges(
        seqnames = pair_dt$motif_chrom[idx],
        ranges   = IRanges(start = pair_dt$motif_start[idx],
                           end   = pair_dt$motif_end[idx]),
        strand   = pair_dt$motif_strand[idx]
    )

    # Extract REF sequences (strand-aware: reverse-complement for "-")
    ref_seqs <- as.character(getSeq(fa, motif_gr))

    # Build ALT sequences
    # snp_position_in_motif is 1-based offset in transcription direction.
    # snp_alt is on genomic + strand; for "-" motifs we need to RC the alt
    # base so it's correct in transcription direction.
    pos_in_motif <- pair_dt$snp_position_in_motif[idx]
    alt_base_genomic <- pair_dt$snp_alt[idx]
    rc_table <- c(A = "T", T = "A", C = "G", G = "C", N = "N")
    alt_base_tx <- ifelse(
        pair_dt$motif_strand[idx] == "+",
        alt_base_genomic,
        rc_table[alt_base_genomic]
    )

    # Apply substitution character-by-character (vectorized via substr<-)
    alt_seqs <- ref_seqs
    for (j in seq_along(alt_seqs)) {
        substr(alt_seqs[j], pos_in_motif[j], pos_in_motif[j]) <- alt_base_tx[j]
    }

    pair_dt[idx, motif_sequence_ref := ref_seqs]
    pair_dt[idx, motif_sequence_alt := alt_seqs]

    # Score: walk the position columns once
    score_seqs <- function(seqs, log_odds_mat) {
        # log_odds_mat: W rows x 4 cols (A C G T)
        n <- length(seqs)
        scores <- numeric(n)
        for (j in seq_len(n)) {
            chars <- charToRaw(seqs[j])
            # Map A=1, C=2, G=3, T=4 (others -> NA)
            idx_letter <- match(rawToChar(chars, multiple = TRUE),
                                 c("A", "C", "G", "T"))
            if (any(is.na(idx_letter))) {
                scores[j] <- NA_real_
            } else {
                # Sum diagonal: log_odds_mat[1, idx[1]] + log_odds_mat[2, idx[2]] + ...
                scores[j] <- sum(log_odds_mat[cbind(seq_along(idx_letter),
                                                    idx_letter)])
            }
        }
        scores
    }

    cat("  Computing scores...\n")
    pair_dt[idx, pwm_score_ref := score_seqs(ref_seqs, log_odds)]
    pair_dt[idx, pwm_score_alt := score_seqs(alt_seqs, log_odds)]
    pair_dt[idx, delta_pwm_score := pwm_score_alt - pwm_score_ref]

    cat("  Score summaries:\n")
    cat("    pwm_score_ref:    min=",
        round(min(pair_dt[idx]$pwm_score_ref, na.rm = TRUE), 2),
        ", max=", round(max(pair_dt[idx]$pwm_score_ref, na.rm = TRUE), 2),
        ", mean=", round(mean(pair_dt[idx]$pwm_score_ref, na.rm = TRUE), 2), "\n")
    cat("    delta (alt-ref):  min=",
        round(min(pair_dt[idx]$delta_pwm_score, na.rm = TRUE), 2),
        ", max=", round(max(pair_dt[idx]$delta_pwm_score, na.rm = TRUE), 2),
        ", median |delta|=",
        round(median(abs(pair_dt[idx]$delta_pwm_score), na.rm = TRUE), 2), "\n")
    cat("    n |delta| > 1.0:", sum(abs(pair_dt[idx]$delta_pwm_score) > 1.0,
                                      na.rm = TRUE), "\n")
}

# =============================================================================
# Write output: only the new PWM columns + pair_id (for joining)
# =============================================================================
out_dt <- pair_dt[, .(
    pair_id, motif_hit_id,
    pwm_scoring_applicable,
    motif_sequence_ref, motif_sequence_alt,
    pwm_score_ref, pwm_score_alt, delta_pwm_score
)]

out_path <- file.path(opt$outdir, paste0(opt$motif_id, "_pwm_scores.tsv.gz"))
fwrite(out_dt, out_path, sep = "\t", na = "NA", compress = "gzip")
cat("\nWrote:", out_path, "\n")
cat("Dimensions:", nrow(out_dt), "rows x", ncol(out_dt), "cols\n")
