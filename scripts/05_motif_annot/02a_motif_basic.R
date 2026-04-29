# =============================================================================
# 02a_motif_basic.R — Sequence extraction + palindromicity
# =============================================================================
#
# OUTPUTS (added to motif_dt):
#   motif_sequence         : Extracted DNA sequence in motif's strand orientation
#   palindromicity_score   : Fraction of positions matching reverse complement
#
# STRAND HANDLING:
#   getSeq() respects the strand of the input GRanges and returns the
#   reverse-complemented sequence for "-" strand motifs. So motif_sequence
#   is always the "as RNAPII would read it" version, which is what we want.
#
# PALINDROMICITY:
#   For a sequence S, palindromicity = mean(S[i] == revcomp(S)[i] over i).
#   1.0 = fully palindromic, 0.0 = no palindromic positions.
#
# THIS MODULE IS UNSTRANDED in the sense that there's no separate
# _sense / _antisense version — the motif sequence is intrinsic to the
# (motif, orientation) pair. If you flip orientation, motif_sequence
# becomes the reverse complement of the current value, which means you
# must re-extract — so the column tracks the "as recorded" orientation only.
# =============================================================================

annotate_basic <- function(motif_dt, motif_gr, resources) {

    # getSeq is strand-aware: for "-" strand motifs, returns reverse-complement
    seqs <- getSeq(resources$genome_fa, motif_gr)
    motif_dt[, motif_sequence := as.character(seqs)]

    # Palindromicity: fraction of positions matching their reverse complement
    motif_dt[, palindromicity_score := compute_palindromicity(motif_sequence)]

    cat("    motif_sequence: extracted (strand-aware via getSeq)\n")
    cat("    palindromicity_score: mean=",
        round(mean(motif_dt$palindromicity_score, na.rm = TRUE), 3),
        " n_perfect_pal=", sum(motif_dt$palindromicity_score == 1, na.rm = TRUE),
        "/", nrow(motif_dt), "\n", sep = "")

    invisible(motif_dt)
}

compute_palindromicity <- function(seq_strings) {
    # Vectorized over a character vector of motif sequences
    if (length(seq_strings) == 0) return(numeric(0))
    dna <- DNAStringSet(seq_strings)
    rc  <- reverseComplement(dna)
    vapply(seq_along(dna), function(i) {
        s <- strsplit(as.character(dna[[i]]), "")[[1]]
        r <- strsplit(as.character(rc[[i]]),  "")[[1]]
        sum(s == r) / length(s)
    }, numeric(1))
}
