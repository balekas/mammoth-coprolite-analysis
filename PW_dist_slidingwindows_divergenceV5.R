#!/usr/bin/env Rscript

###############################################################################
# Sliding-window % divergence analysis for aligned DNA sequences
#
# DESCRIPTION
# -----------
# This script calculates pairwise percent divergence (100 − percent identity)
# between a reference sequence and all other sequences in a multiple-sequence
# DNA alignment (FASTA format), using overlapping sliding windows.
#
# For each window:
#   - Only alignment positions where BOTH sequences are not gaps ('-') are used
#   - Windows with too many gaps are excluded (reported as NA)
#   - Percent divergence is calculated as:
#         % divergence = 100 − (% identity)
#
# The output is written as a tidy CSV file suitable for downstream analysis
# and plotting.
#
#
# HOW TO RUN
# ----------
# From the command line:
#
#   Rscript PW_dist_slidingwindows_divergenceV5.R alignment.fasta
#
#
# INPUT
# -----
# A multiple-sequence DNA alignment in FASTA format.
# All sequences must already be aligned and of equal length.
#
#
# OUTPUT
# ------
# sliding_window_divergence.csv
#
# Columns:
#   start      - window start position (1-based)
#   end        - window end position
#   midpoint   - midpoint of the window
#   reference  - reference sequence ID
#   sequence   - query sequence ID
#   divergence - percent divergence for that window
#
#
# DEPENDENCIES
# ------------
# Biostrings   (for FASTA alignment handling)
# tidyverse    (for data manipulation and CSV output)
#
###############################################################################

# ---------------- LOAD LIBRARIES ----------------
library(Biostrings)
library(tidyverse)

# ---------------- READ INPUT ARGUMENTS ----------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript PW_dist_slidingwindows_divergence.R alignment.fasta")
}

infile <- args[1]

# Output file
outfile <- "sliding_window_divergence.csv"

# ---------------- USER-DEFINED PARAMETERS ----------------

# Sliding window settings
window_size <- 25   # window length in base pairs
step_size   <- 25    # step size between windows (controls overlap)

# Gap handling
max_gap_fraction   <- 0.3   # skip window if >30% of sites are gaps
allow_partial_last <- TRUE  # include a shorter final window?

# Reference sequence selection
# Option 1: specify the reference sequence by name
reference_id <- "Krestovka"

# Option 2: use the first sequence in the FASTA file
# reference_id <- NULL

# -------------------------------------------------------

# ---------------- READ ALIGNMENT ----------------
aln <- readDNAStringSet(infile)

# Alignment length (assumes all sequences have equal length)
aln_length <- width(aln)[1]

# Determine reference sequence
if (is.null(reference_id)) {
  reference_id <- names(aln)[1]
}

if (!reference_id %in% names(aln)) {
  stop("ERROR: Reference sequence not found in alignment.")
}

message("Reference sequence: ", reference_id)

# All other sequences will be compared to the reference
seqIDs <- setdiff(names(aln), reference_id)

# ---------------- SLIDING WINDOW ANALYSIS ----------------

# Store results as a list, then combine into a data frame
results <- list()

# Generate window start positions
starts <- seq(1, aln_length, by = step_size)

for (start in starts) {

  # Calculate window end position
  end <- start + window_size - 1

  # Handle partial last window
  if (end > aln_length) {
    if (!allow_partial_last) break
    end <- aln_length
  }

  # Extract alignment slice for this window
  slice <- subseq(aln, start = start, end = end)

  # Extract reference sequence for this window
  ref_seq <- strsplit(as.character(slice[reference_id]), "")[[1]]

  # Compare each query sequence to the reference
  for (id in seqIDs) {

    qry_seq <- strsplit(as.character(slice[id]), "")[[1]]

    # Identify positions without gaps in either sequence
    valid <- ref_seq != "-" & qry_seq != "-"

    # Fraction of positions excluded due to gaps
    gap_fraction <- 1 - (sum(valid) / length(ref_seq))

    if (gap_fraction > max_gap_fraction || sum(valid) == 0) {
      # Too many gaps or no valid sites
      divergence <- NA
    } else {
      # Calculate percent divergence
      matches <- sum(ref_seq[valid] == qry_seq[valid])
      identity <- (matches / sum(valid)) * 100
      divergence <- round(100 - identity, 2)
    }

    # Store result
    results[[length(results) + 1]] <- tibble(
      start      = start,
      end        = end,
      midpoint   = (start + end) / 2,
      reference  = reference_id,
      sequence   = id,
      divergence = divergence
    )
  }
}

# Combine results into a single data frame
df <- bind_rows(results)

# ---------------- WRITE OUTPUT ----------------
write_csv(df, outfile)

message("Analysis complete.")
message("Results written to: ", outfile)
