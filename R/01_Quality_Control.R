# SETUP ####

# cutadapt must be installed and in your $PATH

## Packages ####
library(readxl)
library(tidyverse) 
library(dada2) 
library(purrr)
library(Biostrings) 
library(ShortRead) 
library(parallel)
sessionInfo()

## Functions ####

# build all orientations of dna seqs
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

# Discover primer matches, regardless of orientation
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rc <- function(s) as.character(Biostrings::reverseComplement(Biostrings::DNAString(s)))

# Build adapter patterns for "within last N bases" at 3'
head_slack_patterns <- function(primer, slack = 2) {
  sapply(0:slack, function(k) paste0("^", strrep("N", k), primer))
}
tail_window_patterns <- function(primer, lastN = 10) {
  base <- rc(primer)
  sapply(0:lastN, function(k) paste0(base, strrep("N", k), "$"))
}

run_cutadapt <- function(R1_in, R2_in, R1_out, R2_out, FWD, REV,
                         lead_slack_5p = 2,              # common default
                         lead_slack_5p_R1 = NULL,        # override if needed
                         lead_slack_5p_R2 = NULL,        # override if needed
                         lastN_3p = 10,
                         err = 0.07, min_overlap = 15, min_len = 50,
                         threads = max(1, parallel::detectCores()-1)) {
  
  slack_R1 <- if (is.null(lead_slack_5p_R1)) lead_slack_5p else lead_slack_5p_R1
  slack_R2 <- if (is.null(lead_slack_5p_R2)) lead_slack_5p else lead_slack_5p_R2
  
  R1_5p <- head_slack_patterns(FWD, slack = slack_R1)  # R1 expects FWD at start
  R2_5p <- head_slack_patterns(REV, slack = slack_R2)  # R2 expects REV at start
  
  R1_3p <- tail_window_patterns(REV, lastN = lastN_3p) # R1 trims REV.rc near tail
  R2_3p <- tail_window_patterns(FWD, lastN = lastN_3p) # R2 trims FWD.rc near tail
  
  dirs <- unique(dirname(c(R1_out, R2_out)))
  invisible(lapply(dirs, function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE)))
  
  args <- c(
    "-j", as.character(threads),
    "--match-read-wildcards",
    "--error-rate", as.character(err),
    "--overlap", as.character(min_overlap),
    "--no-indels",
    "--minimum-length", as.character(min_len)
  )
  for (p in R1_5p) args <- c(args, "-g", p)
  for (p in R2_5p) args <- c(args, "-G", p)
  for (p in R1_3p) args <- c(args, "-a", p)
  for (p in R2_3p) args <- c(args, "-A", p)
  
  args <- c(args, "-o", R1_out, "-p", R2_out, R1_in, R2_in)
  system2("cutadapt", args = args)
}
# Pass-through for extra args like lastN_3p and anchor5
run_domain <- function(df, r1_in, r2_in, r1_out, r2_out, FWD, REV,
                       skip_if_exists = TRUE, ...) {
  idx <- if (skip_if_exists) {
    !(file.exists(df[[r1_out]]) & file.exists(df[[r2_out]]))
  } else rep(TRUE, nrow(df))
  if (!any(idx)) {
    message("All outputs present, nothing to do.")
    return(invisible(NULL))
  }
  extra <- list(...)
  mapply(
    FUN    = run_cutadapt,
    R1_in  = df[[r1_in]][idx],
    R2_in  = df[[r2_in]][idx],
    R1_out = df[[r1_out]][idx],
    R2_out = df[[r2_out]][idx],
    MoreArgs = c(list(FWD = FWD, REV = REV), extra)
  )
}



## LOAD DATA ####
meta <- read_xlsx("./data/metadata/clean_metadata.xlsx")

# check for missing or misnamed files
missing_inputs <- meta %>%
  transmute(sample_name,
            fung_R1_missing = !file.exists(fp_raw_fwd_fung),
            fung_R2_missing = !file.exists(fp_raw_rev_fung),
            bact_R1_missing = !file.exists(fp_raw_fwd_bact),
            bact_R2_missing = !file.exists(fp_raw_rev_bact)) %>%
  rowwise() %>%
  filter(any(c_across(-sample_name)))

if (nrow(missing_inputs) > 0) {
  print(missing_inputs)
  stop("Some raw input FASTQs are missing. Fix paths before proceeding.")
}


## Primer sequences ####
bact_F <- "GTGYCAGCMGCCGCGGTAA"   # 515F-Y
bact_R <- "GGACTACNVGGGTWTCTAAT"  # 806R-B
fung_F <- "GCATCGATGAAGAACGCAGC"  # ITS3
fung_R <- "TCCTCCGCTTATTGATATGC"  # ITS4


# all possible orientations
bact.FWD.orients <- allOrients(bact_F)
bact.REV.orients <- allOrients(bact_R)
fung.FWD.orients <- allOrients(fung_F)
fung.REV.orients <- allOrients(fung_R)




# REMOVE PRIMERS ####
# bacteria
run_domain(
  meta,
  r1_in  = "fp_raw_fwd_bact",
  r2_in  = "fp_raw_rev_bact",
  r1_out = "fp_cutadapt_fwd_bact",
  r2_out = "fp_cutadapt_rev_bact",
  FWD = bact_F, REV = bact_R,
  lead_slack_5p = 3,     # SAME slack on R1 and R2
  lastN_3p = 10
)

# fungi
run_domain(
  meta,
  r1_in  = "fp_raw_fwd_fung",
  r2_in  = "fp_raw_rev_fung",
  r1_out = "fp_cutadapt_fwd_fung",
  r2_out = "fp_cutadapt_rev_fung",
  FWD = fung_F, REV = fung_R,
  lead_slack_5p = 3,
  lastN_3p = 10
)