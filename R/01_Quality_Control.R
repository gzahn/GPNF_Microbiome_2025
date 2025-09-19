# SETUP ####

# cutadapt, itsxpress, and vsearch must be installed and in your $PATH

#see DOI: 10.21769/BioProtoc.4395 for pipeline information

## Packages ####
library(readxl)
library(tidyverse) 
library(dada2) 
library(purrr)
library(Biostrings) 
library(ShortRead) 
library(parallel)
sessionInfo()
set.seed(666)

nthreads <- max(1, parallel::detectCores()-1)

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
                         threads = max(1, nthreads)) {
  
  slack_R1 <- if (is.null(lead_slack_5p_R1)) lead_slack_5p else lead_slack_5p_R1
  slack_R2 <- if (is.null(lead_slack_5p_R2)) lead_slack_5p else lead_slack_5p_R2
  
  R1_5p <- head_slack_patterns(FWD, slack = slack_R1)  # R1 expects FWD at start
  R2_5p <- head_slack_patterns(REV, slack = slack_R2)  # R2 expects REV at start
  
  R1_3p <- tail_window_patterns(REV, lastN = lastN_3p) # R1 trims REV.rc near tail
  R2_3p <- tail_window_patterns(FWD, lastN = lastN_3p) # R2 trims FWD.rc near tail
  
  dirs <- unique(dirname(c(R1_out, R2_out)))
  invisible(lapply(dirs, function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE)))
  
  args <- c(
    "-j", as.character(nthreads),
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



## Load Data ####
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
bact_F <- "GTGYCAGCMGCCGCGGTAA"   # 515F-B
bact_R <- "GGACTACNVGGGTWTCTAAT"  # 806R-B
fung_F <- "GTGAATCATCGAATCTTTGAA"  # ITS86
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

# EXTRACT ITS2 REGION FROM FUNGAL READS ####

# skipping reverse reads due to low quality

# pick path vars for ease
fung_in <- meta$fp_cutadapt_fwd_fung
fung_out <- meta$fp_itsx_fwd_fung


# set itsxpress system path
###### may need to update if itsxpress not in your $PATH
# example:
# itsxpress_path <- "/home/gzahn/.local/bin/itsxpress"
itsxpress_path <- "itsxpress"



# run for-loop, printing commands to file and keeping log files
for(i in seq_along(fung_in)){
  itsxpress <- paste0(itsxpress_path,
                      " --fastq ",fung_in[i],
                      " --outfile ",fung_out[i],
                      " --region ITS2",
                      " --taxa Fungi",
                      " --threads ",nthreads,
                      " --log ",fung_out[i],".log",
                      " --single_end")
  # write commands to file
  sink("./output/itsxpress_commands.sh",append = TRUE)
  cat(itsxpress,"\n")
  sink(NULL)
  
  system(command = itsxpress)
}


# FILTER AND TRIM FOR DADA2 ####

# Assign filepaths for bact and fung
bact_in_f <- meta$fp_cutadapt_fwd_bact
bact_in_r <- meta$fp_cutadapt_rev_bact
bact_out_f <- meta$fp_clean_fwd_bact
bact_out_r <- meta$fp_clean_rev_bact

fung_in <- meta$fp_itsx_fwd_fung
fung_out <- meta$fp_clean_fwd_fung


bact_ft_out <- filterAndTrim(bact_in_f, bact_out_f, bact_in_r, bact_out_r, # input and output file names as denoted above
                          maxN = 0, # uncalled bases are currently not supported in dada2
                          maxEE=c(3,3), # refers to the maximum expected errors allowed
                          truncQ=2, # special value denoting "end of good quality sequence" (optional)
                          rm.phix=TRUE, # automatically remove PhiX spike-in reads from sequencing center
                          compress=TRUE, # compress output files with gzip
                          multithread=nthreads) # On Windows set multithread=FALSE


fung_ft_out <- filterAndTrim(fung_in, fung_out, # input and output file names as denoted above
                          maxN=0, # uncalled bases are currently not supported in dada2
                          maxEE=3, # refers to the maximum expected errors allowed
                          truncQ=2, # special value denoting "end of good quality sequence" (optional)
                          rm.phix=TRUE, # automatically remove PhiX spike-in reads from sequencing center
                          compress=TRUE, # compress output files with gzip
                          multithread=nthreads) # On Windows set multithread=FALSE

# some samples may have no reads pass the filter. Be aware.

saveRDS(bact_ft_out,"./output/bact_filt_out.RDS")
saveRDS(fung_ft_out,"./output/fung_filt_out.RDS")


# reload point
bact_ft_out <- readRDS("./output/bact_filt_out.RDS")
fung_ft_out <- readRDS("./output/fung_filt_out.RDS")


# quickly inspect filtration outcomes
bact_ft_out; fung_ft_out

bact_ft_out %>% colSums()

fung_ft_out %>% 
  as.data.frame() %>% 
  mutate(domain="Fungi") %>% 
  full_join(
    bact_ft_out %>% 
      as.data.frame() %>% 
      mutate(domain="Bacteria")
  ) %>% 
  pivot_longer(-domain,names_prefix = "reads.") %>% 
    ggplot(aes(y=value,color=name)) +
  geom_density(linewidth=2) +
  facet_wrap(~domain) +
  theme_bw() +
  labs(y="Read counts")

