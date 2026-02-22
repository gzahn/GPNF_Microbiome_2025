# SETUP ####
.libPaths(c(.libPaths(), "weka/data/lab/cheeke/R_libs"))

## Packages ####
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(readxl); packageVersion("readxl")

## Options ####
set.seed(666)
nthreads <- parallel::detectCores()-1

## Functions ####
rm_missing_samples <- function(x){
  if(!any(file.exists(x))){
    x <- x[file.exists(x)]
  }
  return(x)
}

## Load Data ####
meta <- read_xlsx("./data/metadata/clean_metadata_Palouse.xlsx")

# Clean seq filepaths (removing any samples with missing files)
bact_fwds <- meta$fp_clean_fwd_bact %>% rm_missing_samples()
bact_revs <- meta$fp_clean_rev_bact %>% rm_missing_samples()
fung_fwds <- meta$fp_clean_fwd_fung %>% rm_missing_samples() # reverse files not used for fungi
fung_fwds_raw <- meta$fp_raw_fwd_fung %>% rm_missing_samples() # reverse files not used for fungi

# stop if sequences are empty RB
stopifnot(
  length(bact_fwds) > 0,
  length(bact_revs) > 0,
  length(fung_fwds) > 0
)

## Check seq data quality ####
#rand_samples <- sample(seq_along(bact_fwds),2)
#plotQualityProfile(bact_fwds[rand_samples])
#plotQualityProfile(bact_revs[rand_samples])
#plotQualityProfile(fung_fwds[rand_samples])
#plotQualityProfile(fung_fwds_raw[rand_samples])



# ASV INFERRENCE ####

## Learn error profiles ####
bact_errF <- learnErrors(bact_fwds, 
                         multithread=nthreads, 
                         MAX_CONSIST = 10,
                         verbose = 1,
                         randomize = TRUE)
saveRDS(bact_errF,"./output/physeq/bact_errF.RDS")

bact_errR <- learnErrors(bact_revs, 
                         multithread=nthreads, 
                         MAX_CONSIST = 10,
                         verbose = 1,
                         randomize = TRUE)
saveRDS(bact_errR,"./output/physeq/bact_errR.RDS")

fung_errF <- learnErrors(fung_fwds, 
                         multithread=nthreads, 
                         MAX_CONSIST = 10,
                         verbose = 1,
                         randomize = TRUE)
saveRDS(fung_errF,"./output/physeq/fung_errF.RDS")

# plot error profiles
plotErrors(bact_errF, nominalQ=FALSE); ggsave("./output/physeq/bact_errF.png")
plotErrors(bact_errR, nominalQ=FALSE); ggsave("./output/physeq/bact_errR.png")
plotErrors(fung_errF, nominalQ=FALSE); ggsave("./output/physeq/fung_errF.png")


## Sample inference ####

### bacteria fwd ####
# dereplicate
bact_derepF <- derepFastq(bact_fwds, verbose=TRUE)
saveRDS(bact_derepF,"./output/physeq/bact_derepF.RDS")
bact_derepF <- readRDS("./output/physeq/bact_derepF.RDS")
# sample inference
bact_dadaF <- dada(bact_derepF, err=bact_errF, multithread=TRUE, selfConsist = FALSE, verbose=TRUE, pool = "pseudo")
# name and save
names(bact_dadaF) <- meta$sample_name
saveRDS(bact_dadaF,"./output/physeq/bact_dadaF.RDS")
# cleanup if it worked
if(file.exists("./output/physeq/bact_dadaF.RDS") & (length(bact_derepF) > 1)){
  rm(bact_derepF)
  file.remove("./output/physeq/bact_derepF.RDS")
}

### bacteria rev ####
# dereplicate
bact_derepR <- derepFastq(bact_revs, verbose=TRUE)
saveRDS(bact_derepR,"./output/physeq/bact_derepR.RDS")
bact_derepR <- readRDS("./output/physeq/bact_derepR.RDS")
# sample inference
bact_dadaR <- dada(bact_derepR, err=bact_errR, multithread=TRUE, selfConsist = FALSE, verbose=TRUE, pool = "pseudo")
# name and save
names(bact_dadaR) <- meta$sample_name
saveRDS(bact_dadaR,"./output/physeq/bact_dadaR.RDS")
# cleanup if it worked
if(file.exists("./output/physeq/bact_dadaR.RDS") & (length(bact_derepR) > 1)){
  rm(bact_derepR)
  file.remove("./output/physeq/bact_derepR.RDS")
}

### fungal fwd ####
# depreplicate
fung_derepF <- derepFastq(fung_fwds, verbose=TRUE)
saveRDS(fung_derepF,"./output/physeq/fung_derepF.RDS")
fung_derepF <- readRDS("./output/physeq/fung_derepF.RDS")
# sample inference
fung_dadaF <- dada(fung_derepF, err=fung_errF, multithread=TRUE, selfConsist = FALSE, verbose=TRUE, pool = "pseudo")
# name and save
names(fung_dadaF) <- meta$sample_name
saveRDS(fung_dadaF,"./output/physeq/fung_dadaF.RDS")
# cleanup if it worked
if(file.exists("./output/physeq/fung_dadaF.RDS") & (length(fung_derepF) > 1)){
  rm(fung_derepF)
  file.remove("./output/physeq/fung_derepF.RDS")
}


# MERGE READS (bact) ####
mergers <- mergePairs(bact_dadaF, bact_fwds, bact_dadaR, bact_revs, verbose=TRUE)


# BUILD ASV TABLES ####
# for bacteria
bact_seqtab <- makeSequenceTable(mergers)
# for fungi
fung_seqtab <- makeSequenceTable(fung_dadaF)

# REMOVE CHIMERAS ####
bact_seqtab.nochim <- removeBimeraDenovo(bact_seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(bact_seqtab.nochim,"./output/physeq/bact_seqtab.nochim.RDS") #save object file

fung_seqtab.nochim <- removeBimeraDenovo(fung_seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(fung_seqtab.nochim,"./output/physeq/fung_seqtab.nochim.RDS") #save object file

# save chimera stats in file
sink("./output/stats/chimera_counts.txt")
cat("Bacteria")
cat("\n")
dim(bact_seqtab.nochim)
sum(bact_seqtab.nochim)/sum(bact_seqtab)
cat("Fungi")
cat("\n")
dim(fung_seqtab.nochim)
sum(fung_seqtab.nochim)/sum(fung_seqtab)
sink(NULL)

# CLEAN UP ####

# Remove all seqs with fewer than 100 nucleotides (if any) ####
bact_keeper_esvs <- nchar(names(as.data.frame(bact_seqtab.nochim))) > 99
bact_seqtab.nochim <- bact_seqtab.nochim[,bact_keeper_esvs]
fung_keeper_esvs <- nchar(names(as.data.frame(fung_seqtab.nochim))) > 99
fung_seqtab.nochim <- fung_seqtab.nochim[,fung_keeper_esvs]

# remove newly singleton taxa
bact_seqtab.nochim <- bact_seqtab.nochim[,colSums(bact_seqtab.nochim) > 1]
fung_seqtab.nochim <- fung_seqtab.nochim[,colSums(fung_seqtab.nochim) > 1]

# REMOVE CONTAMINANTS ####
# find contam seqs
bact_contams = isContaminant(bact_seqtab.nochim, neg = as.logical(meta$neg_control), normalize = TRUE)
table(bact_contams$contaminant) # how many taxa are contaminants?
write_csv(bact_contams, file = "./output/physeq/bact_likely_contaminants.csv")
fung_contams = isContaminant(fung_seqtab.nochim, neg = as.logical(meta$neg_control), normalize = TRUE)
table(fung_contams$contaminant) # how many taxa are contaminants?
write_csv(fung_contams, file = "./output/physeq/fung_likely_contaminants.csv")

# remove contaminant sequences and control samples from both tables, respectively
bact_seqtab.nochim = bact_seqtab.nochim[,(which(bact_contams$contaminant != TRUE))]
bact_seqtab.nochim = bact_seqtab.nochim[!as.logical(meta$neg_control),]
fung_seqtab.nochim = fung_seqtab.nochim[,(which(fung_contams$contaminant != TRUE))]
fung_seqtab.nochim = fung_seqtab.nochim[!as.logical(meta$neg_control),]

# remove samples from metadata too
meta <- meta[!as.logical(meta$neg_control),]

dim(bact_seqtab.nochim)
dim(fung_seqtab.nochim)
dim(meta)

# check for empty samples
which(rowSums(bact_seqtab.nochim) == 0)
which(rowSums(fung_seqtab.nochim) == 0)

# Save "clean" sequence tables
saveRDS(bact_seqtab.nochim,"./output/physeq/bact_seqtab_for_taxonomy.RDS")
saveRDS(fung_seqtab.nochim,"./output/physeq/fung_seqtab_for_taxonomy.RDS")

