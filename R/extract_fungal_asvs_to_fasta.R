library(Biostrings)
library(phyloseq)
library(tidyverse)

fp <- "./output/physeq/fung_seqtab_for_taxonomy.RDS"
fung <- readRDS(fp)
fung_seqs <- colnames(fung)
names(fung_seqs) <- fung_seqs

writeXStringSet(DNAStringSet(fung_seqs),"./output/fungal_asvs.fasta")
