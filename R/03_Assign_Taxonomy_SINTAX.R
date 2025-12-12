# SETUP ####

# run sintax assignment bash script (this takes a few minutes on 16 cores)
system("./bash/sintax_tax_assignment.sh")

## Packages ####
library(tidyverse)
library(phyloseq)

## Functions ####
source("./R/functions.R")

## Metadata ####
meta <- readRDS("./output/physeq/fung_raw_physeq.RDS") %>% sample_data()
sample_names(meta)

# READ TAXONOmeta# READ TAXONOMY OUTPUT FROM VSEARCH ####
fung_otus  <- readRDS("./output/physeq/fung_seqtab_for_taxonomy.RDS")

fung_reads <- colnames(fung_otus)

fung_tax <- parse_sintax_vsearch("./output/fung.sintax.tsv", asv_order = fung_reads)

# lock in that the rows exactly match your ASV vectors
stopifnot(identical(rownames(fung_tax), fung_reads))

# rebuild phyloseq object
otu <- otu_table(fung_otus,taxa_are_rows = FALSE)
tax <- tax_table(fung_tax)

sample_names(meta)
sample_names(tax)
sample_names(otu)


fung_ps_sintax <- phyloseq(otu, tax, meta)

# export
saveRDS(fung_ps_sintax,"./output/physeq/fungal_physeq_with_sintax_taxonomy.RDS")
