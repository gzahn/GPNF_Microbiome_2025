library(tidyverse)
library(phyloseq)
library(readxl)
source("./R/functions.R")


# import the data

# otu tables
fung_otu <- readRDS("./output/physeq/fung_seqtab_for_taxonomy.RDS")
bact_otu <- readRDS("./output/physeq/bact_seqtab_for_taxonomy.RDS")

# taxonomy tables
bact_tax <- readRDS("./output/bact_tax_table.RDS")
fung_tax <- readRDS("./output/fung_tax_table.RDS")

# metadata
meta <- read_xlsx("./data/metadata/clean_metadata.xlsx")


# build physeq for fungi
fung_otu_table <- otu_table(fung_otu,taxa_are_rows = FALSE)
sample_names(fung_otu_table)
fung_tax_table <- tax_table(fung_tax)
sample_names(fung_tax_table) <- sample_names(fung_otu_table)

meta <- meta %>% 
  dplyr::filter(meta$sample_name %in% sample_names(fung_otu_table))
meta_table <- sample_data(meta)

sample_names(meta_table) <- sample_names(fung_otu_table)

fung <- phyloseq(fung_otu_table,
                 fung_tax_table,
                 meta_table)

# physeq for bacteria
bact_otu_table <- otu_table(bact_otu,taxa_are_rows = FALSE)
sample_names(bact_otu_table)
bact_tax_table <- tax_table(bact_tax)
sample_names(bact_tax_table) <- sample_names(bact_otu_table)

bact <- phyloseq(bact_otu_table,
                 bact_tax_table,
                 meta_table)


# we now have phyloseq objects
# Save raw forms:
saveRDS(fung,"./output/physeq/fung_raw_physeq.RDS")
saveRDS(bact,"./output/physeq/bact_raw_physeq.RDS")

# import sintax version of fungal physeq
fung_sintax <- readRDS("./output/physeq/fungal_physeq_with_sintax_taxonomy.RDS")

# clean up non-fungi and suspect (NA) phyla
fung <- 
fung %>% 
  subset_taxa(Kingdom == "Fungi" & !is.na(Phylum))

fung_sintax <- 
  fung_sintax %>% 
  subset_taxa(Kingdom == "Fungi" & !is.na(Phylum))

# clean up non-bacteria
bact <- 
bact %>% 
  subset_taxa(Kingdom == "Bacteria" & Family != "Mitochondria" & Order != "Chloroplast")



fung %>% 
  transform_sample_counts(ra) %>% 
  plot_bar2(fill="Phylum")

fung_sintax %>% 
  transform_sample_counts(ra) %>% 
  plot_bar2(fill="Phylum")

bact %>% 
  plot_bar2(fill="Phylum")

# clean up metadata to remove useless clutter
fung@sam_data <- 
fung@sam_data %>% 
  as("data.frame") %>% 
  dplyr::select(-contains("fp_"), -starts_with("env_"), -starts_with("library_"),
                -notes,-contains("filename"),-contains("sra")) %>% 
  sample_data()

fung_sintax@sam_data <- 
  fung_sintax@sam_data %>% 
  as("data.frame") %>% 
  dplyr::select(-contains("fp_"), -starts_with("env_"), -starts_with("library_"),
                -notes,-contains("filename"),-contains("sra")) %>% 
  sample_data()

bact@sam_data <- 
  bact@sam_data %>% 
  as("data.frame") %>% 
  dplyr::select(-contains("fp_"), -starts_with("env_"), -starts_with("library_"),
                -notes,-contains("filename"),-contains("sra")) %>% 
  sample_data()

# save cleaned up forms
saveRDS(fung,"./output/physeq/fung_clean_physeq.RDS")
saveRDS(fung_sintax,"./output/physeq/fung_sintax_clean_physeq.RDS")
saveRDS(bact,"./output/physeq/bact_clean_physeq.RDS")


# check positive controls

zymo <- 
fung %>% 
  subset_samples(sample_names(fung) == "Zymo")
zymo@tax_table[1,]


zymo %>% 
  subset_taxa(taxa_sums(zymo) > 0) %>% 
  plot_bar2(fill = "Order")
zymo %>% 
  subset_taxa(taxa_sums(zymo) > 0) %>% 
  otu_table() %>% 
  taxa_names()
# dig into taxonomy for fungi... these look weird


fung@sam_data$vigor

fung@sam_data$mergevar <- paste(fung@sam_data$vigor,fung@sam_data$unit,sep="|")
fung_m <- fung %>% 
  merge_samples("mergevar")

fung_m@sam_data$unit <- sample_names(fung_m) %>% str_split("\\|") %>% map_chr(2)



fung_m %>% 
  transform_sample_counts(ra) %>% 
  plot_bar2(fill = "Phylum") +
  facet_wrap(~vigor,scales="free")



zymo <- 
  bact %>% 
  subset_samples(sample_names(bact) == "Zymo")

zymo %>% 
  subset_taxa(taxa_sums(zymo) > 0) %>% 
  plot_bar2(fill = "Order")

zymo %>% 
  subset_taxa(taxa_sums(zymo) > 0) %>% 
  otu_table() %>% 
  taxa_names()
