# SETUP ####

## packages ####
# SETUP ####
.libPaths(c("/data/lab/cheeke/R_libs", .libPaths()))
## packages ####
library(dada2, lib.loc = "/data/lab/cheeke/R_libs")
library(ShortRead, lib.loc = "/data/lab/cheeke/R_libs")
library(Biostrings, lib.loc = "/data/lab/cheeke/R_libs")
library(archive, lib.loc = "/data/lab/cheeke/R_libs")

# options
options(timeout=600) # longer download timeout....change to more seconds if timing out
nthreads=parallel::detectCores()-1

# functions
source("./R/functions.R")

# dir setup
if(!dir.exists("./data/databases")){
  dir.create("./data/databases")
}

## download taxonomy db files (if they don't exist already)
# silva 138.2
if(!file.exists("./data/databases/silva_nr99_v138.2_toSpecies_trainset.fa.gz")){
  download.file("https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toSpecies_trainset.fa.gz?download=1",
                destfile = "./data/databases/silva_nr99_v138.2_toSpecies_trainset.fa.gz")
}
if(!file.exists("./data/databases/silva_v138.2_assignSpecies.fa.gz")){
  download.file("https://zenodo.org/records/14169026/files/silva_v138.2_assignSpecies.fa.gz?download=1",
                destfile = "./data/databases/silva_v138.2_assignSpecies.fa.gz")
}
# eukaryome 1.9.4
if(!file.exists("./data/databases/DADA2_EUK_ITS_v1.9.4.fasta.gz")){
  download.file("https://sisu.ut.ee/wp-content/uploads/sites/643/DADA2_EUK_ITS_v1.9.4.zip",
                destfile = "./data/databases/DADA2_EUK_ITS_v1.9.4.zip")
}
# extract compressed eukaryome database and prep for dada2
if(!file.exists("./data/databases/DADA2_EUK_ITS_v1.9.4.fasta.gz")){
  unzip(zipfile = "./data/databases/DADA2_EUK_ITS_v1.9.4.zip",junkpaths = TRUE)
  archive_extract("./DADA2_EUK_ITS_v1.9.4.7z",dir = "./data/databases")
  system(command = "gzip ./data/databases/DADA2_EUK_ITS_v1.9.4.fasta")
  file.remove("./DADA2_EUK_ITS_v1.9.4.7z")
  file.remove("./data/databases/DADA2_EUK_ITS_v1.9.4.zip")
}

# set taxonomy database objects
bact_species_db <- "./data/databases/silva_nr99_v138.2_toSpecies_trainset.fa.gz"
bact_assign_species_db <- "./data/databases/silva_v138.2_assignSpecies.fa.gz"
unite_db <- "./data/databases/UNITE_public_19.02.2025.fasta.gz"  # CHANGE ME to location on your machine

if(any(!file.exists(bact_species_db,bact_assign_species_db, unite_db))){
  stop("Make sure taxonomy databases are downloaded and in ./data/databases/")
}


## Load ASV tables

bact <- readRDS("./output/physeq/bact_seqtab_for_taxonomy.RDS")
fung <- readRDS("./output/physeq/fung_seqtab_for_taxonomy.RDS")

# ASSIGN TAXONOMY ####

## bacteria ####
bact_tax <- 
  assignTaxonomy(bact,
                 refFasta = bact_species_db,
                 multithread = nthreads,
                 tryRC = TRUE,
                 minBoot = 70)
# run addSpecies algorithm and add to previous tax table
bact_tax_spp <-addSpecies(taxtab = bact_tax,refFasta = bact_assign_species_db,tryRC = TRUE)
# save progress
saveRDS(bact_tax_spp,"./output/bact_tax_table.RDS")

## fungi ####
fung_tax <- 
  assignTaxonomy(fung,
                 refFasta = unite_bd,
                 multithread = nthreads,
                 tryRC = TRUE,
                 minBoot = 70)
# save progress
saveRDS(fung_tax,"./output/fung_tax_table_unite.RDS")