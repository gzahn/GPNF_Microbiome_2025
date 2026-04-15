# BUILD PHYLOGENIES ####

# SETUP ####
set.seed(666)

# setup for Kamiak ###CHANGE ME######
.libPaths(c("/home/zachary.shortt/R/lib",
            "/data/lab/cheeke/R_libs",
            .libPaths()))

## packages ####
library(tidyverse)
library(phyloseq)
library(DECIPHER)
library(janitor)
library(phylogram)
library(Biostrings)

# PROJECT SETTINGS ####
project_name <- "tree_for_fava"
marker <- "16s"

# directories
base_dir      <- "sandbox/zac/analysis/outputs"
physeq_dir    <- file.path(base_dir, "phyloseq_objects")
phylogeny_dir <- file.path(base_dir, project_name, "phylogenies")

dir.create(phylogeny_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(physeq_dir, recursive = TRUE, showWarnings = FALSE)



# file paths
ps_path <- file.path(physeq_dir, paste0(project_name, "_", marker, "_clean_phyloseq_object.RDS"))

if (!file.exists(ps_path)) {
  ps_path <- file.path(physeq_dir, paste0(marker, "_clean_phyloseq_object.RDS"))
}

alignment_rds_path <- file.path(phylogeny_dir, paste0(project_name, "_", marker, "_dna_alignment.RDS"))
aligned_fasta_path <- file.path(phylogeny_dir, paste0(project_name, "_", marker, "_ASV_aligned.fasta"))
tree_nwk_path      <- file.path(phylogeny_dir, paste0(project_name, "_", marker, "_fasttree.nwk"))
ps_tree_path       <- file.path(phylogeny_dir, paste0(project_name, "_", marker, "_Physeq_cleaned_w_tree.RDS"))

## load phyloseq object ####
ps_16s <- readRDS(ps_path)

## extract sequences ####
seqs_16s <- rownames(tax_table(ps_16s))
names(seqs_16s) <- paste0("ASV_", seq_along(seqs_16s))  # propagates to tree tip labels

# ALIGNMENT ####
alignment_16s <- DECIPHER::AlignSeqs(DNAStringSet(seqs_16s), processors = NULL)
saveRDS(alignment_16s, alignment_rds_path)

Biostrings::writeXStringSet(alignment_16s, aligned_fasta_path)

# Build tree with FastTree ####
cmd <- paste(
  "FastTree -nt -gtr -gamma",
  shQuote(aligned_fasta_path),
  ">",
  shQuote(tree_nwk_path)
)

status <- system(cmd)
stopifnot(status == 0, file.exists(tree_nwk_path))

tree_16s <- ape::read.tree(tree_nwk_path)

# relabel tips from ASV_# -> phyloseq taxa_names
map_16s <- setNames(taxa_names(ps_16s), names(seqs_16s))
tree_16s$tip.label <- map_16s[tree_16s$tip.label]

# sanity checks ####
stopifnot(length(tree_16s$tip.label) == ntaxa(ps_16s))
stopifnot(setequal(tree_16s$tip.label, taxa_names(ps_16s)))
stopifnot(!is.null(tree_16s$edge.length))

# attach tree to phyloseq ####
ps_16s_w_tree <- merge_phyloseq(ps_16s, phy_tree(tree_16s))

# export ####
saveRDS(ps_16s_w_tree, ps_tree_path)