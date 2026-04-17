#!/usr/bin/env Rscript

rm(list = ls())
graphics.off()

.libPaths(c("/home/zachary.shortt/R/lib",
            "/data/lab/cheeke/R_libs",
            .libPaths()))

suppressPackageStartupMessages({
  library(phyloseq)
  library(ape)
  library(picante)
  library(vegan)
  library(ecodist)
  library(permute)
  library(gee)
})

set.seed(666)

args <- commandArgs(trailingOnly = TRUE)
mode <- ifelse(length(args) >= 1, args[1], "prep")

# settings ----------------------------------------------------------------
rare.depth    <- 5146
data.set.name <- "snowbrush_april_2026"
project_name  <- "tree_for_fava"
marker        <- "16s"
no.reps       <- 999

# paths -------------------------------------------------------------------
analysis_dir  <- "sandbox/zac/analysis"
outputs_dir   <- file.path(analysis_dir, "outputs")
scripts_dir   <- file.path(analysis_dir, "R", "assembly_model")

phylogeny_dir <- file.path(outputs_dir, project_name, "phylogenies")
ps_tree_path  <- file.path(
  phylogeny_dir,
  paste0(project_name, "_", marker, "_Physeq_cleaned_w_tree.RDS")
)

assembly_dir  <- file.path(outputs_dir, "assembly_model")
null_dir      <- file.path(assembly_dir, "null_reps")
logs_dir      <- file.path(assembly_dir, "logs")

dir.create(assembly_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(null_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

otu_rds_path        <- file.path(assembly_dir, paste0(data.set.name, "_OTU_table.RDS"))
phy_rds_path        <- file.path(assembly_dir, paste0(data.set.name, "_phy.rds"))
comm_rds_path       <- file.path(assembly_dir, paste0(data.set.name, "_comm.rds"))
phylo_dist_rds_path <- file.path(assembly_dir, paste0(data.set.name, "_phylo_dist.rds"))
bmntd_obs_path      <- file.path(assembly_dir, paste0(data.set.name, "_bMNTD_weighted.csv"))
bnti_out_path       <- file.path(assembly_dir, paste0(data.set.name, "_bNTI_weighted.csv"))

source(file.path(scripts_dir, "beta.nti.R"))

# -------------------------------------------------------------------------
# prep mode: rarefy, match tree/otu, save objects, calculate observed bMNTD
# -------------------------------------------------------------------------
if (mode == "prep") {
  
  cat("Running mode = prep\n")
  cat("Input file:", ps_tree_path, "\n")
  
  if (!file.exists(ps_tree_path)) {
    stop("Missing input phyloseq file: ", ps_tree_path)
  }
  
  ps <- readRDS(ps_tree_path)
  
  ps_rare <- prune_samples(sample_sums(ps) >= rare.depth, ps)
  
  if (nsamples(ps_rare) == 0) {
    stop("No samples remain after pruning at rare.depth = ", rare.depth)
  }
  
  ps_rare <- rarefy_even_depth(
    ps_rare,
    sample.size = rare.depth,
    rngseed = 1,
    replace = FALSE,
    verbose = TRUE
  )
  
  otu <- as(otu_table(ps_rare), "matrix")
  if (taxa_are_rows(ps_rare)) otu <- t(otu)
  
  phylo <- phy_tree(ps_rare)
  phylo$tip.label <- unname(gsub("'", "", phylo$tip.label))
  
  if (!all(phylo$tip.label %in% colnames(otu))) {
    stop("Some tree tip labels are missing from OTU table.")
  }
  
  otu <- otu[, phylo$tip.label, drop = FALSE]
  saveRDS(otu, otu_rds_path)
  
  match.phylo.otu <- match.phylo.data(phylo, as.data.frame(t(otu)))
  
  comm <- t(match.phylo.otu$data)
  phy  <- match.phylo.otu$phy
  phy$tip.label <- unname(phy$tip.label)
  
  comm <- comm[, phy$tip.label, drop = FALSE]
  phylo_dist <- cophenetic(phy)
  
  saveRDS(phy, phy_rds_path)
  saveRDS(comm, comm_rds_path)
  saveRDS(phylo_dist, phylo_dist_rds_path)
  
  cat("all(colnames(comm) == phy$tip.label): ",
      all(colnames(comm) == phy$tip.label), "\n")
  cat("NA branch lengths: ", sum(is.na(phy$edge.length)), "\n")
  
  if (sum(is.na(phy$edge.length)) > 0) {
    stop("Tree has NA branch lengths; stop before betaMNTD/betaNTI.")
  }
  
  cat("Calculating observed bMNTD\n")
  print(date())
  
  beta.mntd.weighted <- as.matrix(
    comdistnt(
      t(match.phylo.otu$data),
      cophenetic(match.phylo.otu$phy),
      abundance.weighted = TRUE
    )
  )
  
  print(date())
  write.csv(beta.mntd.weighted, bmntd_obs_path, quote = FALSE)
  
  cat("Done prep mode\n")
  quit(save = "no")
}

# -------------------------------------------------------------------------
# null mode: run a chunk of null reps
# usage: Rscript run_bMNTD_pipeline.R null 1 25
# -------------------------------------------------------------------------
if (mode == "null") {
  
  cat("Running mode = null\n")
  
  if (length(args) < 3) {
    stop("For mode='null', provide start_rep and end_rep.")
  }
  
  start_rep <- as.integer(args[2])
  end_rep   <- as.integer(args[3])
  
  if (is.na(start_rep) || is.na(end_rep)) {
    stop("start_rep and end_rep must be integers.")
  }
  
  otu <- readRDS(otu_rds_path)
  phy <- readRDS(phy_rds_path)
  phy$tip.label <- unname(phy$tip.label)
  
  if (!all(colnames(otu) == phy$tip.label)) {
    stop("OTU columns do not match phy tip labels.")
  }
  
  dist_mat <- cophenetic(phy)
  
  for (i in start_rep:end_rep) {
    if (i > no.reps) break
    
    out_file <- file.path(null_dir, paste0(data.set.name, "_bMNTD_weighted_rep_", i, ".csv"))
    
    if (file.exists(out_file)) {
      cat("Skipping existing rep ", i, "\n", sep = "")
      next
    }
    
    cat("Running rep ", i, "\n", sep = "")
    print(date())
    
    rand.weighted.beta.mntd <- as.matrix(
      comdistnt(
        comm = otu,
        dis = taxaShuffle(dist_mat),
        abundance.weighted = TRUE
      )
    )
    
    write.csv(rand.weighted.beta.mntd, out_file, quote = FALSE)
    rm(rand.weighted.beta.mntd)
    gc()
  }
  
  cat("Done null mode\n")
  quit(save = "no")
}

# -------------------------------------------------------------------------
# finalize mode: check all reps and calculate betaNTI
# -------------------------------------------------------------------------
if (mode == "finalize") {
  
  cat("Running mode = finalize\n")
  
  otu <- readRDS(otu_rds_path)
  
  beta.mntd.weighted <- read.csv(bmntd_obs_path, row.names = 1, check.names = FALSE)
  
  all.files <- sort(as.numeric(
    gsub(
      ".*_rep_|\\.csv",
      "",
      list.files(
        null_dir,
        pattern = paste0(data.set.name, "_bMNTD_weighted_rep_\\d+\\.csv$")
      )
    )
  ))
  
  reps.to.do <- (1:no.reps)[!(1:no.reps %in% all.files)]
  cat("Missing reps:\n")
  print(reps.to.do)
  
  if (length(reps.to.do) > 0) {
    stop("Some null reps are still missing.")
  }
  
  beta.nti.weighted <- beta.nti.calc.stegen(
    samp = otu,
    reps = all.files,
    path.to.reps = file.path(null_dir, paste0(data.set.name, "_bMNTD_weighted_rep_")),
    beta.mntd.obs = beta.mntd.weighted
  )
  
  write.csv(beta.nti.weighted, bnti_out_path, quote = FALSE)
  
  cat("Wrote betaNTI to: ", bnti_out_path, "\n")
  cat("Done finalize mode\n")
  quit(save = "no")
}

stop("Unknown mode. Use one of: prep, null, finalize")