# RUN FAVA
set.seed(666)

suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
  library(FAVA)
  library(ape)
  library(viridis)
})

source("./R/functions.R")

theme_set(
  theme_bw() +
    theme(
      axis.title = element_text(face = "bold", size = 18),
      axis.text = element_text(face = "bold", size = 14),
      legend.title = element_text(face = "bold", size = 18),
      legend.text = element_text(face = "bold", size = 14),
      strip.text = element_text(face = "bold", size = 18),
      plot.title = element_text(face = "bold", size = 18),
      strip.background = element_rect(fill = "white")
    )
)

# -------------------------------------------------------------------
# paths
# -------------------------------------------------------------------
project_name  <- "tree_for_fava"
marker        <- "16s"
base_dir      <- "sandbox/zac/analysis/outputs"
phylogeny_dir <- file.path(base_dir, project_name, "phylogenies")
ps_tree_path  <- file.path(
  phylogeny_dir,
  paste0(project_name, "_", marker, "_Physeq_cleaned_w_tree.RDS")
)

# output dirs
out_dir  <- "./output"
fig_dir  <- file.path(out_dir, "figs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# load phyloseq objects
# -------------------------------------------------------------------
fung <- readRDS("./output/physeq/fung_clean_physeq.RDS") %>%
  subset_samples(!grepl("zymo", sample_name, ignore.case = TRUE))

bact <- readRDS(ps_tree_path) %>%
  prune_samples(!grepl("zymo", sample_names(.), ignore.case = TRUE), .)

# -------------------------------------------------------------------
# add coordinates to sample data
# -------------------------------------------------------------------
coords_df <- read.csv("data/metadata/meta_data_with coords.csv") %>%
  select(sample_id, lattitude, longitude, elevation)

fung_meta <- data.frame(sample_data(fung)) %>%
  left_join(coords_df, by = c("sample_name" = "sample_id")) %>%
  column_to_rownames("sample_name")

sample_data(fung) <- sample_data(fung_meta)

bact_meta <- data.frame(sample_data(bact)) %>%
  left_join(coords_df, by = c("sample_name" = "sample_id")) %>%
  column_to_rownames("sample_name")

sample_data(bact) <- sample_data(bact_meta)

# -------------------------------------------------------------------
# convert to FAVA relative-abundance tables
# -------------------------------------------------------------------
fung_dat <- FAVA::relab_phyloseq(fung)
bact_dat <- FAVA::relab_phyloseq(bact)

# number of metadata columns at front of table
fung_nmeta <- 1 + ncol(data.frame(sample_data(fung)))
bact_nmeta <- 1 + ncol(data.frame(sample_data(bact)))

fung_K <- ncol(fung_dat) - fung_nmeta + 1
bact_K <- ncol(bact_dat) - bact_nmeta + 1

# -------------------------------------------------------------------
# bacterial tree + similarity matrix for weighted FAVA
# -------------------------------------------------------------------
bact_tree <- phy_tree(bact)
stopifnot(!is.null(bact_tree))
stopifnot(length(bact_tree$tip.label) == ntaxa(bact))
stopifnot(setequal(bact_tree$tip.label, taxa_names(bact)))
stopifnot(!is.null(bact_tree$edge.length))
stopifnot(sum(is.na(bact_tree$edge.length)) == 0)
stopifnot(sum(bact_tree$edge.length < 0, na.rm = TRUE) == 0)

bact_dist <- ape::cophenetic.phylo(bact_tree)
bact_order <- colnames(bact_dat)[bact_nmeta:ncol(bact_dat)]

if (!all(
  bact_order == colnames(bact_dist),
  bact_order == rownames(bact_dist)
)) {
  stop("Error: bacterial taxa order does not match distance matrix.")
}

bact_simil_mat <- exp(-bact_dist)

# -------------------------------------------------------------------
# optional relabundance plots
# -------------------------------------------------------------------
bact_palette <- viridis::turbo(nrow(bact_simil_mat))[
  sample(seq_len(nrow(bact_simil_mat)))
]
names(bact_palette) <- colnames(bact_dat)[bact_nmeta:ncol(bact_dat)]

plot_relabund(
  bact_dat,
  group = "unit",
  arrange = "both",
  K = bact_K
) +
  ggplot2::scale_color_manual(values = bact_palette) +
  ggplot2::scale_fill_manual(values = bact_palette)

# fungi plot optional too, but can be huge / messy
# plot_relabund(
#   fung_dat,
#   group = "unit",
#   arrange = "both",
#   K = fung_K
# )

# -------------------------------------------------------------------
# compute FAVA
# -------------------------------------------------------------------

# fungi: unweighted only
fung_fava_uw <- fava(
  relab_matrix = fung_dat,
  group = "unit",
  K = fung_K
)

# bacteria: unweighted
bact_fava_uw <- fava(
  relab_matrix = bact_dat,
  group = "unit",
  K = bact_K
)

# bacteria: weighted by phylogenetic similarity
bact_fava_w <- fava(
  relab_matrix = bact_dat,
  group = "unit",
  K = bact_K,
  S = bact_simil_mat
)

# -------------------------------------------------------------------
# bootstrap
# -------------------------------------------------------------------

# fungi: unweighted bootstrap
fung_fava_bs <- bootstrap_fava(
  relab_matrix = fung_dat,
  n_replicates = 1000,
  seed = 666,
  group = "unit",
  K = fung_K
)

# bacteria: weighted bootstrap
bact_fava_bs <- bootstrap_fava(
  relab_matrix = bact_dat,
  n_replicates = 1000,
  seed = 666,
  group = "unit",
  K = bact_K,
  S = bact_simil_mat
)

saveRDS(fung_fava_bs, "./sandbox/zac/analysis/outputs/bootstrap_fava_fungi.RDS")
saveRDS(bact_fava_bs, "./sandbox/zac/analysis/outputs/bootstrap_fava_bacteria.RDS")

# -------------------------------------------------------------------
# export pairwise bootstrap stats
# -------------------------------------------------------------------
fung_pairwise <- fung_fava_bs$observed_difference %>%
  mutate(Comparison = str_replace_all(Comparison, "\n", " ")) %>%
  full_join(fung_fava_bs$P_values, by = "Comparison")

bact_pairwise <- bact_fava_bs$observed_difference %>%
  mutate(Comparison = str_replace_all(Comparison, "\n", " ")) %>%
  full_join(bact_fava_bs$P_values, by = "Comparison")

write_csv(fung_pairwise, file.path(out_dir, "fava_stats_fungi.csv"))
write_csv(bact_pairwise, file.path(out_dir, "fava_stats_bacteria.csv"))

# -------------------------------------------------------------------
# export per-group FAVA values
# -------------------------------------------------------------------
fung_values <- fung_fava_uw %>%
  mutate(amplicon = "Fungi", weighting = "Unweighted")

bact_values_uw <- bact_fava_uw %>%
  mutate(amplicon = "Bacteria", weighting = "Unweighted")

bact_values_w <- bact_fava_w %>%
  mutate(amplicon = "Bacteria", weighting = "Weighted")

fava_values <- bind_rows(
  fung_values,
  bact_values_uw,
  bact_values_w
)

write_csv(fava_values, file.path(out_dir, "fava_values_all.csv"))

# -------------------------------------------------------------------
# plots of per-group FAVA values
# -------------------------------------------------------------------
p_fung <- fung_fava_uw %>%
  ggplot(aes(x = unit, y = FAVA)) +
  geom_point(size = 4) +
  labs(
    title = "Fungal FAVA by unit",
    y = "Unweighted FAVA",
    x = "Unit"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_fung

p_bact_uw <- bact_fava_uw %>%
  ggplot(aes(x = unit, y = FAVA)) +
  geom_point(size = 4) +
  labs(
    title = "Bacterial unweighted FAVA by unit",
    y = "Unweighted FAVA",
    x = "Unit"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_bact_w <- bact_fava_w %>%
  ggplot(aes(x = unit, y = FAVA)) +
  geom_point(size = 4) +
  labs(
    title = "Bacterial weighted FAVA by unit",
    y = "Weighted FAVA",
    x = "Unit"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "fung_fava_values.png"), p_fung, dpi = 400, height = 8, width = 12)
ggsave(file.path(fig_dir, "bact_fava_values_unweighted.png"), p_bact_uw, dpi = 400, height = 8, width = 12)
ggsave(file.path(fig_dir, "bact_fava_values_weighted.png"), p_bact_w, dpi = 400, height = 8, width = 12)

# -------------------------------------------------------------------
# bootstrap distribution plots
# -------------------------------------------------------------------
ggsave(
  plot = fung_fava_bs$bootstrap_distribution_plot,
  filename = file.path(fig_dir, "fava_bootstrap_plot_fungi.png"),
  dpi = 400,
  width = 16,
  height = 8
)

ggsave(
  plot = bact_fava_bs$bootstrap_distribution_plot,
  filename = file.path(fig_dir, "fava_bootstrap_plot_bacteria.png"),
  dpi = 400,
  width = 16,
  height = 8
)