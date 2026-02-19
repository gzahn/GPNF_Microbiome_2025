
library(tidyverse)
library(vegan)
library(phyloseq)
library(broom)
library(ggthemes)
library(RColorBrewer)
library(fungaltraits)
library(janitor)
library(corncob)
library(ranger)
library(vip)
library(openxlsx)



phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange", "#DA5724", "#508578", "#CD9BCD",
  "#8a592f", "#673770", "#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861",
  
  # Added 5 new distinct colors
  "#1FA187",  # teal
  "#F6E8C3",  # light sand
  "#B8E186",  # soft green (not overlapping your darker greens)
  "#E66101",  # warm orange-brown, distinct from pure orange
  "#4D4D4D"   # neutral charcoal for grounding
)

#a <- as.data.frame(sample_data(ps_guild))

source("./R/functions/funguild.R")
# load phyloseq object -----------------------------------------------------
ps <- readRDS("output/physeq/fungal_phyloseq.rds") %>% 
  prune_samples(!grepl("zymo", sample_names(.), ignore.case = TRUE), .)

ps_ra <- transform_sample_counts(ps, function(x) x / sum(x))


# permanova stuff  --------------------------------------------------------

## ------------------------------------------------------------
## PERMANOVA (adonis2) assumption checks for a continuous trait
## using the equivalent dbRDA/capscale model + dispersion tests
## ------------------------------------------------------------

library(phyloseq)
library(vegan)
library(permute)

## --- 0) Prep metadata (make sure types are correct) ----------
meta <- data.frame(sample_data(ps_ra), check.names = FALSE)

# continuous predictor MUST be numeric
meta$dry_mass <- as.numeric(as.character(meta$dry_mass))

# grouping / blocking factor
meta$unit <- factor(meta$unit)

# sanity checks
stopifnot(is.numeric(meta$dry_mass))
stopifnot(is.factor(meta$unit))
stopifnot(nrow(meta) == attr(bc_dist, "Size"))

## --- 1) Run PERMANOVA (main effects; marginal tests) ---------
perm_mod <- adonis2(
  bc_dist ~ dry_mass + unit,
  data = meta,
  permutations = 999,
  by = "margin"
)
print(perm_mod)

## Optional: test interaction (does slope differ by unit?)
perm_int <- adonis2(
  bc_dist ~ dry_mass * unit,
  data = meta,
  permutations = 999,
  by = "margin"
)
print(perm_int)

## --- 2) Key PERMANOVA assumption: homogeneity of dispersion ---
# For group factor (unit). If significant, interpret PERMANOVA for unit cautiously.
bd <- betadisper(bc_dist, meta$unit)
print(anova(bd))
print(permutest(bd, permutations = 999))

# Quick visualization
plot(bd)

## --- 3) dbRDA model equivalent to adonis2 for diagnostics -----
# Controls for unit, tests dry_mass (matches the idea of "dry_mass + unit" with unit partialled out)
mdb <- capscale(bc_dist ~ dry_mass + Condition(unit), data = meta)

# Model tests (should broadly align with adonis2 results)
print(anova(mdb))                 # overall
print(anova(mdb, by = "terms"))   # term-wise (dry_mass after conditioning)

## --- 4) Residual diagnostics (like lm() but in distance space)
# Residuals vs fitted (pattern = nonlinearity / leverage)
ordiresids(mdb, display = "working")
# Scale-location (heteroscedasticity / variance patterns)
ordiresids(mdb, display = "scale", residuals = "working")
# Normal Q-Q (outliers / heavy tails)
ordiresids(mdb, display = "qq", residuals = "working")

## --- 5) Optional: test nonlinearity in dry_mass --------------
mdb_poly <- capscale(bc_dist ~ poly(dry_mass, 2) + Condition(unit), data = meta)
print(anova(mdb, mdb_poly))  # compare linear vs quadratic

## --- 6) Optional: influence / leverage diagnostics -----------
inf <- influence(mdb)

# Hat values (leverage)
plot(inf$hat, type = "h",
     ylab = "Leverage (hat values)", xlab = "Sample index")
abline(h = 2 * mean(inf$hat, na.rm = TRUE), col = "red", lty = 2)

# Optionally: list top leverage samples (by rowname)
top_hat <- order(inf$hat, decreasing = TRUE)[1:10]
print(data.frame(
  SampleID = rownames(meta)[top_hat],
  hat = inf$hat[top_hat],
  dry_mass = meta$dry_mass[top_hat],
  unit = meta$unit[top_hat]
))

# Bray-Curtis distance matrix
bc_dist <- phyloseq::distance(ps_ra, method = "bray")

# Metadata
meta <- data.frame(sample_data(ps_ra))
meta$dry_mass <- as.numeric(as.character(meta$dry_mass))


# PERMANOVA
a <- adonis2(bc_dist ~ vigor, data = meta, permutations = 999)
a

adonis2(bc_dist ~ dry_mass + unit, data = meta,
        permutations = 999, by = "margin")


nmds <- ordinate(ps_ra, method = "NMDS", distance = "bray")

fit <- envfit(nmds, meta$dry_mass, permutations = 999)

plot(nmds)
plot(fit, col = "red")

capscale(bc_dist ~ dry_mass + Condition(unit), data = meta)



ps_b <- readRDS("output/physeq/bact_phyloseq.rds") %>% 
  prune_samples(!grepl("zymo", sample_names(.), ignore.case = TRUE), .)

alpha_div <- estimate_richness(ps)

sd <- data.frame(sample_data(ps), check.names = FALSE)
sd$SampleID <- rownames(sd)

# 2) Alpha diversity (add whatever measures you want)
alpha_div <- estimate_richness(ps, measures = c("Observed","Chao1","Shannon","Simpson","InvSimpson"))
alpha_div$SampleID <- rownames(alpha_div)

# 3) Join, then restore rownames and drop SampleID
sd2 <- dplyr::left_join(sd, alpha_div, by = "SampleID")
rownames(sd2) <- sd2$SampleID
sd2$SampleID <- NULL

# 4) Put back into phyloseq
sample_data(ps) <- phyloseq::sample_data(sd2)



guilds <- assign_funguild_to_phyloseq(ps)

ps_guild <- guilds [[1]]

ps_ecm <- subset_taxa(ps_guild, guild_fg == "Ectomycorrhizal")


# stacked bar graphs  -----------------------------------------------------
ps_genus_20 <- ps%>%
  tax_glom(taxrank = "Family") %>%                    
  transform_sample_counts(function(x) {x / sum(x)}) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.05) %>%                        
  mutate(Family = ifelse(Abundance < 0.15, "Other", Family))

ggplot(ps_genus_20)+
  aes(x = vigor, y = Abundance, fill = fct_rev(fct_infreq(Family))) +
  geom_col(position = "fill") +  
 # scale_fill_manual(values = phylum_colors) +
  theme(legend.position = "right",                 
        legend.title = element_text(size = 8),      
        legend.text = element_text(size = 6),      
        legend.key.size = unit(0.5, "lines"),       
        legend.spacing = unit(0.5, "lines")) +  
  #facet_wrap(~colonized)+
  guides(alpha = "none")+
  labs(
    x = "Vigor",
    y = "Proportion of bact Families",
    fill = "Family"
  )+
  theme_few()



ps_genus_20 <- ps%>%
  tax_glom(taxrank = "Genus") %>%                    
  transform_sample_counts(function(x) {x / sum(x)}) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.05) %>%                        
  mutate(Family = ifelse(Abundance < 0.1, "Other", Genus))

ggplot(ps_genus_20)+
  aes(x = dry_mass, y = Abundance, fill = fct_rev(fct_infreq(Genus))) +
  geom_col(position = "fill") +  
  #scale_fill_manual(values = phylum_colors) +
  theme(legend.position = "right",                 
        legend.title = element_text(size = 8),      
        legend.text = element_text(size = 6),      
        legend.key.size = unit(0.5, "lines"),       
        legend.spacing = unit(0.5, "lines")) +  
  #facet_wrap(~colonized)+
  guides(alpha = "none")+
  labs(
    x = "Vigor",
    y = "Proportion of Families",
    fill = "Family"
  )+
  theme_few()


ps_genus_20 <- ps_ecm %>%
  tax_glom(taxrank = "Genus") %>%                    
  transform_sample_counts(function(x) {x / sum(x)}) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.05) %>%                        
  mutate(Family = ifelse(Abundance < 0.1, "Other", Genus))

ggplot(ps_genus_20)+
  aes(x = dry_mass, y = Abundance, fill = fct_rev(fct_infreq(Genus))) +
  geom_col(position = "fill") +  
  scale_fill_manual(values = phylum_colors) +
  theme(legend.position = "right",                 
        legend.title = element_text(size = 8),      
        legend.text = element_text(size = 6),      
        legend.key.size = unit(0.5, "lines"),       
        legend.spacing = unit(0.5, "lines")) +  
  #facet_wrap(~colonized)+
  guides(alpha = "none")+
  labs(
    x = "Vigor",
    y = "Proportion of Families",
    fill = "Family"
  )+
  theme_few()

ggplot(ps_genus_20)+
  aes(x = dry_mass , y = Abundance, fill = fct_rev(fct_infreq(Genus))) +
  geom_col(position = "fill") +  
  scale_fill_manual(values = phylum_colors) +
  theme(legend.position = "right",                 
        legend.title = element_text(size = 8),      
        legend.text = element_text(size = 6),      
        legend.key.size = unit(0.5, "lines"),       
        legend.spacing = unit(0.5, "lines")) +  
  #facet_wrap(~colonized)+
  guides(alpha = "none")+
  labs(
    x = "Vigor",
    y = "Proportion of Families",
    fill = "Family"
  )+
  theme_few()

a <- as.data.frame(tax_data(ps_em))


# proportion mutualist ----------------------------------------------------

guild_df <- guilds [[2]] %>% 
  mutate(
    dry_mass = as.numeric(dry_mass),
    log_mutualist = log10(proportion_mutualist + 1e-6),
    log_pathogen = log10(proportion_pathogen + 1e-6),
  )

ggplot(guild_df, aes(x = log_mutualist, y = dry_mass, ))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_few()

ggplot(guild_df, aes(x = proportion_pathogen, y = dry_mass, ))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_few()


# new analysis ------------------------------------------------------------

#div metrics

sam_dat <- data.frame(sample_data(ps), check.names = FALSE) %>% 
  mutate(
    dry_mass = as.numeric(dry_mass)
  )

ggplot(sam_dat, aes(x = vigor, y = InvSimpson, fill = vigor))+
  geom_boxplot()+
    theme_few()

for (a in c("Observed","Chao1","Shannon","Simpson","InvSimpson")) {
  p <- ggplot(sam_dat, aes(x = vigor, y = .data[[a]], fill = vigor)) +
    geom_point(position = position_jitter(width = 0.15), alpha = 0.6) +
    theme_few() +
    labs(title = "bact", x = "Vigor", y = a)
  
  print(p)
}


#sam_dat <- data.frame(sample_data(ps_b), check.names = FALSE)

ggplot(sam_dat, aes(x = vigor, y = InvSimpson, fill = vigor))+
  geom_boxplot()+
  theme_few()

for (a in  c("Observed","Chao1","Shannon","Simpson","InvSimpson")){
  p <- ggplot(sam_dat, aes_string(x = "dry_mass", y = a, fill = "unit"))+
    geom_point()+
    geom_smooth(method = "lm")+
    theme_few()+
    labs(
      x = "dry mass",
      y = a
    )
  
  print(p)
}

for (a in  c("Observed","Chao1","Shannon","Simpson","InvSimpson")){
  p <- ggplot(sam_dat, aes(x = "dry_mass", y = a, fill = vigor))+
    geom_point()+
    geom_smooth(method="lm")+
    theme_few()
  
  print(p)
}

ggplot(sam_dat, aes(x = , y = dry_mass, fill = vigor))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_few()

# running random forest ---------------------------------------------------

library(randomForest)

run_rf_growth <- function(ps,
                          growth_predictors = c("shoot_dm", "leaf_number","final_root_dm"),
                          gene,
                          tax_rank          = "Genus",
                          num_trees         = 999,
                          seed              = 2,
                          num_features_vip = 20,
                          out_dir           = NULL) {
  
  # ---- packages ----
  pkgs <- c("phyloseq","janitor","corncob","ranger","tidyverse")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "))
  
  # ---- checks ----
  if (!inherits(ps, "phyloseq")) stop("ps must be a phyloseq object.")
  meta0 <- as(phyloseq::sample_data(ps), "data.frame")
  
  missing_cols <- setdiff(growth_predictors, colnames(meta0))
  if (length(missing_cols)) stop("Missing columns in sample_data(ps): ", paste(missing_cols, collapse = ", "))
  
  # ---- tax glom (optional) ----
  ps_rank <- ps
  if (!is.null(tax_rank) && nzchar(tax_rank)) {
    ps_rank <- phyloseq::tax_glom(ps, taxrank = tax_rank)
  }
  
  # ---- RA table (samples x taxa) ----
  ps_ra <- phyloseq::transform_sample_counts(ps_rank, function(x) x / sum(x))
  otu_mat <- as(phyloseq::otu_table(ps_ra), "matrix")
  ra_table <- as.data.frame(otu_mat, check.names = FALSE)
  
  # ---- rename taxa columns to taxonomy labels ----
  taxa_print <- corncob::otu_to_taxonomy(phyloseq::taxa_names(ps_ra), ps_ra)
  colnames(ra_table) <- janitor::make_clean_names(taxa_print)
  
  # drop all-zero taxa columns
  ra_table <- ra_table[, colSums(ra_table) > 0, drop = FALSE]
  
  # ---- metadata (plain data.frame) ----
  meta <- data.frame(phyloseq::sample_data(ps_ra), check.names = FALSE, stringsAsFactors = FALSE)
  
  # ---- modeling base df: growth + taxa ----
  df <- meta %>%
    dplyr::select(dplyr::all_of(growth_predictors)) %>%
    dplyr::bind_cols(ra_table)
  
  # ---- fit RF per growth var ----
  top_taxa_dfs <- list()
  models <- list()
  vip_plots <- list()
  
  for (growth_var in growth_predictors) {
    
    df_1 <- df %>%
      dplyr::mutate(
        "{growth_var}" := as.numeric(as.character(.data[[growth_var]]))
      ) %>%
      dplyr::filter(!is.na(.data[[growth_var]])) %>% 
      select(-dplyr::all_of(growth_predictors[growth_predictors != growth_var]))
    
    set.seed(seed)
    rf_mod <- ranger::ranger(
      formula = stats::as.formula(paste0(growth_var, " ~ .")),
      data = df_1,
      importance = "permutation",
      num.trees = num_trees
    )
    
    models[[growth_var]] <- rf_mod
    
    p_vip <- vip::vip(rf_mod, num_features = num_features_vip) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(face = "bold.italic")) +
      ggplot2::labs(
        title = "Top taxa predicting plant growth (RF regression)",
        subtitle = paste0("Response: ", growth_var)
      )
    
    vip_plots[[growth_var]] <- p_vip
    
    if (!is.null(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      ggplot2::ggsave(
        filename = file.path(out_dir, paste0(gene, "/"),paste0(growth_var, ".png")),
        plot = p_vip,
        width = 10, height = 5, units = "in", dpi = 300
      )
    }
    
    top_20 <- vip::vi_model(rf_mod) %>%
      dplyr::arrange(dplyr::desc(Importance)) %>%
      dplyr::slice_head(n = 20)
    
    top_taxa_dfs[[growth_var]] <- top_20
    print (paste0("Completed for ", growth_var))
    
  }
  
  # ---- combine top taxa tables ----
  top_taxa_long <- dplyr::bind_rows(top_taxa_dfs, .id = "growth_metric") %>%
    dplyr::mutate(
      taxa = stringr::str_extract(Variable, "[^_]+_[^_]+$") # keep last two underscore chunks
    )
  
  return(list(
    df = df,
    ra_table = ra_table,
    models = models,
    vip_plots = vip_plots,
    top_taxa_dfs = top_taxa_dfs,
    top_taxa_long = top_taxa_long
  ))
}
# running rf for bact -----------------------------------------------------



# running rf for fung -----------------------------------------------------

res_fung <- run_rf_growth(
  ps = fung,
  gene = "ITS",
  growth_predictors = c("vigor","dry_mass"),
  tax_rank = "Genus",
  out_dir = "./R/Zac mess around/figures/random_forest"
)

a <- res_fung$top_taxa_long
exported_table <- a %>%
  select(growth_metric, taxa, Importance)

write.xlsx(exported_table, "./R/Zac mess around/figures/random_forest/ITS_toptaxa.xlsx")


# pdp --------------------------------------------------------------------

library(pdp)

ordered_tax <- res_fung$top_taxa_long %>%  
  dplyr::arrange(dplyr::desc(Importance))

tax <- ordered_tax$Variable[3]
pdp::partial(
  object = res_fung$models$dry_mass,
  pred.var = tax,
  train = res_fung$df,
  plot = TRUE
)

a <- as.data.frame(sample_data(fung))


