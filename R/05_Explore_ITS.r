library(tidyverse)
library(vegan)
library(phyloseq)
library(broom)
library(ggthemes)


ps <- readRDS("output/physeq/fungal_phyloseq.rds") %>% 
  subset_samples(!grepl("zymo", sample_names(ps), ignore.case = TRUE))


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

# guild assignment  -------------------------------------------------------
###download datasheet from https://pmc.ncbi.nlm.nih.gov/articles/PMC9958157/#sec21##

fungal_traits <- read_csv("data/fungaltraits/fungal_traits.csv") %>% 
  rename (Genus = GENUS)

ps_genus <- subset_taxa(ps, !is.na(Genus) & Genus != "")

tax_df <- as.data.frame(tax_table(ps_genus), stringsAsFactors = FALSE) %>%
  rownames_to_column("ASV")

trait_cols <- c("Genus", "primary_lifestyle", 
                "Secondary_lifestyle",
                "Ectomycorrhiza_lineage_template", 
                "Ectomycorrhiza_exploration_type_template"
                )  

fungal_traits_sub <- fungal_traits %>%
  select(all_of(trait_cols)) %>%     
  distinct(Genus, .keep_all = TRUE)  


tax_with_traits <- tax_df %>%
  left_join(fungal_traits_sub, by = "Genus")

tax_mat <- tax_with_traits %>%
  mutate(across(everything(), as.character)) %>% 
  column_to_rownames("ASV") %>%
  as.matrix()

tax_table(ps_genus) <- tax_table(tax_mat)
#The above code added fungal traits to our phyloseq like fungal lifestyle, meaning we can use 
#commands like subset_taxa now with fungal traits!

# -#####EcM subset##------------------------------------------------------------------------

ps_ecm <- subset_taxa(ps_genus, str_detect(tolower(primary_lifestyle), "ectomycorrhizal"))





# rarefaction curve, and rarefying data -----------------------------------


otu <- as(otu_table(ps), "matrix")

rarecurve(otu, step = 100, sample = min(rowSums(otu)), col = "gray", cex = 0.6, label = FALSE)
# # FALSE
# set.seed(711)
# rarefied_ps <-rarefy_even_depth(ps, sample.size = min(sample_sums(ps)),
#                                 rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

#####Stacked Bar######
ps_family_20 <- ps %>%
  tax_glom(taxrank = "Family") %>%                    
  transform_sample_counts(function(x) {x / sum(x)}) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.01) %>%                        
  mutate(Family = ifelse(Abundance < 0.2, "Other", Family))

ps_family_05 <- ps %>%
  tax_glom(taxrank = "Family") %>%                    
  transform_sample_counts(function(x) {x / sum(x)}) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.01) %>%                        
  mutate(Family = ifelse(Abundance < 0.05, "Other", Family))

ps_genus_05 <- ps %>%
  tax_glom(taxrank = "Genus") %>%                    
  transform_sample_counts(function(x) {x / sum(x)}) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.01) %>%                        
  mutate(Family = ifelse(Abundance < 0.05, "Other", Genus))

ps_genus_20 <- ps %>%
  tax_glom(taxrank = "Genus") %>%                    
  transform_sample_counts(function(x) {x / sum(x)}) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.2) %>%                        
  mutate(Family = ifelse(Abundance < 0.3, "Other", Genus))

ggplot(ps_family_20)+
  aes(x = vigor, y = Abundance, fill = fct_rev(fct_infreq(Family))) +
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
  aes(x = vigor, y = Abundance, fill = fct_rev(fct_infreq(Genus))) +
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
    fill = "Genus"
  )+
  theme_few()

ps_ecm_family_20 <- ps_ecm %>%
  tax_glom(taxrank = "Family") %>%                    
  transform_sample_counts(function(x) {x / sum(x)}) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.01) %>%                        
  mutate(Family = ifelse(Abundance < 0.2, "Other", Family))

ps_ecm_family_05 <- ps_ecm %>%
  tax_glom(taxrank = "Family") %>%                    
  transform_sample_counts(function(x) {x / sum(x)}) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.01) %>%                        
  mutate(Family = ifelse(Abundance < 0.05, "Other", Family))

ps_ecm_genus_05 <- ps_ecm %>%
  tax_glom(taxrank = "Genus") %>%                    
  transform_sample_counts(function(x) {x / sum(x)}) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.01) %>%                        
  mutate(Genus = ifelse(Abundance < 0.05, "Other", Genus))

ps_ecm_genus_20 <- ps_ecm %>%
  tax_glom(taxrank = "Genus") %>%                    
  transform_sample_counts(function(x) {x / sum(x)}) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.2) %>%                        
  mutate(Genus = ifelse(Abundance < 0.3, "Other", Genus))

ggplot(ps_ecm_family_20)+
  aes(x = vigor, y = Abundance, fill = fct_rev(fct_infreq(Family))) +
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

ggplot(ps_ecm_genus_05)+
  aes(x = height, y = Abundance, fill = fct_rev(fct_infreq(Genus))) +
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
    y = "Proportion of Genera",
    fill = "Genus"
  )+
  theme_few()

######ordination#####
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

ord_bray <- ordinate(ps_rel, method = "NMDS", distance = "bray")

df_ord <- plot_ordination(ps_rel, ord_bray, justDF = TRUE)


ggplot(df_ord, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 1) + 
  coord_fixed(ratio = 1)+
  labs(
    title = "fig"
  )+
  theme_few()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position    = "bottom",
        legend.direction   = "vertical")


plot_ordination(ps_rel, ord_bray, color = "vigor") +   # change color variable
  geom_point(size = 2, alpha = 0.9) +
  stat_ellipse(level = 0.8, alpha = 0.2, geom = "polygon") +
  coord_fixed() +
  labs(title = "title") +
  theme_few() +
  theme(
    plot.title      = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

plot_ordination(ps_rel, ord_bray, color = "height") +   # change color variable
  geom_point(size = 2, alpha = 0.9) +
  #stat_ellipse(level = 0.8, alpha = 0.2, geom = "polygon") +
  coord_fixed() +
  labs(title = "title") +
  theme_few() +
  theme(
    plot.title      = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

