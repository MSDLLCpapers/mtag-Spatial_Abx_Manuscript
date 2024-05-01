library("readr")
library(tidyverse)
library(vegan)
library(ggplot2)
theme_set(theme_bw(base_size = 20))
library(RColorBrewer)
library(ape)
library(ggpubr)
library(rstatix)
library(cowplot)
library(pheatmap)
library(umap)

source("Functions.R")
source("Image_Counts_by_FoV.R") # outputs Counts.summary and raw.data.wide


# Data prep
Counts.df <- raw.data.wide %>% group_by(ASV_Cluster_ID, Mouse, Tissue, Section, FoV) %>% summarize(Count = n()) %>% group_by(Mouse, Tissue, Section, FoV) %>%  mutate(Total = sum(Count))

#join it with metadata 
Counts.df <- left_join(Counts.df, md.images, by = "Mouse")
Counts.df$ASV_Cluster_ID <- as.character(Counts.df$ASV_Cluster_ID)

ASV <- read.csv("../data/asv_lineage_abundance.csv")
ASV$ASV_Cluster_ID <- as.character(ASV$ASV_Cluster_ID)
Counts.df <- left_join(Counts.df, select(ASV, 1:2), by = "ASV_Cluster_ID" )

Counts.df <- Counts.df %>% mutate_at("ASV_Cluster_Scientific_Name", ~replace_na(.,"Bacteria"))

# Species counts per tx
Counts.df<- Counts.df %>% group_by(Mouse, Tissue, Section, FoV) %>% mutate(SpeciesCount = n_distinct(ASV_Cluster_ID))

# Fix factor levels and change H2O to water
Counts.df$Treatment <- gsub("H2O", "Water", Counts.df$Treatment) #changes H20 to water with pattern matching
Counts.df$Treatment <- as.factor(Counts.df$Treatment) # make them as factors
Counts.df$Treatment <- factor(Counts.df$Treatment, levels=c('Water', 'Ampicillin', 'Vancomycin')) # change their order
write.csv(Counts.df, "../data/Counts_df.csv")
