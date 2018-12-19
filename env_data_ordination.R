## Environmental covariate ordinations from Phormidium metagenomics
## collections summer 2017.
## Keith Bouma-Gregson


#### Libraries #################################################################
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggrepel)
################################################################################


#### FILE PATHS ################################################################
dir_input <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data", "Scripts")
dir_out_fig <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data", "Output_figures", "EnvData")
dir_out_table <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data", "EnvData")
################################################################################

#### RUN FORMATTING SCRIPT TO GENERATE ORDINATION MATRIX
# Variables have already been normalized by mean and variance
source(file.path(dir_input, "env_data_format.R"))
head(env.ord)

#### PRINCIPAL COMPONENTS ANALYSIS

## ORDINATE
env.pca <- rda(env.ord)
summary(env.pca)
saveRDS(env.pca, file= file.path(dir_out_table, "env.data.pca"))

## EXTRACT ORDINATION SCORES AND LOADINGS
pca.site.scores <- as.data.frame(scores(env.pca, choices= c(1, 2, 3,4))$sites)
row.names(pca.site.scores) <- str_replace(row.names(pca.site.scores), "PH2017_", "")
pca.covariate.scores <- as.data.frame(scores(env.pca, choices= c(1, 2, 3, 4))$species)

ord.sites <- row.names(pca.site.scores)
ord.covariates <- row.names(pca.covariate.scores)
ord.covariates2 <- c("Temp", "pH", "DO", "Cond", "Alk", "Depth", "Flow", "Canopy", "TDP", "NH4", "NO3", "DOC", "TDN", "NP")


#### MAKE FIGURES
source("/Users/KeithBG/Documents/UC Berkeley/CyanoMeta_NSF/Metagenomics/Data/Scripts/ggplot_themes.R")

## PCA 1 & PCA 2
ggplot() +
  geom_hline(yintercept = 0, color= "gray") +
  geom_vline(xintercept= 0, color= "gray") +
  geom_point(data= pca.site.scores, aes(x= PC1, y= PC2), size= 2, shape= 15) +
  geom_text_repel(data= pca.site.scores, aes(x= PC1, y= PC2, label= ord.sites), color= "black", size= 3) +
  geom_segment(data= pca.covariate.scores, aes(x= 0, xend= PC1, y=0, yend= PC2),
               arrow= arrow(length = unit(0.2, "cm")), size= 0.75,
               color= "black", alpha= 0.3) +
  geom_text_repel(data= pca.covariate.scores, aes(x= PC1, y= PC2, label= ord.covariates2), color= "Tomato", size= 5) +
  labs(x= "PC1 (29.2%)", y= "PC2 (20.5%)") +
  coord_equal() +
  theme_pca
#ggsave(last_plot(), file= "env_data_PCA_1&2.pdf", height= 6, width= 6, units= "in", path= dir_out_fig, device = cairo_pdf)

## PCA 2 & PCA 3
ggplot() +
  geom_hline(yintercept = 0, color= "gray") +
  geom_vline(xintercept= 0, color= "gray") +
  geom_point(data= pca.site.scores, aes(x= PC2, y= PC3), size= 2, shape= 15) +
  geom_text_repel(data= pca.site.scores, aes(x= PC2, y= PC3, label= ord.sites), color= "black", size= 3) +
  geom_segment(data= pca.covariate.scores, aes(x= 0, xend= PC2, y=0, yend= PC3),
               arrow= arrow(length = unit(0.2, "cm")), size= 0.75,
               color= "black", alpha= 0.3) +
  geom_text_repel(data= pca.covariate.scores, aes(x= PC2, y= PC3, label= ord.covariates2), color= "Tomato", size= 5) +
  labs(x= "PC2 (20.5%)", y= "PC3 (15.6%)") +
  coord_equal() +
  theme_pca
#ggsave(last_plot(), file= "env_data_PCA_2&3.pdf", height= 6, width= 6, units= "in", path= dir_out_fig, device = cairo_pdf)

