## Script to make figures of water chemistry data from Phormidium metagenomics
## collections summer 2017


#### Libraries #################################################################
library(ggplot2)
library(scales)
library(tidyverse)
################################################################################


#### FILE PATHS ################################################################
dir_input <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","Scripts")
dir_out_fig <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","Output_figures", "EnvData")
################################################################################

#### PLOTTING VARIABLES ################################################################
watershed.color <- scale_color_discrete(name= "Watershed", labels= c("Creek", "Middle Fork Eel", "Mainstem Eel", "Russian River", "South Fork Eel"))
source("/Users/KeithBG/Documents/UC Berkeley/CyanoMeta_NSF/Metagenomics/Data/Scripts/ggplot_themes.R")
################################################################################


#### SOURCE FORMATTING SCRIPT
source(file.path(dir_input, "env_data_format.R"))

#### MAKE PLOTS

## All analytes faceted together
water %>%
  gather(key= analyte, value= ugL, TDP_ugL:TDN_ugL) %>%
  ggplot(aes(x= site, y= ugL), size= 3) +
  geom_point(aes(color= fork)) +
  facet_grid(analyte~., scales= "free_y") +
  watershed.color +
  theme_water


## TDP
summary(water$TDP_ugL)
ggplot(water, aes(x= NA, y= TDP_ugL)) +
  x.int.line +
  geom_boxplot() +
  geom_point(alpha= 0.3, size= 3) +
  ylim(0, 50) +
  labs(x= "", y= expression(paste("Total dissolved phosphorus (", "\U00B5", "g / L)"))) +
  scale_x_discrete(breaks= NULL) +
  theme_water

ggplot(water, aes(x= site, y= TDP_ugL)) +
  x.int.line +
  geom_point(aes(color= fork), size= 3) +
  ylim(0, 50) +
  labs(x= "Site", y= expression(paste("Total dissolved phosphorus (", "\U00B5", "g / L)"))) +
  scale_x_discrete(limits= water$site) +
  watershed.color +
  theme_water
ggsave(last_plot(), file= "water_TDP.pdf", height= 6, width= 8, units= "in", path= dir_out_fig, device = cairo_pdf)


## TDN
summary(water$TDN_ugL)
ggplot(water, aes(x= NA, y= TDN_ugL)) +
  geom_boxplot() +
  geom_point(alpha= 0.3, size= 3) +
  ylim(0, 200) +
  labs(x= "", y= expression(paste("Total dissolved nitrogen (", "\U00B5", "g / L)"))) +
  scale_x_discrete(breaks= NULL) +
  theme_water

ggplot(water, aes(x= site, y= TDN_ugL)) +
  geom_point(aes(color= fork), size= 3) +
  ylim(0, 200) +
  labs(x= "Site", y= expression(paste("Total dissolved nitrogen (", "\U00B5", "g / L)"))) +
  scale_x_discrete(limits= water$site) +
  watershed.color +
  theme_water
ggsave(last_plot(), file= "water_TDN.pdf", height= 6, width= 8, units= "in", path= dir_out_fig, device = cairo_pdf)

## NO3
summary(water$NO3_ugL)
ggplot(water, aes(x= NA, y= NO3_ugL)) +
  geom_boxplot() +
  geom_point(alpha= 0.3, size= 3) +
  ylim(0, 70) +
  labs(x= "", y= expression(paste("Nitrate (", "\U00B5", "g / L)"))) +
  scale_x_discrete(breaks= NULL) +
  theme_water

ggplot(water, aes(x= site, y= NO3_ugL)) +
  geom_point(aes(color= fork), size= 3) +
  ylim(0, 70) +
  labs(x= "Site", y= expression(paste("Nitrate (", "\U00B5", "g / L)"))) +
  scale_x_discrete(limits= water$site) +
  watershed.color +
  theme_water
ggsave(last_plot(), file= "water_NO3.pdf", height= 6, width= 8, units= "in", path= dir_out_fig, device = cairo_pdf)

## NH4
summary(water$NH4_ugL)
ggplot(water, aes(x= NA, y= NH4_ugL)) +
  geom_boxplot() +
  geom_point(alpha= 0.3, size= 3) +
  ylim(0, 10) +
  labs(x= "", y= expression(paste("Ammonium (", "\U00B5", "g / L)"))) +
  scale_x_discrete(breaks= NULL) +
  theme_water

ggplot(water, aes(x= site, y= NH4_ugL)) +
  geom_point(aes(color= fork), size= 3) +
  ylim(0, 10) +
  labs(x= "Site", y= expression(paste("Ammonium (", "\U00B5", "g / L)"))) +
  scale_x_discrete(limits= water$site) +
  watershed.color +
  theme_water
ggsave(last_plot(), file= "water_NH4.pdf", height= 6, width= 8, units= "in", path= dir_out_fig, device = cairo_pdf)


## N:P molar ration
summary(water$NP)
ggplot(water, aes(x= NA, y= NP)) +
  geom_boxplot() +
  geom_point(alpha= 0.3, size= 3) +
  ylim(0, 25) +
  labs(x= "", y= "Molar N:P") +
  scale_x_discrete(breaks= NULL) +
  theme_water

ggplot(water, aes(x= site, y= NP)) +
  geom_point(aes(color= fork), size= 3) +
  ylim(0, 25) +
  labs(x= "Site", y= "Molar N:P") +
  scale_x_discrete(limits= water$site) +
  watershed.color +
  theme_water
ggsave(last_plot(), file= "water_NP.pdf", height= 6, width= 8, units= "in", path= dir_out_fig, device = cairo_pdf)





