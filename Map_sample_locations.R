## Make a map of the sampling locations for the PH2017 Phormidium sampling


#### Libraries #################################################################
library(tidyverse)
library(ggmap)
library(ggrepel)
library(ggplot2)
source("/Users/KeithBG/R_functions/ggplot_scalebar_north_arrow.R")
################################################################################


#### FILE PATHS ################################################################
dir_input <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","EnvData")
dir_input_script <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data", "Scripts")
dir_output_fig <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data", "Output_figures", "EnvData")

################################################################################


## Read in data
latlong <- read_csv(file.path(dir_input, "PhormMeta17_LatLong_combined.csv")) %>%
            mutate(year= as.character(year))
#latlong17 <- read_csv(file.path(dir_input, "PhormMeta17_LatLong.csv"))

## Remove rows with duplicate sites to reduce the number of points
latlong.map <- latlong %>%
  distinct(acronym, .keep_all = TRUE)

## Make base map of Eel and Russian River watersheds
source(file.path(dir_input_script, "Map_eel_russian.R"))
# Returns R object: PH2017_eel_russian_base_map


PH2017_map_theme <- theme(text= element_text(size= 14),
                          panel.background = element_rect(fill= "light blue"),
                          panel.border = element_rect(color= "black", fill= NA),
                          legend.key= element_rect(fill= "transparent"),
                          plot.background = element_rect(fill= "transparent", color= NA),
                          panel.grid.minor = element_blank(),
                          panel.grid.major = element_blank()
)


## Add information to base map
PH2017_eel_russian_base_map +
  geom_point(data= latlong.map, aes(x= long, y= lat, fill= year, shape= year), color= "black", size= 4) +
  #geom_label_repel(data= latlong.map, aes(x= long, y= lat, label= acronym)) +
  annotate("text", x= -123.3, y= 40.35, label= "Eel River", angle= 310, size= 4) +
  annotate("text", x= -122.7, y= 38.9, label= "Russian River", angle= 310, size= 4) +
  scale_shape_manual(values= c(21, 23, 24)) +
  scale_fill_manual(values= c("lawn green", "tomato", "gold2")) +
  ggtitle("Sampling Sites") +
  PH2017_map_theme

ggsave(last_plot(), filename= "Map_sample_locations.pdf", width= 6, height= 8, units= "in",
                path= dir_output_fig)




