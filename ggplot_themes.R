## ggplot themes for 2017 PH2017 Metagenomics collections

## Load libraries
library(ggplot2)
library(grid)

#### COLOR/FILL VARIABLES


#### AXIS LABELS


#### OTHER
x.int.line <- geom_hline(yintercept = 0, size= 0.25, color= "black")



## White background, no lines, white strip background
theme_water <- theme(panel.grid = element_blank(),
                   plot.margin = unit(c(1, 1, 1, 1), "cm"),
                   text = element_text(size= 14),
                   plot.background = element_rect(fill = "transparent"), # bg of the plot
                   panel.background = element_rect(fill= "transparent", color="black"),
                   axis.text = element_text(colour="black"),
                   axis.text.x = element_text(angle= 45, vjust= 1, hjust= 1),
                   axis.title.x = element_text(margin = margin(t= 5)),
                   axis.title.y = element_text(margin = margin(r= 15)),
                   legend.background = element_rect(size=0.25, color="black", fill= "transparent"),
                   legend.key = element_blank(),
                   strip.background=element_rect(fill= "transparent", color=NA))

## PCA ORDINATION
theme_pca <- theme(panel.grid = element_blank(),
                     plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                     text = element_text(size= 16),
                     plot.background = element_rect(fill = "transparent"), # bg of the plot
                     panel.background = element_rect(fill= "transparent", color="black"),
                     axis.text = element_text(colour="black"),
                     axis.title.x = element_text(margin = margin(t= 5)),
                     axis.title.y = element_text(margin = margin(r= 15)),
                     legend.background = element_rect(size=0.25, color="black", fill= "transparent"),
                     legend.key = element_blank(),
                     strip.background = element_rect(fill= "transparent", color= NA))

