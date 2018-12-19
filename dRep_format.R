## Script to input and format dRep analysis of PH2017 genomes


#### Libraries #################################################################
library(tidyverse)
library(stringr)
library(vegan)
library(pheatmap)
################################################################################


#### FILE PATHS ################################################################
dir_input <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "dRep","output_20180927")
dir_out_fig <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "dRep","output_20180927")
################################################################################



## READ dRep OUTPUT
df <- read_csv(file.path(dir_input, "Ndb.csv")) %>%
        select(querry, reference, ani)


## Transform into a matrix
ani.mat <- df %>%
             spread(key= querry, value= ani) %>%
             select(-reference) %>%
             as.matrix(.)
colnames(ani.mat) <- str_replace(colnames(ani.mat), ".contigs.fa", "")
rownames(ani.mat) <- colnames(ani.mat)


## Remove NAs (these are generated because only samples within the same primary cluster have an NA value calculated)
ani.mat[is.na(ani.mat) == TRUE] <- 0.8

# Calculate Bray-Curtis distance
ani.dist.row <- vegdist(ani.mat, method= "bray")
ani.dist.col <- vegdist(t(ani.mat), method= "bray")

## Identify the 4 species
ani.clusters <- hclust(ani.dist.row, method = "ward.D")
ani.species <- data.frame(ani_species= cutree(ani.clusters, k= 4))
ani.species$genome <- str_replace(row.names(ani.species), ".fa", "")

## Change the species numbers to match the PH2015 publication and write a TSV file
ani.species %>%
  mutate(ani_species= ifelse(ani_species == 1, 4,
                             ifelse(ani_species == 2, 1,
                                    ifelse(ani_species == 3, 2, 3)))) %>%
  select(genome, ani_species) %>%
  write_tsv(path= file.path(dir_out_fig, "ani_species.tsv"))


## Color ramp
ani.breaks <- c(0.8, 0.85, seq(0.90, 1, by= 0.005))
legend.breaks <- c(0.8, 0.85, seq(0.90, 1, by= 0.02))

ani.colors <- c("gray", "lightblue",  colorRampPalette(c("cornsilk", "snow", "goldenrod", "salmon", "tomato", "firebrick"))(21))

## Saving heatmap
fig.filename <- file.path(dir_out_fig, "ani_heatmap_20181213.pdf")


pheatmap(ani.mat,
         color= ani.colors,
         breaks= ani.breaks,
         legend_breaks= legend.breaks,
         scale= "none",
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         clustering_distance_rows= ani.dist.row,
         clustering_distance_cols= ani.dist.col,
         clustering_method= "ward.D",
         margins= (c(10, 10)),
         show_rownames= TRUE,
         show_colnames=  FALSE,
         #display_numbers = TRUE,
         #fontsize_number = 4,
         filename= fig.filename,
         height= 10,
         width= 15
)


?pheatmap




