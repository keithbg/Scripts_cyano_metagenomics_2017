## Script to input and format dRep analysis of PH2017 genomes


#### Libraries #################################################################
library(tidyverse)
library(stringr)
library(vegan)
library(pheatmap)
library(ggplot2)
################################################################################


#### FILE PATHS ################################################################
dir_input <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Microcoleus_Analysis","GenomesData", "dRep","output_20180927")
dir_out_fig <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Microcoleus_Analysis","GenomesData", "dRep","output_20180927")
dir_out_tables <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Microcoleus_Analysis","GenomesData", "dRep","Output_tables")

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
ani.dist.row <- vegdist(ani.mat, method= "euclidean")
ani.dist.col <- vegdist(t(ani.mat), method= "euclidean")

## Identify the 4 species
ani.clusters <- hclust(ani.dist.row, method = "ward.D2")
ani.species <- data.frame(ani_species= cutree(ani.clusters, k= 4))
ani.species$genome <- str_replace(row.names(ani.species), ".fa", "")

## Change the species numbers to match the PH2015 publication and write a TSV file
ani.species <- ani.species %>%
  mutate(ani_species= ifelse(ani_species == 1, 4,
                             ifelse(ani_species == 2, 1,
                                    ifelse(ani_species == 3, 2, 3)))) %>%
  select(genome, ani_species)
write_tsv(ani.species, path= file.path(dir_out_tables, "ani_species.tsv"))


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

## SPECIES 1 SUB-CLUSTERS
rownames(ani.mat) 
colnames(ani.mat)

rows.sp1 <- rownames(ani.mat) %in% str_c(filter(ani.species, ani_species == "1")$genome, ".fa")
cols.sp1 <- colnames(ani.mat) %in% str_c(filter(ani.species, ani_species == "1")$genome, ".fa")

ani.mat.sp1 <- ani.mat[rows.sp1, cols.sp1]

ani.mat.sp1 %>% 
  as_tibble() %>% 
  mutate(row_name= colnames(.)) %>% 
  pivot_longer(names_to = "col_name", values_to = "ani", -row_name) %>% 
  write_tsv(., path= file.path(dir_out_tables, "sp1_ani.tsv"))

  

# Calculate Euclidean distance
ani.dist.sp1.row <- vegdist(ani.mat.sp1, method= "euclidean")

## Identify the sub-clusters
ani.sp1.clusters <- hclust(ani.dist.sp1.row, method = "ward.D2")


k_cut <- 9
plot(ani.sp1.clusters)
rect.hclust(ani.sp1.clusters, k= k_cut)

sp1.clusters <- data.frame(sp1_clustID= cutree(ani.sp1.clusters, k= k_cut)) %>% 
  mutate(genome= str_replace(row.names(.), ".fa", ""),
         site= str_replace(genome, "_Oscill.*$", "")) %>% 
  mutate(site= str_replace(site, "_s25", ""))

write_tsv(sp1.clusters, path= file.path(dir_out_tables, "sp1_clusters.tsv"))



pheatmap(ani.mat.sp1,
         color= viridis::viridis_pal(option= "magma")(30),
         #breaks= c(0.995),
         #legend_breaks= sp1.legend.breaks,
         scale= "none",
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         clustering_distance_rows= ani.dist.sp1.row,
         clustering_distance_cols= ani.dist.sp1.row,
         clustering_method= "ward.D2",
         cutree_cols = 10,
         cutree_rows = 10,
         margins= (c(10, 10)),
         show_rownames= TRUE,
         show_colnames=  FALSE,
         display_numbers = TRUE,
         number_format = "%.3f",
         fontsize_number = 4,
         #filename= fig.filename,
         #height= 10,
         #width= 15
)



ani.mat.sp1 %>% 
  as_tibble() %>% 
  mutate(row_name= colnames(.)) %>% 
  pivot_longer(names_to = "col_name", values_to = "ani", -row_name) %>% 
  ggplot() +
  geom_histogram(aes(x= ani), boundary= 0, binwidth = 0.001, fill= "black", color= "gray50") +
  theme_bw()
  