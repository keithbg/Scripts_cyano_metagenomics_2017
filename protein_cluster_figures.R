## Move protein cluster fasta files into their own folder for each cluster
## Output summary table of genomes and associated clusters

#### Libraries #################################################################
library(tidyverse)
library(stringr)
library(ggplot2)
library(RColorBrewer)
################################################################################


#### FILE PATHS ################################################################
dir_input <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "protein_clusters", "PHall_20181008")
dir_output <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "protein_clusters", "PHall_20181008", "figures")
source("/Users/KeithBG/R_functions/ggplot_themes.R")
################################################################################


#### PHEATMAP FIGURE ###########################################################
df <- suppressMessages(read_tsv(file.path(dir_input, "family2genome_matrix.tsv"))) %>%
  rename(genome=X1) %>%
  as.data.frame(.)
rownames(df) <- df$genome
prot.mat <- as.matrix(df[, -1])

## Calculate  distances among rows and columns
dist.rows <- vegan::vegdist(prot.mat, method= "jaccard")
dist.cols <- vegan::vegdist(t(prot.mat), method= "jaccard")

## Cut column clusters
col_hclust <- hclust(dist.cols, method= "ward.D")
col_clusters <- cutree(col_hclust, k= 10) %>%
                 data.frame(family= names(.), cluster= ., stringsAsFactors = FALSE) %>%
                 mutate(cluster= as.character(cluster))
row.names(col_clusters) <- col_clusters$family
col_clusters <- col_clusters %>%
  select(-family)

## Cluster annotation colors
ann_colors <- list(cluster= brewer.pal(n= 10, name= "Paired"))
names(ann_colors[[1]]) <- c("1", "10", "2", "3", "4", "5", "6", "7", "8", "9")

fig.filename <- file.path(dir_output, "heatmap_annotated.pdf")

pheatmap::pheatmap(prot.mat,
         color= c("light gray", "black"),
         breaks= c(0, 0.5, 1),
         legend_breaks= c(0, 1),
         legend_labels= c("Absent", "Present"),
         #scale= "none",
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         clustering_distance_rows= dist.rows,
         clustering_distance_cols= dist.cols,
         clustering_method= "ward.D",
         cutree_cols = 10,
         cutree_rows= 5,
         annotation_col = col_clusters,
         annotation_colors = ann_colors,
         annotation_names_col = TRUE,
         #margins= (c(10, 10)),
         #display_numbers = TRUE,
         show_rownames= TRUE,
         show_colnames=  FALSE,
         filename= fig.filename,
         height= 25,
         width= 40
)

#### ANNOTATION SUMMARIES ####

ani.species <- read_tsv("/Users/KeithBG/Documents/UC Berkeley/CyanoMeta_NSF/Metagenomics/Data/GenomesData/genome_ani_species.tsv")

genome.df <- left_join(read_tsv(file.path(dir_input, "orf2subfamily.tsv")),
          read_tsv(file.path(dir_input, "orf2genome.tsv"))) %>%
  rename(family= subfamily) %>%
  left_join(., read_tsv(file.path(dir_input, "heatmap_output_20181008", "family_clusters.tsv"))) %>%
  left_join(., ani.species) %>%
  rename(cluster.fam= cluster) %>%
  mutate(cluster.fam = as.character(cluster.fam)) %>%
  arrange(family) %>%
  select(-orf) %>%
  distinct()

## Count the number of families in each cluster
cluster.count.totals <- genome.df %>%
  select(-genome, -ani_species) %>%
  distinct() %>%
  count(cluster.fam) %>%
  rename(fam_total= n) %>%
  mutate(cluster.fam= ifelse(is.na(cluster.fam), "unassigned", cluster.fam))

## Calculate the proportion of each family cluster is in each genome
cluster.counts <- genome.df %>%
  mutate(cluster.fam= ifelse(is.na(cluster.fam), "unassigned", cluster.fam)) %>%
  count(genome, cluster.fam) %>%
  left_join(., cluster.count.totals) %>%
  mutate(prop_cluster= round((n / fam_total), 4))

ggplot(data= cluster.counts, aes(x= cluster.fam, y= prop_cluster)) +
  geom_bar(aes(fill= cluster.fam), stat= "identity") +
  labs(x= "Cluster Number", y= "Proportion", title= "Proportion of protein families in each genome") +
  scale_fill_discrete(breaks= c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "unassigned"),
                      labels= c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "U-A"),
                      name= "Protein\nFamily\nCluster") +
  scale_x_discrete(limits= c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "unassigned"),
                   labels= c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "U-A")) +
  scale_y_continuous(expand= c(0, 0)) +
  facet_wrap(~genome, nrow= 11) +
  theme_kbg
ggsave(last_plot(), filename= "family_clusters_per_genome2.pdf", width= 30, height= 30, units= "in", path= dir_output, device= cairo_pdf)

## Calculate the proportion of each family cluster is in each ANI species

cluster.counts.ani <- genome.df %>%
  mutate(cluster.fam= ifelse(is.na(cluster.fam), "unassigned", cluster.fam)) %>%
  filter(cluster.fam != "unassigned") %>%
  count(genome, cluster.fam) %>%
  left_join(., cluster.count.totals) %>%
  mutate(prop_cluster= round((n / fam_total), 4)) %>%
  left_join(., ani.species) %>%
  group_by(ani_species, cluster.fam) %>%
  summarize(
    n= length(prop_cluster),
    mean_prop_cluster= mean(prop_cluster),
    sd_prop_cluster= sd(prop_cluster)
  )



ggplot(data= cluster.counts.ani, aes(x= ani_species,
                                     y= mean_prop_cluster,
                                     ymin= mean_prop_cluster - 0.01,
                                     ymax= mean_prop_cluster + sd_prop_cluster)) +
  geom_errorbar(width = 0.5) +
  geom_bar(aes(fill= cluster.fam), stat= "identity") +
  labs(x= "ANI Species", y= "Mean proportion (± sd)", title= "Proportion of protein families in each ANI species") +
  scale_fill_discrete(breaks= c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "unassigned"),
                      labels= c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "U-A"),
                      name= "Protein\nFamily\nCluster") +
  #scale_x_discrete(limits= c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "unassigned"),
   #                labels= c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "U-A")) +
  scale_y_continuous(expand= c(0, 0)) +
  facet_grid(.~cluster.fam) +
  theme_kbg
ggsave(last_plot(), filename= "family_clusters_per_ANIspecies.pdf", width= 10, height= 8, units= "in", path= dir_output, device= cairo_pdf)

?scale_fill_discrete
