## Script to make heatmap plot from protein clustering script from
## /home/meheurap/proteinClusteringPipeline/scripts/subfamilies.py script
## Also cuts and extracts protein cluster information

#### Libraries #################################################################
library(tidyverse)
library(vegan)

library(ape)
library(gplots)
library(phytools)
library(pheatmap)
################################################################################


#### FILE PATHS ################################################################
dir_input <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "protein_clusters", "PHall")
dir_out <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "protein_clusters", "PHall", "output")
################################################################################


## READ IN MATRIX
df <- read_tsv(file.path(dir_input, "family2genome_matrix.tsv")) %>%
       rename(genome=X1) %>%
       as.data.frame(.)
rownames(df) <- df$genome
prot.mat <- as.matrix(df[, -1])


## Make heatmap
clustering.method <- "ward.D"
prot.heatmap <- heatmap.2(prot.mat,
          Rowv = TRUE,
          Colv=TRUE,
          distfun=function(x) dist(x,method = 'binary'),
          hclustfun=function(x) hclust(x, method = clustering.method),
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(6,12),     # widens margins around plot
          breaks= c(0, 0.5, 1),
          col=c("gray", "black"),       # use on color palette defined earlier
          key = FALSE,
          #RowSideColors = liste_row_color,
          #ColSideColors = liste_col_color,
          xlab = "Families",
          ylab = "Genomes",
          cexRow = 1,
          labCol= FALSE,
          dendrogram="both")


save.heatmap <- function(data, filename, path){
  pdf(file= file.path(path, filename), bg= "transparent", height= 20, width= 25)

  heatmap.2(data,
            Rowv = TRUE,
            Colv=TRUE,
            distfun=function(x) dist(x,method = 'binary'),
            hclustfun=function(x) hclust(x, method = clustering.method),
            notecol="black",      # change font color of cell labels to black
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(6,12),     # widens margins around plot
            breaks= c(0, 0.5, 1),
            col=c("gray", "black"),       # use on color palette defined earlier
            key = FALSE,
            #RowSideColors = liste_row_color,
            #ColSideColors = liste_col_color,
            xlab = "Families",
            ylab = "Genomes",
            cexRow = 1,
            labCol= FALSE,
            dendrogram="both")
  dev.off()

}
save.heatmap(data= prot.mat, filename= "heatmap_ward.pdf", path= dir_out)


# Extract dendrograms for rows and columns
prot.heatmap$rowDendrogram %>%
  as.hclust(.) %>%
  as.phylo(.) %>%
  write.tree(., file.path(dir_out, 'dendrogram_genome_ward.nwk'))

# COLUMN EXTRACTION NOT CORRECT BUG IN THE PHYLO DATA FRAME
# SO IT CANT BE READ BY WRITE.TREE
# SEE ATTEMPT TO FIX BELOW
# prot.heatmap$colDendrogram %>%
#   as.hclust(.) %>%
#   as.phylo(.) %>%
#   write.tree(., file.path(dir_out, 'dendrogram_subfamily_ward5.nwk'), digits= 20)
#
# col.phylo <- prot.heatmap$colDendrogram %>%
#   as.hclust(.) %>%
#   as.phylo(.)
# write.tree(col.phylo, file.path(dir_out, 'dendrogram_subfamily_ward6.nwk'), digits= 20)


# Write row order and column order
df.sorted <- df[prot.heatmap$rowInd, prot.heatmap$colInd]
write(rownames(df.sorted), file = file.path(dir_out, "genome_order.txt"), sep='\n')
write(colnames(df.sorted), file = file.path(dir_out, "subfamily_order.txt"), sep='\n')

## Cut subfamily dendrogram
col.cut <- prot.heatmap$colDendrogram %>%
  as.hclust(.) %>%
  cutree(., k= 6) %>%
  as.data.frame(.)
colnames(col.cut)[1] <- "cluster_cut"
col.cut$cluster_cut <- paste0("clust_", col.cut$cluster_cut)
col.cut$leaf_tip <- rownames(col.cut)
clust.num <- length(unique(col.cut$cluster_cut))

str(col.cut)
head(col.cut)
table(col.cut)
col.cut[col.cut ==7]
col.cut[col.cut$cluster_cut == "7", ]

# Set colors for cluster groups
color.match <- data.frame(color= rainbow(clust.num), cluster_cut= paste0("clust_", seq(1:clust.num)), stringsAsFactors = FALSE)

# Write a file
col.cut %>%
  left_join(., color.match) %>%
  select(leaf_tip, color, cluster_cut) %>%
  write_delim(path= file.path(dir_out, "subfamily_dendro_cut.txt"), delim= " ")

#### OLD CODE BELOW


#### PHEATMAP FIGURE ###########################################################

## Calculate  distances among samples
dist.rows <- vegdist(prot.mat, method= "jaccard")
dist.cols <- vegdist(t(prot.mat), method= "jaccard")
?vegdist
# #str(dist.rows)
#
# ## Color ramp
 fig.breaks <- c(0, 0.5, 1)
# #length(ani.breaks)
# #ani.colors <- c(rep("gray", 3), colorRampPalette(c("snow", "tomato", "firebrick"))(7))
#
# ## Saving heatmap
# library(pheatmap)
fig.filename <- file.path(dir_out, "heatmap_cut_5x8.pdf")
#
#
pheatmap(prot.mat,
         color= c("light gray", "black"),
         breaks= fig.breaks,
         legend_breaks= c(0, 1),
         legend_labels= c("Absent", "Present"),
         #scale= "none",
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         clustering_distance_rows= dist.rows,
         clustering_distance_cols= dist.cols,
         clustering_method= "ward.D",
         cutree_cols = 8,
         cutree_rows= 5,
         #margins= (c(10, 10)),
         #display_numbers = TRUE,
         show_rownames= TRUE,
         show_colnames=  FALSE,
         filename= fig.filename,
         height= 25,
         width= 40
)

test <- hclust(dist.rows, method= "ward.D")
plot(test)

cutree(test, k= 5)
#### DEBUGING WRITE.TREE BELOW #################################################
################################################################################



# my.row <- prot.heatmap$rowDendrogram %>%
#   as.hclust(.) %>%
#   as.phylo(.)
#
# plotTree(my.phylo, ftype="i",fsize=0.6,lwd=1)
#
#
# str(my.row)
# summary(my.row$edge[, 1])
# summary(my.row$edge[, 2])
#
# summary(my.phylo$edge[, 1])
# my.phylo$edge[my.phylo$edge[, 1] ==0, ]
#
# which(my.phylo$edge[, 1] ==0)
# head(my.phylo$edge)
# summary(my.phylo$edge[, 2])
#
#
#
#
#
# my.phylo <- prot.heatmap$colDendrogram %>%
#   as.hclust(.) %>%
#   as.phylo(.)
#
# write.tree(my.phylo, file.path(dir_out, 'dendrogram_subfamily_ward5.nwk'), digits= 20)
# str(my.phylo)
#
# table(my.phylo$edge < 0.01)
# table(children < 0.01)
#
# parent[which(parent < 0.1)]
#
# parent[which(parent < 0.1)]
# head(parent)
#
# head(my.phylo$edge)
# length(my.phylo$Nnode)
#
# length(my.phylo$edge[, 1])
#
# #my.phylo$edge2 <- my.phylo$edge
# #my.phylo$edge[, 1] <- ifelse(my.phylo$edge[, 1] == 0, 0.00001, my.phylo$edge[, 1])
# #my.phylo$edge[, 2] <- ifelse(my.phylo$edge[, 2] == 0, 0.00001, my.phylo$edge[, 2])
#
# n <- length(my.phylo$tip.label)
# parent <- my.phylo$edge[, 1]
# parent <- parent + 1
# children <- my.phylo$edge[, 2]
# children <- children + 1
# kids <- vector("list", n + my.phylo$Nnode)
# #for (i in 1:length(parent)){
#
# i <- length(parent)-7
# #for (i in 7500:(length(parent)-7)){
#
# for (i in 1:length(parent)){
#   #if(parent[i] != 0){
#   kids[[parent[i]]] <- c(kids[[parent[i]]], children[i])
#   #}
# }
# max(parent)
# max(children)
# str(kids[[parent[i]]])
# dim(my.phylo$edge)
# my.phylo$edge[n, ]
# parent[n]
#
# children[n]
#
#
# n <- seq(length(parent)-10, length(parent))
#
# write.tree(my.phylo, file.path(dir_out, 'test.nwk'))
#
# head(parent)
# head(children)
#
#


#Error in kids[[parent[i]]] :
#  attempt to select less than one element in integerOneIndex



