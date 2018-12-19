#!/bin/Rscript

## Script to make heatmap plot from protein clustering script from
## /home/meheurap/proteinClusteringPipeline/scripts/subfamilies.py script
## Also cuts and extracts protein cluster information

#### LIBRARIES #################################################################
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(vegan))
suppressPackageStartupMessages(require(pheatmap))
################################################################################



#### OPTIONS #################################################################
options(warn = -1) # suppress warnings

option_list = list(
  make_option(c("-f", "--file"), action="store", default=NA, type='character',
              help="family2genome.tsv input file (must be tab separated)"),
  make_option(c("-o", "--output"), action="store", default="output", type='character',
              help="Output folder name"),
  make_option(c("-p", "--path"), action="store", default=".",
              help="Path to files (default= current directory)"),
  make_option(c("-r", "--row_cut"), action="store", default= 0,
              help="Number of groups to cut rows into (Default = 1)"),
  make_option(c("-c", "--col_cut"), action="store", default= 0,
              help="Number of groups to cut columns into (Default = 1)"),
  make_option(c("-m", "--move_fasta"), action="store", default= FALSE,
              help="Move protein family fasta files to unique folders for each cluster (Default = FALSE)"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

################################################################################


#### FILE PATHS ################################################################
#dir_input <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "protein_clusters")
#dir_out <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "protein_clusters", "PHall", "output")

dir_input <- file.path(opt$path)
dir_out <- file.path(opt$path, opt$output)
dir.create(dir_out, showWarnings = FALSE)
################################################################################


## READ IN MATRIX ##############################################################
df <- suppressMessages(read_tsv(file.path(dir_input, opt$file))) %>%
       rename(genome=X1) %>%
       as.data.frame(.)
rownames(df) <- df$genome
prot.mat <- as.matrix(df[, -1])


#### PHEATMAP FIGURE ###########################################################
## Calculate  distances among rows and columns
dist.rows <- vegdist(prot.mat, method= "jaccard")
dist.cols <- vegdist(t(prot.mat), method= "jaccard")

## Make heatmap
fig.filename <- file.path(dir_out, "family_genome_heatmap.pdf")

pheatmap(prot.mat,
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
         cutree_cols = opt$col_cut,
         cutree_rows= opt$row_cut,
         #margins= (c(10, 10)),
         #display_numbers = TRUE,
         show_rownames= TRUE,
         show_colnames=  FALSE,
         filename= fig.filename,
         height= 25,
         width= 40
)
pheatmap

## Cut rows and column clusters
row_clusters <- hclust(dist.rows, method= "ward.D") %>%
  cutree(., k= opt$row_cut) %>%
  data.frame(genome= names(.), cluster= ., stringsAsFactors = FALSE)
write.table(row_clusters, file= file.path(dir_out, "genome_clusters.tsv"), sep= "\t", quote= FALSE, row.names= FALSE)

col_clusters <- hclust(dist.cols, method= "ward.D") %>%
  cutree(., k= opt$col_cut) %>%
  data.frame(family= names(.), cluster= ., stringsAsFactors = FALSE)
write.table(col_clusters, file= file.path(dir_out, "family_clusters.tsv"), sep= "\t", quote= FALSE, row.names= FALSE)


#### MOVE FASTA FILES ##########################################################
move_fasta_files <- function(){
  #clust.cut.df <- read_delim(file.path(dir_out, "family_clusters.tsv"), delim= " ")
  fasta.files <- list.files(file.path(dir_input, "subfamiliesFasta"), "*.fa")

  for(clust in unique(col_clusters$cluster)){
    cluster.output.dir <- file.path(dir_out, "family_clusters")
    dir.create(cluster.output.dir)

    # Create folder for the files for each cluster
    individual.clust.output.dir <- file.path(cluster.output.dir, paste0("cluster_", clust))
    dir.create(individual.clust.output.dir)

    # Extract the subfamilies in each cluster
    clust.family.list <- col_clusters %>%
      filter(cluster == clust) %>%
      select(family) %>%
      pull()

    # Generate pathways for files in the specific cluster
    clust.fasta.list <- str_replace(fasta.files, "\\.fa", "") %in% clust.family.list %>%
      fasta.files[.] %>%
      paste(file.path(dir_input, "subfamiliesFasta"), ., sep ="/")

    # Copy files to a new folder
    file.copy(clust.fasta.list, individual.clust.output.dir)
  }
}

if(opt$move_fasta == TRUE){
  move_fasta_files()
}





