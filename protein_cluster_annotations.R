## Extract KEGG, UNIREF, and UNIPROT annotations for each ORF in the subFamily
## Uses the "Feature Table" download from ggkbase

#### Libraries #################################################################
library(tidyverse)
library(stringr)
library(ggplot2)
library(Hmisc)
################################################################################


#### FILE PATHS ################################################################
dir_input <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "protein_clusters", "output")
source("/Users/KeithBG/R_functions/ggplot_themes.R")
################################################################################


#### READ IN FILES

orf.df <- left_join(read_tsv(file.path(dir_input, "orf2subfamily.tsv")),
                       read_tsv(file.path(dir_input, "orf2genome.tsv")))

clust.cut.df <- read_delim(file.path(dir_input, "subfamily_dendro_cut.txt"), delim= " ") %>%
                 rename(subfamily= leaf_tip) %>%
                 select(subfamily, cluster_cut)

master.df <- left_join(orf.df, clust.cut.df) %>%
               mutate(cluster_cut= ifelse(is.na(cluster_cut), "unassigned", cluster_cut))


#### Feature table from ggkbase
#file.name.feature <- "PH2017_22_RUC_O_B_Oscillatoriales_46_93.ql"
file.name.feature <- "PH2017_11_PCY_U_A_Oscillatoriales_46_44.ql"

genome_ID <- str_replace(file.name.feature, "_Oscill.*", "")

feature.table <- read_delim(file.path(dir_input, file.name.feature), delim= "\t", col_names = FALSE)
colnames.feature.table <- c("orf", "contig_ID", "feature_number", "contig_length", "GC", "coverage", "codon_table", "winning_taxonomy", "begin_position", "end_position", "complement" , "annotation", "orf_taxonomy", "db_info")
colnames(feature.table) <- colnames.feature.table

# Separate annotation column into the 3 different databases
clean.up.annotation.regex <- " Tax.*$| evalue.*$| \\{.*$| \\(.*$|,.*$| bin.*$"
feature.table <- feature.table %>%
                   separate(annotation, into= c("uniref_all", "uniprot_all", "kegg_all"), sep= "__ ") %>%
                   mutate(uniref_anno= capitalize(str_replace(.$uniref_all, clean.up.annotation.regex, "")),
                          uniprot_anno= capitalize(str_replace(.$uniprot_all, clean.up.annotation.regex, "")),
                          kegg_anno= capitalize(str_replace(.$kegg_all, clean.up.annotation.regex, "")),
                          orf_length= end_position - begin_position + 1) %>%
                   select(orf:end_position, orf_length, everything())


#### EXPORT ANNOTATION TABLE FOR EACH CLUSTER
for(clust_ID in unique(master.df$cluster_cut)){

  if(exists("cluster.annotation.output.dir") == FALSE){
    cluster.annotation.output.dir <- file.path(dir_input, "cluster_annotation")
    dir.create(cluster.annotation.output.dir)
  }

master.df %>%
    filter(genome == genome_ID) %>%
    left_join(., feature.table, by= "orf") %>%
    select(-feature_number, -contig_length, -coverage, -GC, -codon_table, -winning_taxonomy, -complement, -db_info) %>%
    filter(cluster_cut == clust_ID) %>%
    write_tsv(path= file.path(cluster.annotation.output.dir, paste(genome_ID, clust_ID, "annotation.txt", sep= "_")))

}

#### WRITE A TABLE OF ANNOTATION COUNTER PER GENOME PER CLUSTER
# Value of 3 means that possibly each database picked up the same annotation
master.df %>%
  filter(genome == genome_ID) %>%
  left_join(., feature.table, by= "orf") %>%
  select(orf, subfamily, genome, cluster_cut, contig_ID, begin_position, end_position, orf_length, uniref_anno, uniprot_anno, kegg_anno) %>%
  gather(key= db, value= anno, uniref_anno:kegg_anno) %>%
  group_by(genome, cluster_cut, anno) %>%
  dplyr::summarize(anno_count= n()) %>%
  write_tsv(path= file.path(dir_input, paste(genome_ID, "annotation_counts.txt", sep= "_")))



