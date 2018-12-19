## Script to format rp16.py output table and match scaffolds to bins
## using the scaf2bin file downloade from ggkbase

## KBG June, 2018


#### LIBRARIES ################################################################
library(tidyverse)
library(stringr)
library(seqinr)
library(reshape)
###############################################################################


#### FILEPATHS ################################################################
dir_input <- file.path("/Users", "KeithBG", "Documents", "UC\ Berkeley", "CyanoMeta_NSF", "Metagenomics", "Data", "GenomesData")
dir_output <- file.path("/Users", "KeithBG", "Documents", "UC Berkeley", "CyanoMeta_NSF", "Metagenomics", "Data", "GenomesData", "rp16")
###############################################################################

# Remove * at end of each sequence with sed in Terminal
#sed 's/\*//' PH2017_genomes_concat.faa > PH2017_genomes_concat.clean.faa

#### READ IN FASTA FILE
sequences <- read.fasta(file.path(dir_input, "rp16", "PH2017_genomes_concat.clean.faa"), seqtype = "AA", as.string = TRUE, set.attributes = FALSE)
names <- read_tsv(file.path(dir_input, "rp16", "rp16_out_table_20180608_scaf2bin.txt"))


### Remove Genomes with < --minrp option
minrp <- 8
names <- names[(16 - rowSums(is.na(names[,2:17]))) >= minrp,]


## TRANSFORM TO LONG FORMAT
names_melt <- gather(names, key= "RP", value= "ID", L15:S10) %>%
                mutate(ID= ifelse(.$ID == "-", NA, .$ID))
names_melt <- names_melt[complete.cases(names_melt),]

## WRITE A FILE WITH ALL THE AA SEQUENCES FOR EACH RP IN A SINGLE FASTA
for (i in as.vector(unique(names_melt$RP), mode = "character")){
  tmp_df <- subset(names_melt, RP == i)
  tmp_seq <- sequences[names(sequences) %in% tmp_df$ID]
  names(tmp_seq) <- pull(tmp_df[match(names(tmp_seq), tmp_df$ID), 1])
  write.fasta(tmp_seq, names = names(tmp_seq), file.out = paste0(dir_output,"/",i,".faa"))
}





