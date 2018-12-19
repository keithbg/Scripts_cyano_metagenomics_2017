## Script to format rp16.py output table and match scaffolds to bins
## using the scaf2bin file downloade from ggkbase

## KBG June, 2018


#### LIBRARIES ################################################################
library(tidyverse)
library(stringr)
###############################################################################


#### FILEPATHS ################################################################
dir_input <- file.path("/Users", "KeithBG", "Documents", "UC Berkeley", "CyanoMeta_NSF", "Metagenomics", "Data", "GenomesData")
dir_output <- file.path("/Users", "KeithBG", "Documents", "UC Berkeley", "CyanoMeta_NSF", "Metagenomics", "Data", "GenomesData", "rp16")
###############################################################################


#### READ IN OUTPUT TABLE
out.table <- read_tsv(file.path(dir_input, "rp16", "rp16_out_table_20180608.txt")) %>%
               rename("scaffold"= `# scaffold`)
scaf2bin <- read_tsv(file.path(dir_input, "rp16", "PH2017_scaf2bin_all.txt"))


#### ADD BIN COLUMN
out.table$Bin <- rep(NA, nrow(out.table))
for(i in 1:nrow(out.table)){
  out.table$Bin[i] <- as.character(scaf2bin[(scaf2bin$scaffold_name == out.table$scaffold[i]), "bin"])
}
out.table <- out.table %>%
               select(Bin, L15:S10)
#write_tsv(out.table, path= file.path(dir_output, "rp16_out_table_20180608_scaf2bin.txt"))
