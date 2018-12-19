## Script to concatenate parsed sequences from rp16.py
## First run rp16_Parser_Keith.R

## KBG June, 2018


#### LIBRARIES ################################################################
library(tidyverse)
library(stringr)
library(seqinr)
library(reshape)
###############################################################################


#### FILEPATHS ################################################################
dir_input <- file.path("/Users", "KeithBG", "Documents", "UC\ Berkeley", "CyanoMeta_NSF", "Metagenomics", "Data", "GenomesData")
dir_output <- file.path("/Users", "KeithBG", "Documents", "UC Berkeley", "CyanoMeta_NSF", "Metagenomics", "Data", "GenomesData", "rp16", "rp16_sequences")
###############################################################################

rp.names <- gsub(".faa", "", list.files(file.path(dir_input, "rp16", "rp16_sequences"), "*.faa"))
fasta.files <- list.files(file.path(dir_input, "rp16", "rp16_sequences"), "*.faa", full.names=TRUE)
#fasta.list <- lapply(fasta.files, function(x) read.fasta(x, seqtype = "AA", as.string = TRUE, set.attributes = FALSE))
#names(fasta.list) <- rp.names

## LOOP TO EXTRACT EACH RP16 SEQUENCE FROM EACH FILE AND COMBINE INTO A FILE OF ALL RP16 SEQ FOR EACH SAMPLE
for(samp in 1:39){
  for(RP in 1:length(fasta.files)){
    sequence <- read.fasta(fasta.files[RP], seqtype = "AA", as.string = TRUE, set.attributes = FALSE)
    samp.list[RP] <- sequence[samp]
  }
  write.fasta(samp.list, names = rp.names, file.out = file.path(dir_output, paste0(names(sequence)[samp], ".rp16.faa")))
}

## CONCATENATE ALL RP16 TOGETHER
sample.files <- list.files(file.path(dir_input, "rp16", "rp16_sequences"), "*.rp16.faa", full.names=TRUE)

for(i in 1:length(sample.files)){
  sample.name <- str_extract(sample.files[i], "PH.*[^.rp16.faa]")

  concat.seq <- read_delim(sample.files[i], delim= "\t", col_names= FALSE) %>%
    filter(str_detect(.$X1, ">") == FALSE)
  colnames(concat.seq)[1] <- paste0(">", sample.name)

  write_tsv(concat.seq, path= file.path(dir_output, paste0(sample.name, ".concat.rp16.faa")))
}
œ
