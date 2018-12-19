## Script to manipulate and summarize xmfa output from mauve alignment
## KBG July 2018

## Input: mauve .xmfa file

## Output:
## 1) fasta files for each localized collinear blocks (LCBs)
## 2) summary of length of each LCB, gaps per lcb, and number of genomes with an LCB
## 3) table with rows for each LCB and column for each genome. Cells are populated with
##    1's if the genome contains the LCB and 0's if the LCB was absent in a genome

summarize_LCBs <- function(xmfa, path, write_fasta= FALSE){

#### Libraries #################################################################
  library(plyr)
  library(tidyverse)
  library(stringr)
  ################################################################################


  #### FILE PATHS ################################################################
  dir_input <- path
  file.name <- xmfa
  if(dir.exists(file.path(dir_input, "lcb_summary_output")) == TRUE){
    dir_output <- file.path(dir_input, "lcb_summary_output")
  } else {
    dir.create(file.path(dir_input, "lcb_summary_output"))
    dir_output <- file.path(dir_input, "lcb_summary_output")
  }
  if(is.logical(dir_output)){
    stop("output directory is logical")
  }
  ################################################################################


  #### READ IN XMFA FILE
  xmfa <- readLines(file.path(dir_input, file.name))
  xmfa <- xmfa[-grep("#", xmfa)] # remove header information

  # Extract rows where each LCB begins and ends
  xmfa <- c("=", xmfa) # add equals sign to first row to indicate start of LCB
  lcb.rows <- grep("=", xmfa)

  # Make a list where each element is a unique LCB
  # lcb.list <- vector("list", length(lcb.rows)-1) # initialize empty list
  #
  # for(i in 1:(length(lcb.rows)-1)){
  #   lcb.list[[i]] <- xmfa[(lcb.rows[i] + 1):(lcb.rows[i + 1] - 1)]
  # }
  # names(lcb.list) <- paste0("lcb_", seq(1:(length(lcb.rows)-1)))
  #

  ## Create list of LCBs
  lcb.list <- vector("list", length(lcb.rows)-1) # initialize empty list
  for(lcb in 1:(length(lcb.rows)-1)){
    lcb.seq <- xmfa[(lcb.rows[lcb] + 1):(lcb.rows[lcb + 1] - 1)]

    fasta.header <- c(grep(">", lcb.seq), (length(lcb.seq) + 1))

    fasta <- as.character(NULL)
    for(i in 1:(length(fasta.header)-1)){
      seq.rows <- seq(fasta.header[i], fasta.header[i+1])
      seq.rows <- seq.rows[-c(1, length(seq.rows))]

      seq <- as.character(NULL)
      for(j in 1:length(seq.rows)){
        seq <- paste0(seq, lcb.seq[(seq.rows[j])])
      }
      fasta <- rbind(fasta, lcb.seq[fasta.header[i]], seq)

    }
    lcb.list[[lcb]] <- as.character(fasta)
  }
  names(lcb.list) <- paste0("lcb_", seq(1:(length(lcb.rows)-1)))

  #### CALCULATE SUMMARY INFORMATION FOR LCBs

  ## Number of genomes that contain an LCB
  genome.core <- do.call("rbind", lapply(lcb.list, function(x) length(grep(">", x)))) %>%
    as.data.frame(.) %>%
    dplyr::rename("genome_count"=V1) %>%
    mutate(lcb= rownames(.),
           core= ifelse(genome_count == 18, "core", "not_core")) %>%
    select(lcb, genome_count, core)

  ## Nucleotide length and gaps of each lcb
  summarize_lcb <- function(){
    length.per.lcb <- do.call("rbind", lapply(lcb.list, function(x) mean(nchar(x[-grep("[>]", x)]))))
    gaps.mean.per.lcb <- do.call("rbind", lapply(lcb.list, function(x) mean(str_count(x[-grep("[>]", x)], "-")))) %>%
      round(., 1)
    gaps.sd.per.lcb <- do.call("rbind", lapply(lcb.list, function(x) sd(str_count(x[-grep("[>]", x)], "-")))) %>%
      round(., 1)
    gaps.min.per.lcb <- do.call("rbind", lapply(lcb.list, function(x) min(str_count(x[-grep("[>]", x)], "-")))) %>%
      round(., 1)
    gaps.max.per.lcb <- do.call("rbind", lapply(lcb.list, function(x) max(str_count(x[-grep("[>]", x)], "-")))) %>%
      round(., 1)

    # percent.gaps <- ((gaps.mean.per.lcb / length.mean.per.lcb)*100) %>%
    #                   round(., 1)
    lcb.summary <- data.frame(genome.core, length.per.lcb, gaps.mean.per.lcb, gaps.sd.per.lcb, gaps.min.per.lcb, gaps.max.per.lcb) %>%
      as_tibble()
  }
  lcb.summary <- summarize_lcb()

  write.table(lcb.summary,
              file= file.path(dir_output, paste0(file.name, "_lcb_summary.txt")),
              sep="\t",
              row.names= FALSE,
              col.names = TRUE,
              quote= FALSE)

  #### CALCULATE WHICH GENOMES ARE REPRESENTED IN EACH LCB
  ## Read in genomes in each LCB
  genome.df.raw <- lapply(names(lcb.list), function(x) str_replace(lcb.list[[x]][grep(">", lcb.list[[x]])],
                                                                   ">.*/",
                                                                   "")) %>%
    plyr::ldply(., rbind) %>%
    as_tibble() %>%
    mutate(lcb= names(lcb.list)) %>%
    dplyr::mutate_all(funs(as.character))

  genome.id <- genome.df.raw %>%
    as_tibble() %>%
    gather(key= "genome", value= "id", -lcb) %>%
    distinct(id) %>%
    filter(complete.cases(id)) %>%
    pull()

  # genome.id <- sort(unique(pull(genome.df.raw[, 1])))
  # genome.num <-  length(unique(pull(genome.df.raw[, 1])))
  genome.num <- length(genome.id)
  genome.count <- formatC(1:genome.num, width = nchar(genome.num), format = "d", flag = "0")
  colnames(genome.df.raw)[1:genome.num] <- paste0(rep("g", genome.num), genome.count)

  ## Loop to create matrix of presence and absence of each genome in each LCB
  genome.df.clean <- as.data.frame(matrix(rep(NA, nrow(genome.df.raw)*genome.num),
                                          nrow= nrow(genome.df.raw),
                                          ncol= genome.num))
  colnames(genome.df.clean) <- str_replace(genome.id, ".contigs.fa", "")
  rownames(genome.df.clean) <- names(lcb.list)

  for(row in 1:nrow(genome.df.raw)){
    genome.df.clean[row, ] <- as.numeric(genome.id %in% genome.df.raw[row, ])
  }

  genome.df.clean$lcb <- row.names(genome.df.clean)
  genome.df.clean <- select(genome.df.clean, lcb, everything())

  write.table(genome.df.clean,
              file= file.path(dir_output, paste0(file.name, "_genome_table.txt")),
              sep="\t",
              row.names= FALSE,
              col.names = TRUE,
              quote= FALSE)
  print(paste0("OUTPUT: lcb.summary.txt and genome_table.txt written to ", dir_output))

  # genome.df.clean.long <- genome.df.clean %>%
  #   mutate(lcb= row.names(genome.df.clean)) %>%
  #   as_tibble() %>%
  #   gather(key= "genome", value= "Y_N", -lcb) %>%
  #   mutate(Y_N= ifelse(Y_N == 1, "Y", "N")) %>%
  #   arrange(lcb)

  # write.table(genome.df.clean.long,
  #             file= file.path(dir_output, "genome_table_long.txt"),
  #             sep="\t",
  #             row.names= FALSE,
  #             col.names = TRUE,
  #             quote= FALSE)


  #### EXPORT FASTA FILES
  if(write_fasta == TRUE){
    ## Fasta file for each LCB
    dir.create(file.path(dir_output, "lcb_fasta_files"))
    lapply(names(lcb.list), function(x) write.table(lcb.list[[x]],
                                                    file= file.path(dir_output, "lcb_fasta_files", paste0(x, ".fasta")),
                                                    sep=" ",
                                                    row.names= FALSE,
                                                    col.names = FALSE,
                                                    quote= FALSE))
    print(paste0("OUTPUT: fasta files for each LCB written to ", file.path(dir_output, "lcb_fasta_files")))
  } else {
    print("No LCB fasta files exported")
  }
}

