#### Install missing packages
list.of.packages <- c("reshape", "seqinr", "optparse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cloud.r-project.org")

suppressWarnings(library(reshape))
suppressWarnings(library(seqinr))
suppressWarnings(library(optparse))

#### Get I/O Options ###

PWD <- system("pwd", intern = TRUE)

option_list = list(
  make_option(c("-n", "--names"), type="character", default=NULL,
              help="File containing rp16 sequence names and bin correspondence", metavar="character"),

  make_option(c("-s", "--sequences"), type="character", default = NULL,
              help="File containing all concatenated protein sequences from bins", metavar="character"),

  make_option(c("-m", "--minrp"), type="integer", default = 8,
              help="Min number of rp16 proteins for inclusion in set [default = 8]", metavar="numeric"),

  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$names)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (names file).n", call.=FALSE)
}

 if (is.null(opt$sequences)){
   print_help(opt_parser)
   stop("At least one argument must be supplied (sequences file).n", call.=FALSE)
 }



### Load gene name table
print("Loading Files...", quote = FALSE)
names <- read.table(opt$names, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
names[names == "" | names == "-"] <- NA


### Remove Genomes with < --minrp option
names <- names[(16 - rowSums(is.na(names[,2:17]))) >= opt$minrp,]


names_melt <- melt(names, id.vars = "Bin")
names_melt <- names_melt[complete.cases(names_melt),]

print("Parsing rp16 sequences...", quote = FALSE)
system(paste0("pullseq -i ", opt$sequences," -N > all_rp16.faa"), input = as.vector(names_melt$value, mode = "character"))

print("Cleaning Parsed Sequences...", quote = FALSE)
system('sed "s/#.*//" all_rp16.faa > all_rp16_clean.faa')
system('rm all_rp16.faa')


sequences <- read.fasta(file = paste0(PWD,"/all_rp16_clean.faa"),
                        seqtype = "AA",as.string = TRUE, set.attributes = FALSE)


for (i in as.vector(unique(names_melt$variable), mode = "character")){
  tmp_df <- subset(names_melt, variable == i)
  tmp_seq <- sequences[names(sequences) %in% tmp_df$value]
  names(tmp_seq) <- tmp_df[match(names(tmp_seq), tmp_df$value), 1]
  write.fasta(tmp_seq, names = names(tmp_seq), file.out = paste0(PWD,"/",i,".faa"))
}






