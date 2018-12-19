## # Generalized Dissimilarity Matrices of PH2017 Phormidium samples

# Generalized Dissimilarity Matrices
# Fitzpatrick and Keller 2015, Ecology Lettes

## Euclidean distances calculated in Spatial_distances_format.R
## Distances in meters

## Average nucleotide identity (ANI) calculated with dRep

## ani.species.tsv generated in dRep_format.R


#### Libraries #################################################################
library(tidyverse)
# library(gdm) DON"T LOAD THIS PACKAGE AS IT CONFLICTS WITH TIDYVERSE
#library(vegan)
#library(maptools)
#library(ggplot2)
#library(RColorBrewer)
#source("/Users/KeithBG/R_functions/ggplot_themes.R")
################################################################################


#### FILE PATHS ################################################################
dir_input_scripts <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data", "Scripts")
dir_input_ani <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data", "GenomesData", "dRep")
dir_input_species <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data", "GenomesData")
dir_input_env <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data", "EnvData")
dir_output <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data", "dRep")
################################################################################


#### READ AND FORMAT DATA ######################################################

#### PH2015 ENVIRONMENTAL DATA ####

env.2015.df <- read_tsv(file.path(dir_input_env, "PH2015_env_data.tsv"))
env.2015.df <- env.2015.df %>%
  select(-Alk)

## Rename sample ID to match ISME publication
id_info <-  read_tsv("/Users/KeithBG/Documents/UC Berkeley/2015_SummerResearch/Metagenomics/new_sample_IDs2.txt", col_names = FALSE)
for(id in 1:length(env.2015.df$site)){
  env.2015.df$site[id] <- id_info$X2[id_info$X1 %in% env.2015.df$site[id]]
}
rm(id_info)

env.2015.df <- env.2015.df %>%
  rename(temp= temp_c, do_mgL= DO_mgL, cond_ms= cond, canopy_cover_percent= canopy_avg, depth_cm= depth_cont, DOC_ugL= NPOC_mgL) %>%
  select(-watershed_km2) %>%
  mutate(DOC_ugL= DOC_ugL*1000, # convert from mg/L to ug/L
         flow_cont= flow_cont/100, # convert from cm/s to m/s
         NP= (TDN_ugL/14.0067)/(TDP_ugL/30.9737),
         canopy_cover_percent= canopy_cover_percent*100,
         cond_ms= cond_ms/1000)

#### PH2017 ENVIRONMENTAL DATA ####
##(variables have already been normalized by mean and variance)
source(file.path(dir_input_scripts, "env_data_format.R"))
rm(water, env.ord)


env.2017.df <- env %>%
  as_tibble() %>%
  filter(uniqueID != "bear_us") %>%
  select(-site, -date, -uniqueID, -fork, -us_ds, -time, -mbars, -do_perL, -dist_cob,  -alk, -alk_conv, -flow_depth_dft, -biotite_acronym) %>%
  rename(site= biotite_ID, flow_cont= mean_vX)

#### MERGE 2015 and 2017 DATA TOGETHER ####
env.all.df <- full_join(env.2015.df, env.2017.df) %>%
                mutate_at(vars(temp:NP), funs(scale(.))) %>% # scale variables by mean and std. deviation
                select(-pH, -DOC_ugL, -depth_cm) # Remove variables with correlations > 0.5 (see below)

## Test for correlations between variables
   # pairs(env.all.df[, -1])
   # cor(env.all.df[, -1]) > 0.5
   # Removed: pH, DOC, -depth_cm


#### LAT/LONG DATA
lat.long.info <-  read_csv(file.path(dir_input_env, "PhormMeta17_LatLong_combined.csv")) %>%
  select(ggkbase_id, year, fork, lat, long) %>%
  rename(site= ggkbase_id)


## Add lat/long data to environmental data
env.all.df <- lat.long.info %>%
               select(site, lat, long) %>%
               right_join(., env.all.df)


#### ADD GENOMES TO ENVIRONMENTAL DATA ####

## List of ANI species 1 genomes
ani.sp1.genomes <- read_tsv(file.path(dir_input_species, "genome_ani_species.tsv")) %>%
  filter(ani_species == 1) %>%
  select(genome)

genome.site <- read_tsv(file.path(dir_input_env, "site_genome_table.tsv")) %>%
                mutate(site2= str_replace(.$site, "[A|B]$", "")) %>% # Adjust PH2017 site id to account for A and B samples with the same env. data
                mutate(site3= str_replace(.$site2, "PH2017_[0-9].", "")) %>%
                select(-site, -site2)

env.all.df <- env.all.df %>%
                mutate(site2= str_replace(.$site, "[A|B]$", "")) %>%
                mutate(site3= str_replace(.$site2, "PH2017_[0-9].", "")) %>%
                select(-site, -site2) %>%
                full_join(genome.site, .) %>%
                select(-site3)

# Filter by only ANI species 1
env.sp1.df <- env.all.df %>%
                 filter(env.all.df$genome %in% ani.sp1.genomes$genome) %>%
                 as.data.frame()

#### ANI DATA ####
ani.df <- read_csv(file.path(dir_input_ani, "output_20180927",  "Ndb.csv")) %>%
  select(querry, reference, ani) %>%
  spread(key= querry, value= ani) %>%
  select(-reference) %>%
  mutate(genome= str_replace(colnames(ani.mat), ".fa", "")) %>% # Clean up genome names
  mutate(genome= str_replace(.$genome, "_s25", ""))
colnames(ani.df) <- c(ani.df$genome, "genome")

ani.df[is.na(ani.df)] <- 0.5 # Add a low ANI value to all NAs

## Subset only species 1
ani.sp1 <- left_join(ani.sp1.genomes, ani.df)
ani.sp1 <- ani.sp1[, colnames(ani.sp1) %in% c(ani.sp1$genome, "genome")] # select only species 1 columns

## Calculate ANI distances (method Euclidean)
ani.sp1.dist <- as.data.frame(as.matrix(dist(ani.sp1[, -1], method= "euclidean")))
row.names(ani.sp1.dist) <- ani.sp1$genome
colnames(ani.sp1.dist) <- row.names(ani.sp1.dist)

## Add genome as first column
ani.sp1.dist.mat <- as.matrix(data.frame(genome= row.names(ani.sp1.dist), ani.sp1.dist))
ani.sp1.dist.DF <- data.frame(genome= row.names(ani.sp1.dist), ani.sp1.dist, stringsAsFactors = FALSE)


#### RUN GDM MODEL ############################################################


## Format data for model
gdm.format <- gdm::formatsitepair(bioData = ani.sp1.dist.DF,
                             bioFormat = 3,
                             predData=   env.sp1.df,
                             siteColumn = "genome",
                             XColumn = "long",
                             YColumn = "lat")

##  Run GDM
gdm.sp1.fit <- gdm::gdm(gdm.format, geo= TRUE)
summary(gdm.sp1.fit)
str(gdm.sp1.fit)

## Plots
  #plot(gdm.sp1.fit , plot.layout= c(3, 3), include.rug= TRUE, rug.sitepair= gdm.format)


#### OLD CODE ###################################################################
#################################################################################

#### FORMAT FOR GGPLOT
  ## ENVIRONMENTAL PREDICTORS

  ## Extract splines
  gdm.fit.splines <- isplineExtract(gdm.fit)
  str(gdm.fit.splines)

  x.values <- gdm.fit.splines[[1]] %>%
    as_tibble() %>%
    select(pH, DO_mgL, TDN_ugL, TDP_ugL, flow_cont) %>%
    gather(key= "predictor", value= "x", everything())

  y.values <- gdm.fit.splines[[2]] %>%
    as_tibble() %>%
    select(pH, DO_mgL, TDN_ugL, TDP_ugL, flow_cont) %>%
    gather(key= "predictor", value= "y", everything())

  env.pred.df <- data.frame(x.values, y= y.values$y)


  ## OVERALL MODEL
  mod.df <- data.frame(observed= gdm.fit[["observed"]], predicted= gdm.fit[["predicted"]], ecological= gdm.fit[["ecological"]])
  mod.df$eco.fit <- 1 - exp(-mod.df$ecological)
  head(mod.df)


#### GGPLOTS

dir_output_figs <- file.path("/Users", "KeithBG", "Documents", "UC Berkeley", "2015_SummerResearch", "Metagenomics", "Genome_Analyses", "PYANI", "pyani_results_20180907")

## ENVIRONMENTAL PREDICTORS
  legend.name <- "Environmental \npredictors"
  legend.order <- c("TDN_ugL", "DO_mgL", "pH", "TDP_ugL", "flow_cont")
  legend.labels <- c("TDN", "DO", "pH", "TDP", "Velocity")
  predictor.labels <- env.pred.df %>%
                        group_by(predictor) %>%
                        dplyr::summarise(ymax= max(y),
                                  xmax= max(x)) %>%
                        mutate(pred.label= c("DO", "water velocity", "pH", "TDN", "TDP"))
  predictor.colors <- brewer.pal(7, "Dark2")[c(2, 7, 3, 4, 1)]

ggplot(data= env.pred.df, aes(x= x, y= y)) +
  geom_hline(yintercept= 0, size= 0.25, color= "black") +
  geom_line(aes(linetype= predictor, color= predictor), size= 1.5) +
  geom_text(data= predictor.labels, aes(x= xmax, y= ymax, label= pred.label), nudge_x= 0.05, hjust= "left", size= 6) +
  labs(x= "Environmental predictor dissimilarity", y= "ANI dissimilarity") +
  #annotate("text", label= "TDN", x= TDN.label$xmax, y= TDN.label$ymax) +
  scale_color_manual(values= predictor.colors, breaks= legend.order, labels= legend.labels, name= legend.name) +
  scale_linetype_manual(values= c(2, 4, 3, 1, 5), breaks= legend.order, labels= legend.labels, name= legend.name) +
  scale_x_continuous(limits= c(-2.5, 3.6), breaks= seq(-2, 4, by= 1), expand= c(0, 0)) +
  scale_y_continuous(limits= c(0, 0.55), breaks= c(0, 0.25, 0.5)) +
  theme_kbg + theme(legend.position = "none")
ggsave(last_plot(), file= "gdm.env.predictor_20180910.eps", height= 6, width= 8, units= "in", device= cairo_pdf(), path= dir_output_figs)

## OVERALL PREDICTIONS

ggplot(data= mod.df, aes(x= predicted, y= observed)) +
  geom_point() +
  geom_line(aes(x= predicted, y= predicted), color= "red") +
  theme_kbg

ggplot(data= mod.df, aes(x= ecological, y= observed)) +
  geom_hline(yintercept= 0, size= 0.25, color= "black") +
  geom_point(color= "dimgray", size= 2) +
  geom_line(aes(x= ecological, y= eco.fit), linetype= 5, size= 0.75) +
  labs(x= "Predicted ecological distance", y= "Observed ANI distance") +
  scale_x_continuous(limits= c(min(mod.df$ecological), max(mod.df$ecological)), breaks= c(0.1, 0.3, 0.6, 0.9)) +
  theme_kbg
ggsave(last_plot(), file= "gdm.overall.model.eps", height= 6, width= 8, units= "in", device= cairo_pdf(), path= dir_output_figs)





plot(gdm.fit, plot.layout= c(3, 3), include.rug= TRUE, rug.sitepair= gdm.format)

predict(gdm.fit, gdm.format) == mod.df$predicted

##### ATTEMPT 2 JUNE 2018 ######################################################
## MODEL DID NOT WORK ON PRESENCE ABSENCE DATA
# dir_input_gdm <- file.path("/Users", "KeithBG", "Documents", "UC Berkeley", "2015_SummerResearch", "Metagenomics", "Genome_Analyses", "PYANI")
#
# ## PRESENCES ABSENCE OF ANI GROUPS TABLE
# bioTable <- read.table(file.path(dir_input_gdm, "gdm_bio_table.txt"), sep= "\t", header= TRUE, stringsAsFactors = FALSE)
# bioTable$site <- as.character(biotTable$site)
# bioTable_jaccard <- as.data.frame(as.matrix(vegdist(bioTable[, 2:5], method= "jaccard", diag= TRUE, upper= TRUE)))
# bioTable_jaccard$site <- as.character(bioTable$site)
#
# bioTable_jaccard <- bioTable_jaccard[ , c(23, 1:22)]
#
# ## ENVIRONMENTAL VARIABLES
# source("/Users/KeithBG/Documents/UC Berkeley/2015_SummerResearch/Metagenomics/Covariates/Covariates_Format.R")
# str(env.cov.ord.c) # centered variables
# envTable <-  subset(env.cov.ord.c, select= -c(Alk, watershed_km2))
# envTable$site <- bioTable$site
# envTable <- data.frame(envTable[, c(13, 1:12)], long= bioTable$long, lat= bioTable$lat)
#
#
# gdm.formatTable <- formatsitepair(bioData = bioTable,
#                              bioFormat = 1,
#                              predData=   envTable,
#                              siteColumn = "site",
#                              XColumn = "long",
#                              YColumn = "lat")
#
# gdm.formatTable.jaccard <- formatsitepair(bioData = bioTable_jaccard,
#                                   bioFormat = 3,
#                                   predData=   envTable,
#                                   siteColumn = "site",
#                                   XColumn = "long",
#                                   YColumn = "lat")
#
#
# ##  Run GDM
# gdm.fit2 <- gdm(gdm.formatTable, geo= FALSE)
# summary(gdm.fit2)
# str(gdm.fit2)
#
# gdm.fit3 <- gdm(gdm.formatTable.jaccard, geo= TRUE)
# summary(gdm.fit3)
#
# ## Plots
# plot(gdm.fit3, plot.layout= c(3, 3), include.rug= TRUE, rug.sitepair= gdm.format)
#
