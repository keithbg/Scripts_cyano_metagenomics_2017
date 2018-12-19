## Make a map of the sampling locations for the PH2017 Phormidium sampling
## Converts lat/long points to distances between samples (units = meters)


#### Libraries #################################################################
library(tidyverse)
################################################################################
?raster::select

#### FILE PATHS ################################################################
dir_input <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","EnvData")
dir_output <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data", "EnvData")
################################################################################


## Read in data
latlong <- read_csv(file.path(dir_input, "PhormMeta17_LatLong_combined_ggmap.csv")) %>%
            mutate(year= as.character(year))
#latlong17 <- read_csv(file.path(dir_input, "PhormMeta17_LatLong.csv"))

## Calculate Euclidean distance between sample sites
#euc.dist <- raster::pointDistance(latlong[, c(1,2)], latlong[, c(1,2)], lonlat= TRUE, allpairs= TRUE)
euc.dist <- raster::pointDistance(latlong[, c(1,2)], lonlat= TRUE)
rownames(euc.dist) <- latlong$ggkbase_id
colnames(euc.dist) <- latlong$ggkbase_id

euc.df <- as.data.frame(euc.dist) %>%
            mutate(sample.2= row.names(.)) %>%
            as_tibble() %>%
            gather(key= "sample.1", value= "euc_dist", -sample.2) %>%
            select(sample.1, sample.2, euc_dist)
write_tsv(euc.df, path= file.path(dir_output, "Euclidean_distance_samples_meters.tsv"))
