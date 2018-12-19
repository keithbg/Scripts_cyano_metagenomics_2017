## Statistics of ANI and Euclidean distances from PH2017 Phormidium samples

## Euclidean distances calculated in Spatial_distances_format.R
## Distances in meters

## Average nucleotide identity (ANI) calculated with dRep

## ani.species.tsv generated in dRep_format.R


#### Libraries #################################################################
library(tidyverse)
################################################################################


#### FILE PATHS ################################################################
dir_input_dist <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","EnvData")
dir_input_ani <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "dRep")
dir_input_site <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","EnvData")
dir_output <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data")
################################################################################


## Read in data
euclidean.distance <- read_tsv(file.path(dir_input_dist, "Euclidean_distance_samples_meters.tsv"))
site.info <-  read_csv(file.path(dir_input_site, "PhormMeta17_LatLong_combined.csv")) %>%
  select(ggkbase_id, fork)
ani.species <- read_tsv(file.path(dir_input_ani, "output_20180927", "ani_species.tsv")) %>%
  mutate(genome= str_replace(genome, "_s25", ""))

ani.dist <- read_csv(file.path(dir_input_ani, "output_20180927", "Ndb.csv")) %>%
  mutate(querry= str_replace(querry, ".fa", ""),
         reference= str_replace(reference, ".fa", "")) %>%
  mutate(querry= str_replace(querry, "_s25", ""),
         reference= str_replace(reference, "_s25", ""),
         sample.1= str_replace(querry, "(_Oscillatoriales.*)", ""),
         sample.2= str_replace(reference, "(_Oscillatoriales.*)", "")) %>%
  mutate(year.querry= str_replace(querry, "_.*$", ""),
         year.reference= str_replace(reference, "_.*$", "")) %>%
  left_join(., euclidean.distance) %>% # combine euclidean distances with ani distances
  rename(sample.querry = sample.1, sample.reference = sample.2) %>%
  select(querry, reference, sample.querry, sample.reference, year.querry, year.reference, primary_cluster, ani, euc_dist) %>%
  filter(!is.na(euc_dist), querry != reference)


## Add ani species and site information
ani.dist <- ani.species %>%
  rename(querry= genome, species.querry= ani_species) %>%
  right_join(ani.dist)
ani.dist <- ani.species %>%
  rename(reference= genome, species.reference= ani_species) %>%
  right_join(ani.dist)
ani.dist <- site.info %>%
  rename(sample.querry= ggkbase_id, fork.querry= fork) %>%
  right_join(ani.dist)
ani.dist <- site.info %>%
  rename(sample.reference= ggkbase_id, fork.reference= fork) %>%
  right_join(ani.dist)
ani.dist <- ani.dist %>%
  mutate(species.comparison= str_c(species.reference, species.querry, sep= "-"),
         species.match= ifelse(species.querry == species.reference, "Y", "N"),
         year.match= ifelse(year.querry == year.reference, "Y", "N"),
         fork.match= ifelse(fork.querry == fork.reference, "Y", "N"),
         fork.comparison= str_c(fork.reference, fork.querry, sep= "-"))




#### GAM Model for ANI Species 1 ###############################################

## Run GAM model
library(visreg)
library(mgcv)
gam.fit <- gam(ani ~ s(euc_dist, bs= "cr"), data= subset(ani.dist,
                                                         species.match == "Y" &
                                                           species.comparison == "1-1"))
#summary(gam.fit)
#plot(gam.fit)
#gam.check(gam.fit) # check to make sure k value is not too low
#anova(gam.fit)
visualize.gam <- visreg(gam.fit, plot= FALSE)




## Investigating the large variance at distance <10km for species 1
test <- subset(ani.dist,
               species.match == "Y" & species.comparison == "1-1" & euc_dist < dist.limit & ani < 0.985)


dist.limit <- 15000

species1.subset <- subset(ani.dist,
                          species.match == "Y" & species.comparison == "1-1" & euc_dist < dist.limit)

ani.euc.plot3 <- ggplot(data= species1.subset,
                        aes(x= euc_dist, y= ani))

ani.euc.plot3 +
  geom_point(aes(color= fork.comparison, shape= fork.match), size= 3) +
  geom_line(aes(x= euc_dist, y= predict(species1.fit2)))
labs(x= "Euclidean distance (m)", y= "Average nucleotide identity (%)") +
  stat_smooth(method= "lm") +
  #scale_fill_manual(values= fill.year.match) +
  #scale_shape_manual(values= shape.year.match) +
  scale_x_continuous(limits= c(0, dist.limit)) +
  y_axis_format +
  theme_distance


gam.species1.fit1 <- gam(ani ~ s(euc_dist, bs= "cr"), data= species1.subset)
gam.species1.fit2 <- gam(ani ~ fork.match + s(euc_dist, bs= "cr"), data= species1.subset)
summary(gam.species1.fit1)
anova.gam(gam.species1.fit1, gam.species1.fit2, test= "F", freq= TRUE)

species1.fit1 <- lm(ani ~ euc_dist, data= species1.subset)
species1.fit2 <- lm(ani ~ fork.match*euc_dist, data= species1.subset)
predict(species1.fit2)
plot(species1.fit2)
anova(species1.fit1, species1.fit2)


species1.subset.all <- subset(ani.dist,
                          species.match == "Y" & species.comparison == "1-1")



species1.fit.red <- lm(ani ~ 1, data= species1.subset.all)

species1.fit1 <- lm(ani ~ euc_dist, data= species1.subset.all)
species1.fit2 <- lm(ani ~ fork.comparison + euc_dist, data= species1.subset.all)
species1.fit3 <- lm(ani ~ fork.comparison + year.match + euc_dist, data= species1.subset.all)
species1.fit4 <- lm(ani ~ fork.match+euc_dist, data= species1.subset.all)
species1.fit5 <- lm(ani ~ fork.match*euc_dist, data= species1.subset.all)
summary(species1.fit4)
anova(species1.fit4)

anova(species1.fit.red, species1.fit1, species1.fit2, species1.fit3)

length(unique(species1.subset.all$querry))
(46*46)/2-46
