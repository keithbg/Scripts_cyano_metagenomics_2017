## Combine ANI and Euclidean distances from PH2017 Phormidium samples

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
dir_output <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data", "GenomesData", "dRep", "figures_R")
################################################################################


## Read in data
euclidean.distance <- read_tsv(file.path(dir_input_dist, "Euclidean_distance_samples_meters.tsv"))
site.info <-  read_csv(file.path(dir_input_site, "PhormMeta17_LatLong_combined.csv")) %>%
               select(ggkbase_id, fork) # fork indicates which part of the watershed the sample was collected
                                        # cr= creek, sf= south fork, etc.
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


##### PLOTTING PARAMETERS ######################################################

add.gam.mean <- geom_line(data= visualize.gam$fit, aes(x= euc_dist, y=visregFit), size= 0.5)
add.gam.ribbon <- geom_ribbon(data= visualize.gam$fit, aes(x= euc_dist, ymin= visregLwr, ymax= visregUpr), alpha= 0.25, fill= "dodgerblue")


x_axis_format <- scale_x_continuous(breaks= seq(0, 175, by= 25),
                                    labels= c("0", "", "50", "", "100", "", "150", ""),
                                    expand= c(0, 5))

y_axis_format <- scale_y_continuous(limits= c(0.975, 1),
                                    breaks= seq(0.975, 1.000, by= 0.005),
                                    labels= c("", "98", "", "99", "", "100"))
fill.year.match <- c("white", "gray75")
species.facet.labels <- as_labeller(c(`1-1` = "Species 1", `2-2` = "Species 2", `3-3` = "Species 3"))
facet.by.species <- facet_grid(species.comparison~., labeller= labeller(species.comparison= species.facet.labels))


## ggplot theme for response ratio plots
theme_distance <- theme(panel.grid = element_blank(),
                  plot.margin = unit(c(1, 1, 1, 1), "cm"),
                  text = element_text(size= 14),
                  plot.background = element_rect(fill = "transparent"), # bg of the plot
                  panel.background = element_rect(fill= "transparent", color="black"),
                  axis.text = element_text(colour="black"),
                  axis.title.x = element_text(vjust = -0.75),
                  axis.title.y = element_text(vjust = 1.5),
                  legend.background = element_rect(size=0.25, color="black", fill= "transparent"),
                  legend.key = element_blank(),
                  strip.background=element_rect(fill="transparent", color="transparent"),
                  legend.position = "top")

#### PLOT THE DATA #############################################################

ani.euc.plot1 <- ggplot(data= subset(ani.dist,
                                     species.match == "Y" & species.comparison != "4-4"),
                        aes(x= euc_dist/1000, y= ani))

ani.euc.plot1 +
  geom_point(aes(fill= year.match), color= "black", size= 3, shape= 21) +
  labs(x= "Euclidean distance (km)", y= "Average nucleotide identity (%)") +
  scale_fill_manual(values= fill.year.match) +
  scale_shape_manual(values= shape.year.match) +
  x_axis_format +
  y_axis_format +
  facet.by.species +
  theme_distance


## ANI SPECIES 1
ani.euc.plot2 <- ggplot(data= subset(ani.dist,
                                     species.match == "Y" & species.comparison == "1-1"),
                        aes(x= euc_dist, y= ani))

ani.euc.plot2 +
  geom_point(aes(fill= fork.match), color= "black", size= 3, shape= 21) +
  add.gam.ribbon +
  add.gam.mean +
  labs(x= "Euclidean distance (km)", y= "Average nucleotide identity (%)", title= "PH2017 ANI species 1") +
  scale_fill_manual(values= fill.year.match) +
  scale_shape_manual(values= shape.year.match) +
  scale_x_continuous(limits= c(0, 150000),
                     breaks= seq(0, 150000, by= 25000),
                     labels= c("0", "", "50", "", "100", "", "150"),
                     expand= c(0, 2000)) +
  y_axis_format +
  theme_distance
ggsave(last_plot(), filename= "ANI_sp1_euc_km.pdf", width= 8, height= 6.4, units= "in", path= dir_output, device= cairo_pdf)



## Investigating the large variance at distance <10km for species 1
dist.limit <- 15000 # in meters

species1.subset <- subset(ani.dist,
                          species.match == "Y" & species.comparison == "1-1" & euc_dist < dist.limit)

mytest2 <- species1.subset %>%
  filter(ani < 0.985, euc_dist < 5000, fork.match == "N")

ani.euc.plot3 <- ggplot(data= species1.subset,
                        aes(x= euc_dist, y= ani))

ani.euc.plot3 +
  geom_point(aes(color= fork.comparison, shape= fork.match), size= 3) +
  #geom_line(aes(x= euc_dist, y= predict(species1.fit2)))
  labs(x= "Euclidean distance (m)", y= "Average nucleotide identity (%)") +
  stat_smooth(method= "lm") +
  #scale_fill_manual(values= fill.year.match) +
  #scale_shape_manual(values= shape.year.match) +
  scale_x_continuous(limits= c(0, dist.limit)) +
  y_axis_format +
  theme_distance
ggsave(last_plot(), filename= "ANI_sp1_fork_match.pdf", width= 8, height= 6.4, units= "in", path= dir_output, device= cairo_pdf)


