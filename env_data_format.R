## Script to input and format environmental covariate data from Phormidium metagenomics
## collections summer 2017.


#### Libraries #################################################################
library(tidyverse)
library(stringr)
library(lubridate)
################################################################################


#### FILE PATHS ################################################################
dir_input <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","EnvData")
dir_out_fig <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","Output_figures", "EnvData")
################################################################################



## READ AND FORMAT WATER CHEMISTRY DATA
water <- read_tsv(file.path(dir_input, "PhormMeta17_WaterChem.txt")) %>%
  mutate(date= dmy(date),
         NH4_ugL= NH4_uM*14.0067,
         NO3_ugL= NO3_uM*14.0067,
         DOC_ugL= (DOC_mM*12.0107)*1000,
         TDN_ugL= (TDN_mM*14.0067)*1000,
         NP= (TDN_mM*1000)/(TDP_ugL/30.9737)) %>%
  arrange(fork) %>%
  select(date, site, biotite_acronym, fork, contains("_ugL"), NP)


#### READ AND FORMAT ENVIRONMENTAL COVARIATE DATA
env.1 <- read_tsv(file.path(dir_input, "PhormMeta17_EnvData.txt")) %>%
        ## input NA when the vZ probe was exposed (exp) above the water
          mutate(vZ= as.numeric(ifelse(.$vZ == "exp", NA, .$vZ))) %>%
          group_by(date, site, us_ds) %>%
          mutate(mean_vX= mean(vX, na.rm= TRUE)) %>%
          filter(row_number() == 1) %>%
          ungroup() %>%
          select(-starts_with("v"), -flow_prop) %>%
          fill(-dist_cob)


#### CALCULATE AVERAGE CANOPY COVER
canopy.cover <- env.1 %>%
   select(date, site, us_ds, starts_with("can")) %>%
   gather(key= direction, value= dots, can_us:can_rl) %>%
   group_by(date, site, us_ds, can_light_shade) %>%
   dplyr::summarize(mean_dots= mean(dots, na.rm= TRUE)) %>%
   mutate(mean_dots_cover= ifelse(can_light_shade == "l", 100-mean_dots, mean_dots),
          canopy_cover_percent= mean_dots_cover*1.04) %>%
  select(date, site, us_ds, canopy_cover_percent)

#### MERGE CANOPY COVERAGE BACK INTO ENV
env.2 <- env.1 %>%
  select(-starts_with("can"), -notes) %>%
  left_join(., canopy.cover) %>%
  mutate(uniqueID= str_c(site, us_ds, sep= "_"),
         date= dmy(date))

#### MERGE WITH WATER CHEMISTRY DATA
env <- left_join(env.2, water) %>%
  select(date, site, us_ds, uniqueID, biotite_ID, fork, time, everything())


#### PREPARE ORDINATION MATRIX
# variables are normalized by mean and variance with scale() function
env.ord <- env %>%
             filter(uniqueID != "bear_us") %>%
             select(temp, pH, do_mgL, cond_ms, alk_conv, depth_cm, mean_vX, canopy_cover_percent, contains("ugL"), NP) %>%
           as.matrix(.) %>%
           scale(.)
row.names(env.ord) <- subset(env, uniqueID != "bear_us")$biotite_ID

#### REMOVE UNNECESSARY DATA FRAMES
rm(env.1, env.2, canopy.cover)
