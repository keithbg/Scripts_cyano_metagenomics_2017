## SNV linkage format script
## output from Alex CC's ~/programs/strains_analysis/strainRep2.py program


#### Libraries #################################################################
library(tidyverse)
library(ggplot2)
################################################################################

#### SOURCE FUNCTIONS ##########################################################
source("/Users/KeithBG/Documents/UC Berkeley/CyanoMeta_NSF/Metagenomics/Data/Scripts/snv_linkage_functions.R")

#### FILE PATHS ################################################################
dir_input <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "snv_linkage")
dir_input_sp1 <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "snv_linkage", "species_1")
dir_input_sp2 <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "snv_linkage", "species_2")
dir_input_sp3 <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "snv_linkage", "species_3")
dir_input_nostocales <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "snv_linkage", "nostocales")

dir_output_fig <- file.path("/Users","KeithBG","Documents","UC\ Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "snv_linkage", "Figures")
dir_output_table <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "snv_linkage", "Output_tables")
################################################################################


##### PLOTTING PARAMETERS ######################################################
x_axis_format_distance <- scale_x_continuous(expand= c(0, 5))
x_axis_format_window <- scale_x_discrete(labels= NULL, expand= c(0.02, 0))
y_axis_format <- scale_y_continuous(limits= c(0, 1),
                                    breaks= seq(0, 1, by= 0.25),
                                    expand= c(0.02, 0))



pid.facet.labels <- as_labeller(c(`98` = "PID= 0.98", `96` = "PID= 0.96" ))
facet.by.pid <- facet_grid(pid~., labeller= labeller(pid= pid.facet.labels, scales= "free_y"))


## ggplot theme for snv linkage
theme_snv <- theme(panel.grid = element_blank(),
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
################################################################################

## SAMPLE METADATA

samp.md <- read_tsv(file.path(dir_input, "sample.metadata.snv.linkage.tsv"))

## LOG FILES DATA
log_sp1 <- summarize_log_files(path= file.path(dir_input_sp1, "log_files"))
log_sp2 <- summarize_log_files(path= file.path(dir_input_sp2, "log_files"))
log_sp3 <- summarize_log_files(path= file.path(dir_input_sp3, "log_files"))
log_nostocales <- summarize_log_files(path= file.path(dir_input_nostocales, "log_files"))


log_master <- do.call(rbind, list(log_sp1, log_sp2, log_sp3)) %>%
               left_join(., samp.md, by= "sample")
temp_species <- str_extract_all(log_master$species_present, "[0-9]", simplify = TRUE)
log_master <- log_master %>%
  mutate(species_match= ifelse(str_replace_all(log_master$species, "[^0-9]", "") == temp_species[, 1] | str_replace_all(log_master$species, "[^0-9]", "") == temp_species[, 2], "Y", "N"))


## LINKAGE AND FREQ FILES
breadth.minimum <- 0.6
coef.table.sp1.pid99 <- analyse_linkage_data(path= dir_input_sp1, pid= 99, species= 1, output= dir_output_fig, min.breadth= breadth.minimum)
coef.table.sp1.pid98 <- analyse_linkage_data(path= dir_input_sp1, pid= 98, species= 1, output= dir_output_fig, min.breadth= breadth.minimum)

coef.table.sp2.pid99 <- analyse_linkage_data(path= dir_input_sp2, pid= 99, species= 2, output= dir_output_fig, min.breadth= breadth.minimum)
coef.table.sp2.pid98 <- analyse_linkage_data(path= dir_input_sp2, pid= 98, species= 2, output= dir_output_fig, min.breadth= breadth.minimum)

test2 <- analyse_linkage_data2(path= dir_input_sp2, pid= 98, ref_id= "species_2", output= dir_output_fig, min.breadth= breadth.minimum, linkage_plot= TRUE, freq_plot= TRUE)


coef.table.sp3.pid99 <- analyse_linkage_data(path= dir_input_sp3, pid= 99, species= 3, output= dir_output_fig, min.breadth= breadth.minimum)
coef.table.sp3.pid98 <- analyse_linkage_data(path= dir_input_sp3, pid= 98, species= 3, output= dir_output_fig, min.breadth= breadth.minimum)

coef.table.master <- do.call(rbind, list(coef.table.sp1.pid99, coef.table.sp1.pid98, coef.table.sp2.pid99, coef.table.sp2.pid98, coef.table.sp3.pid99, coef.table.sp3.pid98)) %>%
                       mutate(df= as.numeric(df)) %>%
                       left_join(log_master, by = c("sample", "species", "pid"))


## PLOT SIGNIFICANT DISTANCES SLOPES
sig.slopes <- coef.table.master %>%
  filter(p_value < 0.05)

coef.table.master %>%
  count(p_value < 0.05)

coef.table.master %>%
  group_by(p_value < 0.05, slope < 0) %>%
  count()

sig.slope.dist.plots <- with(sig.slopes, str_c(sample, species, str_c("PID=0.", pid, "_dist.jpg"), sep= "-"))
make_multi_panel_fig(path= dir_output_fig, file.list= sig.slope.dist.plots, output.name = "multi_fig_sig_slopes")



## MULTI-PANEL FIGURES OF THE LINKAGE AND FREQ DATA

## Species 1
make_multi_panel_fig(path= dir_output_fig, file.pattern= "species_1-PID=0.99_freq", output.name= "multi_panel-species_1-pid99-freq")
make_multi_panel_fig(path= dir_output_fig, file.pattern= "species_1-PID=0.99_dist", output.name= "multi_panel-species_1-pid99-dist")
make_multi_panel_fig(path= dir_output_fig, file.pattern= "species_1-PID=0.98_freq", output.name= "multi_panel-species_1-pid98-freq")
make_multi_panel_fig(path= dir_output_fig, file.pattern= "species_1-PID=0.98_dist", output.name= "multi_panel-species_1-pid98-dist")

## Species 2
make_multi_panel_fig(path= dir_output_fig, file.pattern= "species_2-PID=0.99_freq", output.name= "multi_panel-species_2-pid99-freq")
make_multi_panel_fig(path= dir_output_fig, file.pattern= "species_2-PID=0.99_dist", output.name= "multi_panel-species_2-pid99-dist")
make_multi_panel_fig(path= dir_output_fig, file.pattern= "species_2-PID=0.98_freq", output.name= "multi_panel-species_2-pid98-freq")
make_multi_panel_fig(path= dir_output_fig, file.pattern= "species_2-PID=0.98_dist", output.name= "multi_panel-species_2-pid98-dist")

## Species 3
make_multi_panel_fig(path= dir_output_fig, file.pattern= "species_3-PID=0.99_freq", output.name= "multi_panel-species_3-pid99-freq")
make_multi_panel_fig(path= dir_output_fig, file.pattern= "species_3-PID=0.99_dist", output.name= "multi_panel-species_3-pid99-dist")
make_multi_panel_fig(path= dir_output_fig, file.pattern= "species_3-PID=0.98_freq", output.name= "multi_panel-species_3-pid98-freq")
make_multi_panel_fig(path= dir_output_fig, file.pattern= "species_3-PID=0.98_dist", output.name= "multi_panel-species_3-pid98-dist")




## FILTER SNVS BY FREQUENCY
snv.minimum <- 100
#sp1.filter.99 <- filter_snv_freq(path= dir_input_sp1, output= dir_output_fig, species= 1, pid= 99, filter.values = c(0.2, 0.1, 0.05), direction= "<", min.snv.sites= snv.minimum)
#sp1.filter.98 <- filter_snv_freq(path= dir_input_sp1, output= dir_output_fig, species= 1, pid= 98, filter.values = c(0.2, 0.1, 0.05), direction= "<", min.snv.sites= snv.minimum)
sp1.filter.98.low <- filter_snv_freq(path= dir_input_sp1, output= dir_output_fig, species= 1, pid= 98, filter.values = c(0.2, 0.1), direction= "<", min.snv.sites= snv.minimum)
sp1.filter.98.high <- filter_snv_freq(path= dir_input_sp1, output= dir_output_fig, species= 1, pid= 98, filter.values = c(0.2, 0.3, 0.4), direction= ">", min.snv.sites= snv.minimum)

sp1.filter.98.window <- filter_snv_freq_window(path= dir_input_sp1, output= dir_output_fig, species= 1, pid= 98, min.freq = c(0.05, 0.1, 0.2, 0.3, 0.4), max.freq = c(0.1, 0.2, 0.3, 0.4, 0.5), min.snv.sites= snv.minimum)


#sp2.filter.99 <- filter_snv_freq(path= dir_input_sp2, output= dir_output_fig, species= 2, pid= 99, filter.values = c(0.2, 0.1, 0.05), direction= "<", min.snv.sites= snv.minimum)
sp2.filter.98.low <- filter_snv_freq(path= dir_input_sp2, output= dir_output_fig, species= 2, pid= 98, filter.values = c(0.2, 0.1), direction= "<", min.snv.sites= snv.minimum)
sp2.filter.98.high <- filter_snv_freq(path= dir_input_sp2, output= dir_output_fig, species= 2, pid= 98, filter.values = c(0.2, 0.3, 0.4), direction= ">", min.snv.sites= snv.minimum)

sp2.filter.98.window <- filter_snv_freq_window(path= dir_input_sp2, output= dir_output_fig, species= 2, pid= 98, min.freq = c(0.05, 0.1, 0.2, 0.3, 0.4), max.freq = c(0.1, 0.2, 0.3, 0.4, 0.5), min.snv.sites= snv.minimum)


test2 <- filter_snv_freq_window(path= dir_input_sp2, output= dir_output_fig, ref_id= "species_2", pid= 98, min.freq = c(0.1, 0.2), max.freq = c(0.2, 0.3), min.snv.sites= snv.minimum)


#sp3.filter.99 <- filter_snv_freq(path= dir_input_sp3, output= dir_output_fig, species= 3, pid= 99, filter.values = c(0.2, 0.1, 0.05), direction= "<", min.snv.sites= snv.minimum)
#sp3.filter.98 <- filter_snv_freq(path= dir_input_sp3, output= dir_output_fig, species= 3, pid= 98, filter.values = c(0.2, 0.1, 0.05), direction= "<", min.snv.sites= snv.minimum)
sp3.filter.98.low <- filter_snv_freq(path= dir_input_sp3, output= dir_output_fig, species= 3, pid= 98, filter.values = c(0.2, 0.1), direction= "<", min.snv.sites= snv.minimum)
sp3.filter.98.high <- filter_snv_freq(path= dir_input_sp3, output= dir_output_fig, species= 3, pid= 98, filter.values = c(0.2, 0.3, 0.4), direction= ">", min.snv.sites= snv.minimum)

sp3.filter.98.window <- filter_snv_freq_window(path= dir_input_sp3, output= dir_output_fig, species= 3, pid= 98, min.freq = c(0.05, 0.1, 0.2, 0.3, 0.4), max.freq = c(0.1, 0.2, 0.3, 0.4, 0.5), min.snv.sites= snv.minimum)




## Combined filtered results into a master data frame
freq.table.files <- list.files(file.path(dir_output_fig), pattern= "filter_summary-pid98")

freq.table.list <- list(NULL)
for(i in 1:length(freq.table.files)){
  freq.table.list[i] <- list(read_tsv(file.path(dir_output_fig, freq.table.files[i])))
}

freq.filter.master <- do.call(rbind, freq.table.list) %>%
  left_join(., samp.md, by= "sample") %>%
  mutate(pid= as.character(pid))

freq.unfiltered.master <- subset(freq.filter.master, direction == "unfiltered") %>%
                            left_join(log_master, by = c("sample", "species", "pid"))

freq.low.master <- subset(freq.filter.master, direction == "<" | direction == "unfiltered") %>%
  left_join(log_master, by = c("sample", "species", "pid"))


freq.high.master <- subset(freq.filter.master, direction == ">" | direction == "unfiltered") %>%
  left_join(log_master, by = c("sample", "species", "pid")) %>%
  filter(complete.cases(.))


#write_tsv(freq.unfiltered.master, path= file.path(dir_output_table, "freq_unfiltered_master.tsv") )

#### MAKE PLOTS ################################################################


#### LOG DATA #####

## BREADTH
ggplot(data= log_master, aes(x= species, y= breadth)) +
  geom_point(aes(color= species_match), position = "jitter", size= 3) +
  facet_grid(species_match~pid, scales= "free_y") +
  theme_snv


## NOSTOCALES
ggplot(data= log_master, aes(x= species, y= breadth)) +
  geom_point(aes(color= species_match), position = "jitter", size= 3) +
  facet_grid(species_match~pid, scales= "free_y") +
  theme_snv



#### LINKAGE DATA #####

ph2015_10s_99 <- read_tsv(file.path(dir_input_sp1, "PH2015_10S_0.99.linkage")) %>%
  mutate(pid= "99") %>%
  filter(complete.cases(.)) %>%
  rename(row_num= X1) %>%
  mutate(r2= ifelse(r2 > 1, 1, r2), # some r2 values seem to be just above 1
         dist.bin= cut(Distance, breaks= 10),
         scaffold= str_replace(Window, ":.*$", ""))
         #freq= colSums(c(count_AB:count_ab)))

ggplot(data= ph2015_10s_99, aes(x= Distance, y= 1/r2)) +
  geom_point() +
  theme_snv


ggplot(data= coef.table.master, aes(x= species, y= slope)) +
  geom_boxplot(aes(fill= species), alpha= 0.3) +
  geom_point(position = "jitter", size= 2) +
  facet_grid(pid~p_value<0.05) +
  theme_snv
ggsave(last_plot(), filename= "linkage_slope1.jpg", width= 8, height= 6, units= "in", dpi= 320, path= file.path(dir_output_fig, "Figs_summary"))


ggplot(data= coef.table.master, aes(x= species, y= exp(intercept))) +
  geom_boxplot(aes(fill= species), alpha= 0.3) +
  geom_point(position = "jitter", size= 2) +
  facet_grid(pid~p_value<0.05) +
  theme_snv
ggsave(last_plot(), filename= "linkage_intercept1.jpg", width= 8, height= 6, units= "in", dpi= 320, path= file.path(dir_output_fig, "Figs_summary"))


ggplot(data= coef.table.master, aes(x= exp(intercept), y= slope)) +
  geom_point(aes(color= species)) +
  facet_grid(pid~p_value<0.05, scales= "free_y") +
  theme_snv
ggsave(last_plot(), filename= "linkage_intercept_slope.jpg", width= 8, height= 6, units= "in", dpi= 320, path= file.path(dir_output_fig, "Figs_summary"))



#### SNV FREQUENCY DATA ####
#freq.unfiltered.master <- read_tsv(file.path(dir_output_table, "freq_unfiltered_master.tsv"))

## ALL ALLELES


## LOW FREQUENCY ALLELES

# Amount of significant slopes
freq.low.master %>%
  group_by(species, p_value < 0.05) %>%
  count()


ggplot(data= freq.low.master, aes(x= species, y= slope)) +
  geom_boxplot(aes(fill= filter_threshold), alpha= 0.3) +
  #geom_point(aes(), position = "jitter", size= 2) +
  facet_grid(p_value<0.05~.) +
  theme_snv
ggsave(last_plot(), filename= "freq_low_slope1.jpg", width= 8, height= 6, units= "in", dpi= 320, path= file.path(dir_output_fig, "Figs_summary"))

ggplot(data= freq.low.master, aes(x= filter_threshold, y= slope)) +
  geom_boxplot() +
  #geom_point() +
  #geom_line() +
  facet_grid(species~., scales= "free_y") +
  theme_snv
ggsave(last_plot(), filename= "freq_low_slope2.jpg", width= 8, height= 6, units= "in", dpi= 320, path= file.path(dir_output_fig, "Figs_summary"))

ggplot(data= freq.low.master, aes(x= filter_threshold, y= slope, group= sample)) +
  geom_point() +
  geom_line() +
  facet_grid(species~., scales= "free_y") +
  theme_snv
ggsave(last_plot(), filename= "freq_low_slope3.jpg", width= 8, height= 6, units= "in", dpi= 320, path= file.path(dir_output_fig, "Figs_summary"))


ggplot(data= freq.low.master, aes(x= filter_threshold, y= prop_r2_equals_1)) +
  geom_boxplot() +
  #geom_point() +
  #geom_line() +
  facet_grid(species~., scales= "free_y") +
  theme_snv




## HIGH FREQUENCY ALLELES

# Amount of significant slopes
freq.high.master %>%
  group_by(species, p_value < 0.05) %>%
  count()

ggplot(data= freq.high.master, aes(x= species, y= slope)) +
  geom_boxplot(aes(fill= filter_threshold), alpha= 0.3) +
  #geom_point(aes(), position = "jitter", size= 2) +
  facet_grid(p_value<0.05~.) +
  theme_snv
ggsave(last_plot(), filename= "freq_high_slope1.jpg", width= 8, height= 6, units= "in", dpi= 320, path= file.path(dir_output_fig, "Figs_summary"))

ggplot(data= freq.high.master, aes(x= filter_threshold, y= slope)) +
  geom_boxplot() +
  #geom_point() +
  #geom_line() +
  facet_grid(species~., scales= "free_y") +
  theme_snv
ggsave(last_plot(), filename= "freq_high_slope2.jpg", width= 8, height= 6, units= "in", dpi= 320, path= file.path(dir_output_fig, "Figs_summary"))

ggplot(data= freq.high.master, aes(x= filter_threshold, y= slope, group= sample)) +
  geom_point() +
  geom_line() +
  facet_grid(species~., scales= "free_y") +
  theme_snv
ggsave(last_plot(), filename= "freq_high_slope3.jpg", width= 8, height= 6, units= "in", dpi= 320, path= file.path(dir_output_fig, "Figs_summary"))

ggplot(data= freq.high.master, aes(x= filter_threshold, y= prop_r2_equals_1)) +
  geom_boxplot() +
  #geom_point() +
  #geom_line() +
  facet_grid(species~., scales= "free_y") +
  theme_snv



## CLONALITY
ggplot(data= subset(freq.unfiltered.master, p_value<0.05), aes(x= clonality, y= slope)) +
  #geom_boxplot(aes(fill= as.character(pid)), alpha= 0.3) +
  geom_point(aes(color= species), size= 3) +
  scale_x_continuous(limits= c(0.9965, 1), expand= c(0.01, 0)) +
  theme_snv
ggsave(last_plot(), filename= "clonality_slope.jpg", width= 8, height= 6, units= "in", dpi= 320, path= file.path(dir_output_fig, "Figs_summary"))

ggplot(data= log_master, aes(x= species, y= clonality)) +
  geom_boxplot(aes(fill= species), alpha= 0.3) +
  geom_point(position = "jitter", size= 2) +
  facet_grid(.~pid) +
  theme_snv
ggsave(last_plot(), filename= "clonality_species.jpg", width= 8, height= 6, units= "in", dpi= 320, path= file.path(dir_output_fig, "Figs_summary"))







  ggplot(data= freq.unfiltered.master, aes(x= species, y= slope)) +
  geom_boxplot(aes(fill= pid), alpha= 0.3) +
  #geom_point(aes(color= as.character(pid)), position = "jitter", size= 3) +
  #facet_grid(.~pid) +
  theme_snv



ggplot(data= freq.unfiltered.master, aes(x= species, y= intercept)) +
  geom_boxplot(aes(fill= as.character(pid)), alpha= 0.3) +
  geom_point(aes(color= as.character(pid)), position = "jitter", size= 3) +
  facet_grid(.~pid) +
  theme_snv

## Samples
ggplot(data= freq.unfiltered.master, aes(x= df, y= slope)) +
  geom_point(aes(color= species_present), size= 3) +
  scale_x_log10() +
  facet_grid(species~pid) +
  theme_snv

ggplot(data= freq.unfiltered.master, aes(x= df, y= prop_r2_equals_1)) +
  geom_point(aes(color= species_present), size= 3) +
  scale_x_log10() +
  facet_grid(species~pid) +
  theme_snv




ggplot(data= freq.unfiltered.master, aes(x= species, y= slope)) +
  geom_boxplot(aes(fill= species), alpha= 0.3) +
  geom_point(position = "jitter", size= 3) +
  facet_grid(.~pid) +
  theme_snv

ggplot(data= freq.unfiltered.master, aes(x= species, y= prop_r2_equals_1)) +
  geom_boxplot(aes(fill= species), alpha= 0.3) +
  geom_point(position = "jitter", size= 3) +
  facet_grid(.~pid) +
  theme_snv



ggplot(data= sp3.filter, aes(x= filter_threshold, y= prop_r2_equals_1, group= sample)) +
  geom_point(aes(color= sample), size= 3) +
  geom_line(aes(color= sample), size= 1) +
  theme_snv

ggplot(data= sp3.filter, aes(x= filter_threshold, y= slope, group= sample)) +
  geom_point(aes(color= sample), size= 3) +
  geom_line(aes(color= sample), size= 1) +
  theme_snv

