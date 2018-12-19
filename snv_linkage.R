## SNV linkage format script
## output from Alex CC's ~/programs/strains_analysis/strainRep2.py program


#### Libraries #################################################################
  library(tidyverse)
  library(progress)
  library(ggplot2)
  library(RColorBrewer)
  library(multipanelfigure)
################################################################################


#### FILE PATHS ################################################################
  dir_input_sp1 <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "snv_linkage", "species_1")
  dir_input_sp2 <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "snv_linkage", "species_1")
  dir_output <- file.path("/Users","KeithBG","Documents","UC\ Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "snv_linkage", "Figures")
  dir_output_table <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "snv_linkage", "Output_tables")
################################################################################

analyse_linkage_data <- function(freq_threshold= 0.05, species, path, linkage_plot= TRUE, freq_plot= TRUE){

## DEFINE VARIABLES
  dir_input <- path
  species.num <- paste0("species_", species)
  linkage.files <- list.files(file.path(dir_input), pattern= "linkage")
  linkage.sample <- str_replace(linkage.files, "_0\\..*$", "")
  freq.files <- list.files(file.path(dir_input), pattern= "freq")
  freq.sample <- str_replace(freq.files, "_0\\..*$", "")

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


  ##### LOOP THROUGH LINKAGE FILES AND MAKE PLOTS #########################################

  if(linkage_plot == TRUE){

    ## progress bar
    pb.linkage <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = length(linkage.files))
    pb.linkage$tick(0)

    rm(counter)
    mod.coef <- data.frame(sample= as.character(NULL), species= as.character(NULL), pid= as.character(NULL), df= as.character(NULL), intercept= as.numeric(NULL), slope= as.numeric(NULL), r2_adj= as.numeric(NULL), stringsAsFactors= FALSE)


    message("Looping over .linkage files")
    for(file in 1:length(linkage.files)){
      pb.linkage$tick() ## Progress Bar


      #### READ IN THE DATA
      snvs.98.coverage <- try(suppressMessages(read_tsv(file.path(dir_input, "log_files", str_c(linkage.sample[file], ".", ref_id, ".pid98", ".log")),
                                                        col_names = FALSE,
                                                        progress= FALSE)) %>%
                                pull() %>%
                                .[str_detect(., "Mean coverage")] %>%
                                str_replace("Mean coverage: ", "") %>%
                                as.numeric(.) %>%
                                round(., 1))

      snvs.96.coverage <- try(suppressMessages(read_tsv(file.path(dir_input, "log_files", str_c(linkage.sample[file], ".", species.num, ".pid96", ".log")),
                                                        col_names = FALSE,
                                                        progress= FALSE)) %>%
                                pull() %>%
                                .[str_detect(., "Mean coverage")] %>%
                                str_replace("Mean coverage: ", "") %>%
                                as.numeric(.) %>%
                                round(., 1))

      if(all(class(snvs.98.coverage) != "try-error", class(snvs.96.coverage) != "try-error")){
        snvs.98 <- suppressMessages(read_tsv(file.path(dir_input, str_c(linkage.sample[file], "_0.", 98, ".linkage")), progress= FALSE)) %>%
          mutate(pid= "98") %>%
          filter(total > snvs.98.coverage*freq_threshold, complete.cases(.))

        snvs.96 <- suppressMessages(read_tsv(file.path(dir_input, str_c(linkage.sample[file], "_0.", 96, ".linkage")), progress= FALSE)) %>%
          mutate(pid= "96") %>%
          filter(total > snvs.96.coverage*freq_threshold, complete.cases(.))

        snvs <- rbind(snvs.98, snvs.96) %>%
          rename(row_num= X1) %>%
          mutate(r2= ifelse(r2 > 1, 1, r2), # some r2 values seem to be just above 1
                 dist.bin= cut(Distance, breaks= 10),
                 scaffold= str_replace(Window, ":.*$", ""))
        #filter(total > 50, complete.cases(.))
        rm(snvs.98, snvs.96)


        #### SUMMARIZE BY DISTANCE BIN
        # snvs.sum <- snvs %>%
        #   group_by(pid, scaffold, dist.bin) %>%
        #   summarize(
        #     N= length(r2),
        #     r2.mean= mean(r2, na.rm= TRUE),
        #     r2.sd= sd(r2, na.rm= TRUE)) %>%
        #   ungroup() %>%
        #   filter(complete.cases(.))

        #### MAKE WIDE FORMAT AND MATRIX FOR HEATMAP
        ## Average r^2 value per scaffold
        #snvs.sum.wide <- snvs.sum %>%
        #                select(-N, -r2.sd) %>%
        #                spread(dist.bin, value= r2.mean)
        #snvs.sum.96.mat <- as.matrix(snvs.sum.wide[snvs.sum.wide$pid == "96", ncol(snvs.sum.wide):3])
        #snvs.sum.98.mat <- as.matrix(snvs.sum.wide[snvs.sum.wide$pid == "98", ncol(snvs.sum.wide):3])

        #### SAMPLE ID INFORMATION FOR OUTPUT FILES
        sample.name <- str_replace(linkage.files[file], "_0\\..*$", "")
        #PID <- str_replace(linkage.files[file], "^PH.*[A-Z]_\\.", "") %>%
        #      str_replace(., ".linkage", "")
        #plot.name <- str_c(sample.name, species.num, PID, sep= "-")
        plot.name <- str_c(sample.name, species.num, sep= "-")


        #### EXPONENTIAL DECAY MODEL
        ## use try() to make sure lm function did not error, otherwise loop crashes
        exp.model.98 <- try(lm(log(r2) ~ Distance, data= subset(snvs, pid== "98" & r2 > 0)))
        exp.model.96 <- try(lm(log(r2) ~ Distance, data= subset(snvs, pid== "96" & r2 > 0)))

        ## Make dataframes of predicted data for ggplots below
        if(class(exp.model.98) != "try-error"){
          exp.model.98.df <- data.frame(pid= "98",
                                        dist.values= subset(snvs, pid== "98")$Distance,
                                        pred.values= exp(predict(exp.model.98, list(Distance= subset(snvs, pid== "98")$Distance))))
        }
        if(class(exp.model.96) != "try-error"){
          exp.model.96.df <- data.frame(pid= "96",
                                        dist.values= subset(snvs, pid== "96")$Distance,
                                        pred.values= exp(predict(exp.model.96, list(Distance=subset(snvs, pid== "96")$Distance))))
        }

        ## Export model coefficients to a table
        if(all(exists("exp.model.96.df"), exists("exp.model.98.df"))){

          if(exists("counter") == FALSE){
            counter <- 1
          }

          mod.coef[counter, ] <- c(sample.name, species.num, "98", summary(exp.model.98)$df[2],
                                   round(exp.model.98$coefficients[1], 5), round(exp.model.98$coefficients[2], 5), round(summary(exp.model.98)$adj.r.squared, 4))
          counter <- counter + 1
          mod.coef[counter, ] <- c(sample.name, species.num, "96", summary(exp.model.96)$df[2],
                                   round(exp.model.96$coefficients[1], 5), round(exp.model.96$coefficients[2], 5), round(summary(exp.model.96)$adj.r.squared, 4))
          counter <- counter + 1
        }

        #### MAKE PLOTS ################################################################


        ## DISTANCE x R^2
        ## Only plot regression line if linear model did not error
        if(all(exists("exp.model.96.df"), exists("exp.model.98.df"))){
          ggplot(data= snvs, aes(x= Distance, y= r2)) +
            geom_hline(yintercept = 0, size= 0.25) +
            geom_point() +
            geom_line(data = exp.model.96.df, aes(x= dist.values, y = pred.values), color= "tomato", size= 1) +
            geom_line(data = exp.model.98.df, aes(x= dist.values, y = pred.values), color= "tomato", size= 1) +
            labs(x= "Pairwise distance (bp)", y= expression(r^2), title= plot.name) +
            y_axis_format +
            x_axis_format_distance +
            facet.by.pid +
            theme_snv
          #ggsave(last_plot(), filename= str_c(plot.name, "_dist", ".pdf"), width= 8, height= 6, units= "in", path= dir_output, device= cairo_pdf)
          ggsave(last_plot(), filename= str_c(plot.name, "_dist", ".jpg"), width= 8, height= 6, units= "in", dpi= 320,  path= dir_output)

        } else {

          ggplot(data= snvs, aes(x= Distance, y= r2)) +
            geom_hline(yintercept = 0, size= 0.25) +
            geom_point() +
            #geom_line(data = exp.model.96.df, aes(x= dist.values, y = pred.values), color= "tomato", size= 1) +
            #geom_line(data = exp.model.98.df, aes(x= dist.values, y = pred.values), color= "tomato", size= 1) +
            labs(x= "Pairwise distance (bp)", y= expression(r^2), title= plot.name) +
            y_axis_format +
            x_axis_format_distance +
            facet.by.pid +
            theme_snv
          #ggsave(last_plot(), filename= str_c(plot.name, "_dist", ".pdf"), width= 8, height= 6, units= "in", path= dir_output, device= cairo_pdf)
          ggsave(last_plot(), filename= str_c(plot.name, "_dist", ".jpg"), width= 8, height= 6, units= "in", dpi= 320,  path= dir_output)
        }

        ## Remove model dataframes for next iteration of the loop
        rm(exp.model.98.df, exp.model.96.df)

        # ## WINDOW x R^2
        #   distance.color.ramp <- colorRampPalette(brewer.pal(9, "PuBu"))(length(levels(snvs$dist.bin)))
        #
        #   ggplot(data= snvs, aes(x= Window, y= r2)) +
        #     geom_hline(yintercept = 0, size= 0.25) +
        #     geom_point(aes(fill= dist.bin), color= "black", size= 2, pch= 21) +
        #     #geom_errorbar() +
        #     labs(x= "Window", y= expression(r^2), title= plot.name) +
        #     scale_fill_manual(values= distance.color.ramp) +
        #     y_axis_format +
        #     x_axis_format_window +
        #     facet.by.pid +
        #     theme_snv
        #   #ggsave(last_plot(), filename= str_c(plot.name, "_window", ".pdf"), width= 8, height= 6, units= "in", path= dir_output, device= cairo_pdf)
        #   ggsave(last_plot(), filename= str_c(plot.name, "_window", ".jpg"), width= 8, height= 6, units= "in", dpi= 320,  path= dir_output)



        ## R^2 HISTOGRAM WITH DIST.BIN
        # ggplot(data= snvs, aes(x= r2)) +
        #   geom_histogram(aes(fill= dist.bin), color= "black", binwidth= 0.01) +
        #   #geom_errorbar() +
        #   #labs(x= "Pairwise distance (bp)", y= expression(r^2), title= plot.name) +
        #   scale_fill_manual(values= distance.color.ramp) +
        #   scale_y_log10(expand= c(0.005, 0)) +
        #   scale_x_continuous(expand= c(0.02, 0)) +
        #   facet.by.pid +
        #   theme_snv

        ## MAKE HEATMAP, IF SPECIFIED IN FUNCTION
        #if(heatmap == TRUE){
        # pheatmap(t(snvs.sum.96.mat),
        #          cluster_rows = FALSE,
        #          cluster_cols = FALSE)
        #         filename= file.path(dir_output, str_c(plot.name,"_heatmap.pdf")))
        #}

      }
    }


    ## Clean up mod.coef table
    mod.coef <- mod.coef %>%
      mutate_at(vars(intercept:r2_adj), funs(as.numeric)) %>%
      as_tibble()
    write_tsv(mod.coef, path= file.path(dir_output_table, str_c("model_coef2_", species.num, ".tsv")))


    ## PLOT HISTOGRAM OF MODEL SLOPES
    ggplot(mod.coef, aes(x= slope)) +
      geom_histogram(color= "black", fill= "snow3", binwidth= 0.001 ) +
      labs(x= "Slope coefficient", y= "Count", title= str_c("Slopes-", species.num)) +
      scale_x_continuous(expand= c(0.01, 0)) +
      facet.by.pid +
      theme_snv
    ggsave(last_plot(), filename= str_c("mod_coef_slope-", species.num,".jpg"), width= 8, height= 6, units= "in", dpi= 320,  path= dir_output)

    ## PLOT HISTOGRAM OF MODEL ADJUSTED R^2
    ggplot(coef.table, aes(x= r2_adj)) +
      geom_histogram(color= "black", fill= "snow3", binwidth= 0.01) +
      labs(x= expression(Adjusted~r^2), y= "Count", title= str_c("r^2_adj-", species.num)) +
      scale_x_continuous(expand= c(0.01, 0)) +
      facet.by.pid +
      theme_snv
    ggsave(last_plot(), filename= str_c("mod_coef_r2adj-", species.num,".jpg"), width= 8, height= 6, units= "in", dpi= 320,  path= dir_output)

    return(mod.coef)
  }


  if(freq_plot == TRUE){
    #### FREQ FILES ####
    ## progress bar
    pb.freq <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = length(linkage.files))
    pb.freq$tick(0)

    message("Looping over .freq files")
    for(file in 1:length(freq.files)){
      pb.freq$tick()

      #### READ IN THE DATA
      freq.df <- suppressMessages(read_tsv(file.path(dir_input, freq.files[file]), progress= FALSE)) %>%
        rename(row_num= X1)

      ## Count histogram frequencies
      freq.counts <- freq.df %>%
        mutate(freq.rounded= round(freq, 3)) %>%
        count(freq.rounded)

      #### SAMPLE ID INFORMATION FOR OUTPUT FILES
      sample.name <- str_replace(freq.files[file], "_0\\..*$", "")
      #PID <- str_replace(linkage.files[file], "^PH.*[A-Z]_", "") %>%
      #      str_replace(., ".linkage", "") %>%
      #      str_c(., "_PID")
      #plot.name <- str_c(sample.name, species.num, PID, sep= "-")
      plot.name <- str_c(sample.name, species.num, sep= "-")

      #### MAKE PLOTS
      # ggplot(data= freq.df, aes(x= freq)) +
      #   geom_histogram(binwidth= 0.01, color= "black", fill= "snow2") +
      #   labs(x= "SNV frequency", y= "Count", title= plot.name) +
      #   scale_x_continuous(expand= c(0.01, 0), breaks= seq(0, 1, by= 0.1)) +
      #   theme_snv
      # #ggsave(last_plot(), filename= str_c(plot.name, "_freq", ".pdf"), width= 8, height= 6, units= "in", path= dir_output, device= cairo_pdf)
      # ggsave(last_plot(), filename= str_c(plot.name, "_freq", ".jpg"), width= 8, height= 6, units= "in", dpi= 320, path= dir_output)

      ggplot(freq.counts, aes(x= freq.rounded, y= n)) +
        geom_vline(xintercept= 0.5, size= 0.25, color= "gray") +
        geom_point(size= 2, alpha= 0.5) +
        geom_smooth(method= "loess", span= 0.1, se= FALSE, size= 1, color= "tomato") +
        labs(x= "SNV frequency", y= "Count",  title= plot.name) +
        scale_x_continuous(breaks= seq(0, 1, 0.1), limits= c(0, 1), expand= c(0.01, 0)) +
        scale_y_continuous(breaks= scales::pretty_breaks()) +
        #facet_wrap(~sample, ncol= 3, scales= "free_y") +
        theme_snv
      ggsave(last_plot(), filename= str_c(plot.name, "_freq", ".jpg"), width= 8, height= 6, units= "in", dpi= 320, path= dir_output)
    }
  }



}



coef.table.sp1 <- analyse_linkage_data(path= dir_input_sp1, species= 1, linkage_plot = FALSE)
coef.table.sp2 <- analyse_linkage_data(path= dir_input_sp2, species= 2, linkage_plot = FALSE)

test.fig <- multi_panel_figure(width= 300, height= 300, unit= "mm", rows= 3, columns= 3)

for(file in list.files(dir_output, pattern= "-species_2_dist")){
  if(exists("multi_fig") == FALSE){
    multi_fig <- multi_panel_figure(width= 300, height= 300, unit= "mm", rows= 10, columns= 6, row_spacing= 0, column_spacing= 0, panel_label_type = "none")
    last_file <- list.files(dir_output, pattern= "-species_2_dist")[length(list.files(dir_output, pattern= "-species_2_dist"))]
  }

  multi_fig <- multi_fig %>%
    fill_panel(file.path(dir_output, file), scaling= "shrink")

  if(file == last_file){
    message("saving file")
    save_multi_panel_figure(multi_fig, filename= file.path(dir_output, "species_2_dist.jpg"), dpi= 320)
    rm(multi_fig, last_file)
  }

}


for(file in list.files(dir_output, pattern= "-species_1_dist")){
  if(exists("multi_fig") == FALSE){
    multi_fig <- multi_panel_figure(width= 300, height= 300, unit= "mm", rows= 5, columns= 6, row_spacing= 0, column_spacing= 0, panel_label_type = "none")
    last_file <- list.files(dir_output, pattern= "-species_1_dist")[length(list.files(dir_output, pattern= "-species_1_dist"))]
  }

  multi_fig <- multi_fig %>%
    fill_panel(file.path(dir_output, file), scaling= "shrink")

  if(file == last_file){
    message("saving file")
    save_multi_panel_figure(multi_fig, filename= file.path(dir_output, "species_1_dist.jpg"), dpi= 320)
    rm(multi_fig, last_file)
  }

}


for(file in list.files(dir_output, pattern= "-species_1_freq")){
  if(exists("multi_fig") == FALSE){
    multi_fig <- multi_panel_figure(width= 300, height= 300, unit= "mm", rows= 5, columns= 6, row_spacing= 0, column_spacing= 0, panel_label_type = "none")
    last_file <- list.files(dir_output, pattern= "-species_1_freq")[length(list.files(dir_output, pattern= "-species_1_freq"))]
  }

  multi_fig <- multi_fig %>%
    fill_panel(file.path(dir_output, file), scaling= "shrink")

  if(file == last_file){
    message("saving file")
    save_multi_panel_figure(multi_fig, filename= file.path(dir_output, "species_1_freq.jpg"), dpi= 320)
    rm(multi_fig, last_file)
  }

}

for(file in list.files(dir_output, pattern= "-species_2_freq")){
  if(exists("multi_fig") == FALSE){
    multi_fig <- multi_panel_figure(width= 300, height= 300, unit= "mm", rows= 10, columns= 6, row_spacing= 0, column_spacing= 0, panel_label_type = "none")
    last_file <- list.files(dir_output, pattern= "-species_2_freq")[length(list.files(dir_output, pattern= "-species_2_freq"))]
  }

  multi_fig <- multi_fig %>%
    fill_panel(file.path(dir_output, file), scaling= "shrink")

  if(file == last_file){
    message("saving file")
    save_multi_panel_figure(multi_fig, filename= file.path(dir_output, "species_2_freq.jpg"), dpi= 320)
    rm(multi_fig, last_file)
  }

}



