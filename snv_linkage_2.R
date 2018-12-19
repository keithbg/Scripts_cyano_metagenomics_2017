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
  dir_input_sp2 <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "snv_linkage", "species_2")
  dir_output_fig <- file.path("/Users","KeithBG","Documents","UC\ Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "snv_linkage", "Figures")
  dir_output_table <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "snv_linkage", "Output_tables")
################################################################################

analyse_linkage_data <- function(freq_threshold= 0.05, species, pid, path, output, linkage_plot= TRUE, freq_plot= TRUE){
## CHECK FOR ERRORS ############################################################
  if(str_detect(as.character(pid), "\\.")){
    stop("pid must be an integer with no decimals (e.g. 0.99 = 99)")
  }

## DEFINE VARIABLES ############################################################
  dir_input <- path
  dir_output <- output
  species.num <- paste0("species_", species)
  pid_value <- as.character(pid)
  linkage.files <- list.files(file.path(dir_input), pattern= "linkage")
  linkage.sample <- str_replace(linkage.files, "_0\\..*$", "")
  freq.files <- list.files(file.path(dir_input), pattern= "freq")
  freq.sample <- str_replace(freq.files, "_0\\..*$", "")

  ##### PLOTTING PARAMETERS ####################################################
  x_axis_format_distance <- scale_x_continuous(expand= c(0, 5))
  x_axis_format_window <- scale_x_discrete(labels= NULL, expand= c(0.02, 0))
  y_axis_format <- scale_y_continuous(limits= c(0, 1),
                                      breaks= seq(0, 1, by= 0.25),
                                      expand= c(0.02, 0))



  #pid.facet.labels <- as_labeller(c(`98` = "PID= 0.98", `96` = "PID= 0.96" ))
  #facet.by.pid <- facet_grid(pid~., labeller= labeller(pid= pid.facet.labels, scales= "free_y"))


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


  ##### LOOP THROUGH .LINKAGE FILES AND MAKE PLOTS #########################################

  if(linkage_plot == TRUE){

    ## progress bar
    pb.linkage <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = length(linkage.files))
    pb.linkage$tick(0)

    rm(counter)

    ## Initialize data frame to store model coefficients
    mod.coef <- data.frame(sample= as.character(NULL), species= as.character(NULL), pid= as.character(NULL), df= as.character(NULL), intercept= as.numeric(NULL), slope= as.numeric(NULL), r2_adj= as.numeric(NULL), stringsAsFactors= FALSE)

    #### THE LOOP ####
    message("Looping over .linkage files")

    for(file in 1:length(linkage.files)){
      pb.linkage$tick() ## Progress Bar


      #### READ IN THE DATA

      ## Coverage data from log file
      snvs.coverage <- try(suppressMessages(read_tsv(file.path(dir_input, "log_files", str_c(linkage.sample[file], ".", species.num, ".pid", pid_value, ".log")),
                                                        col_names = FALSE,
                                                        progress= FALSE)) %>%
                                pull() %>%
                                .[str_detect(., "Mean coverage")] %>%
                                str_replace("Mean coverage: ", "") %>%
                                as.numeric(.) %>%
                                round(., 1))

      if(class(snvs.coverage) != "try-error"){
        ## Linkage data
        snvs <- suppressMessages(read_tsv(file.path(dir_input, str_c(linkage.sample[file], "_0.", 98, ".linkage")), progress= FALSE)) %>%
          mutate(pid= pid_value) %>%
          filter(total > snvs.coverage*freq_threshold, complete.cases(.)) %>%
          rename(row_num= X1) %>%
          mutate(r2= ifelse(r2 > 1, 1, r2), # some r2 values seem to be just above 1
                 dist.bin= cut(Distance, breaks= 10),
                 scaffold= str_replace(Window, ":.*$", ""))


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
        plot.name <- str_c(sample.name, species.num, str_c("PID=0.", pid_value), sep= "-")

        #### EXPONENTIAL DECAY MODEL
        ## use try() to make sure lm function did not error, otherwise loop crashes
        exp.model <- try(lm(log(r2) ~ Distance, data= subset(snvs, r2 > 0)))

        ## Make dataframes of predicted data for ggplots below
        if(class(exp.model) != "try-error"){
          exp.model.df <- data.frame(pid= pid_value,
                                        dist.values= snvs$Distance,
                                        pred.values= exp(predict(exp.model, list(Distance= snvs$Distance))))
        }


        ## Export model coefficients to a table
        if(exists("exp.model.df")){

          if(exists("counter") == FALSE){
            counter <- 1
          }

          mod.coef[counter, ] <- c(sample.name, species.num, pid_value, summary(exp.model)$df[2],
                                   round(exp.model$coefficients[1], 5), round(exp.model$coefficients[2], 5), round(summary(exp.model)$adj.r.squared, 4))
          counter <- counter + 1
        }

        #### MAKE PLOTS ################################################################


        ## DISTANCE x R^2
        ## Only plot regression line if linear model did not error
        if(exists("exp.model.df")){
          ggplot(data= snvs, aes(x= Distance, y= r2)) +
            geom_hline(yintercept = 0, size= 0.25) +
            geom_point() +
            geom_line(data = exp.model.df, aes(x= dist.values, y = pred.values), color= "tomato", size= 1) +
            labs(x= "Pairwise distance (bp)", y= expression(r^2), title= plot.name) +
            y_axis_format +
            x_axis_format_distance +
            theme_snv
          #ggsave(last_plot(), filename= str_c(plot.name, "_dist", ".pdf"), width= 8, height= 6, units= "in", path= dir_output, device= cairo_pdf)
          ggsave(last_plot(), filename= str_c(plot.name, "_dist", ".jpg"), width= 8, height= 6, units= "in", dpi= 320,  path= dir_output)

        } else {

          ggplot(data= snvs, aes(x= Distance, y= r2)) +
            geom_hline(yintercept = 0, size= 0.25) +
            geom_point() +
            #geom_line(data = exp.model.df, aes(x= dist.values, y = pred.values), color= "tomato", size= 1) +
            labs(x= "Pairwise distance (bp)", y= expression(r^2), title= plot.name) +
            y_axis_format +
            x_axis_format_distance +
            theme_snv
          #ggsave(last_plot(), filename= str_c(plot.name, "_dist", ".pdf"), width= 8, height= 6, units= "in", path= dir_output, device= cairo_pdf)
          ggsave(last_plot(), filename= str_c(plot.name, "_dist", ".jpg"), width= 8, height= 6, units= "in", dpi= 320,  path= dir_output)
        }

        ## Remove model dataframes for next iteration of the loop
        rm(exp.model.df)

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
    message("Writing model coefficients table")
    mod.coef <- mod.coef %>%
      mutate_at(vars(intercept:r2_adj), funs(as.numeric)) %>%
      as_tibble()
    write_tsv(mod.coef, path= file.path(dir_output_table, str_c("model_coef_", species.num, ".tsv")))


    ## PLOT HISTOGRAM OF MODEL SLOPES
    ggplot(mod.coef, aes(x= slope)) +
      geom_histogram(color= "black", fill= "snow3", binwidth= 0.001 ) +
      labs(x= "Slope coefficient", y= "Count", title= str_c("Slopes-", species.num)) +
      scale_x_continuous(expand= c(0.01, 0)) +
      theme_snv
    ggsave(last_plot(), filename= str_c("mod_coef_slope-", species.num,".jpg"), width= 8, height= 6, units= "in", dpi= 320,  path= dir_output)

    ## PLOT HISTOGRAM OF MODEL ADJUSTED R^2
    ggplot(mod.coef, aes(x= r2_adj)) +
      geom_histogram(color= "black", fill= "snow3", binwidth= 0.01) +
      labs(x= expression(Adjusted~r^2), y= "Count", title= str_c("r^2_adj-", species.num)) +
      scale_x_continuous(expand= c(0.01, 0)) +
      theme_snv
    ggsave(last_plot(), filename= str_c("mod_coef_r2adj-", species.num,".jpg"), width= 8, height= 6, units= "in", dpi= 320,  path= dir_output)

  }

  #### LOOP OVER .FREQ FILES ###################################################
  if(freq_plot == TRUE){

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
      plot.name <- str_c(sample.name, species.num, str_c("PID=0.", pid_value), sep= "-")

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

  return(mod.coef)
}

coef.table.sp2.pid99 <- analyse_linkage_data(path= dir_input_sp2, pid= 99, species= 2, output= dir_output_fig)
coef.table.sp2.pid98 <- analyse_linkage_data(path= dir_input_sp2, pid= 98, species= 2, output= dir_output_fig)

#coef.table.sp1 <- analyse_linkage_data(path= dir_input_sp1, species= 1, linkage_plot = FALSE)



make_multi_panel_fig <- function(path, file.pattern){

  # Get number of figures in folder
  num.figures <- length(list.files(path, pattern= file.pattern))
  last_file <- list.files(path, pattern= file.pattern)[num.figures]

  # Make blank multi-panel figure
  if(exists("multi_fig") == FALSE){
    multi_fig <- multi_panel_figure(width= 300, height= 30*ceiling(num.figures/6), unit= "mm", rows= ceiling(num.figures/6), columns= 6, row_spacing= 0, column_spacing= 0, panel_label_type = "none")

  }

  # Loop to fill panels
  message("Filling panels")
  pb.panels <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = num.figures)
  pb.panels$tick(0)

  for(file in list.files(path, pattern= file.pattern)){
    pb.panels$tick() ## Progress Bar

    multi_fig <- suppressMessages(fill_panel(figure= multi_fig, panel= file.path(path, file), scaling= "shrink"))
  }

  # Save figure
  message("Saving file")
  save_multi_panel_figure(multi_fig, filename= file.path(path, str_c(file.pattern, ".jpg")), dpi= 320)
  rm(multi_fig, last_file)
}


make_multi_panel_fig(path= dir_output_fig, file.pattern= "species_2-PID=0.99_freq")
make_multi_panel_fig(path= dir_output_fig, file.pattern= "species_2-PID=0.99_dist")

make_multi_panel_fig(path= dir_output_fig, file.pattern= "species_2-PID=0.98_freq")
make_multi_panel_fig(path= dir_output_fig, file.pattern= "species_2-PID=0.98_dist")

