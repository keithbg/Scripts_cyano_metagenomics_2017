## Make a map of the Eel and Russian river network

# HUC8 for Russian River: 18010110
# Drainage area column for Russian = TotDASqKM
# Drainage area column for Eel = CUMDRAINAG

#### Libraries #################################################################
library(tidyverse)
library(ggmap)
#library(ggsn)
library(rgdal)
library(broom)
source("/Users/KeithBG/R_functions/ggplot_scalebar_north_arrow.R")

################################################################################


#### FILE PATHS ################################################################
shp_input_russian <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","EnvData", "GIS_files", "Russian_river_network")
shp_input_eel <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","EnvData", "GIS_files", "Eel_river_network_shapefile")
shp_input_CA <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","EnvData", "GIS_files", "CA_boundary_shapefile")
kml_input <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","EnvData", "GIS_files")
dir_output <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","EnvData")
################################################################################


## River networks
river.network.shp.import <- function(shp_pathway, layer_name, drainage_column= NA, min_drainage_area= 0){
  river.network.shp <- readOGR(dsn= shp_pathway, layer= layer_name)
  # Get column names in the attributes table
  print(ogrInfo(dsn= shp_pathway, layer= layer_name))

  # class(russian.network.shp) # SpatialLinesDataFrame
  # proj4string(russian.network.shp) # projection of shapefile
  if(min_drainage_area > 0){
    print(paste("Removing drainage areas <", min_drainage_area))
    col_num <- which(names(river.network.shp@data) %in% drainage_column) # get column number of drainage area data
    river.network.shp <- river.network.shp[river.network.shp@data[, col_num] >= min_drainage_area, ] # keep drainage areas >10 km^2
  }
  river.network.shp <- spTransform(river.network.shp, CRS("+proj=longlat +datum=WGS84")) # reproject to match ggmap
  river.network.shp <- tidy(river.network.shp)
  return(river.network.shp)
}

russian.network <- river.network.shp.import(shp_pathway = shp_input_russian,
                                            layer_name = "Russian_river_network",
                                            drainage_column = "TotDASqKM",
                                            min_drainage_area = 20)
eel.network <- river.network.shp.import(shp_pathway = shp_input_eel,
                                        layer_name = "eel_rivernetwork_utm",
                                        drainage_column = "CUMDRAINAG",
                                        min_drainage_area = 20)

##  Watershed outlines
eel.watershed <- readOGR(file.path(kml_input, "Eel_River_Drainage_Area.kml"), layer= "Eel_River_Drainage_Area") %>%
                   spTransform(., CRS("+proj=longlat +datum=WGS84")) %>%  # reproject to match ggmap
                   tidy(.)

russian.watershed <- readOGR(file.path(kml_input, "Russian_River_Drainage_Area.kml"), layer= "Russian_River_Drainage_Area") %>%
  spTransform(., CRS("+proj=longlat +datum=WGS84")) %>%  # reproject to match ggmap
  tidy(.)


## California Outline
CA.outline <- readOGR(dsn= shp_input_CA, layer= "CA_boundary_kbg") %>%
                       spTransform(., CRS("+proj=longlat +datum=WGS84")) %>%  # reproject to match ggmap
                       tidy(.)

## Make base_map
  PH2017_map_theme <- theme(text= element_text(size= 14),
                           panel.background = element_rect(fill= "light blue"),
                           panel.border = element_rect(color= "black", fill= NA),
                           legend.key= element_rect(fill= "white"),
                           plot.background = element_rect(fill= "transparent", color= NA),
                           panel.grid.minor = element_blank(),
                           panel.grid.major = element_blank()
  )

  PH2017_eel_russian_base_map <- ggplot() +
    geom_polygon(data= CA.outline, aes(x= long, y= lat, group=group), color= "black", fill= "snow1") +
    geom_polygon(data= russian.watershed, aes(x= long, y= lat, group=group), color= "black", fill= NA) +
    geom_polygon(data= eel.watershed, aes(x= long, y= lat, group=group), color= "tomato", fill= NA) +
    geom_path(data= russian.network, aes(x= long, y= lat, group=group), color= "#2A788EFF") +
    geom_path(data= eel.network, aes(x= long, y= lat, group=group), color= "#2A788EFF") +
    scale_bar(lon = -124.4, lat = 38.3,
              distance_lon = 20, distance_lat = 3, distance_legend = 6, dist_unit = "km",
              arrow_length= 10, arrow_distance = 8) +
    coord_map(xlim= c(-124.5, -122.4), ylim= c(38.25, 40.75)) +
    scale_y_continuous(breaks= c(38.5, 39, 39.5, 40, 40.5), labels= c("", "39", "", "40", "")) +
    scale_x_continuous(breaks= c(-124.5, -124, -123.5, -123.0, -122.5), labels= c("", "-124", "", "-123", "")) +
    labs(x= "Longitude", y= "Latitude") +
    PH2017_map_theme




