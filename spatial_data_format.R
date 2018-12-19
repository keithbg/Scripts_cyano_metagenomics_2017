
#### Libraries #################################################################
library(riverdist)
#source("/Users/KeithBG/R_functions/ggplot_scalebar_north_arrow.R")
################################################################################


#### FILE PATHS ################################################################
dir_input <- file.path("/Users","KeithBG","Documents","UC Berkeley","CyanoMeta_NSF","Metagenomics", "Data","EnvData")
################################################################################




#### CALCULATE DISTANCES BETWEEN SAMPLING SITES


 river.net <- line2network(path= file.path(dir_input, "GIS_data"), layer= "eel_rivernetwork_utm", reproject= "+proj=longlat +datum=WGS84")
# took 48 minutes to run
 saveRDS(river.net, file= file.path(dir_input, "eel_river_latlong_rdist.rds"))

 # river.net <- readRDS(file.path(dir_input, "eel_river_rdist.rds"))
#
#
# river.cleanup <- function(){
# river.net.A <- removeduplicates(river.net)
# river.net.B <- dissolve(river.net.A)
# river.net.C <- splitsegments(river.net.B)
# river.net.D <- addverts(river.net.C, 100)
# }
#
# river.net.clean1 <- river.cleanup()
# saveRDS(river.net.clean1, file= file.path(dir_input, "river.net.clean1.rds"))

#river.net.clean1 <- readRDS(file.path(dir_input, "river.net.clean1.rds"))
#river.net.clean2 <- removemicrosegs(river.net.clean1)
#saveRDS(river.net.clean2, file= file.path(dir_input, "river.net.clean2.rds"))

#river.net2 <- cleanup(river.net)
#,

#river.net.clean2 <- readRDS(file.path(dir_input, "river.net.clean2.rds"))






#### CONVERT LAT LONG TO UTM PROJECTION
# site.lat.long <- read_csv(dir_input, "PhormMeta17_LatLong.csv")
#
# line98albers <- project(line98,proj="+proj=aea +lat_1=55 +lat_2=65
#                         +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
#                         +ellps=GRS80 +towgs84=0,0,0")


#### Make base map ####

## Read in Eel River watershed shape file
#   eel.network.shp.import <- function(){
#     eel.network.shp <- readOGR(dsn= file.path(dir_input, "GIS_data"), layer= "eel_rivernetwork_utm")
#     # ogrInfo(dsn= "/Users/KeithBG/Documents/UC Berkeley/Manuscripts/CyanotoxinSurvey/GIS_data", layer= "eel_rivernetwork_utm")
#     # class(eel.network.shp) # SpatialLinesDataFrame
#     # proj4string(eel.network.shp) # projection of shapefile
#     eel.network.shp.subset <- subset(eel.network.shp, CUMDRAINAG > 20) # remove drainage areas <20 km^2
#     eel.network.shp.reproj <- spTransform(eel.network.shp.subset, CRS("+proj=longlat +datum=WGS84")) # reproject to match ggmap
#     eel.network.df <- tidy(eel.network.shp.reproj)
#     return(eel.network.df)
#   }
#   eel.network.df <- eel.network.shp.import()
#
# ## Read in state boundary shape file
#   CA.map.shp.import <- function(){
#     CA <- readOGR(dsn= file.path(dir_input, "GIS_data"), layer= "CA_boundary_kbg")
#     CA.reproj <- spTransform(CA, CRS("+proj=longlat +datum=WGS84")) # reproject to match ggmap
#     CA.df <- tidy(CA.reproj)
#     return(CA.df)
#   }
#   CA.df <- CA.map.shp.import()
#
# ## Make base_map
#   make_base_map <- function(){
#   spatt_map_theme <- theme(text= element_text(size= 18),
#                            panel.background = element_rect(fill= "light blue"),
#                            panel.border = element_rect(color= "black", fill= NA),
#                            legend.key= element_rect(fill= "white"),
#                            plot.background = element_rect(fill= "transparent", color= NA),
#                            panel.grid.minor = element_blank(),
#                            panel.grid.major = element_blank()
#                            )
#
#   base_map <- ggplot() +
#     geom_polygon(data= CA.df, aes(x= long, y= lat, group=group), color= "black", fill= "snow1") +
#     geom_path(data= eel.network.df, aes(x= long, y= lat, group=group), color= "#2A788EFF") +
#     scale_bar(lon = -124.4, lat = 39.3,
#               distance_lon = 20, distance_lat = 3, distance_legend = 6, dist_unit = "km",
#               arrow_length= 10, arrow_distance = 8) +
#     coord_map(xlim= c(-124.5, -122.7), ylim= c(39.25, 40.75)) +
#     labs(x= "Longitude", y= "Latitude") +
#     spatt_map_theme
#   return(base_map)
#   }
#   base_map <- make_base_map()
#
