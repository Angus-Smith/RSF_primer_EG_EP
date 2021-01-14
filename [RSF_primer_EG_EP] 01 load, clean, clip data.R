## Angus Smith MSc thesis

#####################################################################
#### [RSF_primer_EG_EP] 01 load, clean, clip data
#####################################################################

###########
## TO DOs
# - 



##########################
#### Overview
##########################
## Load point data and landscape data. Prep for use in the RSF cours:
## https://terpconnect.umd.edu/~egurarie/teaching/SpatialModelling_AKTWS2018/6_RSF_SSF.html
## 


##########################
#### Settings
##########################

rm(list = ls())



##########################
#### Setup
##########################

setwd("~/Google_Drive/projects/MSc/RSF_primer")

source("RSF_primer_EG_EP/[RSF_primer_EG_EP] functions.R")

# pcks <- list("sp", "sf", "raster", "rgdal", "rgeos", "ggmap", "adehabitatHR", "sjPlot")

# sapply(pcks, require, char = TRUE)

library(tidyverse)
library(lubridate)
library(sf)
library(stars)

# global variables
# spatial reference
SR <- "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" # NAD83 / BC Albers


##########################
#### Load data
##########################

# road data
road <- st_read("data/DRA_DGTL_ROAD_ATLAS_MPAR_SP/DRA_MPAR_line.shp")

# elevation data
elev <- read_stars("data/dataset/DEM.tif")

# historical fire
fire <- st_read("data/PROT_HISTORICAL_FIRE_POLYS_SP/H_FIRE_PLY_polygon.shp")

# point data
elk <- read.csv("data/SPI_TELEMETRY_OBS_NONSENS_SP/SPTLMTRBSN.csv", header = T, sep = ",")



##########################
#### prep data
##########################

# elk to sf
elk <- st_as_sf(elk, coords = c("X", "Y"), crs = st_crs(road))

# reduce elk to elk
elk <- filter(elk, SPCSNGLSHN == "Roosevelt Elk")

# parse datetime for elk
elk$dttm <- ymd_hm(paste(elk$BSRVTNDT, elk$BSRVTNHR, elk$BSRVTNMNT))

# crop elev to elk (plus buffer)
elk_nad83 <- elk %>% st_buffer(dist = 500) %>% st_transform(st_crs(elev))

elev <- st_crop(elev, elk_nad83 %>% st_bbox())

# elev to proper crs
# write_stars(elev, "data/elevnad83.tif")
## done in QGIS manually
elev <- read_stars("data/elevproj.tif")

# reduce road to road
road <- filter(road, FTYPE %in% c("Road", "Bridge"))

# crop road to elev (plus buffer)
road <- st_crop(road, st_bbox(elev) %>% st_as_sfc() %>% st_buffer(dist = 5000))

# get distance to road raster
# st_write(road, "data/road.shp")
## done in QGIS manually
dist.roads <- read_stars("data/road_dist.tif")

# Get time since fire raster
## vague, will differ by animal. Not accounting for that. Taking max year.
fire$ttfire <- max(year(elk$dttm)) - fire$FIRE_YEAR
fire <- filter(fire, ttfire >= 0)
# st_write(fire, "data/fire.shp")
## done in QGIS manually
time.since.fire <- read_stars("data/fire_time.tif")



##########################
#### save data
##########################

write_stars(elev, "data/RSF_layers/elev.tif")
write_stars(time.since.fire, "data/RSF_layers/tsfire.tif")
write_stars(dist.roads, "data/RSF_layers/distroads.tif")

elk <- as_Spatial(elk)

save(elk, file = "data/RSF_layers/elk.rda")

