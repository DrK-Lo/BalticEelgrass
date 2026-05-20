##########################################
## iNaturalist eelgrass occurrence data ##
##########################################

# set wd
setwd("Documents/Github/BalticEelgrass")

# packages
library(readr)
library(tidyverse)
library(dplyr)
library(ggplot2) # for plotting
library(raster) # for environmental pixels
library(sf) # mapping
library(s2) # mapping
library(rnaturalearth) # mapping
library(maps) # mapping
library(ggspatial) # mapping
library(ggrepel) # mapping
library(stars)

# read in data
occurrence <- read_tsv("data/iNat/occurrence.txt")
seascape_sites <- read.csv("data/seascape_data/sampling_sites_coordinates_Baltic_Sea.csv")

# check it
head(occurrence)
colnames(occurrence)

# bounds
map_bounds <- data.frame(xmin = 10, xmax = 20, ymin = 54, ymax = 60)

# filter
occur <- occurrence %>% 
  dplyr::select(gbifID, rightsHolder, datasetName, scientificName, species,
         eventDate, year, month, day, verbatimLocality, decimalLatitude, 
         decimalLongitude, stateProvince, coordinateUncertaintyInMeters) %>% 
  dplyr::filter(datasetName == "iNaturalist research-grade observations" & 
           coordinateUncertaintyInMeters <= 1000 ) %>%
  as.data.frame()
occur

# with seascape
occur_sea <- seascape_sites %>% 
  dplyr::select(region, site_full, site, lat, long) %>%
  dplyr::filter(lat <= map_bounds$ymax & lat >= map_bounds$ymin & long <= map_bounds$xmax & long >= map_bounds$xmin)

# map these distributions
# set spherical geometry
sf::sf_use_s2(FALSE)

# set object
world <- ne_countries(scale = "medium", returnclass = "sf")
world_crop <- sf::st_crop(world, c(xmin = 10, xmax = 20, ymin = 54, ymax = 60))

map_seagrass <- ggplot(data = world) +
  theme_bw()+
  geom_sf(data = world_crop, fill = 'antiquewhite1') +
  geom_point(data = occur, aes(x = decimalLongitude, y = decimalLatitude), fill = "darkolivegreen1", shape = 21, size = 2) +
  geom_point(data = occur_sea, aes(x = long, y = lat), fill = "tomato", shape = 22, size = 2) +
  theme(plot.title = element_text(size = 18), panel.grid.major = element_line(color = "aliceblue"),
        panel.background = element_rect(fill = "aliceblue"), legend.position = 'right')+
  labs(title = "Eelgrass Occurrence", subtitle = "iNaturalist Observation Data", x = "Longitude", y = "Latitude")
map_seagrass

# this dataset is clearly incomplete, likely because the East coast of Sweden is less populated & meadows are generally deeper/sparser

# add sampling sites from seascape?

