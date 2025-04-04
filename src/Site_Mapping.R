###########################
## EELGRASS SITE MAPPING ##
###########################

## setup
########
# set wd
setwd("~/Eelgrass/BalticEelgrass/data")

# load packages
library(BiocManager) # needed to download some packages
library(tidyr) # for data wrangling
library(ggplot2) # plotting
library(s2) # mapping
library(rnaturalearth) # mapping
library(rnaturalearthdata) # mapping
library(maps) # mapping
library(ggspatial) # mapping
library(viridis) # color palettes
########

## data
#######
# read in data
sites <- read.csv("~/Eelgrass/BalticEelgrass/data/experiment/eelgrass_exp_sites.csv")
sites

seascape <- read.csv("~/Eelgrass/BalticEelgrass/data/seascape_data/sampling_sites_coordinates_Baltic_Sea.csv")
#######

## mapping sites
################
# set s2 false
sf::sf_use_s2(FALSE)

# crop sf
world <- ne_countries(scale = "medium", returnclass = "sf")
world_crop <- sf::st_crop(world, c(xmin =.3 , xmax = 33.5, ymin = 53.1, ymax = 66.1))

# map experimental sites
exp_sites <- ggplot(data = world) + 
  geom_sf(data = world_crop, fill = 'white')+
  theme_classic()+
  labs(title = "Map of Experimental Population Source Sites",
       x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(.3, 33.5), ylim = c(53.1, 66.1), expand = FALSE) +
  geom_point(data = sites, aes(x = (long),
                               y = (lat),
                               fill = site_abbrev),
             color = "black", shape = 21, size = 2.5)+
  theme(plot.title = element_text(size = 18), panel.grid.major = element_line(color = "white"), panel.background = element_rect(fill = "white"), legend.position = 'right') +
scale_fill_manual(name = "Population",
                  values = c("#155d27","#25a244","#4ad66d","#b7efc5",
                             "#caf0f8","#90e0ef","#0096c7","#023e8a"),
                  breaks = c("VIK","VAT","HOG","BAR",
                             "KUR","KAL","HOR","BJO"))
exp_sites

# generate palette
seascape_colors <- magma(39)

# map seascape sites
# this isn't a perfect map, but for my current purposes works
training_sites <- ggplot(data = world) + 
  geom_sf(data = world_crop, fill = 'white')+
  theme_classic()+
  labs(title = "Map of Seascape Population for Training",
       x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(.3, 33.5), ylim = c(53.1, 66.1), expand = FALSE) +
  geom_point(data = seascape, aes(x = (long),
                               y = (lat),
                               fill = as.factor(pop_order)),
             color = "black", shape = 21, size = 2.5) +
  theme(plot.title = element_text(size = 18), panel.grid.major = element_line(color = "white"), panel.background = element_rect(fill = "white"), legend.position = 'right') +
  scale_fill_viridis(discrete = T, option = "A")
training_sites
################