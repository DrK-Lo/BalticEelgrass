###########################
## EELGRASS SITE MAPPING ##
###########################

## setup
########
# set wd
setwd("~/Documents/Github/BalticEelgrass")

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
library(rvg) # convert to vector graphic
library(officer) # export to ppt
########

## data
#######
# read in data
sites <- read.csv("data/experiment/eelgrass_exp_sites.csv")
sites

seascape <- read.csv("data/seascape_data/sampling_sites_coordinates_Baltic_Sea.csv")
#######

## mapping sites
################
all_sites <- ggplot(data = world) + 
  geom_sf(data = world_crop, fill = "gray90")+
  theme_classic()+
  labs(title = "Map of Seascape Training and \nExperimental Source Sites",
       x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(3, 32.5), ylim = c(53.1, 66.1), expand = FALSE) +
  geom_point(data = seascape, aes(x = (long),
                                  y = (lat)),
             color = "black", fill = "gray30", 
             shape = 23, size = 1.5) +
  geom_point(data = sites, aes(x = (long),
                               y = (lat),
                               fill = site_abbrev),
             color = "black", shape = 21, size = 4) +
  scale_fill_manual(name = "Experimental \nSource Site",
                    values = pop_colors,
                    breaks = names(pop_orders)) +
  geom_text(label = "North \nSea", 
            x = 5.8, y = 56.5,
            size = 3, fontface = "italic") +
  geom_text(label = "Baltic \nSea", 
#            x = 18.9, y = 56,
            x = 20.5, y = 58.5,
            size = 3, fontface = "italic") +
  theme(plot.title = element_text(size = 18), 
        panel.grid.major = element_line(color = "aliceblue"), 
        panel.background = element_rect(fill = "aliceblue"), 
        legend.position = "right")
all_sites
################

## export
#########
doc <- read_pptx() %>% 
  add_slide(layout = "Title and Content", master = "Office Theme") %>% 
  ph_with(value = dml(ggobj = all_sites), 
          location = ph_location(left = 1, top = 1, width = 6.5, height = 4))
print(doc, "results/figures/site_map.pptx")
#########

