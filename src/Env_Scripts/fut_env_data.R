#####################
## future env data ##
#####################

## setup
########
setwd("~/Documents/GitHub/BalticEelgrass")

# packages
library(LEA)
library(tidyr)
library(dplyr)
library(readxl) # reading xlsx files
library(psych) # environmental data cleaning
library(vegan) # environmental data cleaning
library(ggplot2) # plotting
library(forcats) # order data
library(sdmpredictors) # data download
library(raster) # spatial data
library(s2) # mapping
library(rnaturalearth) # mapping
library(rnaturalearthdata) # mapping

# data
seascape_dat <- read.csv("data/EnvDat/summary_env_data_meadows_GEA.csv") # seascape env data
exp_dat_cop <- read.csv("data/EnvDat/complete_expSites_copernicus_2025-12-01.csv")[,-1]
exp_dat <- read_excel("data/EnvDat/Camille_env_data_experiment.xlsx")
########

## prep data
############
# seascape env data reads in funky, change column names and remove first row
colnames(seascape_dat) <- seascape_dat[1,]
seascape_dat <- seascape_dat[-1,]

# the column names have spaces in them, change names to be more code friendly
colnames(seascape_dat) <- c("region","site_full","pop","long","lat","mean_temp_sea","mean_sal_sea","temp_8.5_sea","sal_8.5_sea","max_temp_sea","min_sal_sea")

# also deal with lat/long
seascape_dat$long <- as.numeric(gsub(",", ".", gsub("\\.", "", seascape_dat$long)))
seascape_dat$lat <- as.numeric(gsub(",", ".", gsub("\\.", "", seascape_dat$lat)))

# and ensure consistency in site abbreviations by getting rid of funky characters
seascape_dat$pop[seascape_dat$pop == "ÅLA"] <- "ALA"
seascape_dat$pop[seascape_dat$pop == "BÅD"] <- "BAD"

# exp dat also has these column name issues
exp_dat2 <- exp_dat[,1:10]
colnames(exp_dat2) <- c("site_full","pop","long","lat","mean_temp_c",
                        "temp_4.5_c","temp_8.5_c","mean_sal_c",
                        "sal_4.5_c","sal_8.5_c")

# copernicus data
exp_dat_cop2 <- exp_dat_cop %>% 
  dplyr::select(pop, sal_med, sal_min, temp_med, temp_max) %>% 
  rename(median_sal_cop = sal_med, min_sal_cop = sal_min,
         median_temp_cop = temp_med, max_temp_cop = temp_max)
############

## extract layers
#################
# sites
points_sea <- vect(seascape_dat[,c("pop", "lat", "long")], geom = c("long", "lat"), crs = "EPSG:4326")

# current
curr_layers <- list_layers(marine = T)
curr <- terra::rast(load_layers(c("BO2_salinitymean_ss","BO2_salinitymin_ss",
                                  "BO2_tempmean_ss","BO2_tempmax_ss")))
curr_environment <- data.frame(name = seascape_dat$site_full, 
                               pop = seascape_dat$pop,
                               lat = seascape_dat$lat,
                               long = seascape_dat$long,
                               order = 1:39,
                               BO_salmean = terra::extract(curr$BO2_salinitymean_ss, 
                                                           points_sea),
                               BO_salmin = terra::extract(curr$BO2_salinitymin_ss, 
                                                          points_sea),
                               BO_tempmean = terra::extract(curr$BO2_tempmean_ss, 
                                                            points_sea),
                               BO_tempmax = terra::extract(curr$BO2_tempmax_ss, 
                                                           points_sea))
curr_environment2 <- curr_environment[, !grepl(".ID", names(curr_environment))]
colnames(curr_environment2) <- c("name", "pop", "lat", "long", "order",
                                   "BO_salmean", "BO_salmin", 
                                   "BO_tempmean", "BO_tempmax")

# future
# check layers
future_layers <- list_layers_future(marine = T)
unique(future_layers$scenario)
future45 <- get_future_layers(c("BO2_salinitymean_ss","BO2_salinitymin_ss",
                                "BO2_tempmean_ss","BO2_tempmax_ss"),
                            scenario = "RCP45", year = 2100)
future85 <- get_future_layers(c("BO2_salinitymean_ss","BO2_salinitymin_ss",
                                "BO2_tempmean_ss","BO2_tempmax_ss"),
                              scenario = "RCP85", year = 2100)

# load layers
future45_load <- terra::rast(load_layers(layercodes = future45$layer_code,
                             rasterstack = T))
future85_load <- terra::rast(load_layers(layercodes = future85$layer_code,
                             rasterstack = T))

# extract
future_environment <- data.frame(name = seascape_dat$site_full, 
                                 pop = seascape_dat$pop,
                                 lat = seascape_dat$lat,
                                 long = seascape_dat$long,
                                 order = 1:39,
                                 BO_salmean45 = terra::extract(future45_load$BO2_RCP45_2100_salinitymean_ss, 
                                                               points_sea),
                                 BO_salmin45 = terra::extract(future45_load$BO2_RCP45_2100_salinitymin_ss, 
                                                              points_sea),
                                 BO_tempmean45 = terra::extract(future45_load$BO2_RCP45_2100_tempmean_ss, 
                                                               points_sea),
                                 BO_tempmax45 = terra::extract(future45_load$BO2_RCP45_2100_tempmax_ss, 
                                                              points_sea),
                                 BO_salmean85 = terra::extract(future85_load$BO2_RCP85_2100_salinitymean_ss, 
                                                               points_sea),
                                 BO_salmin85 = terra::extract(future85_load$BO2_RCP85_2100_salinitymin_ss, 
                                                              points_sea),
                                 BO_tempmean85 = terra::extract(future85_load$BO2_RCP85_2100_tempmean_ss, 
                                                               points_sea),
                                 BO_tempmax85 = terra::extract(future85_load$BO2_RCP85_2100_tempmax_ss, 
                                                              points_sea))
future_environment2 <- future_environment[, !grepl(".ID", names(future_environment))]

# colnames
colnames(future_environment2) <- c("name", "pop", "lat", "long", "order",
                                  "BO_salmean45", "BO_salmin45", 
                                  "BO_tempmean45", "BO_tempmax45",
                                  "BO_salmean85", "BO_salmin85", 
                                  "BO_tempmean85", "BO_tempmax85")
#################

## compare
##########
# bio-oracle data
bo_env <- merge(curr_environment2, future_environment2, 
                by = c("name", "pop", "lat", "long", "order"))

# seascape env from gea
seascape_env <- merge(curr_environment2, seascape_dat, by = c("pop", "lat", "long"))

# how different are these values?
seascape_env$BOgea_salmean_diff <- as.numeric(seascape_env$BO_salmean) - as.numeric(seascape_env$mean_sal_sea)
seascape_env$BOgea_salmin_diff <- as.numeric(seascape_env$BO_salmin) - as.numeric(seascape_env$min_sal_sea)
seascape_env$BOgea_tempmean_diff <- as.numeric(seascape_env$BO_tempmean) - as.numeric(seascape_env$mean_temp_sea)
seascape_env$BOgea_tempax_diff <- as.numeric(seascape_env$BO_tempmax) - as.numeric(seascape_env$max_temp_sea)

# average/sd of differences
diffs_BOgea <- seascape_env %>% summarise(avg_salmean_BOgea = mean(BOgea_salmean_diff, na.rm = T),
                                    sd_salmean_BOgea = sd(BOgea_salmean_diff, na.rm = T),
                                    avg_salmin_BOgea = mean(BOgea_salmin_diff, na.rm = T),
                                    sd_salmin_BOgea = sd(BOgea_salmin_diff, na.rm = T),
                                    avg_tempmean_BOgea = mean(BOgea_tempmean_diff, na.rm = T),
                                    sd_tempmean_BOgea = sd(BOgea_tempmean_diff, na.rm = T),
                                    avg_tempmax_BOgea = mean(BOgea_tempax_diff, na.rm = T),
                                    sd_tempmax_BOgea = sd(BOgea_tempax_diff, na.rm = T))

# compare exp site values
exp_dat_cop2
exp_dat2
seascape_env

# combine
sea_exp <- merge(seascape_env, exp_dat2[,c("pop", "mean_temp_c", "temp_4.5_c", "temp_8.5_c",
                                           "mean_sal_c", "sal_4.5_c", "sal_8.5_c")], by = c("pop"))
sea_exp$BOc_salmean_diff <- as.numeric(sea_exp$BO_salmean) - as.numeric(sea_exp$mean_sal_c)
sea_exp$BOc_tempmean_diff <- as.numeric(sea_exp$BO_tempmean) - as.numeric(sea_exp$mean_temp_c)

diffs_BOc <- sea_exp %>% summarise(avg_salmean_BOc = mean(BOc_salmean_diff, na.rm = T),
                                   sd_salmean_BOc = sd(BOc_salmean_diff, na.rm = T),
                                   avg_tempmean_BOc = mean(BOc_tempmean_diff, na.rm = T),
                                   sd_tempmean_BOc = sd(BOc_tempmean_diff, na.rm = T))

sea_exp_cop <- merge(seascape_env, exp_dat_cop2, by = c("pop"))
sea_exp_cop$BOcop_salmean_diff <- as.numeric(sea_exp_cop$BO_salmean) - as.numeric(sea_exp_cop$median_sal_cop)
sea_exp_cop$BOcop_salmin_diff <- as.numeric(sea_exp_cop$BO_salmin) - as.numeric(sea_exp_cop$min_sal_cop)
sea_exp_cop$BOcop_tempmean_diff <- as.numeric(sea_exp_cop$BO_tempmean) - as.numeric(sea_exp_cop$median_temp_cop)
sea_exp_cop$BOcop_tempmax_diff <- as.numeric(sea_exp_cop$BO_tempmax) - as.numeric(sea_exp_cop$max_temp_cop)


diffs_BOcop <- sea_exp_cop %>% summarise(avg_salmean_BOcop = mean(BOcop_salmean_diff, na.rm = T),
                                          sd_salmean_BOcop = sd(BOcop_salmean_diff, na.rm = T),
                                          avg_salmin_BOcop = mean(BOcop_salmin_diff, na.rm = T),
                                          sd_salmin_BOcop = sd(BOcop_salmin_diff, na.rm = T),
                                          avg_tempmean_BOcop = mean(BOcop_tempmean_diff, na.rm = T),
                                          sd_tempmean_BOcop = sd(BOcop_tempmean_diff, na.rm = T),
                                          avg_tempmax_BOcop = mean(BOcop_tempmax_diff, na.rm = T),
                                          sd_tempmax_BOcop = sd(BOcop_tempmax_diff, na.rm = T))

# all together
diffs_BOgea
diffs_BOc
diffs_BOcop
##########
