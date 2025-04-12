#### Try to get nc files in right format

# setup
setwd("~/BalticEelgrass/data")

library(raster)
library(CopernicusMarine)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(ncdf4)

# tutorial nc files
# https://rpubs.com/boyerag/297592

# read in and check
nc <- nc_open("EnvDat/cmems_mod_bal_phy_my_P1D-m_so-thetao_10.85E-19.35E_55.16N-60.57N_0.50-2.55m_2011-06-01-2021-09-30.nc")
print(nc)

# need to capture metadata information about the file in lat/lon/time/depth dimensions
lon <- ncvar_get(nc, "longitude")
lat <- ncvar_get(nc, "latitude", verbose = F)
t <- ncvar_get(nc, "time")
d <- ncvar_get(nc, "depth")

# grab the salinity and temp data we want
# recall "so" is salinity and "thetao" is (potential) temp
sal <- ncvar_get(nc, "so")
temp <- ncvar_get(nc, "thetao")

# need to associate the metadata
dim(sal) # 307  326    3 3775
dim(temp) # same dims as above
length(lon) # 307 
length(lat) # 326
dim(t) # 3775
dim(d) # 3

# deal with NAs
fillvalue <- ncatt_get(nc, "so", "_FillValue") # -999
fillvalue <- ncatt_get(nc, "thetao", "_FillValue") # -999
sal[sal == fillvalue$value] <- NA
temp[temp == fillvalue$value] <- NA

# associate time with dates
times <- data.frame(t = t, date = seq(as.Date("2011-06-01"), as.Date("2021-09-30"), by = "days"))
times$month <- str_split_fixed(times$date, "-", 3)[,2]

# get dates that fall into "summer"
summer <- c("06", "07", "08")
summer_dates <- times[times$month %in% summer,]
summer_indices <- which(t %in% summer_dates$t)
summer_sal <- sal[,,,summer_indices]
summer_temp <- temp[,,,summer_indices]

# start by collapsing depth so we only work with three dimensions - not working
#library(s2dv)
sal_means_d <- MeanDims(summer_sal, 3, na.rm = T)
temp_means_d <- MeanDims(summer_temp, 3, na.rm = T)

# determine which lat long points i need
env_org <- as.data.frame(read_excel("EnvDat/Camille_envdata_experiment.xlsx"))
colnames(env_org) <- gsub(" ", "_", colnames(env_org))

# do these points show up exactly in our df?
which(env_org$Longitude %in% lon) # no! need to make this into a 

longs_nc <- matrix(nrow = 8, ncol = 2)
for (i in 1:length(env_org$Longitude)) {
  number = env_org$Longitude[i]
  long = lon[which.min(abs(lon - number))]
  longs_nc[i,1] = env_org$Site_code[i]
  longs_nc[i,2] = long
}
longs_nc <- as.data.frame(longs_nc)
colnames(longs_nc) <- c("Site_code","Longitude_cop")

lats_nc <- matrix(nrow = 8, ncol = 2)
for (i in 1:length(env_org$Latitude)) {
  number = env_org$Latitude[i]
  lats = lat[which.min(abs(lat - number))]
  lats_nc[i,1] = env_org$Site_code[i]
  lats_nc[i,2] = lats
}
lats_nc <- as.data.frame(lats_nc)
colnames(lats_nc) <- c("Site_code","Latitude_cop")

# indices for these lats and longs
lat_indices <- which(lat %in% lats_nc$Latitude_cop)
long_indices <- which(lon %in% longs_nc$Longitude_cop)

# keep just these points
latlongCop <- merge(lats_nc, longs_nc, by = c("Site_code"))
summer_sal[long_indices, lat_indices,,]
summer_temp[long_indices, lat_indices,,]

# salinity expanded
expanddims_sal <- function(depth,time) {
  summer_sal_latlon_d_t <- matrix(nrow = 8, ncol = 4)
  for (i in 1:8) {
    summer_sal_latlon_d_t[i,1] = latlongCop$Site_code[i]
    summer_sal_latlon_d_t[i,2] = as.numeric(latlongCop$Latitude_cop[i])
    summer_sal_latlon_d_t[i,3] = as.numeric(latlongCop$Longitude_cop[i])
    summer_sal_latlon_d_t[i,4] = as.numeric(summer_sal[i,i,depth,time])
  }
  summer_sal_latlon_d_t <- as.data.frame(summer_sal_latlon_d_t)
  colnames(summer_sal_latlon_d_t) <- c("Site_code","lat","long","sal")
  return(summer_sal_latlon_d_t)
}
summer_sals_dt <- data.frame(Site_code = NULL, lat = NULL, long = NULL, sal = NULL, depth = NULL, time = NULL)
for (x in 1:length(d)) {
  for(y in 1:nrow(summer_dates)) {
    df = expanddims_sal(x,y)
    df$depth = x
    df$time = y
    summer_sals_dt <- rbind(summer_sals_dt,df)
  }
}
dim(summer_sals_dt)
head(summer_sals_dt)

# temperature expanded
expanddims_temp <- function(depth,time) {
  summer_temp_latlon_d_t <- matrix(nrow = 8, ncol = 4)
  for (i in 1:8) {
    summer_temp_latlon_d_t[i,1] = latlongCop$Site_code[i]
    summer_temp_latlon_d_t[i,2] = as.numeric(latlongCop$Latitude_cop[i])
    summer_temp_latlon_d_t[i,3] = as.numeric(latlongCop$Longitude_cop[i])
    summer_temp_latlon_d_t[i,4] = as.numeric(summer_temp[i,i,depth,time])
  }
  summer_temp_latlon_d_t <- as.data.frame(summer_temp_latlon_d_t)
  colnames(summer_temp_latlon_d_t) <- c("Site_code","lat","long","temp")
  return(summer_temp_latlon_d_t)
}
summer_temps_dt <- data.frame(Site_code = NULL, lat = NULL, long = NULL, temp = NULL, depth = NULL, time = NULL)
for (x in 1:length(d)) {
  for(y in 1:nrow(summer_dates)) {
    df = expanddims_temp(x,y)
    df$depth = x
    df$time = y
    summer_temps_dt <- rbind(summer_temps_dt, df)
  }
}
dim(summer_temps_dt)
head(summer_temps_dt)

# now use to calculate the avg temp by depth at each time point for each site
summer_temps_dt$time <- as.factor(summer_temps_dt$time)
summer_temps_dt$temp <- as.numeric(summer_temps_dt$temp)

sites_times_temps <- summer_temps_dt %>% group_by(time,Site_code) %>%
  summarise(mean_temp = mean(temp),
            min_temp = min(temp),
            max_temp = max(temp),
            med_temp = median(temp),
            sd_temp = sd(temp))
sites_times_temps
dim(sites_times_temps)

# and for salinity
summer_sals_dt$time <- as.factor(summer_sals_dt$time)
summer_sals_dt$sal <- as.numeric(summer_sals_dt$sal)

sites_times_sals <- summer_sals_dt %>% group_by(time,Site_code) %>%
  summarise(mean_sal = mean(sal),
            min_sal = min(sal),
            max_sal = max(sal),
            med_sal = median(sal),
            sd_sal = sd(sal))
sites_times_sals
dim(sites_times_sals)

# stats for sites across all depths and times
sites_sals <- summer_sals_dt %>% group_by(Site_code) %>%
  summarise(mean_sal = mean(sal),
            min_sal = min(sal),
            max_sal = max(sal),
            med_sal = median(sal),
            sd_sal = sd(sal))
sites_sals
dim(sites_sals)

sites_temps <- summer_temps_dt %>% group_by(Site_code) %>%
  summarise(mean_temp = mean(temp),
            min_temp = min(temp),
            max_temp = max(temp),
            med_temp = median(temp),
            sd_temp = sd(temp))
sites_temps
dim(sites_temps)

# put it together and save
eelgrass_sites_copernicus <- merge(sites_sals, sites_temps, by = c("Site_code"))
write.csv(eelgrass_sites_copernicus, "EnvDat/eelgrass_exppops_copernicus.csv")
