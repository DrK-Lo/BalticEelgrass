## For working with netcdf data

library(dplyr)
library(tidyr)
library(ncdf4)
library(raster)
library(terra)

# experimental site data
exp_pops <- read.csv("data/experiment/eelgrass_exp_sites.csv")

# set up objects
sal_complete <- matrix(nrow = 8)
temp_complete <- matrix(nrow = 8)

# daily data for summer months 2011-2021
for (year in c(2011:2021)) {
  print(paste0("Starting year ",year))
  
  name = paste0("data/largeDat/cmems_mod_bal_phy_my_P1D-m_", year, ".nc")
  nc_year <- nc_open(name)
  
  # get attributes
  lon <- ncvar_get(nc_year, "longitude")
  lat <- ncvar_get(nc_year, "latitude", verbose = F)
  t <- ncvar_get(nc_year, "time")
  tunits <- ncatt_get(nc_year,"time","units")
  d <- ncvar_get(nc_year, "depth")
  
  sal_array <- ncvar_get(nc_year, "so")
  sal_lname <- ncatt_get(nc_year, "so", "long_name")
  salunits <- ncatt_get(nc_year, "so", "units")
  fillvalue <- ncatt_get(nc_year, "so", "_FillValue")
  
  temp_array <- ncvar_get(nc_year, "thetao")
  temp_lname <- ncatt_get(nc_year, "thetao", "long_name")
  tempunits <- ncatt_get(nc_year, "thetao", "units")
  fillvalue2 <- ncatt_get(nc_year, "thetao", "_FillValue")
  
  # get global attributes
  title <- ncatt_get(nc_year,0,"title")
  institution <- ncatt_get(nc_year,0,"institution")
  datasource <- ncatt_get(nc_year,0,"source")
  references <- ncatt_get(nc_year,0,"references")
  history <- ncatt_get(nc_year,0,"history")
  Conventions <- ncatt_get(nc_year,0,"Conventions")
  
  # close file
  nc_close(nc_year) 
  
  # aggregate data
  aggregated_sal <- apply(sal_array, c(1, 2, 4), function (x) {mean(x, na.rm = T)})
  aggregated_temp <- apply(temp_array, c(1, 2, 4), function (x) {mean(x, na.rm = T)})
  
  # prep for raster
  aggregated_sal_ls <- list()
  for (i in 1:dim(aggregated_sal)[3]) { # Iterate through the time/band dimension
    r <- raster(t(aggregated_sal[,,i]), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat))
    crs(r) <- "+proj=longlat +datum=WGS84"
    aggregated_sal_ls[[i]] <- r
  }
  
  aggregated_temp_ls <- list()
  for (i in 1:dim(aggregated_temp)[3]) { # Iterate through the time/band dimension
    r <- raster(t(aggregated_temp[,,i]), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat))
    crs(r) <- "+proj=longlat +datum=WGS84"
    aggregated_temp_ls[[i]] <- r
  }
  
  # build rasterstack
  salrast <- stack(aggregated_sal_ls)
  temprast <- stack(aggregated_temp_ls)
  
  # built flipped on y axis, fix this
  salrastfl <- flip(salrast, direction = "y")
  temprastfl <- flip(temprast, direction = "y")
  
  # we also need to make sure our sites won't have missing data by imputing raster with nearest neighbors
  salrastimp <- focal(rast(salrastfl), w = 3, fun = "modal", na.policy = "only")
  temprastimp <- focal(rast(temprastfl), w = 3, fun = "modal", na.policy = "only")
  
  # now can extract necessary data
  coords_exp <- exp_pops[,c("long","lat")]
  sal_sites <- raster::extract(salrastimp, coords_exp)
  temp_sites <- raster::extract(temprastimp, coords_exp)
  
  sal_complete <- cbind(sal_complete, sal_sites)
  temp_complete <- cbind(temp_complete, temp_sites)
}

sal_complete_df <- as.data.frame(sal_complete)
temp_complete_df <- as.data.frame(temp_complete)

# add sites
sal_complete_df$sal_complete <- temp_complete_df$temp_complete <- exp_pops$site_abbrev
saldf <- sal_complete_df[,!names(sal_complete_df) %in% "ID"]
tempdf <- temp_complete_df[,!names(temp_complete_df) %in% "ID"]

# add dates
colnamesdf <- c(seq(as.Date("2011-06-01"), as.Date("2011-08-31"), by="days"), 
                seq(as.Date("2012-06-01"), as.Date("2012-08-31"), by="days"),
                seq(as.Date("2013-06-01"), as.Date("2013-08-31"), by="days"),
                seq(as.Date("2014-06-01"), as.Date("2014-08-31"), by="days"),
                seq(as.Date("2015-06-01"), as.Date("2015-08-31"), by="days"),
                seq(as.Date("2016-06-01"), as.Date("2016-08-31"), by="days"),
                seq(as.Date("2017-06-01"), as.Date("2017-08-31"), by="days"),
                seq(as.Date("2018-06-01"), as.Date("2018-08-31"), by="days"),
                seq(as.Date("2019-06-01"), as.Date("2019-08-31"), by="days"),
                seq(as.Date("2020-06-01"), as.Date("2020-08-31"), by="days"),
                seq(as.Date("2021-06-01"), as.Date("2021-08-31"), by="days"))
colnamesdf <- as.character(colnamesdf)
colnamesdf <- c("site_name", colnamesdf)

colnames(saldf) <- colnames(tempdf) <- colnamesdf

# save these objects
write.csv(saldf, 
          paste0("data/EnvDat/salinity_expSites_copernicus_", Sys.Date(),".csv"), 
          row.names = F)
write.csv(tempdf, 
          paste0("data/EnvDat/temperature_expSites_copernicus_", Sys.Date(),".csv"), 
          row.names = F)

# overall stats
saldf_sum <- data.frame(sal_mean = apply(saldf[,-1],1,mean), 
                        sal_min = apply(saldf[,-1],1,min),
                        sal_max = apply(saldf[,-1],1,max),
                        sal_quant01 = apply(saldf[,-1],1,function (x) {quantile(x,0.01)}),
                        sal_quant99 = apply(saldf[,-1],1,function (x) {quantile(x,0.99)})) 
tempdf_sum <- data.frame(temp_mean = apply(tempdf[,-1],1,mean), 
                        temp_min = apply(tempdf[,-1],1,min),
                        temp_max = apply(tempdf[,-1],1,max),
                        temp_quant01 = apply(tempdf[,-1],1,function (x) {quantile(x,0.01)}),
                        temp_quant99 = apply(tempdf[,-1],1,function (x) {quantile(x,0.99)})) 

# put together df
site_abbrev <- rownames(saldf_sum)
copernicus_df <- cbind(site_abbrev, saldf_sum, tempdf_sum)
write.csv(copernicus_df, paste0("data/EnvDat/complete_expSites_copernicus", Sys.Date(), ".csv"))

# expand to long format df to calculate monthly/annual trends
sal_expanded <- saldf %>% 
  pivot_longer(
    cols = "2011-06-01":"2021-08-31", 
    names_to = "date",
    values_to = "salinity"
  ) %>% as.data.frame()
temp_expanded <- tempdf %>% 
  pivot_longer(
    cols = "2011-06-01":"2021-08-31", 
    names_to = "date",
    values_to = "temperature"
  ) %>% as.data.frame()

# split date info
dates <- str_split_fixed(sal_expanded$date, "-", 3)
colnames(dates) <- c("year", "month", "day")
copernicus_expanded <- cbind(sal_expanded, temp_expanded$temperature, dates)
colnames(copernicus_expanded) <- c("site_name","date","sal","temp","year","month","day")
copernicus_expanded$date <- as.Date(copernicus_expanded$date)

# annual and monthly trends
annual_trends <- copernicus_expanded %>% group_by(site_name,year) %>%
  summarise(mean_temp = mean(temp),
            min_temp = min(temp),
            max_temp = max(temp),
            med_temp = median(temp),
            sd_temp = sd(temp),
            mean_sal = mean(sal),
            min_sal = min(sal),
            max_sal = max(sal),
            med_sal = median(sal),
            sd_sal = sd(sal),
            .groups = "drop")
annual_trends

monthly_trends <- copernicus_expanded %>% group_by(site_name,month) %>%
  summarise(mean_temp = mean(temp),
            min_temp = min(temp),
            max_temp = max(temp),
            med_temp = median(temp),
            sd_temp = sd(temp),
            mean_sal = mean(sal),
            min_sal = min(sal),
            max_sal = max(sal),
            med_sal = median(sal),
            sd_sal = sd(sal),
            .groups = "drop")
monthly_trends

overall_trends <- copernicus_expanded %>% group_by(site_name,month,year) %>%
  summarise(mean_temp = mean(temp),
            min_temp = min(temp),
            max_temp = max(temp),
            med_temp = median(temp),
            sd_temp = sd(temp),
            mean_sal = mean(sal),
            min_sal = min(sal),
            max_sal = max(sal),
            med_sal = median(sal),
            sd_sal = sd(sal),
            .groups = "drop")
overall_trends

# plot the data by group
ggplot(copernicus_expanded, 
       aes(x = date, y = sal, fill = site_name, col = site_name)) +
  geom_smooth() +
  labs(title = "Salinity at Experimental Sites 2011-2021",
       x = "Date", y = "Salinity", color = "Site", fill = "Site") +
  theme_bw()

ggplot(copernicus_expanded, 
       aes(x = date, y = temp, fill = site_name, col = site_name)) +
  geom_smooth() +
  labs(title = "Temperature at Experimental Sites 2011-2021",
       x = "Date", y = "Temperature", color = "Site", fill = "Site") +
  theme_bw()
