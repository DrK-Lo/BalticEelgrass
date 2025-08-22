## For working with netcdf data

library(dplyr)
library(ncdf4)
library(raster)
library(terra)

# experimental site data
exp_pops <- read.csv("data/experiment/eelgrass_exp_sites.csv")

# set up objects
sal_complete <- c()
temp_complete <- c()

# daily data for summer months 2011-2021
for (year in c(2011:2021)) {
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

sal_complete
temp_complete
