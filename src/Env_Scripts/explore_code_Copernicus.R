install.packages("R.matlab")
library(R.matlab)

###salinity
sal <- readMat("data/EnvDat/SalinityDaily.mat")
#rows are the 39 sites in the same order that I got from you, the columns are the data for the 92 days with 1 June in the first column, and then there are 11 slices, one for each year, where the first slice is year 2011.

head(sal)
summary(sal)
str(sal)

mean(sal$SalinityDaily)
9.16166

# Extracting the slice for the first year (2011)
#year_2011_slice <- data$SalinityDaily[, , 1]

# Calculate the mean for each site over all years, regardless of the day
mean_for_sites <- apply(sal$SalinityDaily, MARGIN = 1, mean)
#margin=1 gives mean for each site regardless of day
min_sal_for_sites <- apply(sal$SalinityDaily, MARGIN = 1, min)
sal_percentile_10_for_sites <- apply(sal$SalinityDaily, MARGIN = 1, function(x) quantile(x, 0.1))
sal_percentile_1_for_sites <- apply(sal$SalinityDaily, MARGIN = 1, function(x) quantile(x, 0.01))
sal_percentile_90_for_sites <- apply(sal$SalinityDaily, MARGIN = 1, function(x) quantile(x, 0.9))

###temp
temp <- readMat("data/EnvDat/TempDaily.mat")
#rows are the 39 sites in the same order that I got from you, the columns are the data for the 92 days with 1 June in the first column, and then there are 11 slices, one for each year, where the first slice is year 2011.

head(temp)
summary(temp)
str(temp)


# Calculate the mean for each site over all years, regardless of the day
T_mean_for_sites <- apply(temp$TempDaily, MARGIN = 1, mean)
#margin=1 gives mean for each site regardless of day
T_max_for_sites <- apply(temp$TempDaily, MARGIN = 1, max)

# Calculate the 90th percentile for each site
T_percentile_90_for_sites <- apply(temp$TempDaily, MARGIN = 1, function(x) quantile(x, 0.9))
T_percentile_99_for_sites <- apply(temp$TempDaily, MARGIN = 1, function(x) quantile(x, 0.99))
T_percentile_10_for_sites <- apply(temp$TempDaily, MARGIN = 1, function(x) quantile(x, 0.1))


str(data)

# Assuming your data is a 3D array with rows as sites, columns as days, and the third dimension as years
# Create a function to find the maximum temperature for at least 5 consecutive days over all years
max_temp_for_5_days <- function(site_data) {
  max_temp = -Inf
  current_streak = 0
  
  for (temp in site_data) {
    if (temp > max_temp) {
      max_temp = temp
    } else {
      current_streak = 0
    }
    
    current_streak = current_streak + 1
    
    if (current_streak >= 5) {
      return(max_temp)
    }
  }
  
  return(NA)  # Return NA if no streak of 5 or more days is found
}

# Apply the function to each site
max_temp_for_5_days_for_sites <- apply(temp$TempDaily, MARGIN = 1, max_temp_for_5_days)



