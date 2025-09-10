###analysis MJ

setwd("~/Documents/projects/VR_2022/mesocosm experiment/HOBOloggers_GOeelgrass2024")


########
#first reduce field logger to every 15min

# Load required packages
library(readxl)    # for reading Excel
library(dplyr)     # for data manipulation
library(lubridate) # for working with dates and times

# Step 1: Load the Excel file
# (Update the path if needed)
data <- read_excel("Gasholmen_21444187_GOeelgrass2024.xlsx")

library(dplyr)
library(lubridate)

# First, check the actual column names
names(data)

# Find the column with "Temperature" and rename it
names(data)[grep("Temperature", names(data))] <- "Temp"
names(data)[grep("Date-Time", names(data))] <- "DateTime"
names(data)[grep("Light", names(data))] <- "Lux"

# Now filter for 15-minute intervals
data_filtered <- data %>%
  filter(minute(DateTime) %in% c(0, 15, 30, 45))

# Preview result
head(data_filtered)
write.csv(data_filtered, "Gasholmen_15min_filtered.csv", row.names = FALSE)


########
#read file with all loggers combined

data <- read_excel("Hobo_loggers_combined.xlsx")
#this file is lightly cleaned so far for obvious temperature spikes

#tidy up
library(dplyr)
library(stringr)

# Clean column names
names(data) <- names(data) %>%
  str_replace_all(" \\(.*?\\)", "") %>%                      # Remove units like "(lux)", "(??C)"
  str_replace_all("(?i)L*light+", "light") %>%               # Fix typos: Llight, llight, ight ??? light
  str_replace_all("[^A-Za-z0-9_]", "_") %>%                  # Replace non-alphanumeric with underscores
  str_replace_all("_+", "_") %>%                             # Collapse multiple underscores
  str_replace_all("(_L)?_?light$", "_light") %>%             # Fix suffixes like _L_light or _light
  str_replace_all("_temperature$", "_temperature") %>%       # Ensure proper _temperature suffix
  str_replace_all("_$", "")                                  # Remove trailing underscores

# Rename first column manually if it's still unnamed or has odd name
names(data)[1] <- "Row"

# Check result
names(data)
#View(data)

######
# calculate offset based to calibration
# Define calibration time
calib_time <- as.POSIXct("2024-07-01 18:45:00", tz = "UTC")

# Extract only the temperature columns at that moment
calib_row <- data %>%
  filter(Date_Time == calib_time) %>%
  dplyr::select(ends_with("_Temperature")) %>%
  pivot_longer(cols = everything(), names_to = "Sensor", values_to = "Temp")

# Calculate mean and offset
calib_offsets <- calib_row %>%
  mutate(
    Mean_Temp = mean(Temp, na.rm = TRUE),
    Offset = Mean_Temp - Temp  # Offset needed to bring each sensor to the mean
  ) %>%
  dplyr::select(Sensor, Offset)


#visualise offset distribution
library(ggplot2)
ggplot(calib_offsets, aes(x = Offset)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
  labs(title = "Distribution of Sensor Offsets",
       x = "Offset (??C)",
       y = "Count of Sensors") +
  theme_minimal()

#max offset was 0.32C so no adjustments were made



#######
#plot

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

data_long <- data %>%
  pivot_longer(
    cols = matches("_Temperature$"),  # only temperature columns
    names_to = "Sensor",
    values_to = "Value"
  ) %>%
  mutate(
    Group = str_extract(Sensor, "^[A-Z]\\d"),  # e.g. A1, B2
    Section = str_extract(Sensor, "^[A-Z]")    # e.g. A, B, C
  )

library(tidyverse)
library(lubridate)

# Convert Date_Time to POSIXct if not already
data$Date_Time <- as.POSIXct(data$Date_Time)

# Tanks to plot (from your naming pattern, e.g. A1, A2, ..., B1, B2, ..., C1, etc.)
tanks <- unique(sub("([A-C]\\d).*", "\\1", names(data)[grep("_Temperature", names(data))]))

# Function to gather and plot temperatures for one tank
plot_tank <- function(tank) {
  cols <- grep(paste0("^", tank), names(data), value = TRUE)
  cols <- cols[grepl("_Temperature", cols)]
  
  temp_data <- data %>%
    dplyr::select(Date_Time, all_of(cols)) %>%  # use dplyr::select explicitly
    pivot_longer(cols = all_of(cols), names_to = "Sensor", values_to = "Temperature") %>%
    mutate(Sensor = factor(Sensor))
  
  p <- ggplot(temp_data, aes(x = Date_Time, y = Temperature, color = Sensor)) +
    geom_line() +
    labs(title = paste("Temperature over time for tank", tank),
         x = "Date Time", y = "Temperature (??C)") +
    theme_minimal()
  
  print(p)
}

for (tank in tanks) {
      plot_tank(tank)
}


# Just in case, convert to tibble
data <- as_tibble(data)

start_date <- as.POSIXct("2024-07-06")
end_date   <- as.POSIXct("2024-08-06")

combined_temp <- data %>%
  filter(Date_Time >= start_date & Date_Time <= end_date) %>%    # Filter for experiment window
  dplyr::select(Date_Time, dplyr::ends_with("_Temperature")) %>%
  pivot_longer(cols = -Date_Time, names_to = "Sensor", values_to = "Temperature") %>%
  mutate(
    Tank = sub("([A-C]\\d).*", "\\1", Sensor)
  )

cool_tanks <- c("A1", "A2", "C1", "C2", "C3", "C4")
cool_colors <- c(
  "#084594",  # A1 - dark blue
  "#2171b5",  # A2 - slightly lighter
  "#2b8cbe",  # C1 - blue-cyan
  "#3690c0",  # C2
  "#4eb3d3",  # C3
  "#6baed6"   # C4
)
# Define blueish tones for cool tanks
# Get all tanks from data
all_tanks <- sort(unique(combined_temp$Tank))
# Warm tanks = all others
warm_tanks <- setdiff(all_tanks, cool_tanks)
# Warm reddish tones
warm_colors <- c(
  "#b10026", "#e31a1c", "#fc4e2a", "#fd8d3c", "#feb24c", "#fed976"
)[1:length(warm_tanks)]

# Combine them in the right order
tank_colors <- c(setNames(cool_colors, cool_tanks), setNames(warm_colors, warm_tanks))


ggplot(combined_temp, aes(x = Date_Time, y = Temperature, color = Tank)) +
  geom_line(alpha = 0.8) +
  labs(title = "Temperature over time for all tanks (6 July ??? 6 August 2024)",
       x = "Date Time",
       y = "Temperature (??C)",
       color = "Tank") +
  scale_color_manual(values = tank_colors) +
  scale_x_datetime(
    date_breaks = "1 day",
    date_labels = "%d %b",
    expand = c(0, 0)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )




#####
#get some stats
# Ensure data is in correct format
data <- as_tibble(data)

###
# for experimental period only
# Define the date range
data$Date_Time <- as.POSIXct(data$Date_Time)

# Define the date range
start_date <- as.POSIXct("2024-07-06")
end_date   <- as.POSIXct("2024-08-06")

# Filter and summarize
tank_temp_stats_filtered <- data %>%
  filter(Date_Time >= start_date & Date_Time <= end_date) %>%
  dplyr::select(Date_Time, dplyr::ends_with("_Temperature")) %>%
  pivot_longer(cols = -Date_Time, names_to = "Sensor", values_to = "Temperature") %>%
  mutate(Tank = sub("([A-C]\\d).*", "\\1", Sensor)) %>%
  filter(!is.na(Temperature)) %>%  # remove rows with NA values
  group_by(Tank) %>%
  summarise(
    Mean_Temp = mean(Temperature),
    Max_Temp = max(Temperature),
    Min_Temp = min(Temperature),
    P90_Temp = quantile(Temperature, probs = 0.9, na.rm = TRUE),
   .groups = "drop"
  )

# View the results
print(tank_temp_stats_filtered)
#  Tank                  Mean_Temp Max_Temp Min_Temp
#<chr>                     <dbl>    <dbl>    <dbl>
#  1 A1                         20.0     23.9     16.2
#2 A2                         20.0     24.1     16.0
#9 C1                         20.1     25.9     16.7
#10 C2                         20.1     24.8     15.5
#11 C3                         20.1     23.8     16.3
#12 C4                         20.3     24.9     16.3
#3 A3                         26.1     32.3     20.2
#4 A4                         26.1     31.4     20.1
#5 B1                         26.3     32.1     21.0
#6 B2                         26.2     31.7     21.4
#7 B3                         26.1     32.0     21.2
#8 B4                         25.8     31.4     19.9
#13 Gasholmen_temperature      19.4     23.3     15.8
