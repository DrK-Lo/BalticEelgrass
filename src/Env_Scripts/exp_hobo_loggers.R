#################################
## Hobo logger data GOeelgrass ##
#################################

## setup
########
# set wd
setwd("~/Documents/Github/BalticEelgrass/")

# libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(lubridate)

# data
setup <- as.data.frame(read.csv("data/experiment/GOEEL-eelgrass-Sweden-exp - Setup.csv"))
A1a <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/A1a.xlsx"))[,-1]
A1b <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/A1b.xlsx"))[,-1]
A1c <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/A1c.xlsx"))[,-1]
A2a <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/A2a.xlsx"))[,-1]
A2b <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/A2b.xlsx"))[,-1]
A2c <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/A2c.xlsx"))[,-1]
A3a <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/A3a.xlsx"))[,-1]
A3b <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/A3b.xlsx"))[,-1]
A3c <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/A3c.xlsx"))[,-1]
A4a <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/A4a.xlsx"))[,-1]
A4b <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/A4b.xlsx"))[,-1]
A4c <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/A4c.xlsx"))[,-1]
B1a <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/B1a.xlsx"))[,-1]
B1b <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/B1b.xlsx"))[,-1]
B1c <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/B1c.xlsx"))[,-1]
B2a <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/B2a.xlsx"))[,-1]
B2b <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/B2b.xlsx"))[,-1]
B2c <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/B2c.xlsx"))[,-1]
B3a <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/B3a.xlsx"))[,-1]
B3b <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/B3b.xlsx"))[,-1]
B3c <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/B3c.xlsx"))[,-1]
B4a <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/B4a.xlsx"))[,-1]
B4b <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/B4b.xlsx"))[,-1]
B4c <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/B4c.xlsx"))[,-1]
C1a <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/C1a.xlsx"))[,-1]
C1b <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/C1b.xlsx"))[,-1]
C1c <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/C1c.xlsx"))[,-1]
C2a <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/C2a.xlsx"))[,-1]
C2b <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/C2b.xlsx"))[,-1]
C2c <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/C2c.xlsx"))[,-1]
C3a <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/C3a.xlsx"))[,-1]
C3b <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/C3b.xlsx"))[,-1]
C3c <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/C3c.xlsx"))[,-1]
C4a <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/C4a.xlsx"))[,-1]
C4b <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/C4b.xlsx"))[,-1]
C4c <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/C4c.xlsx"))[,-1]
gas <- as.data.frame(read_excel("data/EnvDat/HOBOloggers_GOeelgrass2024/gas2024.xlsx"))[,-1]

# also salinity data
bags <- read.csv("data/EnvDat/GOEEL-eelgrass-Sweden-exp - Bags.csv")
daily_checks <- read.csv("data/EnvDat/GOEEL-eelgrass-Sweden-exp - DailyChecks.csv")
########

## temp data pre-cleaning
#########################
# change colnames to be more usable
colnames(gas) <- colnames(A1a) <- colnames(A1b) <- colnames(A1c) <- colnames(A2a) <- 
  colnames(A2b) <- colnames(A2c) <- colnames(A3a) <- colnames(A3b) <- colnames(A3c) <- 
  colnames(A4a) <- colnames(A4b) <- colnames(A4c) <- colnames(B1a) <- colnames(B1b) <- 
  colnames(B1c) <- colnames(B2a) <- colnames(B2b) <- colnames(B2c) <- colnames(B3a) <- 
  colnames(B3b) <- colnames(B3c) <- colnames(B4a) <- colnames(B4b) <- colnames(B4c) <- 
  colnames(C1a) <- colnames(C1b) <- colnames(C1c) <- colnames(C2a) <- colnames(C2b) <- 
  colnames(C2c) <- colnames(C3a) <- colnames(C3b) <- colnames(C3c) <- colnames(C4a) <- 
  colnames(C4b) <- colnames(C4c) <- c("datetime","temp","light")

# filter by time for all df
# time range
start <- ymd_hms("2024-07-03 09:00:00")
end <- ymd_hms("2024-08-07 09:00:00")

# List of data frame names
log_names <- c("A1a", "A1b", "A1c", "A2a", "A2b", "A2c", 
              "A3a", "A3b", "A3c", "A4a", "A4b", "A4c",
              "B1a", "B1b", "B1c", "B2a", "B2b", "B2c",
              "B3a", "B3b", "B3c", "B4a", "B4b", "B4c",
              "C1a", "C1b", "C1c", "C2a", "C2b", "C2c", 
              "C3a", "C3b", "C3c", "C4a", "C4b", "C4c",
              "gas")

# filtering
for (logger in log_names) {
  tank <- substr(logger, 1, 2)
  assign(paste0(logger, "_trim"), 
         get(logger) %>% 
           filter(between(datetime, start, end)) %>% # filter times
           mutate(tank_ID = tank, logger = logger)) # add metadata columns
}


all_tanks <- rbind(A1a_trim,A1b_trim,A1c_trim,A2a_trim,A2b_trim,A2c_trim,A3a_trim,A3b_trim,
                   A3c_trim,A4a_trim,A4b_trim,A4c_trim,B1a_trim,B1b_trim,B1c_trim,B2a_trim,
                   B2b_trim,B2c_trim,B3a_trim,B3b_trim,B3c_trim,B4a_trim,B4b_trim,B4c_trim,
                   C1a_trim,C1b_trim,C1c_trim,C2a_trim,C2b_trim,C2c_trim,C3a_trim,C3b_trim,
                   C3c_trim,C4a_trim,C4b_trim,C4c_trim)
#########################

## analyze temp trends
######################
# summarize by tank
tank_temps <- all_tanks %>% group_by(tank_ID) %>%
  summarise(mean_temp = mean(temp),
            min_temp = min(temp),
            max_temp = max(temp),
            sd_temp = sd(temp))

# notice that max temp of logger C3c spiked in the beginning of the experiment in a way that does not fit with the other two loggers, trim these data then re-do this bit
# remove those funky rows
all_tanks_trim <- all_tanks[-(108035:108044),]

# re-summarize by tank
tank_temps <- all_tanks_trim %>% group_by(tank_ID) %>%
  summarise(mean_temp = mean(temp),
            median_temp = median(temp),
            min_temp = min(temp),
            max_temp = max(temp),
            sd_temp = sd(temp))

# add a column for heat/ctrl
control <- c("A1","A2","C1","C2","C3","C4")
all_tanks_trim$trt <- ifelse(all_tanks_trim$tank %in% control, "ctrl", "heat")

# keep this raw data
write.csv(all_tanks_trim, "data/EnvDat/all_tanks_logger_dat.csv")

# summarize by trt
trt_temps <- all_tanks_trim %>% group_by(trt) %>%
  summarise(mean_temp = mean(temp),
            median_temp = median(temp),
            min_temp = min(temp),
            max_temp = max(temp),
            sd_temp = sd(temp))

# save these data
write.csv(tank_temps, paste0("data/EnvDat/temp_tanks_summary_", Sys.Date() ,".csv"), row.names = F)
write.csv(trt_temps, paste0("data/EnvDat/temp_trts_summary_", Sys.Date() ,".csv"), row.names = F)

# compare to meadow
gas_temps <- data.frame(meadow = "gas", min_temp = min(gas_trim$temp), max_temp = max(gas_trim$temp), avg_temp = mean(gas_trim$temp), sd_temp = sd(gas_trim$temp))

# we have a small time series where the logger experienced temps 10 degrees above any other recorded and light was high enough to suggest it was out of the water, remove these
gas_temps <- data.frame(meadow = "gas", min_temp = min(gas_trim[which(gas_trim$temp < 35),]$temp), 
                        max_temp = max(gas_trim[which(gas_trim$temp < 35),]$temp), 
                        avg_temp = mean(gas_trim[which(gas_trim$temp < 35),]$temp), 
                        sd_temp = sd(gas_trim[which(gas_trim$temp < 35),]$temp)) # looks good!

# remove practice
setup_trim <- setup[!setup$bagKey == "F054C074", c("setupKey","setupLabel","bagKey")]

# get labels
exp_labels <- str_split_fixed(setup_trim$setupLabel, "07/", 2)[,1]

# split this into a larger dataframe
exp_df <- as.data.frame(str_split_fixed(exp_labels, "_", 5))
colnames(exp_df) <- c("bagnum","tank","pop","genet","trt")
exp_df$ind <- paste0(exp_df$pop,"_",exp_df$genet)

# now format temp data by treatment
trt_temp_dat <- data.frame(trt = c("TempControl-21psu", "TempControl-7psu", "TempWarm-16psu", "TempWarm-5psu"), 
                          max_temp = c(trt_temps[trt_temps$trt == "ctrl",]$max_temp, trt_temps[trt_temps$trt == "ctrl",]$max_temp, trt_temps[trt_temps$trt == "heat",]$max_temp, trt_temps[trt_temps$trt == "heat",]$max_temp), 
                          avg_temp = c(trt_temps[trt_temps$trt == "ctrl",]$mean_temp, trt_temps[trt_temps$trt == "ctrl",]$mean_temp, trt_temps[trt_temps$trt == "heat",]$mean_temp, trt_temps[trt_temps$trt == "heat",]$mean_temp), 
                          med_temp = c(trt_temps[trt_temps$trt == "ctrl",]$median_temp, trt_temps[trt_temps$trt == "ctrl",]$median_temp, trt_temps[trt_temps$trt == "heat",]$median_temp, trt_temps[trt_temps$trt == "heat",]$median_temp))

# and combine
exp_ind_env <- merge(exp_df, trt_env_dat, by = c("trt"))

# order this by bagnum
exp_env_ord <- exp_ind_env[order(as.numeric(exp_ind_env$bagnum)),]

# also need pop/trt level
pops <- unique(levels(as.factor(exp_env_ord$pop)))
trt <- unique(levels(as.factor(exp_env_ord$trt)))
pop_trt_df <- expand.grid(pops, trt)
colnames(pop_trt_df) <- c("pop","trt")
pop_trt_env <- merge(pop_trt_df, trt_env_dat, by = c("trt"))
######################

## salinity
###########
# check data
summary(daily_checks)

# there are some massive outliers, remove these
daily_checks_clean <- daily_checks %>%
  filter(dailyTemp < 40 & dailySalinity <30)
summary(daily_checks_clean)

# change bag colnames for consistency
colnames(bags)[colnames(bags) == "PopID"] <- "pop"

# put together the salinity data so we can look at trends by tank/treatment
sal_dat <- merge(bags, daily_checks_clean, by = c("bagKey"))

# expand out the date time info from daily checks
library(lubridate)
sal_dat$date <- str_split_fixed(sal_dat$dailyTimestamp, " ", 2)[,1]
sal_dat$date <- mdy(sal_dat$date)
sal_dat$time <- str_split_fixed(sal_dat$dailyTimestamp, " ", 2)[,2]
sal_dat$poptrt <- sal_dat$trt

# it looks like some bag ids got mixed up between current W and future E
sal_dat_trim <- sal_dat[!(sal_dat$trt == "TempControl-21psu" & sal_dat$dailySalinity < 18),]
sal_dat_clean <- sal_dat_trim[!(sal_dat_trim$trt == "TempWarm-5psu" & sal_dat_trim$dailySalinity > 8),]

# sal summary stats
bag_sals <- sal_dat_clean %>% group_by(bagKey) %>%
  summarise(mean_sal = mean(dailySalinity),
            median_sal = median(dailySalinity),
            min_sal = min(dailySalinity),
            max_sal = max(dailySalinity),
            sd_sal = sd(dailySalinity)) %>%
  as.data.frame()

# merge
bag_sals_df <- merge(bags, bag_sals, by = c("bagKey"))

# and treatment level stats
trt_sals <- sal_dat_clean %>% group_by(trt) %>%
  summarise(mean_sal = mean(dailySalinity),
            median_sal = median(dailySalinity),
            min_sal = min(dailySalinity),
            max_sal = max(dailySalinity),
            sd_sal = sd(dailySalinity)) %>%
  as.data.frame()
###########

## complete exp env dat
#######################
# combine exp env data
exp_env_df <- merge(bag_sals_df, tank_temps, by = c("tank_ID"))

# get the xy position for each bag in 384 cell grid as well
# x goes 1-24, y goes 1-16
# try using rep and seq to make this work
rep(c(rep(1,4), rep(2,4), rep(3,4), rep(4,4)), 4)

coords_bags <- data.frame(x_coord = c(rep(c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), 
                                            rep(5,4), rep(6,4), rep(7,4), rep(8,4)), 
                                          4),
                                      rep(c(rep(9,4), rep(10,4), rep(11,4), rep(12,4), 
                                            rep(13,4), rep(14,4), rep(15,4), rep(16,4)), 
                                          4),
                                      rep(c(rep(17,4), rep(18,4), rep(19,4), rep(20,4), 
                                            rep(21,4), rep(22,4), rep(23,4), rep(24,4)), 
                                          4)),
                          y_coord = rep(c(rep(1:4, 8), rep(5:8, 8), rep(9:12, 8), rep(13:16, 8)), 3),
                          bagnum_num = 1:384)
coords_bags$position <- paste0(coords_bags$x_coord, "_", coords_bags$y_coord)

# merge coords
exp_env_full_df <- merge(exp_env_df, coords_bags, by = c("bagnum_num"))

# treatment level
trt_temps
trt_sals

# expand the temp data
alltrt_temps <- data.frame(trt_sals$trt, trt_temps[rep(seq_len(nrow(trt_temps)), each = 2),])
colnames(alltrt_temps) <- c("trt", "temp_trt", "mean_temp", "median_temp", "min_temp", "max_temp", "sd_temp")

trt_env_full <- merge(alltrt_temps, trt_sals, by = "trt")

# save both the indiv and pop level env dfs
write.csv(exp_env_full_df, paste0("data/EnvDat/indiv_trt_env_dat_", Sys.Date(), ".csv"), row.names = F)
write.csv(trt_env_full, paste0("data/EnvDat/trt_env_dat_", Sys.Date(), ".csv"), row.names = F)
#######################
