# setup workspace
setwd("~/Desktop")
#install.packages("stringr")
library(stringr)
library(lubridate)

# make sure to download the most current daily check data before running this script

# read data
bags <- read.csv("GOEEL-eelgrass-Sweden-exp_Bags.csv")
daily <- read.csv("GOEEL-eelgrass-Sweden-exp_DailyChecks.csv")

# check data
bags
daily

# split date and time
date <- str_split_fixed(daily$dailyTimestamp, " ", n = 2)
daily$date <- date[,1]

# reformat date
newdate <- strptime(as.character(daily$date), "%m/%d/%Y")
daily$dateformatted <- format(newdate, "%Y-%m-%d")

# filter data to only today's date
daily_today <- daily[(daily$dateformatted == Sys.Date()),] # can also replace Sys.Date() with date of interest
daily_prev <- daily[!(daily$dateformatted == Sys.Date()),] # can also replace Sys.Date() with date of interest

# which individuals are dead today?
dead <- daily_today[which(daily_today$isAlive == FALSE),]$bagKey
dead_bags <- bags[which(bags$bagKey %in% dead),]$bagnum_num

# which were dead prior to today?
dead_prev <- daily_prev[which(daily_prev$isAlive == FALSE),]$bagKey
dead_prev_bags <- bags[which(bags$bagKey %in% dead_prev),]$bagnum_num

# merge daily and bags
daily_tot <- merge(bags, daily_today, by = c("bagKey"))
#View(daily_tot)

# order by bagnum
daily_ordered <- daily_tot[order(daily_tot$bagnum_num),]
#View(daily_ordered)

# find number of entries - should be 384
length(levels(as.factor(daily_ordered$bagnum_num)))

# which ones are missing? any duplicates?
bags_shouldbe <- c(1:384)
bags_dead_rm <- bags_shouldbe[!bags_shouldbe %in% dead_prev_bags]
missing <- bags_dead_rm[!(bags_dead_rm %in% daily_ordered$bagnum_num)]
missing # these are missing
duplicated <- daily_ordered$bagnum_num[duplicated(daily_ordered$bagnum_num)] 
duplicated # these are duplicated

# double check that the dead bags not recorded are from previous day, but any dead from today do have data
length(bags_dead_rm) + length(missing) + length(dead_prev_bags) - length(duplicated) == 384 

