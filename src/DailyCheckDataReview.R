# setup workspace
setwd("~/Desktop")
library(stringr)

# read data
bags <- read.csv("GOEEL-eelgrass-Sweden-exp_Bags.csv")
daily <- read.csv("GOEEL-eelgrass-Sweden-exp_DailyChecks.csv")

# check data
bags
daily

# which individuals are dead?
dead <- daily[which(daily$isAlive == FALSE),]$bagKey
dead_bags <- bags[which(bags$bagKey %in% dead),]$bagnum_num

# split date and time
date <- str_split_fixed(daily$dailyTimestamp, " ", n = 2)
daily$date <- date[,1]

# filter data to only today's date (TO DO - or can do this in spreadsheet before upload to R)
daily_today <- daily[(daily$date == "7/18/2024"),] # change for today's date

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
bags_dead_rm <- bags_shouldbe[!bags_shouldbe %in% dead_bags]
missing <- bags_dead_rm[!(bags_dead_rm %in% daily_ordered$bagnum_num)]
missing # these are missing
duplicated <- daily_ordered$bagnum_num[duplicated(daily_ordered$bagnum_num)] 
duplicated # these are duplicated

# double check that the dead bags not recorded are from previous day, but any dead from today do have data
length(bags_dead_rm) + length(missing) + length(dead_bags) - length(duplicated) == 384 
