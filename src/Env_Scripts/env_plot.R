##################
## ENV PLOTTING ##
##################

## setup
########
setwd("~/Documents/GitHub/BalticEelgrass")

# libraries
library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(rvg)
library(officer)

# read in data
copernicus_plot <- read.csv("data/EnvDat/copernicus_plotting_data.csv")
tanks_dat <- read.csv("data/EnvDat/all_tanks_logger_dat.csv")[,-1]
bags <- read.csv("data/EnvDat/GOEEL-eelgrass-Sweden-exp - Bags.csv")
daily_checks <- read.csv("data/EnvDat/GOEEL-eelgrass-Sweden-exp - DailyChecks.csv")
########

## prep pop dat
###############
# pop info
order_pops <- c("VAT","VIK","HOG","BAR","KUR","KAL","HOR","BJO","TempControl-21psu","TempControl-7psu","TempWarm-16psu","TempWarm-5psu")

# order pops_in_exp & set explicit colors
pops_in_df <- data.frame(poptrt = order_pops, 
                         order = 1:12, 
                         cols = c("#03045e","#0077b6","#00b4d8","#90e0ef","#6ede8a","#25a244","#155d27","#002800","#FFA203","#FE691E","deeppink","deeppink4"),
                         coast = c(rep("West",4), rep("East",4), rep("trt",4))) %>% 
  mutate(poptrt2 = poptrt) %>% 
  mutate(poptrt2 = ifelse(poptrt2 == "TempControl-21psu", "Current North Sea", 
                          ifelse(poptrt2 == "TempControl-7psu", "Current Baltic Sea", 
                                 ifelse(poptrt2 == "TempWarm-16psu", "Future North Sea", 
                                        ifelse(poptrt2 == "TempWarm-5psu", "Future Baltic Sea", poptrt2)))))
###############

## copernicus data
##################
copernicus_plot$date <- as.Date(copernicus_plot$date)
copernicus_plot$temp_adjust <- copernicus_plot$temp + 4
copernicus_plot$sal_adjust <- ifelse(copernicus_plot$coast == "West", copernicus_plot$sal - 5, copernicus_plot$sal)

# temp thru time at sites
ggplot(copernicus_plot, aes(x = date, y = temp_adjust, 
                            fill = fct_reorder(site_name, order), 
                            col = fct_reorder(site_name, order))) +
  geom_line()+
  #  geom_line(na.rm = T, aes(group = time_block)) +
  facet_grid(~ year, scales = "free_x") +
  scale_fill_manual(name = "Site",
                    values = pops_in_df$cols,
                    labels = pops_in_df$poptrt2) +
  scale_color_manual(name = "Site",
                     values = pops_in_df$cols,
                     labels = pops_in_df$poptrt2) +
  geom_hline(yintercept = 20.07134, col = "#FFA202", linewidth = 1) +
  geom_hline(yintercept = 25.80038, col = "#FF1393", linewidth = 1) +
  labs(title = "Temperature at Experimental Sites 2011-2021",
       x = "Date", y = "Temperature (ºC)", color = "Site", fill = "Site") +
  theme_bw()

# sal thru time at sites
ggplot(copernicus_plot, aes(x = date, y = sal, 
                            fill = fct_reorder(site_name, order), 
                            col = fct_reorder(site_name, order))) +
  geom_line()+
  #  geom_line(na.rm = T, aes(group = time_block)) +
  facet_grid(~ year, scales = "free_x") +
  scale_fill_manual(name = "Site",
                    values = pops_in_df$cols,
                    labels = pops_in_df$poptrt2) +
  scale_color_manual(name = "Site",
                     values = pops_in_df$cols,
                     labels = pops_in_df$poptrt2) +
  geom_hline(yintercept = 21, col = "#FFA202", linewidth = 1) +
  geom_hline(yintercept = 7, col = "#FF1393", linewidth = 1) +
  geom_hline(yintercept = 16, col = "#FE691F", linewidth = 1) +
  geom_hline(yintercept = 5, col = "#8B0A50", linewidth = 1) +
  labs(title = "Salinity at Experimental Sites 2011-2021",
       x = "Date", y = "Salinity (PSU)", color = "Site", fill = "Site") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme_bw()

# salinity smoothed
ggplot(copernicus_plot, 
       aes(x = date, y = sal, 
           fill = fct_reorder(site_name, order), 
           col = fct_reorder(site_name, order))) +
  geom_smooth() +
  scale_fill_manual(name = "Site",
                    values = pops_in_df$cols,
                    labels = pops_in_df$poptrt2) +
  scale_color_manual(name = "Site",
                     values = pops_in_df$cols,
                     labels = pops_in_df$poptrt2) +
  geom_hline(yintercept = 21, col = "#FFA202", linewidth = 1) +
  geom_hline(yintercept = 7, col = "#FF1393", linewidth = 1) +
  geom_hline(yintercept = 16, col = "#FE691F", linewidth = 1) +
  geom_hline(yintercept = 5, col = "#8B0A50", linewidth = 1) +
  labs(title = "Salinity at Experimental Sites 2011-2021",
       x = "Date", y = "Salinity (PSU)", color = "Site", fill = "Site") +
  theme_bw()
##################

## temp data 
############
# expand out the date time info from tanks
tanks_dat$date <- str_split_fixed(tanks_dat$datetime, " ", 2)[,1]
tanks_dat$date <- as.Date(tanks_dat$date)
tanks_dat$time <- str_split_fixed(tanks_dat$datetime, " ", 2)[,2]
tanks_dat$poptrt <- tanks_dat$trt
tanks_dat <- tanks_dat %>% mutate(poptrt = ifelse(poptrt == "ctrl", "Current", "Future"))

# duplicate the current/future data for later plotting
tanks_dat1 <- tanks_dat %>% mutate(poptrt2 = poptrt) %>% 
  mutate(poptrt2 = ifelse(poptrt2 == "Current", "Current North Sea", "Future North Sea"))
tanks_dat2 <- tanks_dat %>% mutate(poptrt2 = poptrt) %>% 
  mutate(poptrt2 = ifelse(poptrt2 == "Current", "Current Baltic Sea", "Future Baltic Sea"))
tanks_dat3 <- rbind(tanks_dat1, tanks_dat2)

# look at just these data
ggplot(tanks_dat) +
  geom_violin(aes(x = poptrt, y = temp, fill = poptrt)) +
  #  geom_boxplot(aes(x = trt, y = temp, fill = trt), width = 0.1) +
  scale_fill_manual(name = "Treatment",
                    values = c("#FFA202","#FF1393"),
                    labels = c("Current", "Future")) +
  theme_bw()

# with the duplicated treatments
ggplot(tanks_dat3) +
  geom_violin(aes(x = poptrt2, y = temp, fill = poptrt2)) +
  #  geom_boxplot(aes(x = trt, y = temp, fill = trt), width = 0.1) +
  scale_fill_manual(name = "Treatment",
                    values = pops_in_df[9:12,]$cols,
                    labels = pops_in_df[9:12,]$poptrt2) +
  labs(title = "(C) Temperature in Experimental Mesocosms", 
       x = NULL, y = "Temperature (ºC)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# also plot by tank
tank_cols <- data.frame(tank_ID = levels(as.factor(tanks_dat$tank_ID)), 
                        col = c("#03045e","#023e8a","#590d22","#800f2f","#c9184a","#ff4d6d","#ff758f","#ffb3c1","#0077b6","#00b4d8","#48cae4","#ade8f4"))
tanks_plotting <- merge(tank_cols, tanks_dat, by = c("tank_ID"))

ggplot(tanks_plotting, aes(x = datetime, y = temp, 
                           col = tank_ID)) +
  scale_color_manual(values = tank_cols$col)+
  geom_line()+
  labs(title = "Temperature in Mesocosm Tanks", 
       x = "Date", y = "Temperature (ºC)")+
  theme(axis.ticks.x = element_blank()) +
  theme_classic()

# combine with site temps
cop_temps <- copernicus_plot[,c("date","poptrt","temp_adjust")]
colnames(cop_temps) <- c("date","poptrt","temp")

tanks_dat_red <- tanks_dat3[,c("date","poptrt2","temp")]
colnames(tanks_dat_red) <- c("date", "poptrt", "temp")

temps_all <- rbind(tanks_dat_red[,c("date","poptrt","temp")], 
                   cop_temps)

temps_all$poptrt2 <- temps_all$poptrt

temps_plotting <- merge(temps_all, pops_in_df[,colnames(pops_in_df) != "poptrt"], by = c("poptrt2"))

# plot
temp_plot <- ggplot(temps_plotting) +
  geom_violin(aes(x = as.factor(fct_reorder(poptrt2, order)), y = temp, fill = fct_reorder(poptrt, order))) +
  #  geom_boxplot(aes(x = trt, y = temp, fill = trt), width = 0.1) +
  scale_fill_manual(name = "Site",
                    values = pops_in_df$cols,
                    labels = pops_in_df$poptrt2) +
  labs(title = "(B) Temperature at Experimental Source\nSites and in Mesocosms", 
       x = NULL, y = "Temperature (ºC)", color = "Site", fill = "Site") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
temp_plot
############

## salinity data
################
# check data
summary(daily_checks)

# there are some massive outliers, remove these
daily_checks_clean <- daily_checks %>%
  filter(dailyTemp < 40 & dailySalinity <30)

# put together the salinity data so we can look at trends by tank/treatment
sal_dat <- merge(bags, daily_checks_clean, by = c("bagKey"))

# expand out date time info from daily checks
sal_dat$date <- str_split_fixed(sal_dat$dailyTimestamp, " ", 2)[,1]
sal_dat$date <- mdy(sal_dat$date)
sal_dat$time <- str_split_fixed(sal_dat$dailyTimestamp, " ", 2)[,2]
sal_dat$poptrt <- sal_dat$trt
sal_dat1 <- sal_dat %>% mutate(poptrt2 = ifelse(poptrt == "TempControl-21psu", "Current North Sea",
                                                ifelse(poptrt == "TempControl-7psu", "Current Baltic Sea",
                                                       ifelse(poptrt == "TempWarm-16psu", "Future North Sea", "Future Baltic Sea"))))

# plot just these data
ggplot(sal_dat1) +
  geom_violin(aes(x = as.factor(poptrt2), y = dailySalinity, fill = poptrt2)) +
  scale_fill_manual(name = "Treatment",
                    values = pops_in_df[9:12,]$cols,
                    labels = pops_in_df[9:12,]$poptrt2) +
  labs(title = "Salinity in Mesocosm Treatments",
       x = "Site", y = "Salinity (PSU)", color = "Treatment", fill = "Treatment") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# it looks like some bag ids got mixed up between current W and future E
sal_dat_trim <- sal_dat1[!(sal_dat1$trt == "TempControl-21psu" & sal_dat1$dailySalinity < 18),]
sal_dat_clean <- sal_dat_trim[!(sal_dat_trim$trt == "TempWarm-5psu" & sal_dat_trim$dailySalinity > 8),]

# plot again
ggplot(sal_dat_clean) +
  geom_violin(aes(x = as.factor(poptrt2), y = dailySalinity, fill = poptrt2)) +
  scale_fill_manual(name = "Site",
                    values = pops_in_df[9:12,]$cols,
                    labels = pops_in_df[9:12,]$poptrt2) +
  labs(title = "Salinity in Mesocosm Treatments",
       x = "Treatment", y = "Salinity (PSU)", color = "Treatment", fill = "Treatment") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# df for all salinity now
cop_sals <- copernicus_plot[,c("date","poptrt","sal")]
colnames(cop_sals) <- c("date","poptrt","dailySalinity")
cop_sals$poptrt2 <- cop_sals$poptrt

sal_all <- rbind(sal_dat_clean[,c("date", "poptrt", "dailySalinity", "poptrt2")],
                 cop_sals)

sal_plotting <- merge(sal_all, pops_in_df, by = c("poptrt", "poptrt2"))

# plot
sal_plot <- ggplot(sal_plotting) +
  geom_violin(aes(x = as.factor(fct_reorder(poptrt2, order)), 
                  y = (dailySalinity), 
                  fill = fct_reorder(poptrt2, order))) +
  scale_fill_manual(name = "Site",
                    values = pops_in_df$cols,
                    labels = pops_in_df$poptrt2) +
  labs(title = "(D) Salinity at Experimental Source\nSites and in Mesocosms",
       x = NULL, y = "Salinity (PSU)", color = "Site", fill = "Site") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
sal_plot
################

## export
#########
doc <- read_pptx() %>% 
  add_slide(layout = "Title and Content", master = "Office Theme") %>% 
  ph_with(value = dml(ggobj = temp_plot), 
          location = ph_location(left = 1, top = 1, width = 8, height = 3.5)) %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>% 
  ph_with(value = dml(ggobj = sal_plot), 
          location = ph_location(left = 1, top = 1, width = 8, height = 3.5))
  
print(doc, "results/figures/env_plots.pptx")
#########
