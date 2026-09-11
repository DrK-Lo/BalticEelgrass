##################
## ENV DISTANCE ##
##################

## setup
########
# set working directory
setwd("~/Documents/GitHub/BalticEelgrass")

# libraries
library(BiocManager) # needed to download specific packages
library(readxl) # for excel files
library(vegan) # euclidean distance
library(tidyr) # data cleaning
library(dplyr) # data cleaning
library(stringr) # data cleaning
library(sdmpredictors) # for downloading environmental data
######## 

## data
#######
# read in env data
seascape_scaled_ind1 <- read.csv("data/EnvDat/full_seascape_scaled_ind_only_2026-09-04.csv")
exp_scaled_ind1 <- read.csv("data/EnvDat/exp_scaled_ind_only_2026-09-04.csv")
trt_scaled_ind1 <- read.csv("data/EnvDat/trt_scaled_ind_only_2026-09-04.csv")
#######

## calc dist
############
# assign relevant data to new objs
trt_scaled_fordist_ind <- trt_scaled_ind1
exp_scaled_fordist_ind <- exp_scaled_ind1

# calc euclidean distance between matched row/trt vals
ind_dists <- sqrt(rowSums((exp_scaled_fordist_ind[,c("median_temp_std", "max_temp_std", "median_sal_std", "min_sal_std")] - trt_scaled_fordist_ind[,c("median_temp_std", "max_temp_std", "median_sal_std", "min_sal_std")])^2))

# add metadata
ind_dists_df <- cbind(trt_scaled_fordist_ind[,c("pop", "trt", "bagnum_text")], ind_dists)

# aggregate to one distance per trt x site (mean across all bags in that trt)
dist_summary <- ind_dists_df %>%
  group_by(trt, pop) %>%
  summarize(dist = mean(ind_dists, na.rm = TRUE))

# keep these dfs
write.csv(dist_summary, paste0("data/EnvDat/euclidian_trt_exp_summary_", Sys.Date(), ".csv"))
write.csv(ind_dists_df, paste0("data/EnvDat/euclidean_trt_exp_indiv_", Sys.Date(), ".csv"))
############

## plotting
###########
# pop and trt orders for plotting
pops <- c("VAT","VIK","HOG","BAR","KUR","KAL","HOR","BJO")
trts <- c("TempControl-21psu", "TempControl-7psu", "TempWarm-16psu", "TempWarm-5psu")
popsdf <- data.frame(pop = pops, pop_order = 1:length(pops))
trtsdf <- data.frame(trt = trts, trt_order = 1:length(trts))

dist_summary1 <- merge(popsdf, dist_summary, by = "pop")
dist_summary2 <- merge(trtsdf, dist_summary1, by = "trt")

# plotting
euc_dist_ind <- ggplot(dist_summary2, aes(x = fct_reorder(pop, pop_order), y = fct_reorder(trt, trt_order), fill = dist)) +
  theme_classic() +
  geom_tile(color = "black") +
  scale_fill_viridis_c(name = "Euclidean \nDistance", limits = c(0, max(dist_summary2$dist))) +
  theme(plot.title = element_text(size = 20),
        legend.title = element_text(size = 14),
        legend.position = "right", legend.justification = "top") +
  ylab("") + xlab("") +
  ggtitle("Environmental Distance between \nSites and Treatments") +
  theme(plot.title = element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.position = "right", legend.justification = "top",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
euc_dist_ind
###########

