##################
## ENV DISTANCE ##
##################

## SETUP
########
# set working directory
setwd("~/Eelgrass/BalticEelgrass/data")

# libraries
library(BiocManager) # needed to download specific packages
library(readxl) # for excel files
library(vegan) # euclidean distance
library(tidyr) # data cleaning
library(dplyr) # data cleaning
library(stringr) # data cleaning
library(sdmpredictors) # for downloading environmental data
######## 

## DATA
#######
# read in site data
sites <- read.csv("seascape_data/sampling_sites_coordinates_Baltic_Sea.csv")
indivs <- read_excel("seascape_data/MLL_per_individual_list.xlsx")
#######

## ENV DATA
###########
# looking for environmental data from 8 experimental sites, 40 seascape sites
###########

