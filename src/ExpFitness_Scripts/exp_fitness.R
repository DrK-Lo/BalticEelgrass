##########################
## Experimental Fitness ##
##########################

## setup
########
# set wd
setwd("~/Eelgrass/BalticEelgrass/data")

# load packages
library(stringr) # for data cleaning
library(tidyr) # for data cleaning & wrangling
library(dplyr) # for adata analysis
library(ggplot2) # for plotting
library(forcats) # for plotting
library(plotrix) # for std err
########

## data
#######
# read in data
setup <- read.csv("experiment/GOEEL-eelgrass-Sweden-exp - Setup.csv")
takedown <- read.csv("experiment/GO_eelgrass_mesocosm_takedown_data.csv")
head(setup)
head(takedown)

# trim
takedown_trim <- takedown[,c("bagnum_num","tank_ID","DNA_label","new_label",
                          "sal_full","temp_full","nbr_leaves","length_longest_leaf",
                          "total_weight","isAlive")]
head(takedown_trim)

# populations & treatments
takedown_trim$pop <- str_split_fixed(takedown_trim$DNA_label, "_", 2)[,1]
takedown_trim$trt <- paste0(takedown_trim$temp_full, "_", takedown_trim$sal_full)
head(takedown_trim)
takedown_trim

## mortality
############
# pops
pops <- levels(as.factor(takedown_trim$pop))

# trts
trts <- levels(as.factor(takedown_trim$trt))

# mortality by pop
# set up df
mort_pop <- matrix(data = NA, nrow = length(pops), ncol = 2)

# calc
for (i in 1:length(pops)) {
  # split pops
  split_mort = takedown_trim[takedown_trim$pop == pops[i],]
  mort_pop[i,1] = pops[i]
  # number that survived
  mort_pop[i,2] = length(split_mort[split_mort$isAlive == "TRUE",]$isAlive)/length(split_mort$isAlive)
}

# clean up
mort_pop <- as.data.frame(mort_pop)
colnames(mort_pop) <- c("pop","surv")
mort_pop$surv <- as.numeric(mort_pop$surv)

# mortality by treatment
mort_trt <- matrix(data = NA, nrow = length(trts), ncol = 2)
for (i in 1:length(trts)) {
  # split pops
  split_mort = mort_trim[mort_trim$trt == trts[i],]
  mort_trt[i,1] = trts[i]
  # number that survived
  mort_trt[i,2] = length(split_mort[split_mort$isAlive == "TRUE",]$isAlive)/length(split_mort$isAlive)
}
mort_trt <- as.data.frame(mort_trt)
colnames(mort_trt) <- c("trt","surv")
mort_trt$surv <- as.numeric(mort_trt$surv)

# mortality by population & treatment
# add trt_pop column
takedown_trim$trt_pop <- paste0(takedown_trim$trt,"-",takedown_trim$pop)
trt_pop <- levels(as.factor(takedown_trim$trt_pop))

# calc
mort_trt_pop <- matrix(data = NA, nrow = length(trt_pop), ncol = 2)
for (i in 1:length(trt_pop)) {
  # split pops
  split_mort = takedown_trim[takedown_trim$trt_pop == trt_pop[i],]
  mort_trt_pop[i,1] = trt_pop[i]
  # number that survived
  mort_trt_pop[i,2] = length(split_mort[split_mort$isAlive == "TRUE",]$isAlive)/length(split_mort$isAlive)
}
mort_trt_pop <- as.data.frame(mort_trt_pop)
colnames(mort_trt_pop) <- c("trt_pop","surv")
mort_trt_pop$surv <- as.numeric(mort_trt_pop$surv)

# add in trt and pop as separate columns
mort_trt_pop[,3:4] <- str_split_fixed(mort_trt_pop$trt_pop, "-", 2)
colnames(mort_trt_pop) <- c("trt_pop","surv","trt","pop")

# population order
pop_order <- data.frame(pop = c("VIK","VAT","HOG","BAR","KUR","KAL","HOR","BJO"), 
                        order = c(1:8))

mort_trt_pop_mg <- merge(mort_trt_pop, pop_order, by = "pop")

# by coast?
mort_trt_pop_mg$coast <- ifelse(mort_trt_pop_mg$order >= 5, "East", "West")

# more informative trt names
trt <- levels(as.factor(mort_trt_pop_mg$trt))
names <- c("Current W Coast", "Current E Coast", "Future W Coast", "Future E Coast")
trt_nm <- data.frame(trt = trt, trt_name = names)
mort_trt_pop_mg <- merge(mort_trt_pop_mg, trt_nm, by = "trt")

# calc by coast
mortality_coasttrt_surv <- mort_trt_pop_mg %>% group_by(coast,trt) %>%
  summarise(mean_surv = mean(as.numeric(surv)), 
            sd_surv = sd(as.numeric(surv)),
            se_surv = std.error(as.numeric(surv)),
            .groups = 'drop')
mortality_coasttrt_surv <- merge(mortality_coasttrt_surv, trt_nm, by = "trt")
mortality_coasttrt_surv$CurrFut <- str_split_fixed(as.character(mortality_coasttrt_surv$trt_name)," ",2)[,1]
############

## plotting
###########
# pop x trt
ggplot(mort_trt_pop_mg, aes(x = fct_reorder(pop, order), y = surv, fill = trt)) +
  theme_classic() +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(name = "Treatment",
                    values = c("#FFA203","#FE691E","deeppink","deeppink4"),
                    labels = c("Current W Coast", "Current E Coast",
                               "Future W Coast", "Future E Coast")) +
  theme(plot.title = element_text(size = 24), 
        plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "right") +
  xlab("Population") + 
  ylab("Survival") +
  ggtitle("Survival by Population & Treatment", subtitle = "Current-Future Comparison")

# trt x pop
trt_order <- c("Current W Coast","Future W Coast","Current E Coast","Future E Coast")
ggplot(mort_trt_pop_mg, aes(x = factor(trt_name, level = trt_order),
                            y = surv, fill = fct_reorder(pop,order))) +
  theme_classic() +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(name = "Population",
                    values = c("#155d27","#25a244","#4ad66d","#b7efc5",
                               "#caf0f8","#90e0ef","#0096c7","#023e8a"),
                    labels = c("VIK","VAT","HOG","BAR","KUR","KAL","HOR","BJO")) +
  theme(plot.title = element_text(size = 24), 
        plot.subtitle = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "right") +
  xlab("Treatment") + 
  ylab("Survival") +
  ggtitle("Survival by Treatment & Population", subtitle = "Common Garden Comparison")

ggplot(mortality_coasttrt_surv, aes(x = fct_reorder(trt_name, coast),
                                    y = mean_surv, fill = fct_rev(coast))) +
  theme_classic() +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  theme(plot.title = element_text(size = 24), 
        legend.title = element_text(size = 15),
        legend.position = "right") +
  geom_errorbar(aes(ymin=mean_surv-se_surv, ymax=mean_surv+se_surv), width=.2,
                position=position_dodge(.9))+
  scale_fill_manual(name = "Coast",
                    values = c("#95cf92","#369acc"),
                    labels = c("West","East")) +
  xlab("Treatment") + 
  ylab("Mean Survival") +
  ggtitle("Survival by Treatment & Coast")
###########


## growth data
##############
# check data
head(setup)
head(takedown)

# need to process setup data
View(setup)

# remove practice rows (bagKey F054C074) & trim dataset
setup_rm <- setup[!setup$bagKey == "F054C074",]
setup_rm2 <- setup_rm[,c("setupKey","setupLabel","bagKey","setupTimestamp","nbr_leaves",
                         "sideShoot","nbr_sideShootLeaves","length_longest_leaf",
                         "leaf_width_hole","weight","IsBlackLesions")]

# get more meaningful labels in setup
dim(setup_rm2)
setup_rm2[,12:15] <- str_split_fixed(setup_rm$setupLabel, "_", 5)[,1:4]
colnames(setup_rm2) <- c("setupKey","setupLabel","bagKey","setupTimestamp","nbr_leaves",
                         "sideShoot","nbr_sideShootLeaves","length_longest_leaf",
                         "leaf_width_hole","weight","IsBlackLesions","bagnum",
                         "tank_ID","pop","genet")

# we also had to get rid of one individual because of planting error - ind 274
# setup_rm3 <- setup_rm2[!setup_rm2$bagnum == 274,]

# match setup and takedown
# can do this using "DNA_label" stuctured like HOR_gen-01
setup_rm2$DNA_label <- paste0(setup_rm2$pop,"_",setup_rm2$genet)

# check this 
head(setup_rm2)

# now put together setup and takedown data
full_dat <- merge(setup_rm3, takedown, by = c("DNA_label"))
dim(full_dat) # we've lost two individuals... one was the guy we had to remove

# which individual was not supposed to be lost?
setup_rm3[which(!(setup_rm3$DNA_label %in% takedown$DNA_label)),]
# VAT_gen-13 is missing from takedown data, confirmed by inspecting data sheet

takedown[which(!(takedown$DNA_label %in% setup_rm3$DNA_label)),]
# bagnum 274 is missing, which it should be since we removed it

# is there a different one that was duplicated in takedown? in setup?
which(duplicated(takedown$DNA_label)) # nope
which(duplicated(setup_rm3$DNA_label)) # nope

# LEAVE IT HERE FOR NOW

# how can I determine overall length growth?
# add together growth of all existing leaves + length of new leaves
# for second shoots, can also add in length of new shoots' leaves

# weight growth - can look at weight before and after - need the before data

##############