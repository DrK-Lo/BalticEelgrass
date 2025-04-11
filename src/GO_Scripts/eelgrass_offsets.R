######################
## EELGRASS OFFSETS ##
######################

## SETUP
########
# set working directory
setwd("~/Eelgrass/BalticEelgrass/data")

# libraries
library(BiocManager) # needed to download specific packages
library(readxl) # reading xlsx files
library(vcfR) # reading vcf files
library(LEA) # landscape genomics
library(psych) # environmental data cleaning
library(vegan) # environmental data cleaning
library(tidyr) # data cleaning
library(dplyr) # data cleaning
library(stringr) # data cleaning
library(sdmpredictors) # for downloading environmental data - NOT LOADING IN CLUSTER
library(qvalue) # for calculating qvalues
library(ggplot2) # plotting
library(s2) # mapping
library(rnaturalearth) # mapping - NOT LOADING IN CLUSTER
library(rnaturalearthdata) # mapping
library(maps) # mapping
library(ggspatial) # mapping - NOT LOADING IN CLUSTER
########

## DATA
#######
# read in data
sites <- read.csv("seascape_data/sampling_sites_coordinates_Baltic_Sea.csv")
indiv_list <- read_excel("seascape_data/MLL_per_individual_list.xlsx")
gen_div <- read_excel("seascape_data/Gen_div_table.xlsx")

# read in vcf data
eelgrass_vcf_MLL <- read.vcfR(
  "seascape_data/zostera_230619_seascape_miss25_sorted_95_MLL.recode.vcf")
eelgrass_vcf_full <- read.vcfR(
  "seascape_data/zostera_230619_seascape_miss25_sorted_95_no_reps.recode.vcf")
#######

## GEN DATA PREP
################
# investigate vcf data - MLL
dim(eelgrass_vcf_MLL) # 3287 variants, 8 fix cols, 242 genotypes
fix_MLL <- as.data.frame(eelgrass_vcf_MLL@fix)
gt_MLL <- as.data.frame(eelgrass_vcf_MLL@gt)

# how many variants pass filtering?
table(fix_MLL$FILTER) # all pass

# quality summary
summary(as.numeric(fix_MLL$QUAL))

#investigate vcf data - FULL
dim(eelgrass_vcf_full) # 4149 variants, 8 fix cols, 696 genotypes
fix_full <- as.data.frame(eelgrass_vcf_full@fix)
gt_full <- as.data.frame(eelgrass_vcf_full@gt)

# how many variants pass filtering?
table(fix_full$FILTER) # all pass

# how many individuals in each population?
head(indiv_list)
pops_split <- stringr::str_split_fixed(indiv_list$`Indv #`, "-", 2) # split Indv # to have population in a seperate column
indiv_list_full <- cbind(indiv_list, pops_split) # add this info back to original file
colnames(indiv_list_full) <- c("MLL_no","Indv_no","Pop","Pop_no")
table(indiv_list_full$Pop) # between 4 and 20 individuals in each pop - most 17+ inds

# need indiv data for MLL dataset
# individuals are stored in column names of vcf - isolate these names
colnames(eelgrass_vcf_MLL@gt)
MLL_split <- stringr::str_split_fixed(colnames(eelgrass_vcf_MLL@gt), "-", 3)
MLL_inds <- MLL_split[-1,3]

# now get subsetted indiv dataframe from indiv list full
indiv_list_MLL <- indiv_list_full[indiv_list_full$Indv_no %in% MLL_inds,]
################

## LFMM & GENO FORMATS
######################
# convert vcf files to lfmm file format
eelgrass_lfmm_MLL <- vcf2lfmm("seascape_data/zostera_230619_seascape_miss25_sorted_95_MLL.recode.vcf")
eelgrass_lfmm_full <- vcf2lfmm("seascape_data/zostera_230619_seascape_miss25_sorted_95_no_reps.recode.vcf")

# read in geno and lfmm files
eelgrass_lfmm_MLL <- read.lfmm("seascape_data/zostera_230619_seascape_miss25_sorted_95_MLL.recode.lfmm")
eelgrass_geno_MLL <- read.geno("seascape_data/zostera_230619_seascape_miss25_sorted_95_MLL.recode.geno")
eelgrass_lfmm_full <- read.lfmm("seascape_data/zostera_230619_seascape_miss25_sorted_95_no_reps.recode.lfmm")
eelgrass_geno_full <- read.geno("seascape_data/zostera_230619_seascape_miss25_sorted_95_no_reps.recode.geno")
eelgrass_lfmm_full_rm <- read.lfmm("seascape_data/eel_full_rm.lfmm")

# check for missing data
sum((eelgrass_lfmm_MLL == 9)) # 16894/174269 missing data (~9.7% missing data)
sum((eelgrass_lfmm_full_rm == 9)) # 51254/615284 missing data (~8.3% missing data)

# impute missing data - this is a super naive imputation method, replace with more robust method using snmf later
gen_full_imp <- apply(eelgrass_lfmm_full_rm, 2, function(x) replace(x, x == 9, as.numeric(names(which.max(table(x))))))
sum((gen_full_imp == 9)) 

gen_MLL_imp <- apply(eelgrass_lfmm_MLL, 2, function(x) replace(x, x == 9, as.numeric(names(which.max(table(x))))))
sum((gen_MLL_imp == 9))

# remove non-variable genomic sites
gen_full_imp_df <- as.data.frame(gen_full_imp)
gen_full_imp_trim <- Filter(var, gen_full_imp_df)
gen_full_imp_mat <- as.matrix(gen_full_imp_trim)

gen_MLL_imp_df <- as.data.frame(gen_MLL_imp)
gen_MLL_imp_trim <- Filter(var, gen_MLL_imp_df)
gen_MLL_imp_mat <- as.matrix(gen_MLL_imp_trim)

# re-write lfmm files
write.lfmm(gen_full_imp_mat, "seascape_data/gen_full_rm_imp_mat.lfmm")
write.lfmm(gen_MLL_imp_mat, "seascape_data/gen_MLL_imp_mat.lfmm")

# also re-write geno files
write.geno(gen_full_imp_mat, "seascape_data/gen_full_rm_imp_mat.geno")
write.geno(gen_MLL_imp_mat, "seascape_data/gen_MLL_imp_mat.geno")

# read in lfmm files
gen_full <- read.lfmm("seascape_data/gen_full_rm_imp_mat.lfmm") # has 3264 SNPs after filtering
gen_MLL <- read.lfmm("seascape_data/gen_MLL_imp_mat.lfmm") # has 3036 SNPs after filtering
######################

## ANCESTRY
###########
# create snmf object to look at ancestry
MLL_snmf <- snmf("seascape_data/zostera_230619_seascape_miss25_sorted_95_MLL.recode.lfmm", 
                 K = 1:40, project = "new", repetitions = 5, tolerance = 0.00001, 
                 entropy = TRUE, ploidy = 2)

full_snmf <- snmf("seascape_data/zostera_230619_seascape_miss25_sorted_95_no_reps.recode.lfmm", 
                  K = 1:40, project = "new", repetitions = 5, tolerance = 0.00001, 
                  entropy = TRUE, ploidy = 2)

full_snmf_1 <- snmf("seascape_data/zostera_230619_seascape_miss25_sorted_95_no_reps.recode.lfmm", 
                    K = 1:250, project = "new", repetitions = 1, tolerance = 0.00001, 
                    entropy = TRUE, ploidy = 2)

# plot cross entropy criterion
plot(MLL_snmf, cex = 1.2, col = "lightblue", pch = 19, ylim = c(0.076,0.078)) # K = 12 (or something around there)
plot(full_snmf, cex = 1.2, col = "lightblue", pch = 19) # K is seemingly infinite

# select run with lowest cross-entropy
entropy_MLL <- cross.entropy(MLL_snmf, K = 12)
best_MLL <- which.min(entropy_MLL)

entropy_full <- cross.entropy(full_snmf, K = 242)
best_full <- which.min(entropy_full)

# generate Q-matrix
qmatrix_MLL <- Q(MLL_snmf, K = 12, run = best_MLL)
head(qmatrix_MLL)

qmatrix_full <- Q(full_snmf, K = 242, run = best_full)
head(qmatrix_full)

# generate plots
ancestry_MLL <- barplot(t(qmatrix_MLL),
                        col = c("brown2","coral3","darkorange","darkgoldenrod1","darkolivegreen2",
                                "darkolivegreen4","cadetblue1","deepskyblue3","darkorchid1",
                                "darkviolet","palevioletred1","violetred1"),
                        border = NA,
                        xlab = "Individuals", ylab = "Admixture coefficients", 
                        main = "Ancestry Matrix")

ancestry_full <- barplot(t(qmatrix_full),
                         col = c("#648FFF","#FFB000","#DC267F"),
                         border = NA, space = .2,
                         xlab = "Individuals", ylab = "Admixture coefficients", 
                         main = "Ancestry Matrix")

# determine cluster assignment for each individual from MLL run

###########

## ENVR DATA PREP
#################
# inspect site data
sites 
# missing 11/39 temp, 18/39 salinity, 4/39 depth

# separate minimum and maximum depth range
sites[c("min_depth_m", "max_depth_m")] <- strsplit(sites$depth_m, "-", 2)

# turn all envr data into numeric
sites$water_temperature_celcius <- as.numeric(sites$water_temperature_celcius)
sites$salinitu_psu <- as.numeric(sites$salinitu_psu)

# will need to get environmental data from another source to fill gaps
# bio-oracle is being made obsolete - what to use now?
datasets <- list_datasets(marine = T)
list_layers(datasets)
envr_layers <- load_layers(layercodes = c("BO_ph","BO_salinity","BO_sstmax","BO_sstmean","BO_sstmin",
                                          "BO_sstrange"),
                           rasterstack = F)

sites_environment <- data.frame(Name = sites$site_full, 
                                pop = sites$site,
                                lat = sites$lat,
                                long = sites$long,
                                BO_ph = raster::extract(envr_layers$BO_ph, sites[,5:6]),
                                BO_salinity = raster::extract(envr_layers$BO_salinity, sites[,5:6]),
                                BO_sstmax = raster::extract(envr_layers$BO_sstmax, sites[,5:6]),
                                BO_sstmean = raster::extract(envr_layers$BO_sstmean, sites[,5:6]),
                                BO_sstmin = raster::extract(envr_layers$BO_sstmin, sites[,5:6]),
                                BO_sstrange = raster::extract(envr_layers$BO_sstrange, sites[,5:6])
)

# save environmental data as csv
write.csv(sites_environment, "seascape_data/sites_environment_021325.csv")

# future data
list_layers_future(marine = T)

future <- get_future_layers(c("BO_salinity","BO_sstmax","BO_sstmean","BO_sstmin",
                              "BO_sstrange"),
                            scenario = "A2", year = 2100)

future_load <- load_layers(layercodes = c("BO_A2_2100_salinity", "BO_A2_2100_sstmax", "BO_A2_2100_sstmean",
                                          "BO_A2_2100_sstmin", "BO_A2_2100_sstrange"),
                           rasterstack = F)

future_environment <- data.frame(Name = sites$site_full, 
                                 pop = sites$site,
                                 lat = sites$lat,
                                 long = sites$long,
                                 BO_salinity = raster::extract(future_load$BO_A2_2100_salinity, sites[,5:6]),
                                 BO_sstmax = raster::extract(future_load$BO_A2_2100_sstmax, sites[,5:6]),
                                 BO_sstmean = raster::extract(future_load$BO_A2_2100_sstmean, sites[,5:6]),
                                 BO_sstmin = raster::extract(future_load$BO_A2_2100_sstmin, sites[,5:6]),
                                 BO_sstrange = raster::extract(future_load$BO_A2_2100_sstrange, sites[,5:6]))

write.csv(future_environment, "seascape_data/future_environment_021325.csv")

# read back in downloaded environmental data
current_environment <- read.csv("seascape_data/sites_environment_021325.csv")
current_environment <- current_environment[c("Name","BO_ph","BO_salinity",
                                             "BO_sstmax","BO_sstmean","BO_sstmin",
                                             "BO_sstrange")]
colnames(current_environment)[1] <- "site_full"

future_environment <- read.csv("seascape_data/future_environment.csv")
future_environment <- future_environment[c("Name","BO_salinity","BO_sstmax",
                                  "BO_sstmean","BO_sstmin","BO_sstrange")]
colnames(future_environment)[1] <- "site_full"

# combine with other site data to compare accuracy
sites_complete_data <- merge(sites, sites_environment, by = "site_full")
# based on comparison of recorded data and bio-oracle data... 
# bio-oracle is a bad source of data for this region - for now use anyway

# prep environmental data for LEA
# inspect data
current_environment
future_environment
indiv_list_full
indiv_list_MLL

# make sure site codes match
sites$site
table(indiv_list_full$Pop)
sites$site[sites$site == "ÅLA"] <- "ALA"
sites$site[sites$site == "BÅD"] <- "BAD"

# need full site name on indiv lists
site_info <- data.frame(site_full = sites$site_full, Pop = sites$site)
indiv_full_merge <- merge(indiv_list_full, site_info, by = "Pop")
indiv_MLL_merge <- merge(indiv_list_MLL, site_info, by = "Pop")

# create dataframe with environment for each individual
current_env_full <- merge(indiv_full_merge, current_environment, by = "site_full")
current_env_MLL <- merge(indiv_MLL_merge, current_environment, by = "site_full")
future_env_full <- merge(indiv_full_merge, future_environment, by = "site_full")
future_env_MLL <- merge(indiv_MLL_merge, future_environment, by = "site_full")

# also check for correlated data
pairs.panels(current_env_full[,-c(1:5)]) # sstmean and max/min are highly correlated (.95/.77)
pairs.panels(current_env_MLL[,-c(1:5)]) # sstmean and max/min(.97/.94)
pairs.panels(future_env_full[,-c(1:5)]) # sstmean and max/min(.87/.83)
pairs.panels(future_env_MLL[,-c(1:5)]) # sstmean and max/min(.9/.91)

# note: based on correlations, may need to get rid of sstmean (V4 current, V3 future)
# also don't have ph in future, so get rid of that
current_full <- current_env_full[,-c(1:6,9)]
current_MLL <- current_env_MLL[,-c(1:6,9)]
future_full <- future_env_full[,-c(1:5,8)]
future_MLL <- future_env_MLL[,-c(1:5,8)]

# scale data
pred_std_full <- decostand(current_full, "standardize")
pred_scaled_full <- scale(pred_std_full)

pred_std_MLL <- decostand(current_MLL, "standardize")
pred_scaled_MLL <- scale(pred_std_MLL)

pred_std_fullF <- decostand(future_full, "standardize")
pred_scaled_fullF <- scale(pred_std_fullF)

pred_std_MLLF <- decostand(future_MLL, "standardize")
pred_scaled_MLLF <- scale(pred_std_MLLF)

# prep as matrix
current_full_mat <- as.matrix(pred_scaled_full)
current_MLL_mat <- as.matrix(pred_scaled_MLL)
future_full_mat <- as.matrix(pred_scaled_fullF)
future_MLL_mat <- as.matrix(pred_scaled_MLLF)

# and convert to env format
write.env(current_full_mat, "seascape_data/current_full_env.env")
write.env(current_MLL_mat, "seascape_data/current_MLL_env.env")
write.env(future_full_mat, "seascape_data/future_full_env.env")
write.env(future_MLL_mat, "seascape_data/future_MLL_env.env")
#################

#########################
## LEA OFFSET ANALYSIS ##
#########################

## SEASCAPE STUDY
#################
# I'm going to use lfmm to do this initially, will add other methods later
# read in current env data
current_full_env <- read.env("seascape_data/current_full_env.env")
current_MLL_env <- read.env("seascape_data/current_MLL_env.env")

# double check genomic data
gen_full
gen_MLL

# fitting model
mod_lfmm_full <- lfmm2(input = gen_full,
                       env = current_full_env, 
                       K = 242) # using number of MLLs - should I keep K the same for both datasets?
mod_lfmm_full_K12 <- lfmm2(input = gen_full,
                           env = current_full_env,
                           K = 12)
mod_lfmm_MLL <- lfmm2(input = gen_MLL,
                      env = current_MLL_env,
                      K = 12) # using snmf estimate - should I use number MLLs?

# perform test
pv_lfmm_full <- lfmm2.test(object = mod_lfmm_full,
                       input = gen_full,
                       env = current_full_env, 
                       full = TRUE)
pv_lfmm_full_K12 <- lfmm2.test(object = mod_lfmm_full_K12,
                               input = gen_full,
                               env = current_full_env, 
                               full = TRUE)
pv_lfmm_MLL <- lfmm2.test(object = mod_lfmm_MLL,
                          input = gen_MLL,
                          env = current_MLL_env, 
                          full = TRUE)

# pvalues
pv_full <- pv_lfmm_full$pvalues
pv_full_K12 <- pv_lfmm_full_K12$pvalues
pv_MLL <- pv_lfmm_MLL$pvalues

# qvalues
qv_full <- qvalue(pv_full, fdr.level = 0.01)
qv_full_K12 <- qvalue(pv_full_K12, fdr.level = 0.01)
qv_MLL <- qvalue(pv_MLL, fdr.level = 0.01)

# candidates
candidates_full <- which(qv_full$significant)
candidates_pvals_full <- (pv_full)[candidates_full]

candidates_full_K12 <- which(qv_full_K12$significant)
candidates_pvals_full_K12 <- (pv_full_K12)[candidates_full_K12]

candidates_MLL <- which(qv_MLL$significant)
candidates_pvals_MLL <- (pv_MLL)[candidates_MLL]

# using FDR to filter out false positives
length(which(qv_full$qvalues < 0.01)) # how many SNPs have an FDR < 1%? 176
length(which(qv_full_K12$qvalues < 0.01)) # how many SNPs have an FDR < 1%? 447 - big diff from K = 242
length(which(qv_MLL$qvalues < 0.01)) # how many SNPs have an FDR < 1%? 131

# identify which SNPs have FDR less than 1%
FDR_cands_full <- colnames(gen_full)[which(qv_full$qvalues < 0.01)]
FDR_cands_index_full <- which(qv_full$qvalues < 0.01)

FDR_cands_full_K12 <- colnames(gen_full)[which(qv_full_K12$qvalues < 0.01)]
FDR_cands_index_full_K12 <- which(qv_full_K12$qvalues < 0.01)

FDR_cands_MLL <- colnames(gen_MLL)[which(qv_MLL$qvalues < 0.01)]
FDR_cands_index_MLL <- which(qv_MLL$qvalues < 0.01)

# df for candidates
cand_df_full <- data.frame(index = colnames(gen_full), pos = 1:length(pv_full),
                           pvals = pv_full, qvals = qv_full$qvalues)
cand_df_full_K12 <- data.frame(index = colnames(gen_full), pos = 1:length(pv_full_K12),
                           pvals = pv_full_K12, qvals = qv_full_K12$qvalues)
cand_df_MLL <- data.frame(index = colnames(gen_MLL), pos = 1:length(pv_MLL),
                          pvals = pv_MLL, qvals = qv_MLL$qvalues)

# overlap in candidates?
# overlap between full with K = 242 vs full with K = 12
which(FDR_cands_full_K12 %in% FDR_cands_full)
length(which(FDR_cands_full_K12 %in% FDR_cands_full)) # 58 SNPs overlap

# overlap between full with K = 12 and MLL
which(FDR_cands_full_K12 %in% FDR_cands_MLL)
length(which(FDR_cands_full_K12 %in% FDR_cands_MLL)) # 12 SNPs overlap

# overlap between full with K = 242 and MLL
which(FDR_cands_full %in% FDR_cands_MLL)
length(which(FDR_cands_full %in% FDR_cands_MLL)) # 3 SNPs overlap

# manhattan plots
# man_full <- ggplot(cand_df_full, aes(x = pos, y = -log10(pvals)))+
#   theme_classic()
plot(cand_df_full$pos, -log10(cand_df_full$pvals), 
     pch = 19, col = "darkslategray4",
     main = "Candidate SNPs Full Dataset",
     xlab = "Position",
     ylab = "-log10(Pvals)")
points(cand_df_full$pos[cand_df_full$qvals < 0.01], 
       -log10(cand_df_full$pvals[cand_df_full$qvals < 0.01]),
       pch = 19,
       col = "hotpink")

plot(cand_df_full_K12$pos, -log10(cand_df_full_K12$pvals), 
     pch = 19, col = "darkslategray4",
     main = "Candidate SNPs Full Dataset",
     xlab = "Position",
     ylab = "-log10(Pvals)")
points(cand_df_full_K12$pos[cand_df_full_K12$qvals < 0.01], 
       -log10(cand_df_full_K12$pvals[cand_df_full_K12$qvals < 0.01]),
       pch = 19,
       col = "hotpink")

plot(cand_df_MLL$pos, -log10(cand_df_MLL$pvals), 
     pch = 19, col = "darkslategray4",
     main = "Candidate SNPs MLL Dataset",
     xlab = "Position",
     ylab = "-log10(Pvals)")
points(cand_df_MLL$pos[cand_df_MLL$qvals < 0.01], 
       -log10(cand_df_MLL$pvals[cand_df_MLL$qvals < 0.01]),
       pch = 19,
       col = "hotpink")
#################

## OFFSETS
##########
# read in future env data
future_full_env <- read.env("seascape_data/future_full_env.env")
future_MLL_env <- read.env("seascape_data/future_MLL_env.env")

# calculate genetic gap
gap_cands_full <- genetic.gap(input = gen_full, 
                         env = current_full_env, 
                         pred.env = future_full_env, 
                         scale = F, 
                         K = 242)

gap_cands_MLL <- genetic.gap(input = gen_MLL, 
                              env = current_MLL_env, 
                              pred.env = future_MLL_env, 
                              scale = F, 
                              K = 12)

# check how well each genetic gap measure correlates with euclidean environmental distance
Delta <- current_full_env - future_full_env
dist_env <- sqrt(rowSums(Delta^2))
plot(dist_env, gap_cands_full$distance,
     main = "Genetic gap vs Environmental Distance - FULL", 
     xlab ="Euclidean distance",  ylab ="sqrt(genetic gap)", 
     cex = .6, col = "darkslategray4", pch = 19)

Delta <- current_MLL_env - future_MLL_env
dist_env <- sqrt(rowSums(Delta^2))
plot(dist_env, gap_cands_MLL$distance,
     main = "Genetic gap vs Environmental Distance - MLL", 
     xlab ="Euclidean distance",  ylab ="sqrt(genetic gap)", 
     cex = .6, col = "darkslategray4", pch = 19)

# Extract genomic offset from genetic gap object
offsets_full <- gap_cands_full$offset
offsets_MLL <- gap_cands_MLL$offset

# summary stats on offsets
summary(offsets_full)
summary(offsets_MLL)

# offset dfs
offset_full_df <- cbind(indiv_full_merge, offsets_full)
offset_MLL_df <- cbind(indiv_MLL_merge, offsets_MLL)

# can now calculate offset by region, site, etc
# average offset by site
# average offset for each population
offset_site_full <- aggregate(offset_full_df$offsets_full, by = list(offset_full_df$Pop), FUN = mean)
colnames(offset_site_full) <- c("site","offset_full")

offset_site_MLL <- aggregate(offset_MLL_df$offsets_MLL, by = list(offset_MLL_df$Pop), FUN = mean)
colnames(offset_site_MLL) <- c("site","offset_MLL")

# order populations by average offset
(ordered_site_full <- offset_site_full[order(offset_site_full$offset_full),])
ordered_site_full$order_full <- 1:39
(ordered_site_MLL <- offset_site_MLL[order(offset_site_MLL$offset_MLL),])
ordered_site_MLL$order_MLL <- 1:39

# put this together into single df to visualize
offset_sites <- merge(ordered_site_full, ordered_site_MLL, by = "site")
offset_sites <- offset_sites[c("site","offset_full","offset_MLL","order_full","order_MLL")]

# save this and raw offsets as .csv
write.csv(offset_full_df, "offset_full_df.csv")
write.csv(offset_MLL_df, "offset_MLL_df.csv")
write.csv(offset_sites, "offset_sites.csv")
##########

## CORRELATING OFFSETS
######################
# need to first ensure datasets are comparable
# use only indivs in both datasets
offset_full_trimmed_df <- offset_full_df[offset_full_df$Indv_no %in% offset_MLL_df$Indv_no,]
dim(offset_MLL_df)
dim(offset_full_trimmed_df)

# quick look at basic correlation
cor(offset_full_trimmed_df$offsets_full, offset_MLL_df$offsets_MLL) # 0.7374169 correlation

# raw data
plot(offset_full_trimmed_df$offsets_full, offset_MLL_df$offsets_MLL, 
     pch = 19,
     main = "Clone-corrected vs Uncorrected Offsets",
     xlab = "Uncorrected Offsets",
     ylab = "Clone-corrected Offsets",
     cex = .6, col = "darkslategray4")
abline(lm(offset_MLL_df$offsets_MLL ~ offset_full_trimmed_df$offsets_full))

# log transform
plot(log(offset_full_trimmed_df$offsets_full), log(offset_MLL_df$offsets_MLL), 
     pch = 19,
     main = "Clone-corrected vs Uncorrected Offsets",
     xlab = "Uncorrected Offsets",
     ylab = "Clone-corrected Offsets",
     cex = .6, col = "darkslategray4")
abline(lm(log(offset_MLL_df$offsets_MLL) ~ log(offset_full_trimmed_df$offsets_full)))

# zoom in on bottom left corner - GRFP PLOT
plot(offset_full_trimmed_df$offsets_full, offset_MLL_df$offsets_MLL, 
     pch = 19,
     main = "Clone-corrected vs Uncorrected Offsets",
     xlab = "Uncorrected Offsets",
     ylab = "Clone-corrected Offsets",
     xlim = c(0,0.003),
     ylim = c(0,0.003),
     cex = 1, col = "darkslategray4")
abline(lm(offset_MLL_df$offsets_MLL[offset_MLL_df$offsets_MLL < .004] ~ offset_full_trimmed_df$offsets_full[offset_MLL_df$offsets_MLL < .004]))

# zoom in on bottom left corner, log transform
plot(log(offset_full_trimmed_df$offsets_full[offset_MLL_df$offsets_MLL < .004]), log(offset_MLL_df$offsets_MLL[offset_MLL_df$offsets_MLL < .004]), 
     pch = 19,
     main = "Full Offset vs MLL Offsets",
     xlab = "Full Offsets",
     ylab = "MLL Offsets",
     cex = .6, col = "darkslategray4")
abline(lm(log(offset_MLL_df$offsets_MLL[offset_MLL_df$offsets_MLL < .004]) ~ log(offset_full_trimmed_df$offsets_full[offset_MLL_df$offsets_MLL < .004])))

# without removing individuals, what is the corr btwn site level offsets?
cor(offset_sites$order_full, offset_sites$order_MLL) # 0.598583 for site orders
cor(offset_sites$offset_full, offset_sites$offset_MLL) # 0.9281348 for site averages

# raw by site
plot(offset_sites$offset_full, offset_sites$offset_MLL, 
     pch = 19,
     main = "Clone-corrected vs Uncorrected Offsets by Site",
     xlab = "Uncorrected Offsets",
     ylab = "Clone-corrected Offsets",
     cex = .6, col = "darkslategray4")
abline(lm(offset_sites$offset_MLL ~ offset_sites$offset_full))

# log transformed by site
plot(log(offset_sites$offset_full), log(offset_sites$offset_MLL), 
     pch = 19,
     main = "Clone-corrected vs Uncorrected Offsets by Site",
     xlab = "Uncorrected Offsets",
     ylab = "Clone-corrected Offsets",
     cex = .6, col = "darkslategray4")
abline(lm(log(offset_sites$offset_MLL) ~ log(offset_sites$offset_full)))

cor(log(offset_sites$offset_full), log(offset_sites$offset_MLL))

# remove outlier sites with exceptionally high GO
offset_sites_trim <- offset_sites[!(offset_sites$site == "TAL"| offset_sites$site == "HAN" | offset_sites$site == "ING"),]

plot(offset_sites_trim$offset_full,
     offset_sites_trim$offset_MLL, 
     pch = 19,
     main = "Clone-corrected vs Uncorrected Offsets by Site",
     xlab = "Uncorrected Offsets",
     ylab = "Clone-corrected Offsets",
     cex = .6, col = "darkslategray4")
abline(lm(offset_sites_trim$offset_MLL ~ offset_sites_trim$offset_full))
######################

# MAP OFFSETS
#############
# THIS DOESN'T WORK ON CLUSTER, DONE ON LOCAL COMPUTER
# put together full df
site_info <- sites[,2:6]
offset_sites
site_complete <- merge(site_info, offset_sites, by = "site")

# set s2 false
sf::sf_use_s2(FALSE) # sf package doesn't work

# set object
world <- ne_countries(scale = "medium", returnclass = "sf")
world_crop <- sf::st_crop(world, c(xmin =.3 , xmax = 33.5, ymin = 53.1, ymax = 66.1))

# full df offset
map_full <- ggplot(data = world) +
  theme_bw()+
  geom_sf(data = world_crop, fill = 'antiquewhite1') +
  geom_point(data = site_complete, aes(x=long, y=lat, 
                                       fill = log(offset_full)), 
             shape = 21, size = 4)+
  scale_colour_gradient(low = "green", high = "red", space = "Lab")+
  theme(plot.title = element_text(size = 24), panel.grid.major = element_line(color = "aliceblue"),
        panel.background = element_rect(fill = "aliceblue"), legend.position = 'right')
(map_full)

map_MLL <- ggplot(data = world) +
  theme_bw()+
  geom_sf(data = world_crop, fill = 'antiquewhite1') +
  geom_point(data = site_complete, aes(x=long, y=lat, 
                                       fill = log(offset_MLL)), 
             shape = 21, size = 4)+
  scale_fill_gradient(low = "green", high = "red", space = "Lab")+
  theme(plot.title = element_text(size = 24), panel.grid.major = element_line(color = "aliceblue"),
        panel.background = element_rect(fill = "aliceblue"), legend.position = 'right')+
  ggtitle("Genomic Offset of Baltic Eelgrass")
(map_MLL)
#############

#########################
## RDA OFFSET ANALYSIS ##
#########################

## RDA seascape
###############
# read in imputed geno format data
geno_imp_full <- read.geno("seascape_data/gen_full_rm_imp_mat.geno")
geno_imp_mll <- read.geno("seascape_data/gen_MLL_imp_mat.geno")

# missing data has already been imputed, trimmed for MAF

# check env data
curr_full <- as.data.frame(current_full_env)
curr_mll <- as.data.frame(current_MLL_env)
fut_full <- as.data.frame(future_full_env)
fut_mll <- as.data.frame(future_MLL_env)

colnames(curr_full) <- colnames(curr_mll) <- colnames(fut_full) <- 
  colnames(fut_mll) <- c("sal","sstmax","sstmin","sstrange")

# approximate pop structure using pca
pca_full <- rda(geno_imp_full, scale=T)
pca_mll <- rda(geno_imp_mll, scale=T)
screeplot(pca_full, type = "barplot", npcs=10, main="PCA Eigenvalues")
screeplot(pca_mll, type = "barplot", npcs=10, main="PCA Eigenvalues")

# keep first three pcs for downstream analysis in both cases


# putting together all the environmental data
current_full <- cbind(pred_scaled_full,)
pred_scaled_full
sites_environment
###############