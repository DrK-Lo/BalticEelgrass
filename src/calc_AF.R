### CALCULATING PER POPULATION ALLELE FREQUENCIES

## SET UP
#########
# set working directory
setwd("~/Eelgrass/BalticEelgrass/data")

# load libraries
library(vcfR)
library(LEA)
library(tidyr)
library(ggplot2)
#########


## PREP DATA
############
# read in vcf data
eelgrass_vcf_MLL <- read.vcfR(
  "seascape_data/zostera_230619_seascape_miss25_sorted_95_MLL.recode.vcf")
eelgrass_vcf_full <- read.vcfR(
  "seascape_data/zostera_230619_seascape_miss25_sorted_95_no_reps.recode.vcf")

# read in lfmm data
eelgrass_lfmm_MLL <- read.lfmm("seascape_data/zostera_230619_seascape_miss25_sorted_95_MLL.recode.lfmm")
eelgrass_lfmm_full <- read.lfmm("seascape_data/zostera_230619_seascape_miss25_sorted_95_no_reps.recode.lfmm")

# convert lfmm data to dataframe
eel_full_df <- as.data.frame(eelgrass_lfmm_full)
eel_MLL_df <- as.data.frame(eelgrass_lfmm_MLL)

# rotate df so rows are now columns (we want individuals as columns, alleles as rows)
eel_full_df <- as.matrix(t(eel_full_df))
eel_MLL_df <- as.matrix(t(eel_MLL_df))

# replace 9/-9 with NA
eel_full_df[eel_full_df == 9] <- NA
eel_MLL_df[eel_MLL_df == 9] <- NA
# check it worked
sum((eel_full_df == 9)) # returns NA - good
sum((eel_MLL_df == 9)) # returns NA - good

# get vector of populations
pop_vector_full <- indiv_list_full$Pop
pop_vector_MLL <- indiv_list_MLL$Pop

# add pop data to genomic data
rownames(eel_full_df) <- pop_vector_full
rownames(eel_MLL_df) <- pop_vector_MLL
############


## EXTRACT POSITION INFO
########################
# quick check of data
eelgrass_vcf_full
eelgrass_vcf_MLL

# use vcfR to extract positional information
full_fix <- as.data.frame(getFIX(eelgrass_vcf_full))
MLL_fix <- as.data.frame(getFIX(eelgrass_vcf_MLL))

# which SNPs are in both datasets?
dim(MLL_fix[MLL_fix$POS %in% full_fix$POS,]) # all MLL SNPs are in full

# filter out full SNPs that aren't in MLL dataset
full_fix_rm <- full_fix[full_fix$POS %in% MLL_fix$POS,]

# are these the same?
which(full_fix_rm != MLL_fix) # yes

# how will I exclude these SNPs from whole vcf file?
rm <- which(!(full_fix$POS %in% MLL_fix$POS)) # these are indices of the SNPs to remove

# reduce eel_full_df file to match eel_MLL_df file
eel_full_df_rm <- eel_full_df[-rm,]

# sanity check
dim(eel_MLL_df)
dim(eel_full_df_rm)
# looks good

# write the trimmed data back to lfmm format
write.lfmm(t(eel_full_df_rm), "seascape_data/eel_full_rm.lfmm")
########################


## CALCULATE AFs
################
# function to calculate allele frequency @ single locus
calc_freq <- function(x){
  a <- sum(x, na.rm=TRUE)
  b <- (2*length(na.omit(x)))
  a/b
}

# function extending this to calculate allele frequencies across whole pop
calcfreq <- function(a, pop){
  tapply(a, pop, calc_freq)
}

# some practice data to check everything works
a <- c(0, 1, 2, 2, NA)
pop1 <- c("A", "A", "B", "B", "C")

amat <- matrix(rep(a,3), ncol=3, byrow=TRUE)
rownames(amat) <- pop1
amat

(freqs_amat <- apply(amat, 2, calcfreq, pop1)) # looks good!)

# let's try it with real data
# full dataset
str(eel_full_df_rm)
freqs_full <- apply(t(eel_full_df_rm), 2, calcfreq, pop_vector_full)
str(freqs_full)
colSums(freqs_full[,1:5])
View(freqs_full[,1:10])

# MLL dataset
freqs_MLL <- apply(t(eel_MLL_df), 2, calcfreq, pop_vector_MLL)
str(eel_MLL_df)
str(freqs_MLL)
colSums(freqs_MLL[,1:5])
View(freqs_MLL[,1:10])

# something weird... freqs_MLL has 392 NA values, while freqs_full has 1
################


## PLOTTING
###########
#line up the loci, make the full matrix the same loci as MLL
plot(freqs_full, freqs_MLL)

# histogram of diff between allele frequencies calculated from MLL vs full
freqs_diff <- freqs_full - freqs_MLL
hist(freqs_diff)
hist(freqs_diff, ylim = c(0,4000))
###########
