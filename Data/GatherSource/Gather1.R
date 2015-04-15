##############
# Gather1 - Create 1st Order data frame
# data is cleaned up
# Earthworm Abundances and Biomass are summed ogether for the respective Functional groups 
# Biodiversity Indices are calculated: Total Abundance, Species Richness, Shannon, Pielou
# Quentin Schorpp
# Updated 15 April 2015
##############


setwd("D:/Quentin_Schorpp/Arbeitsprozess/git_repositories/F2_Earthworms")
# KnitHTML brings an error message if not included; but for repsoducibility the full path should not be given!

# Load data
data <- read.delim("Data/F2_EW_Total.txt", header=TRUE)



## Clean Up Data ####
#________________________________________________

# change variable classes
data$field.ID <- as.factor(data$field.ID)
data$ID <- as.factor(data$ID)
data$hole <- as.factor(data$hole)
data$samcam <- as.factor(data$samcam)
data$age_class <- factor(data$age_class, ordered=TRUE)
data$soil_tillage <- ordered(data$soil_tillage)
data$fertilisation <- ordered(data$fertilisation)
data$weed_control <- ordered(data$weed_control)
data$soil_compaction <- ordered(data$soil_compaction)
data$compaction <- ordered(data$compaction)
data$date <- strptime(data$date, "%d.%m.%Y")

summary(data)
str(data)
dim(data)
#_________________________________________________



## Functional Groups ####
#________________________________________________

# calculate soil functional groups, abundances
data$anc <- rowSums(data[,c("LTR", "ALO")])
data$endo <- rowSums(data[,c("ARO", "ACA", "ACH", "END","OCY","OLA")])
data$epi <- rowSums(data[, c("LRU", "LCA")])


# calculate soil functional groups, biomass
data$anc.bm <- rowSums(data[,c("LTR_BM", "ALO_BM")])
data$endo.bm <- rowSums(data[,c("ARO_BM", "ACA_BM", "ACH_BM", "END_BM","OCY_BM","OLA_BM")])
data$epi.bm <- rowSums(data[, c("LRU_BM", "LCA_BM")])
#_________________________________________________



## Calculate Biodiversity Indices ####
#______________________________________________________________________________
spe <- data[,which(colnames(data)=="ACA"):which(colnames(data)=="UNK")]

# Total Abundances
data$N <- rowSums(spe)

# Species Richness
data$SR <- rowSums(spe[,-c(12,9,5)] >0)

# Shannon entropy
data$H <- vegan::diversity(spe[,-c(12,9,5)])

# Pielou Evenness
data$J <- data$H/log(data$SR)
#_______________________________________________________________________________






