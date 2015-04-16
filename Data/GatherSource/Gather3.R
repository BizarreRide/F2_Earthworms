##############
# Gather3 - Create 3rd Order data frame
# The whole dataset is summarized for "hole", 
# samples taken from all fields at all sampling campaigns are summed up and the average is build 
# Abundance and biomass data then is the mean value per squaremeter over the sampling campaigns [m²];(sqm)
# Biodiversity Indices are calculated new: Total Abundance, Species Richness, Shannon, Pielou
# Quentin Schorpp
# Updated 16 April 2015
##############

# 2nd Order Data

options(digits=3)

# data.sqm <- read.table("data/F2_EW_Total.txt", header=T)
data3 <- data2[,-c(66:76)]
data3$soil_tillage  <- as.numeric(data3$soil_tillage)
data3$fertilisation <- as.numeric(data3$fertilisation)
data3$weed_control <- as.numeric(data3$weed_control)
data3$soil_compaction <- as.numeric(data3$soil_compaction)
data3$compaction <- as.numeric(data3$compaction)

# Extract factors
df1 <- data3[sapply(data3,is.factor)]

# Extract non-factors
dfx <- data3[sapply(data3,is.numeric)]


#library(data.table)
#dfx <- data.table(dfx)
#dfy <- dfx[, c(mean = lapply(.SD, mean), sd = lapply(.SD, sd)), by = df1$field.ID]

# Create Averages over sampling campaign
df2 <- aggregate(dfx, list(df1$field.ID), mean) # better use ddply???

# Create standard deviation over sampling campaign
data3.sd <- aggregate(dfx, list(df1$field.ID), sd) # better use ddply???

# Create Standard Error
se <- function(x) sd(x)/sqrt(length(x))
data3.se <- aggregate(dfx, list(df1$field.ID), se) # better use ddply???

# Melt factors to field.ID
df1 <- df1[order(df1$field.ID),]
df1 <- unique(df1[,-c(which(colnames(df1)=="samcam"),
                     which(colnames(df1)=="season"))])

# Create data frame for means
colnames(df2)[1] <- "field.ID"
data3 <- merge(df1,df2,by="field.ID")
data3$field.ID <- as.numeric(data3$field.ID)
data3[order(data3$field.ID),]
rm(df1,df2,dfx)

## Functional Groups ####
#________________________________________________

# calculate soil functional groups, abundances
data3$anc <- rowSums(data3[,c("LTR", "ALO")])
data3$endo <- rowSums(data3[,c("ARO", "ACA", "ACH", "END","OCY","OLA")])
data3$epi <- rowSums(data3[, c("LRU", "LCA")])


# calculate soil functional groups, biomass
data3$anc.bm <- rowSums(data3[,c("LTR_BM", "ALO_BM")])
data3$endo.bm <- rowSums(data3[,c("ARO_BM", "ACA_BM", "ACH_BM", "END_BM","OCY_BM","OLA_BM")])
data3$epi.bm <- rowSums(data3[, c("LRU_BM", "LCA_BM")])
#_________________________________________________


## Calculate Biodiversity Indices ####
#______________________________________________________________________________
spe <- data3[,which(colnames(data3)=="ACA"):which(colnames(data3)=="UNK")]
bm <- data3[,which(colnames(data3)=="ACA_BM"):which(colnames(data3)=="UNK_BM")]

# Total Abundances
data3$N <- rowSums(spe)

# Total Biomass
data3$N.bm <- rowSums(bm)

# Species Richness
data3$SR <- rowSums(spe[,-c(12,9,5)] >0)

# Shannon entropy
data3$H <- vegan::diversity(spe[,-c(12,9,5)])

# Pielou Evenness
data3$J <- data3$H/log(data3$SR)
#_______________________________________________________________________________

rm(spe,bm)
