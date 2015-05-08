##############
# Gather2 - Create 2nd Order data frame
# The whole dataset is summarized for "hole", 
# samples taken from all holes at one field are summed up and the average is build 
# Abundance and biomass data then is the value per squaremeter [m²];(sqm)
# Biodiversity Indices are calculated new: Total Abundance, Species Richness, Shannon, Pielou
# Quentin Schorpp
# Updated 16 April 2015
##############

# 2nd Order Data

# data2 <- read.table("data/F2_EW_Total.txt", header=T)

data2 <- data[,-c(68:80)] # remove Biodiversity Indices + Functional Groups


# Explanatory variablaes
data21 <- data2[,-c(which(colnames(data2)=="hole"):which(colnames(data2)=="UNK_BM") )]
data21 <- unique(data21[,-1])

# Response variables
data22 <- data2[,c(which(colnames(data2)=="hole"):which(colnames(data2)=="UNK_BM"))]
data22 <- aggregate(.~interaction(data2$field.ID,data2$samcam),data22[,-1], sum)

data2 <- cbind(data21,data22[,-1])


## Functional Groups ####
#________________________________________________

# calculate soil functional groups, abundances
data2$anc <- rowSums(data2[,c("LTR", "ALO")])
data2$endo <- rowSums(data2[,c("ARO", "ACA", "ACH", "END","OCY","OLA")])
data2$epi <- rowSums(data2[, c("LRU", "LCA")])


# calculate soil functional groups, biomass
data2$anc.bm <- rowSums(data2[,c("LTR_BM", "ALO_BM")])
data2$endo.bm <- rowSums(data2[,c("ARO_BM", "ACA_BM", "ACH_BM", "END_BM","OCY_BM","OLA_BM")])
data2$epi.bm <- rowSums(data2[, c("LRU_BM", "LCA_BM")])
#_________________________________________________


## Calculate Biodiversity Indices ####
#______________________________________________________________________________
spe <- data2[,which(colnames(data2)=="ACA"):which(colnames(data2)=="UNK")]
bm <- data2[,which(colnames(data2)=="ACA_BM"):which(colnames(data2)=="UNK_BM")]

# Total Abundances
data2$N <- rowSums(spe)

# Total Biomass
data2$N.bm <- rowSums(bm)

# Species Richness
data2$SR <- rowSums(spe[,-c(12,9,5)] >0)

# Shannon entropy
data2$H <- vegan::diversity(spe[,-c(12,9,5)])

# Pielou Evenness
data2$J <- data2$H/log(data2$SR)
#_______________________________________________________________________________

## Adult vs. Juvenile Endogeics

# Separate abundance data of endogeic earthworms into adult and juveniles
# Separation took place a priori with an Excel Pivot-table

EndoAdult <- read.delim("Data/F2_EW_EndoAdult.txt")

EndoAdult$OCY <- EndoAdult$OCYa + EndoAdult$OCYj 

# Response variables
data23 <- EndoAdult[,c(which(colnames(EndoAdult)=="hole"):which(colnames(EndoAdult)=="OCY"))]
data23 <- aggregate(.~interaction(data$field.ID,data$samcam),data23[,-1], sum)

data2$endad <- rowSums(data23[,c("ACAad", "ACHad", "AROad", "ENDad","OCY","OLA")])

# Separate abundance data of endogeic earthworms into adult and juveniles
# Separation took place a priori with an Excel Pivot-table

EndoAdult_Bm <- read.delim("Data/F2_EW_EndoAdult_Bm.txt")

EndoAdult_Bm$OCY <- EndoAdult_Bm$OCYa + EndoAdult_Bm$OCYj 

# Response variables
data24 <- EndoAdult_Bm[,c(which(colnames(EndoAdult_Bm)=="hole"):which(colnames(EndoAdult_Bm)=="OCY"))]
data24 <- aggregate(.~interaction(data$field.ID,data$samcam),data24[,-1], sum)

data2$endad.bm <- rowSums(data24[,c("ACAad", "ACHad", "AROad", "ENDad","OCY","OLA")])

rm(EndoAdult, data23, EndoAdult_Bm, data24, spe, bm, data21, data22)