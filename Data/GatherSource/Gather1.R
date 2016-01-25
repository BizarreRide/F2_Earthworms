##############
# Gather1 - Create 1st Order data frame
# data is cleaned up
# Earthworm Abundances and Biomass are summed ogether for the respective Functional groups 
# Biodiversity Indices are calculated: Total Abundance, Species Richness, Shannon, Pielou
# Quentin Schorpp
# Updated 15 April 2015
##############


#setwd("D:/Quentin_Schorpp/Arbeitsprozess/git_repositories/F2_Earthworms")
# KnitHTML brings an error message if not included; but for repsoducibility the full path should not be given!

# Load data
#data <- read.delim("Data/F2_EW_Total.txt", header=TRUE)
data <- read.delim("Data/F2_EW_Total_NewCLimate.txt", header=TRUE)



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
data$location <- plyr::revalue(data$location, c("Ronneberg" = "Ronnenberg"))
# Create variable d(delta)samcam as the time difference of each sampling date to the first date, that each field was sampled; first date is Zero for all fields. ####
beginning <- strptime(rep(("01.01.2012"),60),"%d.%m.%Y")
data$date[1:60]
data$dsamcam <- with(data, round(c(date[1:60]-beginning,date[61:120]-beginning,date[121:180]-beginning),0))

# create covariate for offset, area = area of the pit where the soil-core was taken from
data$area <- rep(0.25,180)
summary(data)
str(data)
dim(data)
#_________________________________________________



## Functional Groups ####
#________________________________________________

# calculate soil functional groups, abundances
data$anc <- rowSums(data[,c("LTR", "ALO")])
data$endo <- rowSums(data[,c("ARO", "ACA", "ACH", "END","OCY","OLA")])
#data$endo2 <- rowSums(data[,c("ARO", "ACA", "ACH","OCY","OLA")])
data$epi <- rowSums(data[, c("LRU", "LCA")])


# calculate soil functional groups, biomass
data$anc.bm <- rowSums(data[,c("LTR_BM", "ALO_BM")])
data$endo.bm <- rowSums(data[,c("ARO_BM", "ACA_BM", "ACH_BM", "END_BM","OCY_BM","OLA_BM")])
data$epi.bm <- rowSums(data[, c("LRU_BM", "LCA_BM")])
#_________________________________________________



## Calculate Biodiversity Indices ####
#______________________________________________________________________________
spe <- data[,which(colnames(data)=="ACA"):which(colnames(data)=="UNK")]
bm <- data[,which(colnames(data)=="ACA_BM"):which(colnames(data)=="UNK_BM")]

# Total Abundances
data$N <- rowSums(spe)

# Total Biomass
data$N.bm <- rowSums(bm)

# Species Richness
data$SR <- rowSums(spe[,-c(12,9,5)] >0)

# Shannon entropy
data$H <- vegan::diversity(spe[,-c(12,9,5)])

# Pielou Evenness
data$J <- data$H/log(data$SR)
#_______________________________________________________________________________


# Separate abundance data of endogeic earthworms into adult and juveniles
# Separation took place a priori with an Excel Pivot-table

EndoAdult<- read.delim("Data/F2_EW_EndoAdult.txt")

EndoAdult$OCY <- EndoAdult$OCYa + EndoAdult$OCYj 

data$endad <- rowSums(EndoAdult[,c("ACAad", "ACHad", "AROad", "ENDad","OCY","OLA")])


# Separate biomass data of endogeic earthworms into adult and juveniles
# Separation took place a priori with an Excel Pivot-table

EndoAdult_Bm <- read.delim("Data/F2_EW_EndoAdult_Bm.txt")

EndoAdult_Bm$OCY <- EndoAdult_Bm$OCYa + EndoAdult_Bm$OCYj 

data$endad.bm <- rowSums(EndoAdult_Bm[,c("ACAad", "ACHad", "AROad", "ENDad","OCY","OLA")])

rm(EndoAdult, EndoAdult_Bm)

#_______________________________________________________________________________


# Separate abundance data of endogeic earthworms into adult and juveniles
# Separation took place a priori with an Excel Pivot-table

AncAdult<- read.delim("Data/F2_EW_AncAdult.txt")

data$ancad <- rowSums(AncAdult[,c("ALOad", "LTRad")])


# Separate biomass data of endogeic earthworms into adult and juveniles
# Separation took place a priori with an Excel Pivot-table

AncAdult_Bm <- read.delim("Data/F2_EW_AncAdult_Bm.txt")

data$ancad.bm <- rowSums(AncAdult_Bm[,c("ALOad", "LTRad")])

rm(AncAdult, AncAdult_Bm)

#_______________________________________________________________________________



##################
# Wassergehalte
# Berechnung des gravimetrischen Wassergehalts
# Quentin Schorpp
# 29.04.2015
#################

# load data
moisture <- read.delim("Data/WC_Roh_Gesamt.txt")

# Id to maintain row order
moisture$superID <- rep(1:62,each=5)

# create a factor to be able to correct for double measurements from Ronneberg (2  Samplingdates)
moisture$f.aid  <- factor(c(rep(1,300),rep(2,10)))

# Create factor for building averages
moisture$f.scfield <- interaction(moisture$SamplingCampaign, moisture$field, moisture$f.aid, drop=TRUE)


## Claculation of Water Content ####

# Subtract Tara from dry weights
moisture$TGNew  <- moisture$TG - moisture$Tara

# Calculate Water content as differnce between fresh weight and dry weight
moisture$water <- moisture$FG - moisture$TGNew

# Claculate gravimetric Water content for each sample
moisture$gravimetric <- (moisture$water/moisture$TGNew)*100

#####

# Variance Inspection
  # plot(moisture$f.scfield, moisture$water)
  # scfield.variance <- tapply(moisture$water,moisture$f.scfield, var)
  # points(scfield.variance, col="red", pch=19)

# Due to Sabines lazy measurements, all water contents from spring have high variance

# detect outlier in gravimetric water contents
  # plot(moisture$f.scfield, moisture$gravimetric)

# eliminate outliers
moisture.outlier <- subset(moisture, gravimetric<40)
  # with(moisture.outlier, plot(f.scfield, gravimetric)) # make a new plot

# samples with largest variance
gravimetric.variance <- with(moisture.outlier, tapply(gravimetric,f.scfield, var))
  # plot(gravimetric.variance)
grvar <- which(gravimetric.variance>5)
grvar1 <- names(grvar)
moisture.outlier[moisture.outlier$f.scfield%in%grvar1,c("field", "FG", "TGNew", "water", "gravimetric")] # large varaince seems to be due to real varaince, no bias
# It is asumed to be due to true heterogeneity of the fields, hence not corrected

# Calculate average gravimetric water content for field ad sampling campaign
mc <- aggregate(gravimetric ~ f.scfield + superID, moisture.outlier,mean)

library(splitstackshape)
mc <- cSplit(mc, 1, ".")
setnames(mc,3:5,c("samcam","field","f.aid"))


# create data for earthworms only
mc.ew <- subset(mc,samcam!="F12", drop=TRUE)
mc.ew$samcam <- factor(mc.ew$samcam)
mc.ew <- mc.ew[-c(45,42),]
mc.ew$superID <- c(1:41,43,44,42,45)
setorder(mc.ew,superID)

# repeat each row 4 times to fit to 1st order earthworm data
mc.ew <- mc.ew[rep(seq_len(nrow(mc.ew)), each=4),]
mc.ew$ID <- c(1:180) # new ID, similar to "data"

mc.ew$field <- plyr::revalue(mc.ew$field, c(DO_MaisI = "DO_Mais", RO_MaisI="RO_Mais")) # equal factor levels; optional

# switch the values
data$mc <- mc.ew$gravimetric

# clean up:
rm(mc, mc.ew, moisture, moisture.outlier, grvar, grvar1, gravimetric.variance)


