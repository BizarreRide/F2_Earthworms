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

# create a factor to be able to correct for double measurements from Ronneberg
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

plot(moisture$f.scfield, moisture$water)
scfield.variance <- tapply(moisture$water,moisture$f.scfield, var)
points(scfield.variance, col="red")

# detect outlier in gravimetric water contents
plot(moisture$f.scfield, moisture$gravimetric)

# eliminate outliers
moisture.outlier <- subset(moisture, gravimetric<40)
with(moisture.outlier, plot(f.scfield, gravimetric)) # make a new plot

# samples with largest variance
gravimetric.variance <- with(moisture.outlier, tapply(gravimetric,f.scfield, var))
plot(gravimetric.variance)
grvar <- which(gravimetric.variance>5)
grvar1 <- names(grvar)
moisture.outlier[moisture.outlier$f.scfield%in%grvar1,] # large varaince seems to be due to real varaince, no bias


# Calculate average gravimetric water content for field ad sampling campaign
mc <- aggregate(gravimetric ~ f.scfield + superID, moisture.outlier,mean)

library(splitstackshape)
mc <- cSplit(mc, 1, ".")
colnames(mc)[3:5] <- c("samcam","field","f.aid")


# create data for earthworms only
mc.ew <- subset(mc,samcam!="F12", drop=TRUE)
mc.ew$samcam <- factor(mc.ew$samcam)
mc.ew <- mc.ew[-c(45,42),]
mc.ew$superID <- c(1:41,43,44,42,45)
setorder(mc.ew,superID)

# repeat each row 4 times to fit to 1st order earthworm data
mc.ew <- mc.ew[rep(seq_len(nrow(mc.ew)), each=4),]
mc.ew$ID <- c(1:180) # new ID, similar to "data"

mc.ew$field <- revalue(mc.ew$field, c(DO_MaisI = "DO_Mais", RO_MaisI="RO_Mais")) # equal factor levels; optional

# switch the values
data$mc <- mc.ew$gravimetric

# clean up:
rm(mc, mc.ew, moisture, moisture.outlier, grvar, grvar1, scfield.variance, gravimetric.variance)


