#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# F2_Eearthworms
# Multivariate analyisis of environmental variables
# Quentin Schorpp
# 29.12.2015
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Data/GatherSource/F2_EW_MakeLikeFile.R")
source("Data/GatherSource/Functions/evplot.R")
source("Data/GatherSource/Functions/cleanplot.pca.R")

# Data processing ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Subsetting 
## Response variables are the soil properties of the plots
df.rsp1 <- data[,c("clay","sand","pH","mc","cn")]
df.rsp2 <- data2[,c("clay","sand","pH","mc","cn")]
df.rsp3 <- data3[,c("clay","sand","pH","mc","cn")]

## Explanatory variable are abundances of the functional groups
df.exp1 <- data[,c("ancad", "endad", "epi")]
df.exp2 <- data2[,c("ancad", "endad", "epi")]
df.exp3 <- data3[,c("ancad", "endad", "epi")]

df.groups1 <- data[,c("age_class","samcam","field.ID", "hole")]
df.groups2 <- data2[,c("age_class","samcam","field.ID")]
df.groups3 <- data3[,c("age_class","field.ID")]

# Transformations
df.rsp1.log <- log(df.rsp1)                   # log transformation
df.rsp1.z <- scale(df.rsp1.log)               # z transformation
minShift <- function(x) {x + abs(min(x))}
df.rsp1.z <- apply(df.rsp1.z, 2, minShift) # eliminate negative values

# Plot Aid
## Background color of symbols
p_bg <- grey.colors(5, start = 0.8, end = 0, gamma = 2.2, alpha = NULL)
p_bg <- rep(rep(p_bg, each=12),3)
## shape of symbols
p_shp <- c(rep(rep(c(24,25,23,22,21), each=12),3))
p_shp.legend <- c(2,6,5,0,1)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# df.rsp1 <- df.rsp3
# df.exp1 <- df.exp3

# Test for multivariate Normality ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(mvnormtest)
mnorm.rsp1 <- t(df.rsp1.z)
mshapiro.test(mnorm.rsp1) # No multiple NV. ALso no multivariate normality?

library(MVN)
mardiaTest(df.rsp1)
hzTest(df.rsp1)
mardiaTest(df.rsp1.z)
hzTest(df.rsp1.z)        # Data are not multivariate normal!
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Non metric multidimensional scaling ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Prepare overlay of explanatory variables
dca.rsp1 <- decorana(df.rsp1)
ef <- envfit(dca.rsp1, df.exp1, perm = 1000)
ef
scores(ef, "vectors")
p.dca <- plot(dca.rsp1, type="points")
points(p.dca, "sites", pch=25, bg=df.groups1$age_class, cex=1.0)
identify(p.dca,"sp", labels=names(df.rsp1),cex=1.0)
plot(ef, p.max = 0.05, col="blue") # not working



# NMDS
nmds.rsp1 <- metaMDS(df.rsp1,"euc", k=2, trymax=1000)
nmds.rsp1 <- metaMDS(df.rsp1,"euc", k=2, trymax=1000, previous.best = nmds.rsp1)
nmds.rsp1
envfit(p.nmds.rsp1, df.exp1, perm=999)

# NMDS plot
p.nmds.rsp1 <- plot(nmds.rsp1, type="n")
points(p.nmds.rsp1, "sites", pch = p_shp, bg=p_bg, cex=0.8)
text(nmds.rsp1, display="species")
plot(envfit(p.nmds.rsp1, df.exp1), add = TRUE, col="blue")



# Validate NMDS
par(mfrow=c(1,2))
stressplot(nmds.rsp1, main = "Shepard plot")
gof <- goodness(nmds.rsp1); gof
p.nmds.rsp1.gof <- plot(nmds.rsp1, type="t", main = "Goodness of fit") 
points(nmds.rsp1, display="sites", cex=2*gof/mean(gof))

# Rotate NMDS plot
nmds.rsp1.rtt <- with(df.rsp1, MDSrotate(nmds.rsp1,clay))

# Final NMDS plot
a = 7.3
b = 7.3
c = 9
d = 9.2

p <- length(fam.db.repmes$CA$eig)


# Version 1 mit ordiplot()
par(mfrow=c(1,1),
    mar=c(2,2,2,0.2),
    oma=c(0,0,0,0),
    mgp=c(01,0.2,0),
    cex.lab=0.5,
    cex.axis=0.5,
    cex.main=0.5,
    lwd=0.3,
    tcl=NA,
    pin=c(width=a/2.54, height=b/2.54))


p.nmds.rsp1.rtt <- plot(nmds.rsp1.rtt, type="n", xlim=c(-0.15, 0.15))
points(nmds.rsp1.rtt, display="sites", pch = p_shp, bg = p_bg)
#text(nmds.rsp1.rtt, display="species")
plot(envfit(p.nmds.rsp1.rtt, df.rsp1), add = TRUE, col="black") 
plot(envfit(p.nmds.rsp1.rtt, df.exp1), add = TRUE, col="black") 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# PCA ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Too hard to decide which PCA - scaling definitions are to choose!!!!!!!

pca.rsp1 <- rda(df.rsp1, scale=F) # Inertia is the sum of all species variances.  The eigenvalues sum up to total inertia.  In other
                                  # words,  they  each explain a  certain  proportion  of  total  variance.
pca.rsp2 <- rda(df.rsp1, scale=T) # Now inertia is correlation, and the correlation of a variable with itself is
                                  # one.  Thus the total inertia is equal to the number of variables (species)
summary(pca.rsp1)
# The percentage explained by the rst axis decreased from the previous pca.
# This is natural,  since previously we needed to explain only the abundant  
# species  with  high  variances,  but  now  we  have  to  explain  all
# species equally.  We should not look blindly at percentages, but the result
# we get.

par(mfrow=c(2,3))
p1 <- plot(pca.rsp1, scaling = -1, main="scale=F, scaling=-1")
p1 <- plot(pca.rsp1, scaling = -2, main="scale=F, scaling=-2")
p1 <- plot(pca.rsp1, scaling = -3, main="scale=F, scaling=-3")
p1 <- plot(pca.rsp1, scaling = 1, main="scale=F, scaling=1")
p1 <- plot(pca.rsp1, scaling = 2, main="scale=F, scaling=2")
p1 <- plot(pca.rsp1, scaling = 3, main="scale=F, scaling=3")

par(mfrow=c(2,3))
p1 <- plot(pca.rsp2, scaling = -1, main="scale=T, scaling=-1")
p1 <- plot(pca.rsp2, scaling = -2, main="scale=T, scaling=-2")
p1 <- plot(pca.rsp2, scaling = -3, main="scale=T, scaling=-3")
p1 <- plot(pca.rsp2, scaling = 1, main="scale=T, scaling=1")
p1 <- plot(pca.rsp2, scaling = 2, main="scale=T, scaling=2")
p1 <- plot(pca.rsp2, scaling = 3, main="scale=T, scaling=3")
# For this graph we specifed scaling = -1. The  negative  values  mean  that  species
# scores are divided by the species standard deviations so that abundant and scarce species 
# will be approximately as far away from the origin.
# We already saw an example of scaling = 3 or symmetric scaling in pca.
# The other two integers mean that either species are weighted averages of
# sites (2) or sites  are  weighted  averages  of  species (1). When  we  take
# weighted  averages,  the  range  of  averages  shrinks  from  the  original  values.
# The shrinkage factor is equal to the eigenvalue of CA, which has a theoretical maximum of 1.


# Chossed to scale the environmental soil variables in prior to PCA, 
# since they were all measured on different scales

p1 <- plot(pca.rsp2)
points(p1, "sites", pch=25, bg=df.groups1$age_class, cex=0.7)
ev <- pca.rsp2$CA$eig
x11()
evplot(ev)
screeplot(pca.rsp2, bstick=T, type="lines")
p1 <- cleanplot.pca(pca.rsp2, ahead=0, point=T)
points(p1, "sites", pch=25, bg=df.groups1$age_class, cex=1.0)


# Final PCA plot ####

a = 7.3
b = 7.3
c = 9
d = 9.2

p <- length(fam.db.repmes$CA$eig)


# Version 1 mit ordiplot()
par(mfrow=c(1,2),
    mar=c(2,2,2,0.2),
    oma=c(0,0,0,0),
    mgp=c(01,0.2,0),
    cex.lab=0.5,
    cex.axis=0.5,
    cex.main=0.5,
    lwd=0.3,
    tcl=NA,
    pin=c(width=a/2.54, height=b/2.54))


library(BiodiversityR)
plot1 <- ordiplot(pca.rsp2, choices=c(1,2),type="n", scaling=1, display=c("wa","sp"), main="Scaling 1 wa-scores") #, display=c("sp","lc","cn")
ordiequilibriumcircle(pca.rsp2,plot1)
points(plot1, "sites", pch=p_shp, bg=p_bg, cex=1.0)
identify(plot1,"sp", labels=names(df.rsp1), cex=1.0)

plot2 <- ordiplot(pca.rsp2,choices=c(1,2), type="n", scaling=2, display=c("wa","sp","cn"), main="Scaling 2 wa-scores")
points(plot2, "sites", pch=p_shp, bg=p_bg, cex=1.0)
points(plot2, "species", pch=4, cex=1.0)
legend("bottomleft", legend=c("Sp_O", "SP_I2", "SP_I1", "Sp_Y", "Cm"), text.col="black", pch=pca.p_shp.legend, cex=0.4, pt.cex=0.7)


# Axis 1 + 2
#************
# Scaling 2 
ordiellipse(plot2, df.groups1$age_class, scaling=2, kind="se", conf=0.95)
points(plot2, "centroids", pch=3, bg="black", cex=0.8) 
points(plot2, "sites", pch=shape_rda2, bg="black", cex=0.8) 
#points(plot2, "species", pch=4, cex=0.8) 
text(plot2, "species",label=abv, pch=4, cex=0.4) 
sit.sc2 <- scores(pca.rsp1, display="wa", choices=c(1:p))
spe.sc2 <- scores(pca.rsp1, display="sp", choices=c(1:p))
arrows(0, 0, spe.sc2[,1], spe.sc2[,2], length=0.07, angle=20, col="black")
#identify(plot2, "sp", label=abv,cex=0.5)
#legend("bottomleft", legend=c("Sp_O", "SP_I2", "SP_I1", "Sp_Y", "Cm"), text.col="black", pch=shape_rda1, cex=0.4, pt.cex=0.5)
#dev.copy2pdf(file="Results/dbRDA12sc2.pdf", width=c/2.54, height=d/2.54, useDingbats=F, out.type = "pdf")# 

# Scaling 1
#pcacircle(pca.rsp1)
ordiellipse(plot1, df.groups1$age_class, scaling=2, kind="se", conf=0.95)
points(plot1, "centroids", pch=3, bg="black", cex=0.8) 
#points(plot1, "sites", pch=shape_rda2, bg="black", cex=0.8) 
#points(plot1, "species", pch=4, cex=0.8) 
text(plot1, "species",label=abv, pch=4, cex=0.4) 
#legend("bottomleft", legend=c("Sp_O", "SP_I2", "SP_I1", "Sp_Y", "Cm"), text.col="black", pch=shape_rda1, cex=0.4, pt.cex=0.5)
#identify(plot1,"sp", labels=names(fam.usc), cex=1.0)
#dev.copy2pdf(file="Results/dbRDA12sc1.pdf", width=c/2.54, height=d/2.54, useDingbats=F, out.type = "pdf")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Correspondence Analysis ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Not that the data submitted ot ca must be frequencies or frequency-like, 
# dimensionally homogeneous and non-negative; that is the case of species counts 
# or presence-absence data
# Not for environmental data

ca.rsp1 <- cca(df.rsp1)
plot(ca.rsp1)
chisq.test(df.rsp1/sum(df.rsp1))
plot(ca.rsp1, scaling=1)
plot(dca.rsp1, scaling=1)

# Use DCA when arch effect was observed (in plot?)
dca.rsp1 <- decorana(df.rsp1)
# It is often said that if the axis length is shorter than two units, 
# the data are linear, and pca should be used.
# This is only folklore and not based on research which shows that
# ca is at least as good as pca with short gradients, and usually better.

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




































