#%%%%%%%%%%%%%%%%%%%%
# F2 Earthworms
# Multivariate Analyses
# Exploratory Data Analysis
#%%%%%%%%%%%%%%%%%%%%


# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Data/GatherSource/F2_EW_MakeLikeFile.R")
source("Data/GatherSource/Functions/panelutils.R")
source("Data/GatherSource/Functions/evplot.R")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Data Processing ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Age Class as non-orderd factor
data$age_class <- as.factor(as.character(data$age_class))
data2$age_class <- as.factor(as.character(data2$age_class))
data3$age_class <- as.factor(as.character(data3$age_class))

# Subsetting 
#*************************************************************
## Mulivaraite Response is species abundances 
df.spe1 <- data[,c("ACA", "ARO", "ACH", "OCY", "OLA", "LRU","LCA", "LTR", "ALO")]
df.spe2 <- data2[,c("ACA", "ARO", "ACH", "OCY", "OLA", "LRU","LCA", "LTR", "ALO")]
df.spe3 <- data3[,c("ACA", "ARO", "ACH", "OCY", "OLA", "LRU","LCA", "LTR", "ALO")]
rownames(df.spe2) <- 1:length(rownames(df.spe2))

### Delete all-zero samples
df.spe <- df.spe2[-which(rowSums(df.spe2)==0),]
#df.spe <- df.spe[df.groups$age_class != "A_Cm",]


## Explanatory variables 
#### Categorical variables
df.groups1 <- data[,c("age_class","samcam","field.ID", "hole")]
df.groups2 <- data2[,c("age_class","samcam","field.ID")]
df.groups3 <- data3[,c("age_class","field.ID")]

df.groups <- df.groups2[-which(rowSums(df.spe2)==0),]

#### Soil and CLimate Variables
df.exp1 <- data[,c("age","clay","sand","pH","mc","cn","ats1","hum1", "rad1", "prec1")]
df.exp2 <- data2[,c("age","clay","sand","pH","mc","cn","ats1", "hum1", "rad1", "prec1")]
df.exp3 <- data3[,c("age","clay","sand","pH","mc","cn","ats1", "hum1", "rad1", "prec1")]

df.exp <- df.exp2[-which(rowSums(df.spe2)==0),]

# df.soil2 <- data2[,c("clay","sand","pH","mc","cn")]
# df.climate2 <- data2[,c("ats1", "ata1", "atb1","hum1", "rad1", "prec1")]

#### Biodiversity Variables
df.bdv1 <- data[,c("SR", "H")]
df.bdv2 <- data2[,c("SR", "H")]
df.bdv3 <- data3[,c("SR", "H")]

df.bdv <- data2[-which(rowSums(df.spe2)==0),c("SR", "H")]



# Merged explanatory sets
#data2
df.env <- cbind(df.exp[,c("age","clay","pH","mc","cn","ats1","prec1")],
                df.groups[,c("age_class","samcam")], 
                df.bdv)
#df.env <- subset(df.env, age_class != "A_Cm")

#data3
df.env3 <- cbind(df.exp3[,c("age","clay","pH","mc","cn","ats1","prec1")],
                df.groups3[,c("age_class")], 
                df.bdv3)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Plot Aid
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Background color of symbols
pt_bg <- grey.colors(5, start = 0.8, end = 0, gamma = 2.2, alpha = NULL)
# p_bg1 <- rep(rep(p_bg, each=12),3)
pt_bg2 <- rep(rep(pt_bg, each=3),3)
pt_bg3 <- c(rep(pt_bg, each=3),rep(pt_bg[5],3))
# ## shape of symbols
# p_shp1 <- c(rep(rep(c(24,25,23,22,21), each=12),3))
pt_shp2 <- c(rep(rep(c(24,25,23,22,21), each=3),3))
pt_shp3 <- c(rep(c(24,25,23,22,21), each=3),rep(21,3))
pt_shp.legend <- c(24,25,23,22,21)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Data Inspection #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Minimum and maximum abundances in data set
range(df.spe)

# count cases for each abundance (class)
ab <- table(unlist(df.spe))
ab
 # barplot of the distribution, all species confounded
barplot(ab, las=1, xlab= "abundance", ylab = "Frequency", col=gray(5:0/5))

# Number of absences
sum(df.spe==0)

# Propoortion of zeroes in teh community data set
sum(df.spe==0)/(nrow(df.spe)*ncol(df.spe)) # 60%


require(ggplot2)
ggplot(data2, aes(age_class, OCY+OLA)) + geom_bar(stat="identity")+ geom_point() 
ggplot(data2, aes(age_class, LCA)) + geom_bar(stat="identity")+ geom_point() 
ggplot(data2, aes(age_class, LRU)) + geom_bar(stat="identity")+ geom_point() 
ggplot(data2, aes(age_class, SR)) + geom_bar(stat="identity")+ geom_point() 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Correspondnce Analysis
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

df.spe.pa <- decostand(df.spe, method="pa")


require(vegan)
spe.ca <- cca(df.spe)

summary(spe.ca, display=NULL)
summary(spe.ca)
summary(spe.ca,scaling=1)

# In CA, Eigenvalues over 0.6 indicates a very strong gradient in the data.
# scaling affects the eigenvectors, but not the eigenvalues

# Kaiser Gutman and Broken Stick Criteria
x11()
(ev2 <- spe.ca$CA$eig)
evplot(ev2)

## CA biplots ####
#*************************************************************

a = 7.3 # height
# Figure sizes in cm
b = 7.3 # width
c = 9   # canvas height
d = 9.2 # width

# Plot parameters
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


# scaling 1:
# sites are centroids of species
plot.ca1 <- plot(spe.ca, scaling=1, main="CA earthworm abundances - biplot scaling 1", type="t")
points(plot.ca1, "sites", pch=pt_shp2, bg=pt_bg2, cex=1.0)

# Scaling 2: species are centroids of sites
plot.ca2 <- plot(spe.ca, scaling=2, main="CA fish abundances - biplot scaling 2", type="t")
points(plot.ca2, "sites", pch=pt_shp2, bg=pt_bg2, cex=1.0)
legend("bottomleft", legend=c("Sp_O", "SP_I2", "SP_I1", "Sp_Y", "Cm"), text.col="black",  cex=0.54, pt.cex=0.64, pch=pt_shp.legend, pt.bg=pt_bg)


## Canonical Correspondence Analysis 
#*************************************************************
spe.cca <- cca(df.spe ~ . , df.env[,-1])

plot.cca <- plot(spe.cca, scaling=2, display=c("sp", "lc","cn"), type="p")
points(plot.cca, "constraints", pch=pt_shp2, bg=pt_bg2, cex=1.0)

anova(spe.cca, step=1000)
anova(spe.cca, by="axis",step=1000)
anova(spe.cca, by="term",step=1000)

cca.step.forward <- ordistep(cca(df.spe~1, data=df.env), scope=formula(spe.cca), direction="forward", pstep=1000)

# Most parsimonious CCA
spe.cca <- cca(df.spe ~ SR + age_class + clay , df.env)
plot.cca <- plot(spe.cca, scaling=2, display=c("sp", "lc","cn"), type="p")
points(plot.cca, "constraints", pch=pt_shp2, bg=pt_bg2, cex=1.0)
legend("bottomleft", legend=c("Sp_O", "SP_I2", "SP_I1", "Sp_Y", "Cm"), text.col="black",  cex=0.54, pt.cex=0.64, pch=pt_shp.legend, pt.bg=pt_bg)

anova(spe.cca, step=1000)
anova(spe.cca, by="axis",step=1000)
anova(spe.cca, by="term",step=1000)

# for data3
spe.cca <- cca(df.spe3 ~ . , df.env3)
plot.cca <- plot(spe.cca, scaling=2, display=c("sp", "lc","cn"), type="p")
points(plot.cca, "constraints", pch=pt_shp3, bg=pt_bg3, cex=1.0)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# NMDS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# NMDS
nmds.spe <- metaMDS(df.spe,"bray",binary=F, k=2, trymax=1000)
nmds.spe <- metaMDS(df.spe,"bray",binary=F, k=2, trymax=1000, previous.best = nmds.spe3)
nmds.spe
envfit(nmds.spe, df.env, perm=999)
# Try also Jaccard and binary=T

# NMDS plot
p.nmds.spe <- plot(nmds.spe, type="n")
points(p.nmds.spe, "sites", pch = pt_shp2, bg=pt_bg2, cex=0.8)
text(nmds.spe, display="species")
plot(envfit(nmds.spe, df.env[,c("age","SR","H")]), add = TRUE, col="blue")
legend("bottomleft", legend=c("Sp_O", "SP_I2", "SP_I1", "Sp_Y", "Cm"), text.col="black",  cex=0.54, pt.cex=0.64, pch=pt_shp.legend, pt.bg=pt_bg)


# data3
nmds.spe3 <- metaMDS(df.spe3,"bray",binary=F, k=2, trymax=1000)
nmds.spe3 <- metaMDS(df.spe3,"bray",binary=F, k=2, trymax=1000, previous.best = nmds.spe3)
nmds.spe3
envfit(nmds.spe3, df.env3, perm=999)

p.nmds.spe <- plot(nmds.spe3, type="n")
points(p.nmds.spe, "sites", pch = pt_shp3, bg=pt_bg3, cex=0.8)
text(nmds.spe3, display="species")
plot(envfit(nmds.spe3, df.env3[,c("age","SR","H")]), add = TRUE, col="blue")
legend("bottomleft", legend=c("Sp_O", "SP_I2", "SP_I1", "Sp_Y", "Cm"), text.col="black",  cex=0.54, pt.cex=0.64, pch=pt_shp.legend, pt.bg=pt_bg)





