############
# F2_EW Cluster Analysis
# Quentin Schorpp
# 11.05.2015
############


##############################################################
# Cluster Analysis

# Means for age class
# select species
spe.dbln <- subset(data3, select=c(age_class,ACA:ARO,LCA:LTR,OCY,OLA))
# calculate means
spe.dbln <- aggregate( .~age_class, spe.dbln,mean) 
# set rownames
spe.dbln <- data.frame(spe.dbln[,-1], row.names=spe.dbln[,1])


# Calculate Bray curtis similarity between age classes
spe.dbln <- vegdist(log1p(spe.dbln))*100 # Bray Curtis dissimilarity by default 
                                         # *100 to get values in percentage instead of fractions

# Single linkage agglomerative clustering
spe.dbln.single <- hclust(spe.dbln, method="single")
BC1 <- as.dendrogram(spe.dbln.single)

# Plot Dendrogram
par(mfrow=c(1,1))
plot(BC1, h=-1, xlim=c(100,0), axes=F, center=F)
axis(side=1, at=seq(0,100, by=10))
##############################################################

rm(spe.dbln, spe.dbln.single, BC1)
