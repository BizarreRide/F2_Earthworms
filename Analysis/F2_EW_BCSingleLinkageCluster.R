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
#spe.dbln <- data.frame(spe.dbln[,-1], row.names=spe.dbln[,1])
spe.dbln <- data.frame(spe.dbln[,-1], row.names=c("Cm", "Sp_Y", "Sp_I1", "Sp_I2", "Sp_O"))

# Calculate Bray curtis similarity between age classes
spe.dbln <- vegdist(log1p(spe.dbln))*100 # Bray Curtis dissimilarity by default 
                                         # *100 to get values in percentage instead of fractions

# Single linkage agglomerative clustering
spe.dbln.single <- hclust(spe.dbln, method="single")
BC1 <- as.dendrogram(spe.dbln.single)

# Plot Dendrogram
par(mfrow=c(1,1),
    oma=c(0,0,0,0),
    mar=c(5.5,2,0,5))
plot(BC1, h=-1, xlim=c(100,0), axes=F, center=F, xlab="Bray Curtis similarity [%]")
axis(side=1, at=seq(0,100, by=10))

library(ape)
plot(as.phylo(spe.dbln.single), type = "cladogram", cex = 0.9, label.offset = 1)

library(ggdendro)
BC1 <- dendro_data(BC1)
Figure1.Cluster <- ggplot(segment(BC1)) + 
                      geom_segment(aes(x = y, y = x, xend = yend, yend = xend), size=0.25) + 
                      xlim(100, -20) + 
                      xlab("Bray Curtis similarity [%]") +
                      ylab("") +
                      scale_y_discrete(labels=c("Cm", "Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +
                      geom_text(data = label(BC1), aes(x = y-10, y = x, label = label), cex=3, face="bold", family="Arial") +
                      mytheme +
                      theme(axis.text.y=element_blank(), 
                            axis.ticks.y=element_blank(), 
                            panel.border=element_blank(), 
                            axis.line.y = element_blank(),
                            axis.line.x = element_line(size=0.25))
Figure1.Cluster
#ggsave(Figure1.Cluster,filename="Analysis/Figures/Figure1.pdf", width=7.6, height=5, units="cm", useDingbats=FALSE)

##############################################################

rm(spe.dbln, spe.dbln.single, BC1)
