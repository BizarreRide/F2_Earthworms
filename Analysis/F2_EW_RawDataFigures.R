#################
# F2_Eearthworms
# Produce Graphics of raw data
# Quentin Schorpp
# 07.05.2015
#################


##############################################################
# source data

source("Data/GatherSource/F2_EW_MakeLikeFile.R")

##############################################################


##############################################################
# Figures of raw data for functional groups
# plot averages with 1.96 * standard errors


## Version 1 with "1st Order" data ####


# create subset of real data with age_class, abundance (plus N) and biomass
data1.raw.fig <- data[,c(which(colnames(data)=="age_class"),                      
                        which(colnames(data)=="anc"):which(colnames(data)=="epi"),
                        which(colnames(data)=="endad"),
                        which(colnames(data)=="N"),
                        which(colnames(data)=="anc.bm"):which(colnames(data)=="epi.bm"),
                        which(colnames(data)=="endad.bm"),
                        which(colnames(data)=="N.bm"))]


# Melt dataframe for use in ggplot
# Abundance
data1.rf.melted <- melt(data1.raw.fig, id.vars=c(1,7:11), variable="sfg", value.name="abundance")
data1.rf.melted <- data1.rf.melted[,c(1,7,8)]
data1.rf.melted <- ddply(data1.rf.melted, c("age_class", "sfg"), summarise, 
                        abc.n=length(abundance),
                        abc.sum=sum(abundance),
                        abc.mean=mean(abundance),
                        abc.sd=sd(abundance),
                        abc.se=abc.sd/sqrt(abc.n)) 

# Biomass
data1.rf.melted2 <- melt(data1.raw.fig, id.vars=c(1,2:6), variable="sfg.bm", value.name="biomass")
data1.rf.melted2 <- data1.rf.melted2[,c(1,7,8)]
data1.rf.melted2 <- ddply(data1.rf.melted2, c("age_class", "sfg.bm"), summarise, 
                         bm.n=length(biomass),
                         bm.sum=sum(biomass),
                         bm.mean=mean(biomass),
                         bm.sd=sd(biomass),
                         bm.se=bm.sd/sqrt(bm.n)) 

data1.rf <- cbind(data1.rf.melted,data1.rf.melted2[,c(3:7)])
data1.rf2 <- data1.rf  # second data frame to plot fraction of adult endogeics
data1.rf2[data1.rf2$sfg!="endad",]$abc.mean=0 # set all other groups to zero
data1.rf2[data1.rf2$sfg!="endad",]$bm.mean=0 # set all other groups to zero
data1.rf2$sfg <- factor(data1.rf2$sfg , levels=c("anc", "endo", "endad", "epi", "N")) # set order of levels


# Barplots
# Abundance
rfig1 <- ggplot(data1.rf[data1.rf$sfg!="endad",], aes(x=age_class, y=4*abc.mean, fill=sfg)) +  
            geom_bar(stat="identity", position="dodge") + 
            geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, show_guide=FALSE) + 
            geom_bar(stat="identity", position="dodge", data=data1.rf2[data1.rf2$sfg!="endo",]) +
            geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, show_guide=FALSE, data=data1.rf2[data1.rf2$sfg!="endo",]) +
            geom_errorbar(aes(ymin=4*abc.mean-4*1.96*abc.se, ymax=4*abc.mean+4*1.96*abc.se), position=position_dodge(0.9),width=0.15, size=0.15) +
            xlab("Age Class") + 
            ylab("Abundance") +
            ylim(-10,max(data1.rf$abc.mean+data1.rf$abc.se)) +
            labs(fill="Functional Group") +
            scale_fill_grey(labels=c("anecic","endogeic adult", "endogeic juvenile","epigeic", "total")) +
            scale_y_continuous(breaks=pretty_breaks(n=10)) +
            scale_x_discrete(labels=c("Cm", "Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +  
            mytheme +
            guides(fill=guide_legend(keywidth=0.5, keyheight=0.5)) +
            theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1),
                  legend.title=element_text(size=6),
                  legend.text=element_text(size=7),
                  legend.position=c(0.18,0.68))

#ggsave(rfig1,filename="Analysis/Figures/Figure1.pdf", width=9, height=7, units="cm", useDingbats=FALSE)

# Biomass
rfig2 <- ggplot(data1.rf[data1.rf$sfg!="endad",], aes(x=age_class, y=4*bm.mean, fill=sfg)) +  
            geom_bar(stat="identity", position="dodge") + 
            geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, show_guide=FALSE) + 
            geom_bar(stat="identity", position="dodge", data=data1.rf2[data1.rf2$sfg!="endo",]) +
            geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, show_guide=FALSE, data=data1.rf2[data1.rf2$sfg!="endo",]) +
            geom_errorbar(aes(ymin=4*bm.mean-4*1.96*bm.se, ymax=4*bm.mean+4*1.96*bm.se), position=position_dodge(0.9),width=0.15, size=0.15) +
            xlab("Age Class") + 
            ylab("Abundance") +
            ylim(-10,max(data1.rf$bm.mean+data1.rf$bm.se)) +
            labs(fill="Functional Group", cex=6) +
            scale_fill_grey(labels=c("anecic","endogeic adult", "endogeic juvenile", "epigeic", "total")) +
            scale_y_continuous(breaks=pretty_breaks(n=10)) +
            scale_x_discrete(labels=c("Cm", "Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +  
            mytheme +
            guides(fill=guide_legend(keywidth=0.5, keyheight=0.5)) +
            theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1),
                  legend.title=element_text(size=6),
                  legend.text=element_text(size=7),
                  legend.position=c(0.18,0.68))

#ggsave(rfig2,filename="Analysis/Figures/Figure2.pdf", width=9, height=7, units="cm", useDingbats=FALSE)


## Version 2 with "2nd Order" data ####


# create subset of raw data2 with age_class, abundance (plus N) and biomass
data2.raw.fig <- data2[,c(which(colnames(data2)=="age_class"),                      
                          which(colnames(data2)=="anc"):which(colnames(data2)=="epi"),
                          which(colnames(data2)=="endad"),
                          which(colnames(data2)=="N"),
                          which(colnames(data2)=="anc.bm"):which(colnames(data2)=="epi.bm"),
                          which(colnames(data2)=="endad.bm"),
                          which(colnames(data2)=="N.bm"))]


# melt frame for use in ggplot
# Abundance
data2.rf.melted <- melt(data2.raw.fig, id.vars=c(1,7:11), variable="sfg", value.name="abundance")
data2.rf.melted <- data2.rf.melted[,c(1,7,8)]
data2.rf.melted <- ddply(data2.rf.melted, c("age_class","sfg"), summarise, 
                         abc.n=length(abundance),
                         abc.sum=sum(abundance),
                         abc.mean=mean(abundance),
                         abc.sd=sd(abundance),
                         abc.se=abc.sd/sqrt(abc.n)) 
# Biomass
data2.rf.melted2 <- melt(data2.raw.fig, id.vars=c(1,2:6), variable="sfg.bm", value.name="biomass")
data2.rf.melted2 <- data2.rf.melted2[,c(1,7,8)]
data2.rf.melted2 <- ddply(data2.rf.melted2, c("age_class","sfg.bm"), summarise, 
                          bm.n=length(biomass),
                          bm.sum=sum(biomass),
                          bm.mean=mean(biomass),
                          bm.sd=sd(biomass),
                          bm.se=bm.sd/sqrt(bm.n)) 

data2.rf <- cbind(data2.rf.melted,data2.rf.melted2[,c(2:7)])
data2.rf2 <- data2.rf  # second data frame to plot fraction of adult endogeics
data2.rf2[data1.rf2$sfg!="endad",]$abc.mean=0 # set all other groups to zero
data2.rf2[data1.rf2$sfg!="endad",]$bm.mean=0 # set all other groups to zero
data2.rf2$sfg <- factor(data2.rf2$sfg , levels=c("anc", "endo", "endad", "epi", "N")) # set order of levels


# Barplots
# Abundance
rfig3 <- ggplot(data2.rf[data1.rf$sfg!="endad",], aes(x=age_class, y=abc.mean, fill=sfg)) +  
            geom_bar(stat="identity", position="dodge") + 
            geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, guide=FALSE) + 
            geom_bar(stat="identity", position="dodge",data=data1.rf2[data1.rf2$sfg!="endo",]) + 
            geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, guide=FALSE, data=data2.rf2[data2.rf2$sfg!="endo",]) + 
            geom_errorbar(aes(ymin=abc.mean-1.96*abc.se, ymax=abc.mean+1.96*abc.se), position=position_dodge(0.9),width=0.15, size=0.15) +
            ggtitle("Data 2nd Order") +
            xlab("Age Class") + 
            ylab("Abundance") +
            ylim(-10,max(data2.rf$abc.mean+data2.rf$abc.se)) +
            labs(fill="Functional Group") +
            scale_fill_grey(labels=c("anecic","endogeic", "epigeic", "total")) +
            scale_y_continuous(breaks=pretty_breaks(n=10)) +
            scale_x_discrete(labels=c("Cm", "Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +  
            mytheme +
            guides(fill=guide_legend(keywidth=0.5, keyheight=0.5)) +
            theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1),
                  legend.title=element_text(size=6),
                  legend.text=element_text(size=7),
                  legend.position=c(0.18,0.68))

# Biomass
rfig4 <- ggplot(data2.rf[data1.rf$sfg!="endad",], aes(x=age_class, y=bm.mean, fill=sfg)) +  
            geom_bar(stat="identity", position="dodge") + 
            geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, show_guide=FALSE) + 
            geom_bar(stat="identity", position="dodge",data=data1.rf2[data1.rf2$sfg!="endo",]) + 
            geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, guide=FALSE, data=data2.rf2[data2.rf2$sfg!="endo",]) + 
            geom_errorbar(aes(ymin=bm.mean-1.96*bm.se, ymax=bm.mean+1.96*bm.se), position=position_dodge(0.9),width=0.15, size=0.15) +
            ggtitle("Data 2nd Order") +
            xlab("Age Class") + 
            ylab("Abundance") +
            ylim(-10,max(data2.rf$bm.mean+data2.rf$bm.se)) +
            labs(fill="Functional Group") +
            scale_fill_grey(labels=c("anecic","endogeic", "epigeic", "total")) +
            scale_y_continuous(breaks=pretty_breaks(n=10)) +
            scale_x_discrete(labels=c("Cm", "Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +  
            mytheme +
            guides(fill=guide_legend(keywidth=0.5, keyheight=0.5)) +
            theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1),
                  legend.title=element_text(size=6),
                  legend.text=element_text(size=7),
                  legend.position=c(0.18,0.68))


## Version 3 with "3rd Order" data ####


# create subset of real data2 with age_class, abundance (plus N) and biomass
data3.raw.fig <- data3[,c(which(colnames(data3)=="age_class"),                      
                          which(colnames(data3)=="anc"):which(colnames(data3)=="epi"),
                          which(colnames(data3)=="endad"),
                          which(colnames(data3)=="N"),
                          which(colnames(data3)=="anc.bm"):which(colnames(data3)=="epi.bm"),
                          which(colnames(data3)=="endad.bm"),
                          which(colnames(data3)=="N.bm"))]


# melt frame for use in ggplot
# Abundance
data3.rf.melted <- melt(data3.raw.fig, id.vars=c(1,7:11), variable="sfg", value.name="abundance")
data3.rf.melted <- data3.rf.melted[,c(1,7,8)]
data3.rf.melted <- ddply(data3.rf.melted, c("age_class","sfg"), summarise, 
                         abc.n=length(abundance),
                         abc.sum=sum(abundance),
                         abc.mean=mean(abundance),
                         abc.sd=sd(abundance),
                         abc.se=abc.sd/sqrt(abc.n)) 
# Biomass
data3.rf.melted2 <- melt(data3.raw.fig, id.vars=c(1,2:6), variable="sfg.bm", value.name="biomass")
data3.rf.melted2 <- data3.rf.melted2[,c(1,7,8)]
data3.rf.melted2 <- ddply(data3.rf.melted2, c("age_class","sfg.bm"), summarise, 
                          bm.n=length(biomass),
                          bm.sum=sum(biomass),
                          bm.mean=mean(biomass),
                          bm.sd=sd(biomass),
                          bm.se=bm.sd/sqrt(bm.n)) 

data3.rf <- cbind(data3.rf.melted,data3.rf.melted2[,c(2:7)])
data3.rf2 <- data3.rf  # second data frame to plot fraction of adult endogeics
data3.rf2[data1.rf2$sfg!="endad",]$abc.mean=0 # set all other groups to zero
data3.rf2[data1.rf2$sfg!="endad",]$bm.mean=0 # set all other groups to zero
data3.rf2$sfg <- factor(data3.rf2$sfg , levels=c("anc", "endo", "endad", "epi", "N")) # set order of levels

# Barplots 

# Abundance
rfig5 <- ggplot(data3.rf[data3.rf$sfg!="endad",], aes(x=age_class, y=abc.mean, fill=sfg)) +  
            geom_bar(stat="identity", position="dodge") + 
            geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, show_guide=FALSE) + 
            geom_bar(stat="identity", position="dodge",data=data3.rf2[data3.rf2$sfg!="endo",]) + 
            geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, show_guide=FALSE, data=data3.rf2[data2.rf2$sfg!="endo",]) + 
            geom_errorbar(aes(ymin=abc.mean-abc.se, ymax=abc.mean+abc.se), position=position_dodge(0.9),width=0.15, size=0.15) +
            ggtitle("Data 3rd Order") +
            xlab("Age Class") + 
            ylab("Abundance") +
            ylim(-10,max(data3.rf$abc.mean+data3.rf$abc.se)) +
            labs(fill="Functional Group") +
            scale_fill_grey(labels=c("anecic","endogeic", "epigeic", "total")) +
            scale_y_continuous(breaks=pretty_breaks(n=10)) +
            scale_x_discrete(labels=c("Cm", "Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +  
            mytheme +
            guides(fill=guide_legend(keywidth=0.5, keyheight=0.5)) +
            theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1),
                  legend.title=element_text(size=6),
                  legend.text=element_text(size=7),
                  legend.position=c(0.18,0.68))

# Biomass
rfig6 <- ggplot(data3.rf[data3.rf$sfg!="endad",], aes(x=age_class, y=bm.mean, fill=sfg)) +  
            geom_bar(stat="identity", position="dodge") + 
            geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, show_guide=FALSE) + 
            geom_bar(stat="identity", position="dodge",data=data3.rf2[data3.rf2$sfg!="endo",]) + 
            geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, show_guide=FALSE, data=data3.rf2[data2.rf2$sfg!="endo",]) + 
            geom_errorbar(aes(ymin=bm.mean-bm.se, ymax=bm.mean+bm.se), position=position_dodge(0.9),width=0.15, size=0.15) +
            ggtitle("Data 3rd Order") +
            xlab("Age Class") + 
            ylab("Biomass") +
            ylim(-10,max(data3.rf$bm.mean+data3.rf$bm.se)) +
            labs(fill="Functional Group") +
            scale_fill_grey(labels=c("anecic","endogeic", "epigeic", "total")) +
            scale_y_continuous(breaks=pretty_breaks(n=10)) +
            scale_x_discrete(labels=c("Cm", "Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +  
            mytheme +
            guides(fill=guide_legend(keywidth=0.5, keyheight=0.5)) +
            theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1),
                  legend.title=element_text(size=6),
                  legend.text=element_text(size=7),
                  legend.position=c(0.18,0.68))

##############################################################

# Why are the error bars becoming larger and larger with the Order of the dataset?
  
# One reason is that the standard error is the standard deviation divided by the squareroot of number of items.
# The latter decreases with increasing order of dataset:  
# The squareroot in 1st order is ``r sqrt(4*3*3)``  
# The squareroot in 2nd order is ``r sqrt(3*3)``  
# The squareroot in 3rd order is ``r sqrt(3)``   


# Clean up
rm(data1.raw.fig, data1.rf, data1.rf2, data1.rf.melted, data1.rf.melted2,
   data2.raw.fig, data2.rf, data2.rf2, data2.rf.melted, data2.rf.melted2,
   data3.raw.fig, data3.rf, data3.rf2, data3.rf.melted, data3.rf.melted2)