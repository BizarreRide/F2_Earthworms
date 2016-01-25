
#################
# F2_Eearthworms
# Prediciton Plots
# Quentin Schorpp
# 07.05.2015
#################

# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Predictions
#load(file="Analysis/F2_EW_lsPred.rda")
load(file="Analysis/OutputTables/Delta4_OutlierNewClimate/F2_EW_lsPred_OutlierNewCLimate.rda")

# Raw Data
source("Data/GatherSource/F2_EW_MakeLikeFile.R")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Data processing ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

df.raw.fig <- data[,c(which(colnames(data)=="age_class"),                      
                         which(colnames(data)=="samcam"),
                         which(colnames(data)=="anc"):which(colnames(data)=="epi"),
                         which(colnames(data)=="ancad"),
                         which(colnames(data)=="endad"),
                         which(colnames(data)=="N"),
                         which(colnames(data)=="anc.bm"):which(colnames(data)=="epi.bm"),
                         which(colnames(data)=="ancad.bm"),
                         which(colnames(data)=="endad.bm"),
                         which(colnames(data)=="N.bm"))]
df.raw.fig$samcam <- revalue(df.raw.fig$samcam , c("1" ="autumn 2012",  "2" ="spring 2013", "3"="autumn 2013"))
df.raw.fig$age_class <- revalue(df.raw.fig$age_class , c("A_Cm" ="Cm",  "B_Sp_young" ="Sp_Y", "C_Sp_int1"="Sp_I1", "D_Sp_int2" = "Sp_I2", "E_Sp_old" = "Sp_O"))


# Melt dataframes for use in ggplot

# Abundance ####
#************************************************************

df.raw.melted <- melt(df.raw.fig, id.vars=c(1,2,9:14), variable="sfg", value.name="abundance")
df.raw.melted$sfg <- revalue(df.raw.melted$sfg , c("endo" ="end"))
df.raw.melted$sfg <- factor(df.raw.melted$sfg , levels=c("anc", "ancad","end", "endad", "epi", "N")) # set order of levels
df.raw.melted <- df.raw.melted[,c(1,2,9,10)]
df.raw.melted <- ddply(df.raw.melted, c("age_class", "samcam", "sfg"), summarise, 
                         abc.n=length(abundance),
                         abc.sum=sum(abundance),
                         fit=mean(abundance),
                         abc.sd=sd(abundance),
                         abc.se=abc.sd/sqrt(abc.n),
                         upr =fit + 1.96*abc.se, 
                         lwr =fit - 1.96*abc.se) 

df.raw.abn <- df.raw.melted[,c("age_class","samcam","sfg","fit","lwr","upr")]
df.raw.abn <- df.raw.abn[df.raw.abn$sfg!="N",]
df.raw.abn2 <- df.raw.abn
df.raw.abn2[!df.raw.abn2$sfg %in% c("endad","ancad"),]$fit=0 # set all other groups to zero; only MEAN!!


sfg <- rep(c("ancad", "endad"), each=45)

df.prd.abn <- rbind(ls.pred[[2]], ls.pred[[19]])
df.prd.abn$sfg <- sfg
df.prd.abn$fit <- exp(df.prd.abn$fit)
df.prd.abn$upr <- exp(df.prd.abn$upr)
df.prd.abn$lwr <- exp(df.prd.abn$lwr)
df.prd.abn <- unique(df.prd.abn[,c("age_class","samcam","sfg","fit","lwr","upr")])
df.prd.abn$samcam <- revalue(df.prd.abn$samcam , c("1" ="autumn 2012",  "2" ="spring 2013", "3"="autumn 2013"))


df.plot.abn <- rbind(df.raw.abn, df.prd.abn)
df.plot.abn$measure <- c(rep("observed",nrow(df.raw.abn)), rep("predicted",nrow(df.prd.abn)))
df.plot.abn2 <- df.plot.abn
df.plot.abn2[!df.plot.abn2$sfg %in% c("endad","ancad"),]$fit=0 # set all other groups to zero; only MEAN!!

# Biomass ####
#************************************************************

df.raw.melted2 <- melt(df.raw.fig, id.vars=c(1,2,3:8), variable="sfg", value.name="biomass")
df.raw.melted2$sfg <- revalue(df.raw.melted2$sfg , c("anc.bm" = "anc", "ancad.bm" = "ancad", "endo.bm" ="end", "endad.bm" = "endad", "epi.bm" = "epi", "N.bm" = "N"))
df.raw.melted2$sfg <- factor(df.raw.melted2$sfg , levels=c("anc", "ancad","end", "endad", "epi", "N")) # set order of levels
df.raw.melted2 <- df.raw.melted2[,c(1,2,9,10)]
df.raw.melted2 <- ddply(df.raw.melted2, c("age_class", "samcam","sfg"), summarise, 
                          bm.n=length(biomass),
                          bm.sum=sum(biomass),
                          fit=mean(biomass),
                          bm.sd=sd(biomass),
                          bm.se=bm.sd/sqrt(bm.n),
                          upr =fit + 1.96*bm.se, 
                          lwr =fit - 1.96*bm.se) 

df.raw.bms <- df.raw.melted2[,c("age_class","samcam","sfg","fit","lwr","upr")]
df.raw.bms <- df.raw.bms[df.raw.bms$sfg!="N",]
df.raw.bms2 <- df.raw.bms
df.raw.bms2[!df.raw.bms2$sfg %in% c("endad","ancad"),]$fit=0 # set all other groups to zero; only MEAN!!


sfg <- rep(c("ancad", "endad"), each=45)

df.prd.bms <- rbind(ls.pred[[10]], ls.pred[[12]])
df.prd.bms$sfg <- sfg
df.prd.bms$fit <- exp(df.prd.bms$fit)
df.prd.bms$upr <- exp(df.prd.bms$upr)
df.prd.bms$lwr <- exp(df.prd.bms$lwr)
df.prd.bms <- unique(df.prd.bms[,c("age_class","samcam","sfg","fit","lwr","upr")])
df.prd.bms$samcam <- revalue(df.prd.bms$samcam , c("1" ="autumn 2012",  "2" ="spring 2013", "3"="autumn 2013"))


df.plot.bms <- rbind(df.raw.bms, df.prd.bms)
df.plot.bms$measure <- c(rep("observed",nrow(df.raw.bms)), rep("predicted",nrow(df.prd.bms)))
df.plot.bms2 <- df.plot.bms
df.plot.bms2[!df.plot.bms2$sfg %in% c("endad","ancad"),]$fit=0 # set all other groups to zero; only MEAN!!

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Barplots ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Abundance ####

plot.abn <- ggplot(df.plot.abn[!df.plot.abn$sfg %in% c("endad","ancad"),], aes(x=age_class, y=fit, fill=sfg)) +  
  geom_bar(stat="identity", position="dodge") + 
  geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, show_guide=FALSE) + 
  geom_bar(stat="identity", position="dodge", data=df.plot.abn2[!df.plot.abn2$sfg %in% c("anc","end"),]) +
  geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, show_guide=FALSE, data=df.plot.abn2[!df.plot.abn2$sfg %in% c("anc","end"),]) +
  geom_errorbar(aes(ymin=lwr, ymax=upr), position=position_dodge(0.9),width=0.15, size=0.15) +
  geom_errorbar(aes(ymin=lwr, ymax=upr), position=position_dodge(0.9),width=0.15, size=0.15, data=df.plot.abn[!df.plot.abn$measure %in% "observed",]) +
  facet_grid(measure~samcam, scales="free_y") +
  xlab("Age Class") + 
  ylab("Abundance") +
  #ylim(-10,max(df.plot.abn$fit+df.plot.abn$abc.se)) +
  labs(fill="Functional Group") +
  scale_fill_grey(labels=c("anecic juvenile","anecic adult","endogeic juvenile", "endogeic adult","epigeic", "total")) +
  scale_y_continuous(breaks=pretty_breaks(n=10)) +
  scale_x_discrete(labels=c("Cm", "Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +  
  mytheme +
  guides(fill=guide_legend(keywidth=0.5, keyheight=0.5)) +
  theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1),
        legend.title=element_text(size=6),
        legend.text=element_text(size=7),
        legend.position=c(0.1,0.9))

#ggsave(plot.abn,filename="Analysis/Figures/Figure1.pdf", width=14, height=11, units="cm", useDingbats=FALSE)


# Biomass ####

plot.bms <- ggplot(df.plot.bms[!df.plot.bms$sfg %in% c("endad","ancad"),], aes(x=age_class, y=fit, fill=sfg)) +  
  geom_bar(stat="identity", position="dodge") + 
  geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, show_guide=FALSE) + 
  geom_bar(stat="identity", position="dodge", data=df.plot.bms2[!df.plot.bms2$sfg %in% c("anc","end"),]) +
  geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, show_guide=FALSE, data=df.plot.bms2[!df.plot.bms2$sfg %in% c("anc","end"),]) +
  geom_errorbar(aes(ymin=lwr, ymax=upr), position=position_dodge(0.9),width=0.15, size=0.15) +
  geom_errorbar(aes(ymin=lwr, ymax=upr), position=position_dodge(0.9),width=0.15, size=0.15, data=df.plot.bms[!df.plot.bms$measure %in% "observed",]) +
  facet_grid(measure~samcam, scales="free_y") +
  xlab("Age Class") + 
  ylab("Biomass") +
  #ylim(-10,max(df.plot.bms$fit+df.plot.bms$abc.se)) +
  labs(fill="Functional Group") +
  scale_fill_grey(labels=c("anecic juvenile","anecic adult","endogeic juvenile", "endogeic adult","epigeic", "total")) +
  scale_y_continuous(breaks=pretty_breaks(n=10)) +
  scale_x_discrete(labels=c("Cm", "Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +  
  mytheme +
  guides(fill=guide_legend(keywidth=0.5, keyheight=0.5)) +
  theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1),
        legend.title=element_text(size=6),
        legend.text=element_text(size=7),
        legend.position=c(0.1,0.9))


#ggsave(plot.bms,filename="Analysis/Figures/Figure2.pdf", width=14, height=11, units="cm", useDingbats=FALSE)


# Unimodal responses #####

# Anecic vs Temperature
newdata <- expand.grid(age_class=unique(data$age_class),
                       samcam = unique(data$samcam),               
                       scl.ats1 = seq(min(dt.exp2$scl.ats1),max(dt.exp2$scl.ats1), by=0.2),
                       area = 1)

f1 <- paste("~",formula(ls.bestmodels[[1]])[3])
f1 <- gsub("\\+ \\(1 \\| field.ID\\)", "", f1)
#f1 <- gsub("\\+ offset\\(log\\(area\\)\\)", "", f1)
#f1 <- paste(f1, "+ offset(log(area))")
#f1 <- gsub("\\(1.*\\)", "1",f1)
f1 <- formula(print(f1, quote=F))
X <- model.matrix( f1, newdata)
newdata$fit <- X %*% fixef(ls.bestmodels[[1]])
newdata$SE <- sqrt(  diag(X %*%vcov(ls.bestmodels[[1]]) %*% t(X))  )
newdata$upr=newdata$fit+1.96*newdata$SE
newdata$lwr=newdata$fit-1.96*newdata$SE
pred <- newdata

pred$ats1 <- pred$scl.ats1* sd(data$ats1) + mean(data$ats1)


AnecicTemp <- ggplot(pred, aes(x = ats1, y = exp(fit), ymin = exp(lwr), ymax = exp(upr), col=samcam)) + 
  geom_point() +
  #geom_bar(stat="identity",position = position_dodge(1), col="454545", size=0.15, fill="grey") +
  geom_errorbar(position = position_dodge(1),width=0.15, size=0.15) + 
  facet_grid(.~age_class) +
  geom_hline(xintercept = 1, size=0.15) +
  ylab("Anecic Abundance Ind./m²") +
  xlab(expression(paste("T3",0[surface]))) +
  mytheme +
  theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1)) +
  guides(fill=guide_legend(keywidth=0.5, keyheight=0.5)) +
  theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1),
        legend.title=element_text(size=6),
        legend.text=element_text(size=7),
        legend.position=c(0.1,0.9))

# ggsave(AnecicTemp,filename="Analysis/Figures/FigureAncT.pdf", width=14, height=11, units="cm", useDingbats=FALSE)



# Endogeic vs soil misture
newdata <- expand.grid(age_class=unique(data$age_class),
                      samcam = unique(data$samcam),               
                      scl.pH = mean(dt.exp2$scl.pH),
                      scl.mc = seq(min(dt.exp2$scl.mc),max(dt.exp2$scl.mc), by=0.2),
                      area = 1)

f1 <- paste("~",formula(ls.bestmodels[[2]])[3])
f1 <- gsub("\\+ \\(1 \\| field.ID\\)", "", f1)
#f1 <- gsub("\\+ offset\\(log\\(area\\)\\)", "", f1)
#f1 <- paste(f1, "+ offset(log(area))")
#f1 <- gsub("\\(1.*\\)", "1",f1)
f1 <- formula(print(f1, quote=F))
X <- model.matrix( f1, newdata)
newdata$fit <- X %*% fixef(ls.bestmodels[[2]])
newdata$SE <- sqrt(  diag(X %*%vcov(ls.bestmodels[[2]]) %*% t(X))  )
newdata$upr=newdata$fit+1.96*newdata$SE
newdata$lwr=newdata$fit-1.96*newdata$SE
pred <- newdata

pred$mc <- pred$scl.mc* sd(data$mc) + mean(data$mc)


Endogeic.Moist <- ggplot(pred, aes(x = mc, y = exp(fit), ymin = exp(lwr), ymax = exp(upr), col=samcam)) + 
  geom_point() +
  #geom_bar(stat="identity",position = position_dodge(1), col="454545", size=0.15, fill="grey") +
  geom_errorbar(position = position_dodge(1),width=0.15, size=0.15) + 
  facet_grid(.~age_class) +
  geom_hline(xintercept = 1, size=0.15) +
  ylab("Endogeic Abundance Ind./m²") +
  xlab(expression(paste("T3",0[surface]))) +
  mytheme +
  theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1)) +
  guides(fill=guide_legend(keywidth=0.5, keyheight=0.5)) +
  theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1),
        legend.title=element_text(size=6),
        legend.text=element_text(size=7),
        legend.position=c(0.1,0.9))

# ggsave(Endogeic.Moist,filename="Analysis/Figures/FigureEndMC.pdf", width=14, height=11, units="cm", useDingbats=FALSE)


