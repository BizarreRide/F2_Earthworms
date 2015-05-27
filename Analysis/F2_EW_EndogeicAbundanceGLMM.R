#################
# F2_Eearthworms
# GLMM for Endogeic Abundance
# Quentin Schorpp
# 07.05.2015
#################

# Standardize Covariates ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Data/GatherSource/F2_EW_MakeLikeFile.R")
source("Data/GatherSource/CovariateStandardization.R")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#!!!! I model only adult endogeics, since this decreases the varaince and buffers very large abundance values!!!!!!!!

# plot of anecic abundance ####
# *The data we want to model*
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# !!!! requires data from raw figure plots
source("Analysis/F2_EW_RawDataFigures.R")

endad.raw <-  ggplot(data1.rf[data1.rf$sfg=="endad",], aes(x=age_class, y=abc.mean, fill=sfg)) +   
              geom_bar(stat="identity", position="dodge") + 
              geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, show_guide=FALSE) + 
              #geom_bar(stat="identity", position="dodge", data=data1.rf2[data1.rf2$sfg!=c("anc","end"),]) +
              #geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, show_guide=FALSE, data=data1.rf2[data1.rf2$sfg!=c("anc","end"),]) +
              geom_errorbar(aes(ymin=abc.mean-1.96*abc.se, ymax=abc.mean+1.96*abc.se), position=position_dodge(0.9),width=0.15, size=0.15) +
              facet_grid(.~samcam) +
              xlab("Age Class") + 
              ylab("Abundance") +
              #ylim(-10,max(data1.rf$abc.mean+data1.rf$abc.se)) +
              labs(fill="Functional Group") +
              scale_fill_grey(labels=c("endogeic adults")) +
              scale_y_continuous(breaks=pretty_breaks(n=10)) +
              scale_x_discrete(labels=c("Cm", "Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +  
              mytheme +
              guides(fill=guide_legend(keywidth=0.5, keyheight=0.5)) +
              theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1),
                    legend.title=element_text(size=10),
                    legend.text=element_text(size=10),
                    legend.position=c(0.2,0.92))
endad.raw
#ggsave(endad.raw, filename="Analysis/Figures/Figure4_EndadRaw.pdf", width=15, height=11, units="cm", useDingbats=FALSE)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Assess variability in random effects ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Variability within sites
boxplot(endad~field.ID, data, col="grey", main="Variability within sites", xlab="sites orderd in decreasing age", ylab="endogeic earthworm abundance")

# variability within sampling campaigns
boxplot(endad~samcam, data,col="grey", main="Variability within sampling campaigns", xlab="seasons:\n autumn2012, spring2013, autumn2013", ylab="endogeic earthworm abundance")

# variability within sites and sampling campaigns
boxplot(endad~field.ID+samcam, data, las=2, col="grey", main="Variability within sites \n differing at sampling campaigns", xlab="field+seasons", ylab="endogeic earthworm abundance")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## GLMM using glmer & glmmADMB ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Global model Formulation ####

# I Formulate two different candidate models, 
# since some covariates can't be used in the same model due to colinearity

# Model Formulation
endad.glob1 <- glmer(endad ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))
endad.glob2 <- glmer(endad ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))
# offsset is used due to the personal advice by T. Onkelinx:
# You better use an offset if you want to express the model in terms of m². Just add offset(log(0.25)) to the model. 

summary(endad.glob1)
summary(endad.glob2)

# Overdispersion
E1 <- resid(endad.glob2, type="pearson")
N <- nrow(data)
p <- length(coef(endad.glob2)) #  "+1"  would be used negative in neg binomial distribution for determining the number of parameters (p) due to the "k" 
Dispersion <- sum(E1^2)/(N-p)
Dispersion

## Multimodel averaging ####

load("Analysis/F2_EW_AbundanceDredge_glmer.RData")

#endad.dredge1 <- dredge(endad.glob1)
head(endad.dredge1,10)
endad.avgmod1.d4 <- model.avg(endad.dredge1, subset = delta < 4)
summary(endad.avgmod1.d4)
importance(endad.avgmod.d4) 
endad.avgmod.95p <- model.avg(endad.dredge1, cumsum(weight) <= .95)
# AICc range 832 - 836

#endad.dredge2 <- dredge(endad.glob2)
head(endad.dredge2,10)
endad.avgmod2.d4 <- model.avg(endad.dredge2, subset = delta < 4)
summary(endad.avgmod2.d4)
importance(endad.avgmod.d4) 
endad.avgmod.95p <- model.avg(endad.dredge2, cumsum(weight) <= .95)
#AICc range 845 - 849

write.csv(data.frame(endad.avgmod1.d4$importance), "Analysis/OutputTables/EndadImportance.csv")
write.csv(data.frame(endad.avgmod1.d4$coef.shrinkage), "Analysis/OutputTables/EndadShrinkage.csv")
write.csv(data.frame(endad.avgmod1.d4$msTable), "Analysis/OutputTables/EndadSubsetModels.csv")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Fit the best model ####
# with both glmer() and glmmadmb():
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# endad.best2 <- glmmadmb(endad ~ age_class*samcam + scl.prec1 + scl.mc*scl.pH  + (1|field.ID) + offset(log(area)) ,data=data,family="poisson")
endad.best <- glmer(endad ~ age_class*samcam + scl.prec1 + scl.mc*scl.pH  + (1|field.ID) + offset(log(area)) ,data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))

# **The best model includes an Interaction term!!!**

# Summary output ####
# summary
summary(endad.best)

# anova
summary(aov(endad.best))

write.csv(summary(endad.best)$coefficients, "Analysis/OutputTables/EndadBestCoef.csv")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## Model Validation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Confidence Intervals ####
endad.confint <- confint(endad.best)
coefplot2(endad.best)

# Check Model Assumptions ####

E1 <- resid(endad.best, type="pearson")
F1 <- fitted(endad.best, type="response")
P1 <- predict(endad.best, type="response")

par(mfrow=c(2,2),
    mar=c(4,4.5,1,2))
# Plot fitted vs. residuals
scatter.smooth(F1, E1, cex.lab = 1.5, xlab="Fitted values", ylab=" Residuals")
abline(h = 0, v=0, lty=2)


# plot predicted vs. residuals
scatter.smooth(P1, E1, cex.lab = 1.5, xlab="Predicted values", ylab=" Residuals")
abline(h = 0, v=0, lty=2)

# plot fitted vs. predicted
scatter.smooth(F1, P1, cex.lab = 1.5, xlab="Fitted values", ylab="Predicted")
abline(h = 0, v=0, lty=2)

# Histogram of Residuals
hist(E1, prob=TRUE, main = "", breaks = 20, cex.lab = 1.5, xlab = "Response Residuals", col="PapayaWhip")
lines(density(E1), col="light blue", lwd=3)
lines(density(E1, adjust=2), lty="dotted", col="darkgreen", lwd=2) 

# Normal QQ Plot
qqnorm(y=resid(anc.bm.best))
qqline(y=resid(anc.bm.best))

# Cooks Distances

par(lo)

# plot residuals vs. all covariates ####

# Dataset with all variables IN the model
env1 <- cbind(E1, F1, data[,c("endad", "age_class", "samcam", "ats1")]) # should "ats1" be squared?
# covariates NOT in the model
env0 <- cbind(E1, F1,data[,c("mc", "pH", "cn", "clay", "ats2", "hum1","dgO")])

# IN the model
par(mfrow=c(2,3),
    mar=c(3.8,4,1,3))
for (i in 1:length(env1)) {
  scatter.smooth(env1[,i],env1[,1], ylab="Residuals", xlab=colnames(env1)[i])
  abline(h=0, lty=2, col="red")
}
par(mfrow=c(1,1))

# NOT in the model
par(mfrow=c(3,3),
    mar=c(3.8,4,1,3))
for (i in 1:length(env0)) {
  scatter.smooth(env0[,i],env0[,1], ylab="Residuals", xlab=colnames(env0)[i])
  abline(h=0, lty=2, col="red")
}
par(lo)

# Overdispersion ####
N <- nrow(data)
p <- length(coef(endad.best)) # +1 in case of negbin
Dispersion <- sum(E1^2)/(N-p)
Dispersion
# glmer has less dispersion

overdisp_fun(endad.best)

# Explained variation ####

r.squaredGLMM(endad.best)

## Some more plots of fitted and predicted values ####
# Plots of Predicted Values (- random term)
par(mfrow=c(2,2))
plot(data$age_class,P1)
plot(data$samcam,P1)
plot(data$mc,P1)
plot(data$pH,P1)
plot(data$hum1,P1)

# Plots of fitted Values (+ random term)
par(mfrow=c(2,2))
plot(data$age_class,F1)
plot(data$samcam,F1)
plot(data$mc,F1)
plot(data$pH,F1)
plot(data$hum1,F1)

par(lo)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## Predictionplots ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Predictionplots for Age Class x SamCam ####

# create test data set with all covariates IN the model
endad.td = expand.grid(age_class=unique(data$age_class),
                       samcam = unique(data$samcam),               
                       scl.prec1 = mean(data$scl.prec1),
                       scl.mc = mean(data$scl.mc),
                       scl.pH = mean(data$scl.pH),
                       area = 0.25)


## calculate confidence intervals for predictions from test dataset

# In case of glmmadmb:
#endad.pred <- cbind(endad.td, predict(endad.best2, newdata = endad.td, interval = "confidence")) 

# In case of glmer:
      X <- model.matrix(~ ~ age_class*samcam + scl.prec1 + scl.mc*scl.pH, data = endad.td)
      endad.td$fit <- X %*% fixef(endad.best)
      endad.td$SE <- sqrt(  diag(X %*%vcov(endad.best) %*% t(X))  )
      endad.td$upr=endad.td$fit+1.96*endad.td$SE
      endad.td$lwr=endad.td$fit-1.96*endad.td$SE
      endad.pred <- endad.td


# Rename samcam for facetting
endad.pred$samcam2 <- plyr::revalue(endad.td$samcam,c("1" ="autumn 2012",  "2" ="spring 2013", "3"="autumn 2013"))
endad.pred$age_class <- plyr::revalue(endad.td$age_class,c("A_Cm"="Cm","B_Sp_young" ="Sp_Y","C_Sp_int1" ="Sp_I1","D_Sp_int2" ="Sp_I2","E_Sp_old" ="Sp_O"))

# Reverse scaling of covariate
endad.pred$prec1 <- endad.pred$scl.prec1* sd(data$prec1) + mean(data$prec1)
endad.pred$mc <- endad.pred$scl.mc* sd(data$mc) + mean(data$mc)
endad.pred$pH <- endad.pred$scl.pH* sd(data$pH) + mean(data$pH)

# plot predictions with error bars // confidence intervals???
predfig.endad1 <- ggplot(endad.pred, aes(x = age_class, y = exp(fit), ymin = exp(lwr), ymax = exp(upr))) + 
  geom_bar(stat="identity",position = position_dodge(1), col="454545", size=0.15, fill="grey") +
  geom_errorbar(position = position_dodge(1),col="black",width=0.15, size=0.15) + 
  facet_grid(.~samcam2) +
  geom_hline(xintercept = 1, size=0.15) +
  ylab("Endogeic Abundance Ind./m²") +
  xlab("Age Class") +
  scale_x_discrete(labels=c("Cm", "Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +
  mytheme +
  theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1))
predfig.endad1# to predict for age_class only, take the mean of all continous covariates

#ggsave(predfig.endad1,filename="Analysis/Figures/Figure4_EndadPredGlmer.pdf", width=15, height=11, units="cm", useDingbats=FALSE)

# Prediction plots for average temperature! ####
endad.td = expand.grid(age_class=unique(data$age_class),
                       samcam = unique(data$samcam),               
                       scl.prec1 = mean(data$scl.prec1),#seq(min(data$scl.prec1),max(data$scl.prec1), by=0.2),
                       scl.mc = round(seq(min(data$scl.mc),max(data$scl.mc), length.out=5),2),
                       scl.pH = seq(min(data$scl.pH),max(data$scl.pH), by=0.2),
                       area = 1)

## calculate confidence intervals for predictions from test dataset

# In case of glmmadmb:
      #endad.pred <- cbind(endad.td, predict(endad.best2, newdata = endad.td, interval = "confidence")) 

# In case of glmer:
      X <- model.matrix(~ ~ age_class*samcam + scl.prec1 + scl.mc*scl.pH, data = endad.td)
      endad.td$fit <- X %*% fixef(endad.best)
      endad.td$SE <- sqrt(  diag(X %*%vcov(endad.best) %*% t(X))  )
      endad.td$upr=endad.td$fit+1.96*endad.td$SE
      endad.td$lwr=endad.td$fit-1.96*endad.td$SE
      endad.pred <- endad.td


# Rename samcam for facetting
endad.pred$samcam2 <- plyr::revalue(endad.td$samcam,c("1" ="autumn 2012",  "2" ="spring 2013", "3"="autumn 2013"))
endad.pred$age_class <- plyr::revalue(endad.td$age_class,c("A_Cm"="Cm","B_Sp_young" ="Sp_Y","C_Sp_int1" ="Sp_I1","D_Sp_int2" ="Sp_I2","E_Sp_old" ="Sp_O"))

# Reverse scaling of covariate
endad.pred$prec1 <- endad.pred$scl.prec1* sd(data$prec1) + mean(data$prec1)
endad.pred$mc <- round(endad.pred$scl.mc* sd(data$mc) + mean(data$mc),2)
endad.pred$pH <- endad.pred$scl.pH* sd(data$pH) + mean(data$pH)


predfig.endad2 <- ggplot(endad.pred, aes(x = pH, y = exp(fit), ymin = exp(lwr), ymax = exp(upr), col=samcam2)) + 
  geom_point() +
  #geom_bar(stat="identity",position = position_dodge(1), col="454545", size=0.15, fill="grey") +
  geom_errorbar(position = position_dodge(1),width=0.15, size=0.15) + 
  facet_grid(mc~age_class) +
  geom_hline(xintercept = 1, size=0.15) +
  ylab("Endogeic Abundance Ind./m²") +
  xlab("pH") +
  mytheme +
  theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1))
predfig.endad2
#ggsave(predfig.endad2,filename="Analysis/Figures/Figure4_EndadPredGlmerPHxMc.pdf", width=15, height=11, units="cm", useDingbats=FALSE)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Coefficients and Statistics ####
# Claculate a whole lot of coefficients and statistics
endad.pred <- within(endad.pred, {
  AIC <- AIC(endad.best)
  Rrandom <- summary(lm(fitted(endad.best)~ data$endad))$adj.r.squared
  Rsquared <- summary(lm(predict(endad.best, type="response")~ data$endad))$adj.r.squared  
})
head(endad.pred)

endad.stat = round(summary(endad.best)$coefficients[, c(3,4)],4)
endad.coef = coeftab(endad.best)

endad.env <- glmmadmb(endad ~ samcam + I(ats1^2) + (1|field.ID) + offset(log(area)),data=data,family="poisson")
endad.pvalue <- anova(endad.env, endad.best)$"Pr(>Chi)"

endad.OUT1 = endad.pred
endad.OUT2 = cbind(endad.stat,endad.coef,LogLikP=endad.pvalue[2],fixef=fixef(endad.best))
endad.OUT1
endad.OUT2

#write.table(endad.OUT1, "Predictions+Stats+ConfIntervals.csv", sep=";", append=TRUE)
#write.table(endad.OUT2, "Predictions+Stats+ConfIntervals.csv", sep=";", append=TRUE)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## Post-Hoc Multicomparisons ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Post-Hoc Multicomparisons for age class x samcam, ####
# interaction effect - reduced contrast matrix

# Model with the interaction term
data$ia.acl.smc <- interaction(data$age_class, data$samcam)
endad.tukey <- glmer(endad ~ ia.acl.smc + scl.prec1 + scl.mc*scl.pH + (1|field.ID) + offset(log(area)) ,data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))
# summary(endad.tukey)

# Create contrast matrix
source("Analysis/F2_EW_ContrastMatrix.R")

# Pairwise comparisons (with interaction term)
endad.pairwise <- glht(endad.tukey, linfct=mcp(ia.acl.smc = cm1))
endad.pw.ci <- confint(endad.pairwise)
summary(endad.pairwise, test=adjusted(type="fdr"))

# Confidence intervals including 0
endad.pw.sig <- which(endad.pw.ci$confint[,2]>0)
data.frame(names(endad.pw.sig))

# Plot Errorbars
phfig1 <- ggplot(endad.pw.ci, aes(y = lhs, x = exp(estimate), xmin = exp(lwr), xmax = exp(upr))) + 
  geom_errorbarh() + 
  geom_point() + 
  geom_vline(xintercept = 1) +
  mytheme
phfig1

# plot confidence intervals
par(mar=c(2,15,2,2))
plot(endad.pw.ci) 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Post-Hoc Multicomparisons for age class, ####
# main effect

# Pairwise comparisons (without interaction term)
endad.pairwise <- glht(endad.best, mcp(age_class = "Tukey"))
endad.pw.ci <- confint(endad.pairwise)
summary(endad.pairwise, test=adjusted(type="none"))

# Confidence intervals including 0
endad.pw.sig <- which(endad.pw.ci$confint[,2]>0)
data.frame(names(endad.pw.sig))


# Plot Errorbars
ggplot(endad.pw.ci, aes(y = lhs, x = exp(estimate), xmin = exp(lwr), xmax = exp(upr))) + 
  geom_errorbarh() + 
  geom_point() + 
  geom_vline(xintercept = 1) +
  mytheme

# plot confidence intervals
par(mar=c(2,15,2,2))
plot(endad.pw.ci) 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Post-Hoc Multicomparisons for samcam, ####
# main effect

# Pairwise comparisons (without interaction term)
endad.pairwise <- glht(endad.best, mcp(samcam = "Tukey"))
endad.pw.ci <- confint(endad.pairwise)

# Confidence intervals including 0
endad.pw.sig <- which(endad.pw.ci$confint[,2]>0)
data.frame(names(endad.pw.sig))


# Plot Errorbars
ggplot(endad.pw.ci, aes(y = lhs, x = exp(estimate), xmin = exp(lwr), xmax = exp(upr))) + 
  geom_errorbarh() + 
  geom_point() + 
  geom_vline(xintercept = 1) +
  mytheme

# plot confidence intervals
par(mar=c(2,15,2,2))
plot(endad.pw.ci) 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
