#################
# F2_Eearthworms
# GLMM for Anecic Abundance
# Quentin Schorpp
# 07.05.2015
#################

# Standardize Covariates ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Data/GatherSource/F2_EW_MakeLikeFile.R")
source("Data/GatherSource/CovariateStandardization.R")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# plot of anecic abundance ####
# *The data we want to model*
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# !!!! requires data from raw figure plots
source("Analysis/F2_EW_RawDataFigures.R")

anc.raw <-  ggplot(data1.rf[data1.rf$sfg=="anc",], aes(x=age_class, y=abc.mean, fill=sfg)) +  
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
              scale_fill_grey(labels=c("anecic total")) +
              scale_y_continuous(breaks=pretty_breaks(n=10)) +
              scale_x_discrete(labels=c("Cm", "Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +  
              mytheme +
              guides(fill=guide_legend(keywidth=0.5, keyheight=0.5)) +
              theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1),
                    legend.title=element_text(size=10),
                    legend.text=element_text(size=10),
                    legend.position=c(0.2,0.92))
anc.raw
#ggsave(anc.raw, filename="Analysis/Figures/Figure4_AncRaw.pdf", width=15, height=11, units="cm", useDingbats=FALSE)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Assess variability in random effects ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Variability within sites
boxplot(anc~field.ID, data, col="grey", main="Variability within sites", xlab="sites orderd in decreasing age", ylab="anecic earthworm abundance")

# variability within sampling campaigns
boxplot(anc~samcam, data,col="grey", main="Variability within sampling campaigns", xlab="seasons:\n autumn2012, spring2013, autumn2013", ylab="anecic earthworm abundance")

# variability within sites and sampling campaigns
boxplot(anc~field.ID+samcam, data, las=2, col="grey", main="Variability within sites \n differing at sampling campaigns", xlab="field+seasons", ylab="anecic earthworm abundance")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



## GLMM using glmer ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Global model Formulation ####

# I Formulate two different candidate models, 
# since some covariates can't be used in the same model due to colinearity

# Model Formulation
anc.glob1 <- glmer(anc ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))
anc.glob2 <- glmer(anc ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))
# offsset is used due to the personal advice by T. Onkelinx:
# You better use an offset if you want to express the model in terms of m². Just add offset(log(0.25)) to the model. 

summary(anc.glob1)
summary(anc.glob2)

# Overdispersion
E1 <- resid(anc.glob2, type="pearson")
N <- nrow(data)
p <- length(coef(anc.glob2)) #  "+1"  would be used negative in neg binomial distribution for determining the number of parameters (p) due to the "k" 
Dispersion <- sum(E1^2)/(N-p)
Dispersion

# Multimodel averaging ####

load("Analysis/F2_EW_AbundanceDredge_glmer.RData")

#anc.dredge1 <- dredge(anc.glob1)
head(anc.dredge1,10)
anc.avgmod1.d4 <- model.avg(anc.dredge1, subset = delta < 4)
summary(anc.avgmod1.d4)
data.frame(importance(anc.avgmod1.d4))
# AICc range 741-745

#anc.dredge2 <- dredge(anc.glob2)
head(anc.dredge2,10)
anc.avgmod2.d4 <- model.avg(anc.dredge2, subset = delta < 4)
summary(anc.avgmod2.d4)
data.frame(importance(anc.avgmod2.d4)) 
# AICc range 735-739

write.csv(data.frame(anc.avgmod2.d4$importance), "Analysis/OutputTables/AncImportance.csv")
write.csv(data.frame(anc.avgmod2.d4$coef.shrinkage), "Analysis/OutputTables/AncShrinkage.csv")
write.csv(data.frame(anc.avgmod2.d4$msTable), "Analysis/OutputTables/AncSubsetModels.csv")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Fit the best model ####
# with both glmer() and glmmadmb():
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# anc.best2 <- glmmadmb(anc ~ age_class*samcam + I(scl.ats1^2) + scl.prec1 + (1|field.ID) + offset(log(area)) ,data=data,family="poisson")
anc.best <- glmer(anc ~ age_class*samcam + I(scl.ats1^2) + scl.prec1 + (1|field.ID) + offset(log(area)) ,data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))

# **The best model includes an Interaction term!!!**

# Summary output ####
# summary
summary(anc.best)

# anova
summary(aov(anc.best))

write.csv(summary(anc.best)$coefficients, "Analysis/OutputTables/AncBestCoef.csv")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## Model Validation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Confidence Intervals ####
confint(anc.best)
coefplot2(anc.best)

# Check Model Assumptions ####

E1 <- resid(anc.best, type="pearson")
E2 <- resid(anc.best, type="response")
F1 <- fitted(anc.best, type="response")
P1 <- predict(anc.best, type="response")

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
qqnorm(E2)
qqline(E2)

# Cooks Distances

par(lo)

# plot residuals vs. all covariates ####

# Dataset with all variables IN the model
env1 <- cbind(E1, F1, data[,c("anc", "age_class", "samcam", "ats1", "prec1")]) # should "ats1" be squared?
# covariates NOT in the model
env0 <- cbind(E1, F1,data[,c("mc", "pH", "cn", "clay", "ats2", "hum1","dgO")])

# in the model
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
p <- length(coef(anc.best)) # +1 in case of negbin
Dispersion <- sum(E1^2)/(N-p)
Dispersion
# glmer has less dispersion

overdisp_fun(anc.best)

# Explained variation ####

r.squaredGLMM(anc.best)

# Some more plots of fitted and predicted values ####
# Plots of Predicted Values
par(mfrow=c(2,2))
plot(data$age_class,P1)
plot(data$samcam,P1)
plot(data$ats1^2,P1)
plot(data$prec1^2,P1)

# Plots of fitted Values
par(mfrow=c(2,2))
plot(data$age_class,F1)
plot(data$samcam,F1)
plot(data$ats1^2,F1)
plot(data$prec1,F1)

par(lo)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## Prediction plots ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Prediction plots for Age Clas x SamCam ####

## create test data set with all covariates IN the model
# to predict for age_class only, take the mean of all continous covariates
anc.td = expand.grid(age_class=unique(data$age_class),
                     samcam = unique(data$samcam),               
                     scl.ats1 = mean(data$scl.ats1),
                     scl.prec1 = mean(data$scl.prec1),
                     area = 1)


## calculate confidence intervals for predictions from test dataset

# In case of glmmadmb:
    #anc.pred <- cbind(anc.td, predict(anc.best2, newdata = anc.td, interval = "confidence")) 

# In case of glmer
    X <- model.matrix(~ age_class*samcam + I(scl.ats1^2) + scl.prec1, data = anc.td)
    anc.td$fit <- X %*% fixef(anc.best)
    anc.td$SE <- sqrt(  diag(X %*%vcov(anc.best) %*% t(X))  )
    anc.td$upr=anc.td$fit+1.96*anc.td$SE
    anc.td$lwr=anc.td$fit-1.96*anc.td$SE
    anc.pred <- anc.td

# Rename samcam for facetting
anc.pred$samcam2 <- plyr::revalue(anc.td$samcam,c("1" ="autumn 2012",  "2" ="spring 2013", "3"="autumn 2013"))
anc.pred$age_class <- plyr::revalue(anc.td$age_class,c("A_Cm"="Cm","B_Sp_young" ="Sp_Y","C_Sp_int1" ="Sp_I1","D_Sp_int2" ="Sp_I2","E_Sp_old" ="Sp_O"))

# Reverse scaling of covariate
anc.pred$ats1 <- anc.pred$scl.ats1* sd(data$ats1) + mean(data$ats1)

## plot predictions with error bars // confidence intervals???
predfig.anc1 <- ggplot(anc.pred, aes(x = age_class, y = exp(fit), ymin = exp(lwr), ymax = exp(upr))) + 
  geom_bar(stat="identity",position = position_dodge(1), col="454545", size=0.15, fill="grey") +
  geom_errorbar(position = position_dodge(1),col="black",width=0.15, size=0.15) + 
  facet_grid(.~samcam2) +
  geom_hline(xintercept = 1, size=0.15) +
  ylab("Anecic Abundance Ind./m²") +
  xlab("Age Class") +
  scale_x_discrete(labels=c("Cm", "Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +
  scale_y_continuous(limits=c(0,110)) +
  mytheme +
  theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1))
predfig.anc1

#ggsave(predfig.anc1,filename="Analysis/Figures/Figure3_AncPredGlmer.pdf", width=15, height=11, units="cm", useDingbats=FALSE)

# Prediction plots for average temperature! ####
anc.td = expand.grid(age_class=unique(data$age_class),
                     samcam = unique(data$samcam),               
                     scl.prec1 = mean(data$scl.prec1),
                     scl.ats1 = seq(min(data$scl.ats1),max(data$scl.ats1), by=0.2),
                     area = 1)


## calculate confidence intervals for predictions from test dataset

# In case of glmmadmb:
    #anc.pred <- cbind(anc.td, predict(anc.best, newdata = anc.td, interval = "confidence")) 

# In case of glmer
    X <- model.matrix(~ age_class*samcam + I(scl.ats1^2) + scl.prec1, data = anc.td)
    anc.td$fit <- X %*% fixef(anc.best)
    anc.td$SE <- sqrt(  diag(X %*%vcov(anc.best) %*% t(X))  )
    anc.td$upr=anc.td$fit+1.96*anc.td$SE
    anc.td$lwr=anc.td$fit-1.96*anc.td$SE
    anc.pred <- anc.td

# Rename samcam for facetting
anc.pred$samcam2 <- plyr::revalue(anc.td$samcam,c("1" ="autumn 2012",  "2" ="spring 2013", "3"="autumn 2013"))
anc.pred$age_class <- plyr::revalue(anc.td$age_class,c("A_Cm"="Cm","B_Sp_young" ="Sp_Y","C_Sp_int1" ="Sp_I1","D_Sp_int2" ="Sp_I2","E_Sp_old" ="Sp_O"))

# Reverse scaling of covariate
anc.pred$ats1 <- anc.pred$scl.ats1* sd(data$ats1) + mean(data$ats1)

predfig.anc2 <- ggplot(anc.pred, aes(x = ats1, y = exp(fit), ymin = exp(lwr), ymax = exp(upr), col=samcam2)) + 
  geom_point() +
  #geom_bar(stat="identity",position = position_dodge(1), col="454545", size=0.15, fill="grey") +
  geom_errorbar(position = position_dodge(1),width=0.15, size=0.15) + 
  facet_grid(.~age_class) +
  geom_hline(xintercept = 1, size=0.15) +
  ylab("Anecic Abundance Ind./m²") +
  xlab(expression(paste("T3",0[surface]))) +
  mytheme +
  theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1))
predfig.anc2
# ggsave(predfig.anc2,filename="Analysis/Figures/Figure3_AncPred2Glmer.pdf", width=15, height=11, units="cm", useDingbats=FALSE)

# Averaged T30(surface) response
anc.pred2 <- aggregate(cbind(fit,lwr, upr) ~ ats1, anc.pred, mean)
colnames(anc.pred2)[2:4] <- c("fit","lwr","upr")

predfig.anc3 <- ggplot(anc.pred2, aes(x = ats1, y = exp(fit), ymin = exp(lwr), ymax = exp(upr))) + 
  geom_point() +
  geom_line() +
  #geom_bar(stat="identity",position = position_dodge(1), col="454545", size=0.15, fill="grey") +
  geom_errorbar(position = position_dodge(1),width=0.15, size=0.15) + 
  geom_hline(xintercept = 1, size=0.15) +
  ylab("Anecic Abundance Ind./m²") +
  xlab(expression(paste("T3",0[surface]))) +
  mytheme +
  theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1))
predfig.anc3
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Coefficients and Statistics ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Claculate a whole lot of coefficients and statistics
anc.pred <- within(anc.pred, {
  AIC <- AIC(anc.best)
  Rrandom <- summary(lm(fitted(anc.best)~ data$anc))$adj.r.squared
  Rsquared <- summary(lm(predict(anc.best, type="response")~ data$anc))$adj.r.squared  
})

anc.stat = round(summary(anc.best)$coefficients[, c(3,4)],4)
anc.coef = coeftab(anc.best)

anc.env <- glmmadmb(anc ~ samcam + scl.prec1 + I(scl.ats1^2) + (1|field.ID) + offset(log(area)),data=data,family="poisson")
anc.pvalue <- anova(anc.env, anc.best)$"Pr(>Chi)"

anc.OUT1 = anc.pred
anc.OUT2 = cbind(anc.stat,anc.coef,LogLikP=anc.pvalue[2],fixef=fixef(anc.best))
anc.OUT1
anc.OUT2

#write.table(anc.OUT1, "Predictions+Stats+ConfIntervals.csv", sep=";", append=TRUE)
#write.table(anc.OUT2, "Predictions+Stats+ConfIntervals.csv", sep=";", append=TRUE)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## Post-Hoc Multicomparisons ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Post-Hoc Multicomparisons for age class x samcam, ####
# interaction effect - reduced contrast matrix 

# Model with the interaction term
data$ia.acl.smc <- interaction(data$age_class, data$samcam)
anc.tukey <- glmer(anc ~ ia.acl.smc + I(scl.ats1^2) + scl.prec1 + (1|field.ID) + offset(log(area)) ,data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))
# summary(anc.tukey)

# Create contrast matrix
source("Analysis/F2_EW_ContrastMatrix.R")

# Pairwise comparisons (with interaction term)
anc.pairwise <- glht(anc.tukey, linfct=mcp(ia.acl.smc = cm1))
anc.pw.ci <- confint(anc.pairwise)
summary(anc.pairwise, test=adjusted(type="fdr"))

# Confidence intervals including 0
anc.pw.sig <- which(anc.pw.ci$confint[,2]>0)
data.frame(names(anc.pw.sig))

# Plot Errorbars
phfig1 <- ggplot(anc.pw.ci, aes(y = lhs, x = exp(estimate), xmin = exp(lwr), xmax = exp(upr))) + 
  geom_errorbarh() + 
  geom_point() + 
  geom_vline(xintercept = 1) +
  mytheme
phfig1

# plot confidence intervals
par(mar=c(2,15,2,2))
plot(anc.pw.ci) 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Post-Hoc Multicomparisons for age class, ####
# main effect 

# Pairwise comparisons (without interaction term)
anc.pairwise <- glht(anc.best, mcp(age_class = "Tukey"))
anc.pw.ci <- confint(anc.pairwise)
summary(anc.pairwise, test=adjusted(type="none"))

# Confidence intervals including 0
anc.pw.sig <- which(anc.pw.ci$confint[,2]>0)
data.frame(names(anc.pw.sig))


# Plot Errorbars
ggplot(anc.pw.ci, aes(y = lhs, x = exp(estimate), xmin = exp(lwr), xmax = exp(upr))) + 
  geom_errorbarh() + 
  geom_point() + 
  geom_vline(xintercept = 1) +
  mytheme

# plot confidence intervals
par(mar=c(2,15,2,2))
plot(anc.pw.ci) 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Post-Hoc Multicomparisons for samcam, ####
# main effect 

# Pairwise comparisons (without interaction term)
anc.pairwise <- glht(anc.best, mcp(samcam = "Tukey"))
anc.pw.ci <- confint(anc.pairwise)

# Confidence intervals including 0
anc.pw.sig <- which(anc.pw.ci$confint[,2]>0)
data.frame(names(anc.pw.sig))


# Plot Errorbars
ggplot(anc.pw.ci, aes(y = lhs, x = exp(estimate), xmin = exp(lwr), xmax = exp(upr))) + 
  geom_errorbarh() + 
  geom_point() + 
  geom_vline(xintercept = 1) +
  mytheme

# plot confidence intervals
par(mar=c(2,15,2,2))
plot(anc.pw.ci) 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%