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
#source("Analysis/F2_EW_RawDataFigures.R")

ancad.raw <-  ggplot(data1.rf[data1.rf$sfg=="ancad",], aes(x=age_class, y=abc.mean, fill=sfg)) +  
  geom_bar(stat="identity", position="dodge") + 
  geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, show_guide=FALSE) + 
  #geom_bar(stat="identity", position="dodge", data=data1.rf2[data1.rf2$sfg!=c("ancad","end"),]) +
  #geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, show_guide=FALSE, data=data1.rf2[data1.rf2$sfg!=c("ancad","end"),]) +
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
        legend.position=c(0.1,0.92))
ancad.raw
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Assess variability in random effects ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Variability within sites
boxplot(ancad~field.ID, data, col="grey", main="Variability within sites", xlab="sites orderd in decreasing age", ylab="anecic earthworm abundance")

# variability within sampling campaigns
boxplot(ancad~samcam, data,col="grey", main="Variability within sampling campaigns", xlab="seasons:\n autumn2012, spring2013, autumn2013", ylab="anecic earthworm abundance")

# variability within sites and sampling campaigns
boxplot(ancad~field.ID+samcam, data, las=2, col="grey", main="Variability within sites \n differing at sampling campaigns", xlab="field+seasons", ylab="anecic earthworm abundance")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



## GLMM using glmer ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Global model Formulation ####

# I Formulate two different candidate models, 
# since some covariates can't be used in the same model due to colinearity

# Model Formulation
ancad.glob1 <- glmer(ancad ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))
ancad.glob2 <- glmer(ancad ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))

# offsset is used due to the personal advice by T. Onkelinx:
# You better use an offset if you want to express the model in terms of m². Just add offset(log(0.25)) to the model. 

summary(ancad.glob1)
summary(ancad.glob2)

# Overdispersion
E1 <- resid(ancad.glob2, type="pearson")
N <- nrow(data)
p <- length(coef(ancad.glob2)) #  "+1"  would be used negative in neg binomial distribution for determining the number of parameters (p) due to the "k" 
Dispersion <- sum(E1^2)/(N-p)
Dispersion

# Multimodel averaging ####

load("Analysis/F2_EW_AbundanceDredge_glmer.RData")

options(na.action = "na.fail")
#ancad.dredge1 <- dredge(ancad.glob1)
head(ancad.dredge1,10)
ancad.avgmod1.d4 <- model.avg(ancad.dredge1, subset = delta < 4)
summary(ancad.avgmod1.d4)
data.frame(importance(ancad.avgmod1.d4))
# AICc range 741-745

#ancad.dredge2 <- dredge(ancad.glob2)
head(ancad.dredge2,10)
ancad.avgmod2.d4 <- model.avg(ancad.dredge2, subset = delta < 4)
summary(ancad.avgmod2.d4)
data.frame(importance(ancad.avgmod2.d4)) 
# AICc range 735-739

save(list=c("ancad.dredge1", "ancad.dredge2"), file="Analysis/OutputTables/ancadImportance.rda")
load( file="Analysis/OutputTables/F2_EW_DredgeAll.rda")



write.csv(data.frame(ancad.avgmod2.d4$importance), "Analysis/OutputTables/ancadImportance.csv")
write.csv(data.frame(ancad.avgmod2.d4$coef.shrinkage), "Analysis/OutputTables/ancadShrinkage.csv")
write.csv(data.frame(ancad.avgmod2.d4$msTable), "Analysis/OutputTables/ancadSubsetModels.csv")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Fit the best model ####
# with both glmer() and glmmadmb():
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ancad.best2 <- glmmadmb(ancad ~ age_class*samcam + I(scl.ats1^2) + scl.prec1 + (1|field.ID) + offset(log(area)) ,data=data,family="poisson")
ancad.best <- glmer(ancad ~ age_class*samcam + I(scl.ats1^2) + scl.ats1 + (1|field.ID) + offset(log(area)) ,data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))

# **The best model includes an Interaction term!!!**

# Summary output ####
# summary
summary(ancad.best)

# anova
lmerTest::anova(ancad.best)
car::Anova(ancad.best)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## Model Validation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Confidence Intervals ####
ancad.confint <- confint(ancad.best)
coefplot2(ancad.best)

# Check Model Assumptions ####

E1 <- resid(ancad.best, type="pearson")
E2 <- resid(ancad.best, type="response")
F1 <- fitted(ancad.best, type="response")
P1 <- predict(ancad.best, type="response")

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
env1 <- cbind(E1, F1, data[,c("ancad", "age_class", "samcam", "ats1", "prec1")]) # should "ats1" be squared?
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
p <- length(coef(ancad.best)) # +1 in case of negbin
Dispersion <- sum(E1^2)/(N-p)
Dispersion
# glmer has less dispersion

overdisp_fun(ancad.best)

# Explained variation ####

r.squaredGLMM(ancad.best)

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
ancad.td = expand.grid(age_class=unique(data$age_class),
                     samcam = unique(data$samcam),               
                     scl.ats1 = mean(data$scl.ats1),
                     area = 4)


## calculate confidence intervals for predictions from test dataset

# In case of glmmadmb:
#ancad.pred <- cbind(ancad.td, predict(ancad.best2, newdata = ancad.td, interval = "confidence")) 

# In case of glmer
X <- model.matrix(~ age_class*samcam + I(scl.ats1^2) + scl.ats1, data = ancad.td)
ancad.td$fit <- X %*% fixef(ancad.best)
ancad.td$SE <- sqrt(  diag(X %*%vcov(ancad.best) %*% t(X))  )
ancad.td$upr=ancad.td$fit+1.96*ancad.td$SE
ancad.td$lwr=ancad.td$fit-1.96*ancad.td$SE
ancad.pred <- ancad.td

# Rename samcam for facetting
ancad.pred$samcam2 <- plyr::revalue(ancad.td$samcam,c("1" ="autumn 2012",  "2" ="spring 2013", "3"="autumn 2013"))
ancad.pred$age_class <- plyr::revalue(ancad.td$age_class,c("A_Cm"="Cm","B_Sp_young" ="Sp_Y","C_Sp_int1" ="Sp_I1","D_Sp_int2" ="Sp_I2","E_Sp_old" ="Sp_O"))

# Reverse scaling of covariate
ancad.pred$ats1 <- ancad.pred$scl.ats1* sd(data$ats1) + mean(data$ats1)

## plot predictions with error bars // confidence intervals???
predfig.ancad1 <- ggplot(ancad.pred, aes(x = age_class, y = exp(fit), ymin = exp(lwr), ymax = exp(upr))) + 
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
predfig.ancad1

#ggsave(predfig.ancad1,filename="Analysis/Figures/Figure3_ancadPredGlmer.pdf", width=15, height=11, units="cm", useDingbats=FALSE)

# Prediction plots for average temperature! ####
ancad.td = expand.grid(age_class=unique(data$age_class),
                     samcam = unique(data$samcam),               
                     scl.prec1 = mean(data$scl.prec1),
                     scl.ats1 = seq(min(data$scl.ats1),max(data$scl.ats1), by=0.2),
                     area = 1)


## calculate confidence intervals for predictions from test dataset

# In case of glmmadmb:
#ancad.pred <- cbind(ancad.td, predict(ancad.best, newdata = ancad.td, interval = "confidence")) 

# In case of glmer
X <- model.matrix(~ age_class*samcam + I(scl.ats1^2) + scl.prec1, data = ancad.td)
ancad.td$fit <- X %*% fixef(ancad.best)
ancad.td$SE <- sqrt(  diag(X %*%vcov(ancad.best) %*% t(X))  )
ancad.td$upr=ancad.td$fit+1.96*ancad.td$SE
ancad.td$lwr=ancad.td$fit-1.96*ancad.td$SE
ancad.pred <- ancad.td

# Rename samcam for facetting
ancad.pred$samcam2 <- plyr::revalue(ancad.td$samcam,c("1" ="autumn 2012",  "2" ="spring 2013", "3"="autumn 2013"))
ancad.pred$age_class <- plyr::revalue(ancad.td$age_class,c("A_Cm"="Cm","B_Sp_young" ="Sp_Y","C_Sp_int1" ="Sp_I1","D_Sp_int2" ="Sp_I2","E_Sp_old" ="Sp_O"))

# Reverse scaling of covariate
ancad.pred$ats1 <- ancad.pred$scl.ats1* sd(data$ats1) + mean(data$ats1)


predfig.ancad2 <- ggplot(ancad.pred, aes(x = ats1, y = exp(fit), ymin = exp(lwr), ymax = exp(upr), col=samcam2)) + 
  geom_point() +
  #geom_bar(stat="identity",position = position_dodge(1), col="454545", size=0.15, fill="grey") +
  geom_errorbar(position = position_dodge(1),width=0.15, size=0.15) + 
  facet_grid(.~age_class) +
  geom_hline(xintercept = 1, size=0.15) +
  ylab("Anecic Abundance Ind./m²") +
  xlab(expression(paste("T3",0[surface]))) +
  mytheme +
  theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1))
predfig.ancad2
# ggsave(predfig.ancad2,filename="Analysis/Figures/Figure3_ancadPred2Glmer.pdf", width=15, height=11, units="cm", useDingbats=FALSE)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Coefficients and Statistics ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Claculate a whole lot of coefficients and statistics
ancad.pred <- within(ancad.pred, {
  AIC <- AIC(ancad.best)
  Rrandom <- summary(lm(fitted(ancad.best)~ data$ancad))$adj.r.squared
  Rsquared <- summary(lm(predict(ancad.best, type="response")~ data$ancad))$adj.r.squared  
})

ancad.stat = round(summary(ancad.best)$coefficients[, c(3,4)],4)
ancad.coef = coeftab(ancad.best)

ancad.env <- glmmadmb(ancad ~ samcam + scl.prec1 + I(scl.ats1^2) + (1|field.ID) + offset(log(area)),data=data,family="poisson")
ancad.pvalue <- anova(ancad.env, ancad.best)$"Pr(>Chi)"

ancad.OUT1 = ancad.pred
ancad.OUT2 = cbind(ancad.stat,ancad.coef,LogLikP=ancad.pvalue[2],fixef=fixef(ancad.best))
ancad.OUT1
ancad.OUT2

#write.table(ancad.OUT1, "Predictions+Stats+ConfIntervals.csv", sep=";", append=TRUE)
#write.table(ancad.OUT2, "Predictions+Stats+ConfIntervals.csv", sep=";", append=TRUE)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## Post-Hoc Multicomparisons ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Post-Hoc Multicomparisons for age class x samcam, ####
# interaction effect - reduced contrast matrix 

# Model with the interaction term
data$ia.acl.smc <- interaction(data$age_class, data$samcam)
ancad.tukey <- glmer(ancad ~ ia.acl.smc + I(scl.ats1^2) + scl.prec1 + (1|field.ID) + offset(log(area)) ,data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))
# summary(ancad.tukey)

# Create contrast matrix
source("Analysis/F2_EW_ContrastMatrix.R")

# Pairwise comparisons (with interaction term)
ancad.pairwise <- glht(ancad.tukey, linfct=mcp(ia.acl.smc = cm1))
ancad.pw.ci <- confint(ancad.pairwise)
summary(ancad.pairwise, test=adjusted(type="fdr"))

# Confidence intervals including 0
ancad.pw.sig <- which(ancad.pw.ci$confint[,2]>0)
data.frame(names(ancad.pw.sig))

# Plot Errorbars
phfig1 <- ggplot(ancad.pw.ci, aes(y = lhs, x = exp(estimate), xmin = exp(lwr), xmax = exp(upr))) + 
  geom_errorbarh() + 
  geom_point() + 
  geom_vline(xintercept = 1) +
  mytheme
phfig1

# plot confidence intervals
par(mar=c(2,15,2,2))
plot(ancad.pw.ci) 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Post-Hoc Multicomparisons for age class, ####
# main effect 

# Pairwise comparisons (without interaction term)
ancad.pairwise <- glht(ancad.best, mcp(age_class = "Tukey"))
ancad.pw.ci <- confint(ancad.pairwise)
summary(ancad.pairwise, test=adjusted(type="none"))

# Confidence intervals including 0
ancad.pw.sig <- which(ancad.pw.ci$confint[,2]>0)
data.frame(names(ancad.pw.sig))


# Plot Errorbars
ggplot(ancad.pw.ci, aes(y = lhs, x = exp(estimate), xmin = exp(lwr), xmax = exp(upr))) + 
  geom_errorbarh() + 
  geom_point() + 
  geom_vline(xintercept = 1) +
  mytheme

# plot confidence intervals
par(mar=c(2,15,2,2))
plot(ancad.pw.ci) 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Post-Hoc Multicomparisons for samcam, ####
# main effect 

# Pairwise comparisons (without interaction term)
ancad.pairwise <- glht(ancad.best, mcp(samcam = "Tukey"))
ancad.pw.ci <- confint(ancad.pairwise)

# Confidence intervals including 0
ancad.pw.sig <- which(ancad.pw.ci$confint[,2]>0)
data.frame(names(ancad.pw.sig))


# Plot Errorbars
ggplot(ancad.pw.ci, aes(y = lhs, x = exp(estimate), xmin = exp(lwr), xmax = exp(upr))) + 
  geom_errorbarh() + 
  geom_point() + 
  geom_vline(xintercept = 1) +
  mytheme

# plot confidence intervals
par(mar=c(2,15,2,2))
plot(ancad.pw.ci) 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%