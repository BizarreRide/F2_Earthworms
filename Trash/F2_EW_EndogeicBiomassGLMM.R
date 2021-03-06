#################
# F2_Eearthworms
# GLMM for Anecic Biomass
# Quentin Schorpp
# 07.05.2015
#################


# Standardize Covariates ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Data/GatherSource/F2_EW_MakeLikeFile.R")
source("Data/GatherSource/CovariateStandardization.R")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# plot of anecic biomass ####
# The data we want to model
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# !!!! requires data from raw figure plots
source("Analysis/F2_EW_RawDataFigures.R")
source("Data/GatherSource/CovariateStandardization.R")

endad.bm.raw <- ggplot(data1.rf[data1.rf$sfg.bm=="endad.bm",], aes(x=age_class, y=bm.mean, fill=sfg.bm)) +  
  geom_bar(stat="identity", position="dodge") + 
  geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, show_guide=FALSE) + 
  #geom_bar(stat="identity", position="dodge", data=data1.rf2[data1.rf2$sfg!=c("anc","end"),]) +
  #geom_bar(stat="identity", position="dodge", colour="#454545", size=0.15, show_guide=FALSE, data=data1.rf2[data1.rf2$sfg!=c("anc","end"),]) +
  geom_errorbar(aes(ymin=bm.mean-1.96*bm.se, ymax=bm.mean+1.96*bm.se), position=position_dodge(0.9),width=0.15, size=0.15) +
  facet_grid(.~samcam) +
  xlab("Age Class") + 
  ylab(expression(paste("Biomass \u00B1 CI ","[g x ",0.25,m^-2," ]"))) +
  #ylim(-10,max(data1.rf$abc.mean+data1.rf$abc.se)) +
  labs(fill="Functional Group") +
  scale_fill_grey(labels=c("anecic juvenile","anecic adult","endogeic juvenile", "endogeic adult","epigeic", "total")) +
  scale_y_continuous(breaks=pretty_breaks(n=10)) +
  scale_x_discrete(labels=c("Cm", "Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +  
  mytheme +
  guides(fill=guide_legend(keywidth=0.5, keyheight=0.5)) +
  theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1),
        legend.title=element_text(size=6),
        legend.text=element_text(size=7),
        legend.position=c(0.12,0.92))

endad.bm.raw
#ggsave(endad.bm.raw, filename="Analysis/Figures/Figure7_EndadBmRaw.pdf", width=16.5, height=11, units="cm", useDingbats=FALSE)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Assess variability in random effects ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Variability within sites
boxplot(endad.bm ~ field.ID, data, col="grey", main="Variability within sites", xlab="sites orderd in decreasing age", ylab="endogeic earthworm Biomass")

# variability within sampling campaigns
boxplot(endad.bm ~ samcam, data,col="grey", main="Variability within sampling campaigns", xlab="seasons:\n autumn2012, spring2013, autumn2013", ylab="endogeic earthworm Biomass")

# variability within sites and sampling campaigns
boxplot(endad.bm ~ field.ID+samcam, data, las=2, col="grey", main="Variability within sites \n differing at sampling campaigns", xlab="field+seasons", ylab="endogeic earthworm biomass")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



## GLMM using glmer ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Global model Formulation ####

# I Formulate two different candidate models, 
# since some covariates can't be used in the same model due to colinearity

# Model Formulation
endad.bm.glob1 <- lmer(log1p(endad.bm) ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc:scl.pH + + scl.ats1 + I(scl.ats1^2) + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data)
endad.bm.glob0 <- lmer(log1p(endad.bm) ~ 1 + (1|field.ID) + offset(log(area)),data=data)
# offsset is used due to the personal advice by T. Onkelinx:
# You better use an offset if you want to express the model in terms of m². Just add offset(log(0.25)) to the model. 

summary(endad.bm.glob1)

# Overdispersion
E1 <- resid(endad.bm.glob1, type="pearson")
N <- nrow(data)
p <- length(coef(endad.bm.glob1)) #  "+1"  would be used negative in neg binomial distribution for determining the number of parameters (p) due to the "k" 
Dispersion <- sum(E1^2)/(N-p)
Dispersion

# Multimodel averaging ####

load("Analysis/F2_EW_BiomassDredge.RData")

## Suppose we want to have set of models that exclude combinations of colinear
# variables, that are significantly (p < 0.05) correlated, with Pearson
# correlation coefficient larger than r = 0.5.

std.var2 <- std.var[,c("scl.ats1","scl.cn", "scl.prec1", "scl.clay", "scl.mc", "scl.pH")]
std.var2 <- cbind(std.var2, "I(scl.mc^2)"=std.var$scl.mc^2,"I(scl.ats1^2)" = std.var$scl.ats1^2)#,age_class=data$age_class, age_class=data$samcam)

# Create logical matrix
# Use functions from the MakeLikeFile!
smat <- outer(1:8, 1:8, vCorrelated, data = std.var2)

nm <- colnames(std.var2[1:8])

# Although the squared terms seem not to correlate between mc and ats, we won't fit them within the same model
smat[5,1] <- FALSE
smat[7,1]<- FALSE
smat[8,5]<- FALSE
smat[8,7]<- FALSE

smat

## Compute all possible models 
#endad.bm.dredge1 <- dredge(endad.bm.glob1, subset=smat)
head(endad.bm.dredge1,10)
endad.bm.avgmod1.d4 <- model.avg(endad.bm.dredge1, subset = delta < 4)
summary(endad.bm.avgmod1.d4)

# Null Model:
which(anc.bm.dredge1[2]==NA,)
anc.bm.dredge1[,2:14]

write.csv(data.frame(endad.bm.avgmod1.d4$importance), "Analysis/OutputTables/EndadBmImportance.csv")
write.csv(data.frame(endad.bm.avgmod1.d4$coef.shrinkage), "Analysis/OutputTables/EndadBmShrinkage.csv")
write.csv(data.frame(endad.bm.avgmod1.d4$msTable), "Analysis/OutputTables/EndadBmSubsetModels.csv")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Fit the best model ####
# with glmer()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#endad.bm.best <- lmer(log1p(endad.bm) ~ age_class + samcam + I(scl.mc^2) + scl.pH + scl.prec1 + (1|field.ID) + offset(log(area)) ,data=data)
endad.bm.best <- lmer(log1p(endad.bm) ~ age_class + samcam + I(scl.mc^2) + scl.pH + scl.prec1 + (1|field.ID) ,data=data)

# **The best model includes an Interaction term!!!**

# Summary output ####
# summary
summary(endad.bm.best)

# anova
car::Anova(endad.bm.best, type="III")
summary(aov(endad.bm.best))

write.csv(summary(endad.bm.best)$coefficients, "Analysis/OutputTables/EndadBmBestCoef.csv")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## Model Validation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Confidence Intervals ####
endad.bm.confint <- confint(endad.bm.best)
coefplot2(endad.bm.best)
#write.csv(data.frame(endad.bm.confint), "Analysis/OutputTables/EndadBmConfint.csv")

# Check Model Assumptions ####

E1 <- resid(endad.bm.best, type="response")
F1 <- fitted(endad.bm.best, type="response")
P1 <- predict(endad.bm.best, type="response")

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
qqnorm(E1)
qqline(E1)

# Cooks Distances

par(lo)

# plot residuals vs. all covariates ####

# Dataset with all variables IN the model
env1 <- cbind(E1, F1, data[,c("anc", "age_class", "samcam", "mc", "pH","prec1")]) # should "ats1" be squared?
# covariates NOT in the model
env0 <- cbind(E1, F1,data[,c("cn", "clay", "ats2", "hum1","dgO")])

# in the model
par(mfrow=c(3,3),
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
p <- length(coef(endad.bm.best)) # +1 in case of negbin
Dispersion <- sum(E1^2)/(N-p)
Dispersion

overdisp_fun(endad.bm.best)

# Explained variation ####

r.squaredGLMM(endad.bm.best)

# Some more plots of fitted and predicted values####
# Plots of Predicted Values
par(mfrow=c(2,2))
plot(data$age_class,P1)
plot(data$samcam,P1)
plot(data$ats1^2,P1)
plot(data$prec1,P1)

# Plots of fitted Values
par(mfrow=c(2,2))
plot(data$age_class,F1)
plot(data$samcam,F1)
plot(data$ats1^2,F1)
plot(data$prec1,F1)

par(lo)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## Prediction plots ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Prediction plots for Age Clas x SamCam ####

## create test data set with all covariates IN the model
# to predict for age_class only, take the mean of all continous covariates
endad.bm.td = expand.grid(age_class=unique(data$age_class),
                        samcam = unique(data$samcam),
                        scl.mc = mean(data$scl.ats1),
                        scl.pH = mean(data$scl.ats1),
                        scl.prec1 = mean(data$scl.cn))
                        #area = 1)


## calculate confidence intervals for predictions from test dataset

# In case of glmmadmb:
#endad.bm.pred <- cbind(endad.bm.td, predict(endad.bm.best2, newdata = endad.bm.td, interval = "confidence")) 

# In case of glmer
X <- model.matrix(~ age_class + samcam + scl.mc + scl.pH + scl.prec1, data = endad.bm.td)
endad.bm.td$fit <- X %*% fixef(endad.bm.best)
endad.bm.td$SE <- sqrt(  diag(X %*%vcov(endad.bm.best) %*% t(X))  )
endad.bm.td$upr=endad.bm.td$fit+1.96*endad.bm.td$SE
endad.bm.td$lwr=endad.bm.td$fit-1.96*endad.bm.td$SE
endad.bm.pred <- endad.bm.td

# Rename samcam for facetting
endad.bm.pred$samcam2 <- plyr::revalue(endad.bm.td$samcam,c("1" ="autumn 2012",  "2" ="spring 2013", "3"="autumn 2013"))
endad.bm.pred$age_class <- plyr::revalue(endad.bm.td$age_class,c("A_Cm"="Cm","B_Sp_young" ="Sp_Y","C_Sp_int1" ="Sp_I1","D_Sp_int2" ="Sp_I2","E_Sp_old" ="Sp_O"))

# Reverse scaling of covariate
endad.bm.pred$prec1 <- endad.bm.pred$scl.prec1* sd(data$prec1) + mean(data$prec1)
endad.bm.pred$mc <- endad.bm.pred$scl.mc* sd(data$mc) + mean(data$mc)
endad.bm.pred$pH <- endad.bm.pred$scl.pH* sd(data$pH) + mean(data$pH)

## plot predictions with error bars // confidence intervals???
predfig.endad.bm1 <- ggplot(endad.bm.pred, aes(x = age_class, y = exp(fit), ymin = exp(lwr), ymax = exp(upr))) + 
  geom_bar(stat="identity",position = position_dodge(1), col="454545", size=0.15, fill="grey") +
  geom_errorbar(position = position_dodge(1),col="black",width=0.15, size=0.15) + 
  facet_grid(.~samcam2) +
  geom_hline(xintercept = 1, size=0.15) +
  ylab(expression(paste("Biomass \u00B1 CI ","[g x ",0.25,m^-2," ]"))) +
  xlab("Age Class") +
  scale_x_discrete(labels=c("Cm", "Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +
  #scale_y_log10() +
  scale_y_continuous( limits=c(0,25), breaks=pretty_breaks()) +
  mytheme +
  theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1))
predfig.endad.bm1

ggsave(predfig.endad.bm1,filename="Analysis/Figures/Figure7_EndadBmPredGlmer.pdf", width=16.5, height=11, units="cm", useDingbats=FALSE)

# Prediction plots for average temperature! ####
endad.bm.td = expand.grid(age_class=unique(data$age_class),
                        samcam = unique(data$samcam),
                        scl.prec1 = seq(min(data$scl.prec1),max(data$scl.prec1), length.out=5),
                        scl.pH = seq(min(data$scl.pH),max(data$scl.pH), length.out=5),
                        scl.mc = seq(min(data$scl.mc),max(data$scl.mc), length.out=5),
                        area = 1)


## calculate confidence intervals for predictions from test dataset

# In case of glmmadmb:
#endad.bm.pred <- cbind(endad.bm.td, predict(endad.bm.best, newdata = endad.bm.td, interval = "confidence")) 

# In case of glmer
X <- model.matrix(~ age_class + samcam + scl.mc + scl.pH + scl.prec1, data = endad.bm.td)
endad.bm.td$fit <- X %*% fixef(endad.bm.best)
endad.bm.td$SE <- sqrt(  diag(X %*%vcov(endad.bm.best) %*% t(X))  )
endad.bm.td$upr=endad.bm.td$fit+1.96*endad.bm.td$SE
endad.bm.td$lwr=endad.bm.td$fit-1.96*endad.bm.td$SE
endad.bm.pred <- endad.bm.td

# Rename samcam for facetting
endad.bm.pred$samcam2 <- plyr::revalue(endad.bm.td$samcam,c("1" ="autumn 2012",  "2" ="spring 2013", "3"="autumn 2013"))
endad.bm.pred$age_class <- plyr::revalue(endad.bm.td$age_class,c("A_Cm"="Cm","B_Sp_young" ="Sp_Y","C_Sp_int1" ="Sp_I1","D_Sp_int2" ="Sp_I2","E_Sp_old" ="Sp_O"))

# Reverse scaling of covariate
endad.bm.pred$prec1 <- endad.bm.pred$scl.prec1* sd(data$prec1) + mean(data$prec1)
endad.bm.pred$mc <- endad.bm.pred$scl.mc* sd(data$mc) + mean(data$mc)
endad.bm.pred$pH <- endad.bm.pred$scl.pH* sd(data$pH) + mean(data$pH)

predfig.endad.bm2 <- ggplot(endad.bm.pred, aes(x = pH, y = exp(fit), ymin = exp(lwr), ymax = exp(upr), col=samcam)) + 
  geom_point() +
  #geom_bar(stat="identity",position = position_dodge(1), col="454545", size=0.15, fill="grey") +
  geom_errorbar(position = position_dodge(1),width=0.15, size=0.15) + 
  facet_grid(mc~age_class) +
  geom_hline(xintercept = 1, size=0.15) +
  ylab("Biomass [g]") +
  xlab("Gravimetric Soil Water Content") +
  mytheme +
  theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1))
predfig.endad.bm2
# ggsave(predfig.endad.bm2,filename="Analysis/Figures/Figure3_AncPred2Glmer.pdf", width=15, height=11, units="cm", useDingbats=FALSE)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Coefficients and Statistics ####
# Claculate a whole lot of coefficients and statistics
endad.bm.bm.pred <- within(endad.bm.bm.pred, {
  AIC <- AIC(endad.bm.bm.best)
  Rrandom <- summary(lm(fitted(endad.bm.bm.best)~ data$anc))$adj.r.squared
  Rsquared <- summary(lm(predict(endad.bm.bm.best, type="response")~ data$anc))$adj.r.squared  
})

endad.bm.bm.stat = round(summary(endad.bm.bm.best)$coefficients[, c(3,4)],4)
endad.bm.coef = coeftab(endad.bm.best)

endad.bm.env <- glmmadmb(anc ~ samcam + scl.prec1 + I(scl.ats1^2) + (1|field.ID) + offset(log(area)),data=data,family="poisson")
endad.bm.pvalue <- anova(endad.bm.env, endad.bm.best)$"Pr(>Chi)"

endad.bm.OUT1 = endad.bm.pred
endad.bm.OUT2 = cbind(endad.bm.stat,endad.bm.coef,LogLikP=endad.bm.pvalue[2],fixef=fixef(endad.bm.best))
endad.bm.OUT1
endad.bm.OUT2

#write.table(endad.bm.OUT1, "Predictions+Stats+ConfIntervals.csv", sep=";", append=TRUE)
#write.table(endad.bm.OUT2, "Predictions+Stats+ConfIntervals.csv", sep=";", append=TRUE)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## Post-Hoc Multicomparisons ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Post-Hoc Multicomparisons for age class x samcam, ####
# interaction effect - reduced contrast matrix 

# Model with the interaction term
data$ia.acl.smc <- interaction(data$age_class, data$samcam)
endad.bm.tukey <- glmer(anc ~ ia.acl.smc + I(scl.ats1^2) + scl.prec1 + (1|field.ID) + offset(log(area)) ,data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))
# summary(endad.bm.tukey)

# Create contrast matrix
source("Analysis/F2_EW_ContrastMatrix.R")

# Pairwise comparisons (with interaction term)
endad.bm.pairwise <- glht(endad.bm.tukey, linfct=mcp(ia.acl.smc = cm1))
endad.bm.pw.ci <- confint(endad.bm.pairwise)
summary(endad.bm.pairwise, test=adjusted(type="fdr"))

# Confidence intervals including 0
endad.bm.pw.sig <- which(endad.bm.pw.ci$confint[,2]>0)
data.frame(names(endad.bm.pw.sig))

# Plot Errorbars
phfig1 <- ggplot(endad.bm.pw.ci, aes(y = lhs, x = exp(estimate), xmin = exp(lwr), xmax = exp(upr))) + 
  geom_errorbarh() + 
  geom_point() + 
  geom_vline(xintercept = 1) +
  mytheme
phfig1

# plot confidence intervals
par(mar=c(2,15,2,2))
plot(endad.bm.pw.ci) 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Post-Hoc Multicomparisons for age class, ####
# main effect 

# Pairwise comparisons (without interaction term)
endad.bm.pairwise <- glht(endad.bm.best, mcp(age_class = "Tukey"))
endad.bm.pw.ci <- confint(endad.bm.pairwise)
summary(endad.bm.pairwise, test=adjusted(type="none"))

# Confidence intervals including 0
endad.bm.pw.sig <- which(endad.bm.pw.ci$confint[,2]>0)
data.frame(names(endad.bm.pw.sig))


# Plot Errorbars
ggplot(endad.bm.pw.ci, aes(y = lhs, x = exp(estimate), xmin = exp(lwr), xmax = exp(upr))) + 
  geom_errorbarh() + 
  geom_point() + 
  geom_vline(xintercept = 1) +
  mytheme

# plot confidence intervals
par(mar=c(2,15,2,2))
plot(endad.bm.pw.ci) 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Post-Hoc Multicomparisons for samcam, ####
# main effect 

# Pairwise comparisons (without interaction term)
endad.bm.pairwise <- glht(endad.bm.best, mcp(samcam = "Tukey"))
endad.bm.pw.ci <- confint(endad.bm.pairwise)

# Confidence intervals including 0
endad.bm.pw.sig <- which(endad.bm.pw.ci$confint[,2]>0)
data.frame(names(endad.bm.pw.sig))


# Plot Errorbars
ggplot(endad.bm.pw.ci, aes(y = lhs, x = exp(estimate), xmin = exp(lwr), xmax = exp(upr))) + 
  geom_errorbarh() + 
  geom_point() + 
  geom_vline(xintercept = 1) +
  mytheme

# plot confidence intervals
par(mar=c(2,15,2,2))
plot(endad.bm.pw.ci) 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



