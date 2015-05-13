#################
# F2_Eearthworms
# GLMM for Endogeic Abundance
# Quentin Schorpp
# 07.05.2015
#################


##############################################################
# Standardize Covariates
source("Data/GatherSource/CovariateStandardization.R")
##############################################################

# I model only adult endogeics, since this decreases the varaince and buffers very large abundance values.


##############################################################
# Assess variability in random effects

# Variability within sites
boxplot(endad~field.ID, data, col="grey", main="Variability within sites", xlab="sites orderd in decreasing age", ylab="anecic earthworm abundance")

# variability within sampling campaigns
boxplot(endad~samcam, data,col="grey", main="Variability within sampling campaigns", xlab="seasons:\n autumn2012, spring2013, autumn2013", ylab="anecic earthworm abundance")

# variability within sites and sampling campaigns
boxplot(endad~field.ID+samcam, data, las=2, col="grey", main="Variability within sites \n differing at sampling campaigns", xlab="field+seasons", ylab="anecic earthworm abundance")
##############################################################

## GLMM using glmmADMB ####

##############################################################
# Global model Formulation

# I Formulate two different candidate models, 
# since some covariates can't be used in the same model due to colinearity

# Model Formulation
endad.glob1 <- glmmadmb(endad ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.pH*scl.cn  + scl.hum1 + I(scl.hum1^2) + (1|field.ID) + offset(log(area)),data=data,family="poisson")
endad.glob2 <- glmmadmb(endad ~ age_class*samcam + scl.ats1 + I(scl.ats1^2) + scl.pH*scl.cn  +I(scl.hum1^2) + (1|field.ID) + offset(log(area)),data=data,family="poisson")
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
#endad.dredge1 <- dredge(endad.glob1)
head(endad.dredge1,10)
endad.avgmod1.d4 <- model.avg(endad.dredge1, subset = delta < 4)
summary(endad.avgmod1.d4)
importance(endad.avgmod.d4) 
endad.avgmod.95p <- model.avg(endad.dredge1, cumsum(weight) <= .95)
# AICc begins at 751.23

#endad.dredge2 <- dredge(endad.glob2)
head(endad.dredge2,10)
endad.avgmod2.d4 <- model.avg(endad.dredge2, subset = delta < 4)
summary(endad.avgmod2.d4)
importance(endad.avgmod.d4) 
endad.avgmod.95p <- model.avg(endad.dredge2, cumsum(weight) <= .95)
# AICc ends at 743.64

save(endad.dredge1, endad.dredge2, file="Analysis/AnecicDredge.RData")
rm(endad.dredge1)
load("Analysis/AnecicDredge.RData")
#####
##############################################################


##############################################################
# Fit the best model
# with both glmer() and glmmadmb():

#endad.best <- glmer(endad ~ age_class*samcam + I(scl.ats1^2) + (1|field.ID) + offset(log(area)) ,data=data,family="poisson")
endad.best <- glmmadmb(endad ~ age_class*samcam + I(scl.ats1^2) +  (1|field.ID) + offset(log(area)) ,data=data,family="poisson")

# **The best model includes an Interaction term!!!**


# Summary output 
# summary
summary(endad.best)

# anova
summary(aov(endad.best))
##############################################################


##############################################################
# Model Validation

## Confidence Intervals ####
confint(endad.best)
coefplot2(endad.best)


## Residual plots ####

E1 <- resid(endad.best, type="pearson")
F1 <- fitted(endad.best, type="response")
P1 <- resid(endad.best, type="response")


## Check Model assumptions ####

par(mfrow=c(2,2),
    mar=c(4,4.5,1,2))
# Plot fitted vs. residuals
scatter.smooth(F1, E1, cex.lab = 1.5, ylab=" Residuals", xlab="Fitted values")
abline(h = 0, v=0, lty=2)
data$field[[69]]

# plot fitted vs. predicted
scatter.smooth(F1, P1, cex.lab = 1.5, ylab=" Residuals", xlab="Predicted values")
abline(h = 0, v=0, lty=2)

# Histogram of Residuals
hist(E1, prob=TRUE, main = "", breaks = 20, cex.lab = 1.5, xlab = "Response Residuals", col="PapayaWhip")
lines(density(E1), col="light blue", lwd=3)
lines(density(E1, adjust=2), lty="dotted", col="darkgreen", lwd=2) 

# Normal QQ Plot
qqplot(E1)

# Cooks Distances

par(lo)
#####


# Dataset with all variables IN the model
env1 <- cbind(E1, F1, data[,c("endad", "age_class", "samcam", "ats1")]) # should "ats1" be squared?
# covariates NOT in the model
env0 <- cbind(E1, F1,data[,c("mc", "pH", "cn", "clay", "ats2", "hum1","dgO")])

## plot residual versus all covariates in the model ####

par(mfrow=c(2,3),
    mar=c(3.8,4,1,3))
for (i in 1:length(env1)) {
  scatter.smooth(env1[,i],env1[,1], ylab="Residuals", xlab=colnames(env1)[i])
  abline(h=0, lty=2, col="red")
}
par(mfrow=c(1,1))

# plot residual versus all covariates NOT in the model
par(mfrow=c(3,3),
    mar=c(3.8,4,1,3))
for (i in 1:length(env0)) {
  scatter.smooth(env0[,i],env0[,1], ylab="Residuals", xlab=colnames(env0)[i])
  abline(h=0, lty=2, col="red")
}
par(lo)
##### 

## Overdispersion ####
N <- nrow(data)
p <- length(coef(endad.best)) # +1 in case of negbin
Dispersion <- sum(E1^2)/(N-p)
Dispersion
# glmer has less dispersion

overdisp_fun(endad.best)
#####

## Explained variation ####

r.squaredGLMM(endad.best)

#####

## Some more plots ####
# Plots of Predicted Values
par(mfrow=c(2,2))
plot(data$age_class,predict(endad.best, type="response"))
plot(data$samcam,predict(endad.best, type="response"))
plot(data$ats1^2,predict(endad.best, type="response"))
plot(data$hum1^2,predict(endad.best, type="response"))

# Plots of fitted Values
par(mfrow=c(2,2))
plot(data$age_class,F1)
plot(data$samcam,F1)
plot(data$ats1^2,F1)

par(lo)
##############################################################


##############################################################
#Post-Hoc Multicomparisons for age class

## Model with the interaction term ####
data$ia.acl.smc <- interaction(data$age_class, data$samcam)
endad.tukey <- glmmadmb(endad ~ ia.acl.smc + I(scl.ats1^2) + (1|field.ID) + offset(log(area)) ,data=data,family="poisson")
# summary(endad.tukey)

# Create contrast matrix
k <- 2 # pairwise comparison
n <- 5 # Nr of facot levels to be compared, Interaction factor 1
groups <- 3 # Nr of groups to draw comparisons within Interaction factor 2

x <- groups*choose(n,k) # row number, i.e. all comparisons that should be drawn

cm1 <- matrix(0,x,length(levels(data$ia.acl.smc))+2) # empty contrast matrix with two mor columns for pairwise factor combinations

# Fill in pairwise factor combinations, first create interaction factor!
for(i in 1:2) {
  y <- rbind(t(combn(levels(data$ia.acl.smc)[1:5],k)),
             t(combn(levels(data$ia.acl.smc)[6:10],k)),
             t(combn(levels(data$ia.acl.smc)[11:15],k)))
  cm1[,i] <- y[,i]
}
cm1 <- data.frame(cm1, stringsAsFactors=FALSE) # turn into data frame


colnames(cm1)[3:17] <- levels(data$ia.acl.smc) # fill in column names

# write 1 or -1 if colnames match the names of rows 1 or 2
for (i in 3:17) {
  for (j in 1:30) {
    if(colnames(cm1)[i]==cm1$X1[j]) {cm1[j,i] = -1}
    if(colnames(cm1)[i]==cm1$X2[j]) {cm1[j,i] = 1}
  }
}

rownames(cm1) <- paste(cm1$X1, cm1$X2, sep=" - ") # create rownames
cm1= as.matrix(cm1[,-c(1,2)]) # delete columns 1 and 2
class(cm1) <- "numeric"

# Pairwise comparisons (with interaction term)
endad.pairwise <- glht(endad.tukey, linfct=mcp(ia.acl.smc = cm1))
endad.pw.ci <- confint(endad.pairwise)
summary(endad.pairwise, test=adjusted(type="none"))

# Confidence intervals including 0
endad.pw.sig <- which(endad.pw.ci$confint[,2]>0)
data.frame(names(endad.pw.sig))
#####


## Pairwise comparisons (without interaction term) ####
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
plot(endad.pw.ci) # only maize is significantly different from all SIlphie fields. Within SIlphie there are no differences
##############################################################


##############################################################
#Post-Hoc Multicomparisons for samcam

## Model with the interaction term ####
data$ia.acl.smc <- interaction(data$age_class, data$samcam)
endad.tukey <- glmmadmb(endad ~ ia.acl.smc + I(scl.ats1^2) + (1|field.ID) + offset(log(area)) ,data=data,family="poisson")
# summary(endad.tukey)

# Pairwise comparisons (with interaction term)
endad.pairwise <- glht(endad.tukey, mcp(ia.acl.smc = "Tukey"))
endad.pw.ci <- confint(endad.pairwise)

endad.contrast <- matrix(30,15,0)

# Confidence intervals including 0
endad.pw.sig <- which(endad.pw.ci$confint[,2]>0)
data.frame(names(endad.pw.sig))
#####


## Pairwise comparisons (without interaction term) ####
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
plot(endad.pw.ci) # only maize is significantly different from all SIlphie fields. Within SIlphie there are no differences
##############################################################


##############################################################
# Extract Predictions and plot predictions with error bars

# create test data set with all covariates IN the model
# to predict for age_class only, take the mean of all continous covariates
endad.td = expand.grid(age_class=unique(data$age_class),
                     samcam = unique(data$samcam),               
                     scl.ats1 = mean(data$scl.ats1),
                     area = 1)

# calculate confidence intervals for predictions from test dataset
endad.pred <- cbind(endad.td, predict(endad.best, newdata = endad.td, interval = "confidence")) 

# Rename samcam for facetting
endad.pred$samcam2 <- plyr::revalue(endad.td$samcam,c("1" ="autumn 2012",  "2" ="spring 2013", "3"="autumn 2013"))

# Reverse scaling of covariate
endad.pred$ats1 <- endad.pred$scl.ats1 * attr(endad.pred$scl.ats1, 'scaled:scale') + attr(endad.pred$scl.ats1, 'scaled:center')

# plot predictions with error bars // confidence intervals???
predfig1 <- ggplot(endad.pred, aes(x = age_class, y = exp(fit), ymin = exp(lwr), ymax = exp(upr))) + 
  geom_bar(stat="identity",position = position_dodge(1), col="454545", size=0.15, fill="grey") +
  geom_errorbar(position = position_dodge(1),col="black",width=0.15, size=0.15) + 
  facet_grid(.~samcam2) +
  geom_hline(xintercept = 1, size=0.15) +
  ylab("Anecic Abundance Ind./m²") +
  xlab("Age Class") +
  scale_x_discrete(labels=c("Cm", "Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +
  mytheme +
  theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1))

#ggsave(predfig1,filename="Analysis/Figures/Figure3.pdf", width=15, height=11, units="cm", useDingbats=FALSE)
##############################################################


##############################################################
Coefficients and Statistics
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

##### Prediction plots for average temperature!
endad.td = expand.grid(age_class=unique(data$age_class),
                     samcam = unique(data$samcam),               
                     scl.ats1 = seq(min(data$scl.ats1),max(data$scl.ats1), by=0.2),
                     area = 1)

endad.pred <- cbind(endad.td, predict(endad.best, newdata = endad.td, interval = "confidence")) 

ggplot(endad.pred, aes(x = scl.ats1, y = exp(fit), ymin = exp(lwr), ymax = exp(upr),col=samcam)) + 
  geom_errorbar(position = position_dodge(1)) + 
  geom_point(position = position_dodge(1)) +
  facet_grid(.~age_class) +
  mytheme
