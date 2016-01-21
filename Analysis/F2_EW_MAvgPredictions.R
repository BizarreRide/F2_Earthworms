#§§§§§§§§§§§§§§§§§§§§§§§§§§§
# F2_Eearthworms
# GLMM with Dredge
# Quentin Schorpp
# 29.12.2015
#§§§§§§§§§§§§§§§§§§§§§§§§§§§

# Load Data and packages ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Data/GatherSource/F2_EW_MakeLikeFile.R")
library(MuMIn)
library(lme4)
library(data.table)
library(coefplot2)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Data Processing ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Data/GatherSource/CovariateStandardization.R")

abbreviate(c("response", "explanatory", "abundance", "biomass", "biodiversity"),3)

# Response variables for Abundances
dt.rsp.abn <- as.data.table(data[,c("anc", "ancad", "endo", "endad", "N")])

#Juveniles
dt.rsp.abn[, anc.juv:=anc-ancad]
dt.rsp.abn[, endo.juv:=endo-endad]
dt.rsp.abn[,juv:= anc.juv + endo.juv]


# Response variables for biomass
dt.rsp.bms <- as.data.table(data[,c("anc.bm", "ancad.bm", "endo.bm", "endad.bm", "N.bm")])

# Juveniles
dt.rsp.bms[, anc.juv.bm:=anc.bm-ancad.bm]
dt.rsp.bms[, endo.juv.bm:=endo.bm-endad.bm]
dt.rsp.bms[,juv.bm:= anc.juv.bm + endo.juv.bm]


# Response variables for biodiversity
dt.rsp.bdv <- as.data.table(data[,c("SR", "H", "J")])

dt.rsp <- data.table(cbind(dt.rsp.abn, dt.rsp.bms, dt.rsp.bdv))

dt.exp<- as.data.table(data[,c("age_class", "samcam", "field.ID", "location", "area")])
dt.exp <- cbind(dt.exp, std.var)

responses <-  c(colnames(dt.rsp.abn), colnames(dt.rsp.bms), colnames(dt.rsp.bdv)[1:2])
covariates <- c("age_class", "samcam", "scl.mc", "I(scl.mc^2)", "scl.cn", "scl.pH", "scl.clay", "scl.ats1", "I(scl.ats1^2)", "scl.prec1", "age_class:samcam", "scl.mc:scl.pH", "scl.cn:scl.pH")
coefficients <- c("(Intercept)","age_class.L", "age_class.Q", "age_class.C", "age_class^4", 
                  "samcam2", "samcam3",
                  "scl.mc", "I(scl.mc^2)",
                  "scl.cn", "scl.pH", "scl.clay",
                  "scl.ats1", "I(scl.ats1^2)",
                  "scl.prec1",
                  "age_class.L:samcam2", "age_class.Q:samcam2", "age_class.C:samcam2", "age_class^4:samcam2",
                  "age_class.L:samcam3","age_class.Q:samcam3", "age_class.C:samcam3", "age_class^4:samcam3",
                  "scl.mc:scl.pH","scl.cn:scl.pH")
p <- ncol(dt.rsp)-1
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Data summary ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for(i in 1:p) {
  par(mfrow = c(2,2),
      mar = c(3,3,0,1),
      mgp = c(1.5,0.5,0),
      tck = -0.03,
      oma = c(0,0,2,0))
  dt.exp$y <- dt.rsp[,i, with=F]
  plot(dt.exp$y)
  boxplot(dt.exp$y)
  hist(dt.exp$y, main="")
  plot(dt.exp$y^2)
  title(names(dt.rsp)[i],outer=TRUE)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Detecting Outliers ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
row.names(dt.rsp) <- 1:nrow(dt.rsp)
for(i in 1:p) {
  par(mfrow = c(2,2),
      mar = c(3,3,0,1),
      mgp = c(1.5,0.5,0),
      tck = -0.03,
      oma = c(0,0,2,0))
  dt.exp$y <- dt.rsp[,i, with=F]
  car::Boxplot(dt.exp$y ~ dt.exp$age_class)
  car::Boxplot(dt.exp$y ~ dt.exp$samcam)
  car::Boxplot(dt.exp$y ~ interaction(dt.exp$age_class, dt.exp$samcam))
  title(names(dt.rsp)[i],outer=TRUE)
}

# outliers:
# fungivores: 2
# omnivores: 1


# age_class:
# bacterivores decrease
# herbivores increase
# fungivores have an polynomial relationship

# samcam:
# carnivores increase in the second year, only in Silphie!
# no big differences for the others
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# drop Outliers ####
outlier <- list(abn.anc <- -c(43,65),
                abn.ancad <- -c(35,43,60,110),
                abn.endo <- -c(23,175,151),
                abn.endad <- -c(44,23,11,175,132),
                abn.N <- -c(43,44,151),
                abn.anc.juv <- -c(57,58,118,149,137,138,127),
                abn.endo.juv <- -c(10,151),
                abn.juv <- -c(10,151),
                bms.anc <- -c(43,57,110,102,165),
                bms.ancad <- -c(43,57,110,102,165),
                bms.endo <- -c(19,23,142,175),
                bms.endad <- -c(19,23,26,27,92,175),
                bms.N <- -c(57,43,102,165,151),
                bms.anc.juv <- -c(57,58,118,149,138,137,127),
                bms.endo.juv <- -c(34,10,97,175,174,137),
                bms.juv <- -c(57,97,175,174,137,138),
                bdv.SR <- -c(43,44,173,137,138),
                bdv.H.juv <- -c(53,43,44,81,173,175)
                )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Load dredge objects ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#load( file="Analysis/OutputTables/F2_EW_PdredgeAll.rda") # Parrallelized dredge, probably works with optimizer when rsquared is disabled
load( file="F2_EW_pDredgeAll.rda")  # Non parallelized calculation
load( file="F2_EW_DredgeAll1.rda")  # Non parallelized calculation
load(file="Analysis/F2_EW_Mselect.rda")

ls.abn.dredge[[2]] <- Mselect.ancad
compare_pdredge1 <- ls.abn.dredge[[1]]
compare_pdredge10 <- ls.bms.dredge[[1]]
compare_dredge1 <- ls.abn.dredge[[1]]
compare_dredge10 <- ls.bms.dredge[[1]]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 1. Best Models (lowest AICc) #####

# 1.a get the best models ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ls.bestmodels <- list()
ls.dredge <- c(ls.abn.dredge, ls.bms.dredge, ls.bdv.dredge)
p <- ncol(dt.rsp.abn) + ncol(dt.rsp.abn) + ncol(dt.rsp.bdv) -1

for ( i in 1:p) {
#   dt.exp2 <- dt.exp[outlier[[i]],]
#   dt.exp2$y <- dt.rsp[outlier[[i]],i, with=F]
  dt.exp$y <- dt.rsp[,i, with=F]
  M.best <- get.models(ls.dredge[[i]], 1)[[1]]
  name <- paste("BestModel",i,responses[i], sep = ".")
  assign(name, M.best)
  ls.bestmodels[[i]] <- assign(name, M.best)
  names(ls.bestmodels)[[i]] <- name
}

dt.exp$y <- dt.rsp[,"endad", with=F]
endad.best <- glmer(y ~ age_class*samcam + scl.prec1 + scl.mc*scl.pH  + (1|field.ID),dt.exp,family=poisson, control=glmerControl(optimizer="bobyqa"))

ls.bestmodels[[19]] <- endad.best

dt.exp$y <- dt.rsp[,"ancad.bm", with=F]
ancad.bm.glog <- glmer(y+1 ~ age_class + scl.ats1 + I(scl.ats1^2) + (1 | field.ID), family= gaussian(link="log"), data=dt.exp, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

ls.bestmodels[[20]]  <- ancad.bm.glog

# 1.b Check Convergence ####

p = 19
failconv <- matrix(NA,p,2)
for (i in 1:p) {
  relgrad <- with(ls.bestmodels[[i]]@optinfo$derivs,solve(Hessian,gradient))
  print(max(abs(relgrad)))
  failconv[i,2] <- max(abs(relgrad))
}

# if max(abs(relgrad)) is <0.001 then things might be ok... 
# for details see 
# http://stats.stackexchange.com/questions/110004/how-scared-should-we-be-about-convergence-warnings-in-lme4

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 1.c Model Validation #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for(k in 1:20){ 
  # print(list(summary(ls.bestmodels[[k]]),Anova(ls.bestmodels[[k]], type="II")))
  #corvif(ls.bestmodels[[k]])
  
  E1 <- resid(ls.bestmodels[[k]], type="pearson")
  E2 <- resid(ls.bestmodels[[k]], type="response")
  F1 <- fitted(ls.bestmodels[[k]], type="response")
  P1 <- predict(ls.bestmodels[[k]], type="response")
  
  par(mfrow=c(2,2),
      mar=c(4,4.5,1,2),
      oma=c(0,0,2,0)
  )
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
  
  title(names(ls.bestmodels)[k], outer=TRUE)
  
  # Normal QQ Plots
  qqnorm(E2)
  qqline(E2)
  
  
  title(names(ls.bestmodels)[k], outer=TRUE)
  
}

coefplot2(ls.bestmodels[[2]])


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 1.d Best Model Coefficients ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


df.bstCoef1 <- data.frame(matrix(NA,length(coefficients), (p+1)*3))
rownames(df.bstCoef1) <- coefficients
df.bstCoef1[,1] <- coefficients

df.bstCoef2 <- data.frame(row.names = coefficients, id = 1:25) 

for (i in 1:p) {
  #ERROR HANDLING
  possibleError <- tryCatch(
    coeftab(ls.bestmodels[[i]]),
    error=function(e) e
  )
  
  if(inherits(possibleError, "error")){
    df.bstCoef1[,1+i] <- df.bstCoef1[,1+i] # should be akaike weights of the best model
  }
  
  if(!inherits(possibleError, "error")){
    #REAL WORK
    coefMbst <- data.frame(round(coeftab(ls.bestmodels[[i]])[,c(1,3,6)],2))
    df.help <- merge(df.bstCoef2,coefMbst,by="row.names",all.x=TRUE, sort=FALSE)
    df.help <- df.help[order(df.help$id),]
    #df.bstCoef1[,(1+i*3):(3+i*3)] <- df.help[,3:5]
    df.bstCoef1[,(1+i*3)] <- round(df.help[,3],2)
    df.bstCoef1[,(2+i*3)] <- round(df.help[,4],2)
    df.bstCoef1[,(3+i*3)] <- round(df.help[,5],2)
    #useful(i); fun(i); good(i);
  }
} 
colnames <- rep(c("covariate",responses,"endad.old"), each=3)
numeration <- rep(c(1,2,3), length(colnames)/3)
colnames(df.bstCoef1)  <- paste(colnames, numeration, sep="")


df.add <- matrix(NA,5,(p+1)*3)
row.names(df.add) <- c("Rrandom", "Rsquared", "R2m", "R2c", "AICc")
i=19
for(i in 1:p){
  dt.exp$y <- dt.rsp[,i,with=FALSE]
  df.add[1,(3+i*3)] <- round(summary(lm(fitted(ls.bestmodels[[i]]) ~ dt.exp$y))$adj.r.squared,2)
  df.add[2,(3+i*3)] <- round(summary(lm(predict(ls.bestmodels[[i]], type="response") ~ dt.exp$y))$adj.r.squared,2)
  df.add[3:4,(3+i*3)] <- round(r.squaredGLMM(ls.bestmodels[[i]]),2)
  df.add[5,(3+i*3)] <- round(AIC(ls.bestmodels[[i]]),1)
}

measure <- rep(c("estimate", "2.5%", "97.5%"),p+1)
df.bstCoef1 <- rbind(measure[4:length(measure)],df.bstCoef1[,4:ncol(df.bstCoef1)],df.add[,4:length(measure)])
#df.bstCoefA <- df.bstCoef1

#round(summary(ls.bestmodels[[1]])$coefficients[, c(3,4)],4)

# write.csv(df.bstCoef1, file="Analysis/OutputTables/bstCoef.csv")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 1.e get best model predictions with CI ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Generate Test data-set
newdata <- expand.grid(age_class=rep(unique(dt.exp$age_class), each=3),
                       samcam = unique(dt.exp$samcam))
newdata <- cbind(newdata, field.ID = c(1:12,13,14,15,rep(c(1:12,16,17,18),2)))
df.exp.num <- as.data.frame(lapply(lapply(dt.exp[, 5:33, with=FALSE], mean),rep, nrow(newdata)))
newdata <- cbind(newdata, df.exp.num)
newdata$area <- newdata$area*4

# In case of glmmadmb:
# pred <- cbind(newdata, predict(ls.bestmodels[[1]], newdata, se.fit=TRUE)) 

# In case of glmer
ls.pred <- list()
for(i in 1:18) {
  f1 <- paste("~",formula(ls.bestmodels[[i]])[3])
  f1 <- gsub("\\+ \\(1.*\\)", "", f1)
  f1 <- gsub("\\(1.*\\)", "1",f1)
  f1 <- formula(print(f1, quote=F))
  X <- model.matrix( f1, newdata)
  newdata$fit <- X %*% fixef(ls.bestmodels[[i]])
  newdata$SE <- sqrt(  diag(X %*%vcov(ls.bestmodels[[i]]) %*% t(X))  )
  newdata$upr=newdata$fit+1.96*newdata$SE
  newdata$lwr=newdata$fit-1.96*newdata$SE
  pred <- newdata
  ls.pred[[i]] <- pred
}

X <- model.matrix(~  age_class*samcam + scl.prec1 + scl.mc*scl.pH, data = newdata)
newdata$fit <- X %*% fixef(endad.best)
newdata$SE <- sqrt(  diag(X %*%vcov(endad.best) %*% t(X))  )
newdata$upr=newdata$fit+1.96*newdata$SE
newdata$lwr=newdata$fit-1.96*newdata$SE
endad.pred <- newdata

ls.pred[[19]] <- endad.pred

X <- model.matrix(~ age_class + scl.ats1 + I(scl.ats1^2), data = newdata)
newdata$fit <- X %*% fixef(ancad.bm.glog)
newdata$SE <- sqrt(  diag(X %*%vcov(ancad.bm.glog) %*% t(X))  )
newdata$upr=newdata$fit+1.96*newdata$SE
newdata$lwr=newdata$fit-1.96*newdata$SE
ancad.bm.glog.pred <- newdata

ls.pred[[20]] <- ancad.bm.glog.pred


# Rename samcam for facetting
for (i in 1:20) {
  ls.pred[[i]]$samcam2 <- plyr::revalue(pred$samcam,c("1" ="autumn 2012",  "2" ="spring 2013", "3"="autumn 2013"))
  ls.pred[[i]]$age_class <- plyr::revalue(pred$age_class,c("A_Cm"="Cm","B_Sp_young" ="Sp_Y","C_Sp_int1" ="Sp_I1","D_Sp_int2" ="Sp_I2","E_Sp_old" ="Sp_O"))
  # Reverse scaling of covariate
  ls.pred[[i]]$ats1 <- ls.pred[[i]]$scl.ats1* sd(data$ats1) + mean(data$ats1)
  ls.pred[[i]]$prec1 <- ls.pred[[i]]$scl.prec1* sd(data$prec1) + mean(data$prec1)
  ls.pred[[i]]$mc <- ls.pred[[i]]$scl.mc* sd(data$mc) + mean(data$mc)
  ls.pred[[i]]$cn <- ls.pred[[i]]$scl.cn* sd(data$cn) + mean(data$cn)
  ls.pred[[i]]$clay <- ls.pred[[i]]$scl.clay* sd(data$clay) + mean(data$clay)
  ls.pred[[i]]$pH <- ls.pred[[i]]$scl.pH* sd(data$pH) + mean(data$pH)
  ls.pred[[i]]$area <- ls.pred[[i]]$area*4
}

save(ls.pred, file="Analysis/F2_EW_lsPred.rda")

## plot predictions with error bars // confidence intervals???

# !!!! Attention calculating the exp() is not valid for all models (i.e. Species Richness) !!!!
require(ggplot2)
windows(record=T)
for(i in 1:20){
  predfig <- ggplot(ls.pred[[i]], aes(x = age_class, y = exp(fit), ymin = exp(lwr), ymax = exp(upr))) + 
    geom_bar(stat="identity",position = position_dodge(1), col="454545", size=0.15, fill="grey") +
    geom_errorbar(position = position_dodge(1),col="black",width=0.15, size=0.15) + 
    facet_grid(.~samcam2) +
    geom_hline(xintercept = 1, size=0.15) +
    ylab("Anecic Abundance Ind./m²") +
    xlab("Age Class") +
    scale_x_discrete(labels=c("Cm", "Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +
    #scale_y_continuous(limits=c(0,110)) +
    mytheme +
    theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1))
  print(predfig)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 1.f1 Post Hoc data inspection with lsmeans package ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
detach("package:piecewiseSEM", unload=TRUE)
detach("package:lmerTest", unload=TRUE)
detach("package:afex", unload=TRUE)
detach("package:lsmeans", unload=TRUE)
library(lsmeans)
library(multcompView)

ls.lsm<- list()

df.posthoc <- matrix(NA,15,2+(2*p))

for (i in 1:p) {
  #ERROR HANDLING
  possibleError <- tryCatch(
    lsmeans::lsmeans(ls.bestmodels[[i]],  ~ age_class*samcam, contr= "cld"),
    error=function(e) e
  )
  
  if(inherits(possibleError, "error")){
    df.posthoc[,2+((2*i)-1)] <- NA
    df.posthoc[,2+((2*i)-0)] <- NA
    name <- paste("lsm",i,names(dt.rsp)[i], sep = ".")
    ls.lsm[[i]] <- assign(name, lsm)
    # should be akaike weights of the best model
  }
  
  if(!inherits(possibleError, "error")){
  
  # get the results on a back transformed scale:
  lsm <- lsmeans::lsmeans(ls.bestmodels[[i]],  ~ age_class*samcam, contr= "cld")
  x <- cld(lsm, type = "response", adjust = "bon", sort=FALSE)
  df.posthoc[,2+((2*i)-1)] <- x[,3]
  df.posthoc[,2+((2*i)-0)] <- x$".group"
  print(x)
  name <- paste("lsm",i,names(dt.rsp)[i], sep = ".")
  ls.lsm[[i]] <- assign(name, lsm)
  # to see the results graphically
  p1 <- plot(lsm, by = "samcam", intervals = TRUE, type = "response")
  print(p1)
  #title(names(ncr.biglmer)[i], outer=TRUE)
}
}
df.posthoc[,1] <- paste(x$"age_class")
df.posthoc[,2] <- paste(x$"samcam")

colnames(df.posthoc) <- c("Factor1", "Factor1", rep(c(colnames(dt.rsp)[1:18], "endad.old"),each=2))
df.posthoc <- rbind(c("age_class", "samcam", rep(c("lsmean", "group"), p)), df.posthoc)


for(i in 1:p){
  #ERROR HANDLING
  possibleError <- tryCatch(
    lsmip(ls.lsm[[i]], age_class ~ samcam, type = "response"),
    error=function(e) e
  )
  
#   if(inherits(possibleError, "error")){
#     df.posthoc[,2+((2*i)-1)] <- NA
#     df.posthoc[,2+((2*i)-0)] <- NA
#     name <- paste("lsm",i,names(dt.rsp)[i], sep = ".")
#     ls.lsm[[i]] <- assign(name, lsm)
#     # should be akaike weights of the best model
#   }
  
  if(!inherits(possibleError, "error")){
  lsmip(ls.lsm[[i]], age_class ~ samcam, type = "response")

}
}

#************************************************************************


ls.lsmAC <- list()
df.posthocAC <- matrix(NA,30,2+(2*p))
colnames(df.posthocAC) <- c("contrast", "samcam", rep(c(colnames(dt.rsp)[1:18], "endad.old"),each=2))

for (i in 1:p) {
  #ERROR HANDLING
  possibleError <- tryCatch(
    lsmeans::lsmeans(ls.bestmodels[[i]],  ~ age_class|samcam, at=list(samcam=c("1","2","3"))),
    error=function(e) e
  )
  
  if(inherits(possibleError, "error")){
    df.posthocAC[,2+((2*i)-1)] <- NA
    df.posthocAC[,2+((2*i)-0)] <- NA
    name <- paste("lsm",i,names(dt.rsp)[i], sep = ".")
    ls.lsm[[i]] <- assign(name, lsm)
    # should be akaike weights of the best model
  }
  
  if(!inherits(possibleError, "error")){
  # get the results on a back transformed scale:
  lsm <- lsmeans(ls.bestmodels[[i]],  ~ age_class|samcam, at=list(samcam=c("1","2","3")))
  x <- contrast(lsm, "pairwise" , type = "response")
  print(x)
  xx <- summary(x)
  df.posthocAC[,2+((2*i)-1)] <- round(xx$"estimate",2)
  df.posthocAC[,2+((2*i)-0)] <- round(xx$"p.value",3)
  name <- paste("lsm",i,names(dt.rsp)[i], sep = ".")
  ls.lsmAC[[i]] <- assign(name, lsm)
}
}
df.posthocAC[,1] <- paste(xx$"contrast")
df.posthocAC[,2] <- paste(xx$"samcam")


#************************************************************************


ls.lsmSC <- list()
df.posthocSC <- matrix(NA,15,2+(2*p))
colnames(df.posthocSC) <- c("contrast", "age_class", rep(c(colnames(dt.rsp)[1:18], "endad.old"),each=2))

for (i in 1:p) {
  #ERROR HANDLING
  possibleError <- tryCatch(
    lsmeans(ls.bestmodels[[i]], pairwise ~ samcam|age_class),
    error=function(e) e
  )
  
  if(inherits(possibleError, "error")){
    df.posthocAC[,2+((2*i)-1)] <- NA
    df.posthocAC[,2+((2*i)-0)] <- NA
    name <- paste("lsm",i,names(dt.rsp)[i], sep = ".")
    ls.lsm[[i]] <- assign(name, lsm)
    # should be akaike weights of the best model
  }
  if(!inherits(possibleError, "error")){
  # get the results on a back transformed scale:
  lsm <- lsmeans(ls.bestmodels[[i]], pairwise ~ samcam|age_class)
  x <- contrast(lsm, "pairwise" , type = "response")
  print(x)
  xx <- summary(x)
  df.posthocSC[,2+((2*i)-1)] <- round(xx$"estimate",2)
  df.posthocSC[,2+((2*i)-0)] <- round(xx$"p.value",3)
  name <- paste("lsm",i,names(dt.rsp)[i], sep = ".")
  ls.lsmSC[[i]] <- assign(name, lsm)
}
}

df.posthocSC[,1] <- paste(xx$"contrast")
df.posthocSC[,2] <- paste(xx$"age_class")

# write.csv(df.posthoc, file="Analysis/OutputTables/df.posthoc.csv")
# write.csv(df.posthocAC, file="Analysis/OutputTables/df.posthocAC.csv")
# write.csv(df.posthocSC, file="Analysis/OutputTables/df.posthocSC.csv")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Post Hoc Multicomparisons with function glht
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 1.f2 Post-Hoc Multicomparisons with glht() ####
# interaction effect - reduced contrast matrix 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# # Model with the interaction term
# dt.exp$interaction <- factor(interaction(dt.exp$age_class, dt.exp$samcam))
# dt.exp$y <- dt.rsp$ancad
# ancad.tukey <- glmer(y ~ interaction + I(scl.ats1^2) + scl.prec1 + (1|field.ID) + offset(log(area)) ,data=dt.exp,family=poisson, control=glmerControl(optimizer="bobyqa"))
# # summary(anc.tukey)
# 
# 
# # Create contrast matrix
# source("Analysis/F2_EW_ContrastMatrix.R")
# data2 <- data
# data <- dt.exp
# source("Analysis/F2_EW_ContrastMatrix2.R")
# data <- data2
# 
# # Pairwise comparisons (with interaction term)
# anc.pairwise <- glht(ancad.tukey, linfct=mcp(interaction = cm1))
# anc.pw.ci <- confint(anc.pairwise)
# summary(anc.pairwise, test=adjusted(type="fdr"))
# 
# # Confidence intervals including 0
# anc.pw.sig <- which(anc.pw.ci$confint[,2]>0)
# data.frame(names(anc.pw.sig))
# 
# # Plot Errorbars
# phfig1 <- ggplot(anc.pw.ci, aes(y = lhs, x = exp(estimate), xmin = exp(lwr), xmax = exp(upr))) + 
#   geom_errorbarh() + 
#   geom_point() + 
#   geom_vline(xintercept = 1) +
#   mytheme
# phfig1
# 
# # plot confidence intervals
# par(mar=c(2,15,2,2))
# plot(anc.pw.ci) 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# 2. Multimodel Averaging - Component Models: AICc delta < 4 #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 2.a Extract component Models #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

df.compM2 <- vector()

for(i in 1:18){
  delta4 <- subset(ls.dredge[[i]], delta < 4)
  df.compM <- data.frame(delta4)
  df.compM[,ncol(df.compM)+1] <- rep(colnames(dt.rsp)[i], length(rownames(df.compM)))
  df.compM2 <- rbind(df.compM2, df.compM)
}

# write.csv(df.compM2, file="Analysis/OutputTables/ComponentModels.csv")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 2.b Get Variable importance #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

df.relImportance1 <- matrix(NA,length(covariates), 19)
rownames(df.relImportance1) <- covariates
df.relImportance1[,1] <- covariates
colnames(df.relImportance1) <- c("covariates",responses)
ncol(df.relImportance1) 

df.relImportance2 <- data.frame(row.names = covariates, id = 1:length(covariates))


for (i in 1:18) {
  #ERROR HANDLING
  possibleError <- tryCatch(
    model.avg(ls.dredge[[i]], subset=delta < 4),
    error=function(e) e
  )
  
  if(inherits(possibleError, "error")){
    df.relImportance1[,1+i] <- df.relImportance1[,1+i] # should be akaike weights of the best model
  }
  
  if(!inherits(possibleError, "error")){
    #REAL WORK
    avgMd4 <- model.avg(ls.dredge[[i]], subset=delta < 4) 
    avgMd4Imp <- data.frame(importance(avgMd4))
    df.help <- merge(df.relImportance2,avgMd4Imp,by="row.names",all.x=TRUE, sort=FALSE)
    df.help <- df.help[order(df.help$id),]
    df.relImportance1[,1+i] <- df.help[,3]
    #useful(i); fun(i); good(i);
  }
  
}  #end for

# write.csv(df.relImportance1, file="Analysis/OutputTables/Importance.csv")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 2.c Get averaged coefficients #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

df.avCoef1 <- matrix(NA,length(coefficients), 19)
rownames(df.avCoef1) <- coefficients
df.avCoef1[,1] <- coefficients
colnames(df.avCoef1) <- c("covariates",responses)

df.avCoef2 <- data.frame(row.names = coefficients, id = 1:25) 

for (i in 1:18) {
  #ERROR HANDLING
  possibleError <- tryCatch(
    model.avg(ls.dredge[[i]], subset=delta < 4),
    error=function(e) e
  )
  
  if(inherits(possibleError, "error")){
    df.avCoef1[,1+i] <- df.avCoef1[,1+i] # should be akaike weights of the best model
  }
  
  if(!inherits(possibleError, "error")){
    #REAL WORK
    avgMd4 <- model.avg(ls.dredge[[i]], subset=delta < 4) 
    avgMd4coef <- data.frame(coef(avgMd4,1))
    data.frame(coef(avgMd4,1))
    df.help <- merge(df.avCoef2,avgMd4coef,by="row.names",all.x=TRUE, sort=FALSE)
    df.help <- df.help[order(df.help$id),]
    df.avCoef1[,1+i] <- df.help[,3]
    #useful(i); fun(i); good(i);
  }
  
} 

# write.csv(df.avCoef1, file="Analysis/OutputTables/AVerageCoefficients.csv")


# Weitere Befehle: 
# print(avgMd4)
# modelList(avgMd4)
# confint(avgMd4)
# confset95p <- get.models(ls.abn.dredge[[1]], cumsum(ls.abn.dredge[[1]]$weight) <= .95)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 2.d Get Average Model predictions #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Predictions from each of the models in a set, and with averaged coefficients

ls.delta4 <- list()

for (i in 1:18) {
  dt.exp$y <- dt.rsp[,i, with=F]
  delta4.1 <- get.models(ls.dredge[[i]], delta < 4)
  ls.delta4[[i]] <- delta4.1
}

for(i in 1:18) {
  model.preds = sapply(ls.delta4[[i]], predict, type="response", newdata)
  delta4.2 <- subset(ls.dredge[[i]], delta < 4)
  mod.ave.preds<- model.preds %*% Weights(delta4.2)
  head(mod.ave.preds)
  p1 <- boxplot(mod.ave.preds ~ newdata$age_class)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 3. save all tables ####
# write.csv(df.bstCoef1, file="Analysis/OutputTables/bstCoef.csv")
# save(ls.pred, file="Analysis/F2_EW_lsPred.rda")
# write.csv(df.posthoc, file="Analysis/OutputTables/df.posthoc.csv")
# write.csv(df.posthocAC, file="Analysis/OutputTables/df.posthocAC.csv")
# write.csv(df.posthocSC, file="Analysis/OutputTables/df.posthocSC.csv")
# write.csv(df.compM2, file="Analysis/OutputTables/ComponentModels.csv")
# write.csv(df.relImportance1, file="Analysis/OutputTables/Importance.csv")
# write.csv(df.avCoef1, file="Analysis/OutputTables/AVerageCoefficients.csv")
