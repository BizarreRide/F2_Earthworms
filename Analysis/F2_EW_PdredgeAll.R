#%%%%%%%%%%%%%%%%%%%%%%%%%%%
# F2_Eearthworms
# GLMM with Dredge
# Quentin Schorpp
# 29.12.2015
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Data/GatherSource/F2_EW_MakeLikeFile.R")
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

p <- ncol(dt.rsp)-1
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
outlier.abn <- list(abn.anc <- -c(43,65),
                    abn.ancad <- -c(35,43,60,110),
                    abn.endo <- -c(23,175,151),
                    abn.endad <- -c(44,23,11,175,132),
                    abn.N <- -c(43,44,151),
                    abn.anc.juv <- -c(57,58,118,149,137,138,127),
                    abn.endo.juv <- -c(10,151),
                    abn.juv <- -c(10,151)
)

outlier.bms <-  list(bms.anc <- -c(43,57,110,102,165),
                     bms.ancad <- -c(43,57,110,102,165),
                     bms.endo <- -c(19,23,142,175),
                     bms.endad <- -c(19,23,26,27,92,175),
                     bms.N <- -c(57,43,102,165,151),
                     bms.anc.juv <- -c(57,58,118,149,138,137,127),
                     bms.endo.juv <- -c(34,10,97,175,174,137),
                     bms.juv <- -c(57,97,175,174,137,138)
)

outlier.bdv <-  list(bdv.SR <- -c(43,44,173,137,138),
                     bdv.H.juv <- -c(53,43,44,81,173,175)
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Global Models Abundance ####

 frm.global <- as.formula(y ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.ats1 + I(scl.ats1^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)))
 frm.global.log <- as.formula(log1p(y) ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.ats1 + I(scl.ats1^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)))
#frm.global <- as.formula(y ~ age_class*date + scl.mc + scl.pH + scl.cn + scl.clay + scl.ats1 + I(scl.ats1^2) + scl.prec1 + (1|field.ID) + offset(log(area)))
#frm.global.log <- as.formula(log1p(y) ~ age_class*date + scl.mc + scl.pH + scl.cn + scl.clay + scl.ats1 + I(scl.ats1^2) + scl.prec1 + (1|field.ID) + offset(log(area)))

dt.exp$date <- data$dsamcam
dt.exp2 <- dt.exp


#library(optimx)
#con = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
#con = glmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE))
con = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = FALSE, check.conv.grad="ignore")

p <- ncol(dt.rsp.abn)
ls.abn.Mglobal <- list()

for(i in 1:p){
        #dt.exp <- dt.exp2[outlier.abn[[i]],]
        #dt.exp$y <- dt.rsp.abn[outlier.abn[[i]],i, with=F]
        dt.exp$y <- dt.rsp.abn[,i, with=F]
        Mglobal <- glmer(frm.global,data=dt.exp,family=poisson, control=con)
        print(summary(Mglobal))
        print(overdisp_fun(Mglobal))
        name <- paste("Ew_Model",i,names(dt.rsp.abn)[i], sep = ".")
        assign(name, Mglobal)
        ls.abn.Mglobal[[i]] <- assign(name, Mglobal)
}

# Global Models Biomass ####
p <- ncol(dt.rsp.bms)
ls.bms.Mglobal <- list()
for(i in 1:p){
        #dt.exp <- dt.exp2[outlier.bms[[i]],]
        #dt.exp$y <- dt.rsp.bms[outlier.bms[[i]],i, with=F]
        dt.exp$y <- dt.rsp.bms[,i, with=F]
        Mglobal <- lmer(frm.global.log,data=dt.exp) #, control=glmerControl(optimizer="bobyqa"))
        #Mglobal <- glmer(frm.global.log,data=dt.exp, family=gaussian(link="log"), control=con)
        print(summary(Mglobal))
        print(overdisp_fun(Mglobal))
        name <- paste("Ew_Model",i,names(dt.rsp.bms)[i], sep = ".")
        assign(name, Mglobal)
        ls.bms.Mglobal[[i]] <- assign(name, Mglobal)
}

# Global Models Biodiversity ####
p <- ncol(dt.rsp.bdv)-1
ls.bdv.Mglobal <- list()

for(i in 1:p){
        #dt.exp <- dt.exp2[outlier.bdv[[i]],]
        #dt.exp$y <- dt.rsp.bdv[outlier.bdv[[i]],i, with=F]
        dt.exp$y <- dt.rsp.bdv[,i, with=F]
        Mglobal <- lmer(frm.global,data=dt.exp) #, control=glmerControl(optimizer="bobyqa"))
        print(summary(Mglobal))
        print(overdisp_fun(Mglobal))
        name <- paste("Ew_Model",i,names(dt.rsp.bdv)[i], sep = ".")
        assign(name, Mglobal)
        ls.bdv.Mglobal[[i]] <- assign(name, Mglobal)
}

# 1.c Model Validation #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ls.globalModels <- c(ls.abn.Mglobal, ls.bms.Mglobal,ls.bdv.Mglobal)

for(k in 1:18){ 
  # print(list(summary(ls.globalModels[[k]]),Anova(ls.globalModels[[k]], type="II")))
  #corvif(ls.globalModels[[k]])
  
  E1 <- resid(ls.globalModels[[k]], type="pearson")
  E2 <- resid(ls.globalModels[[k]], type="response")
  F1 <- fitted(ls.globalModels[[k]], type="response")
  P1 <- predict(ls.globalModels[[k]], type="response")
  
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
  
  title(names(ls.globalModels)[k], outer=TRUE)
  
  # Normal QQ Plots
  qqnorm(E2)
  qqline(E2)
  
  
  title(names(ls.globalModels)[k], outer=TRUE)
  
}

coefplot2(ls.globalModels[[2]])


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Dredge Models # demo(dredge.subset) ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Generate Cluster ####
library(parallel)
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 4), type = clusterType))
clusterEvalQ(clust, c(library("lme4")))
#clusterExport(clust, varlist=c("dt.exp", "con"))

# Subsets of models excluded from dredge: ####
opo <- dt.exp[,c("scl.ats1","scl.mc","scl.pH","scl.cn", "scl.prec1", "scl.clay"), with=FALSE]
opo <- as.data.frame(opo)
opo$"I(scl.ats1^2)" <- dt.exp$scl.ats1^2
opo$"I(scl.mc^2)" <- dt.exp$scl.mc^2

smat <- outer(1:8,1:8, vCorrelated, data=opo)
nm <- colnames(opo)
dimnames(smat) <- list(nm, nm)
smat[8,1] <- FALSE
smat[7,2] <- FALSE
smat[8,7] <- FALSE

### A simpler case: exclude only pairs of variables having cor. coefficient
### r > 0.5
# smat <- abs(cor(opo)) <= .5
# smat[!lower.tri(smat)] <- NA


i <- as.vector(smat == FALSE & !is.na(smat))

sexpr <-parse(text = paste("!(", paste("(",
                                       nm[col(smat)[i]], " && ",
                                       nm[row(smat)[i]], ")",
                                       sep = "", collapse = " || "), ")"))
# sexpr <- parse(text = paste(sexpr, "&& dc(scl.ats1, I(scl.ats1^2)) && dc(scl.mc, I(scl.mc^2))"))

options(na.action = na.fail)

# Dredge Abundance ####
p <- ncol(dt.rsp.abn)
ls.abn.dredge <- list()
for(i in 1:p) {
        #dt.exp <- dt.exp2[outlier.abn[[i]],]
        #dt.exp$y <- dt.rsp.abn[outlier.abn[[i]],i, with=F]
        dt.exp$y <- dt.rsp.abn[,i, with=F]
        clusterExport(clust, varlist=c("dt.exp", "con"))
        GM.dredge <- pdredge(ls.abn.Mglobal[[i]],  subset=smat, cluster=clust)
        #fixed=c("age_class", "samcam"),
        name <- paste("Ew_abn.dredge",i,names(dt.rsp.abn)[i], sep = ".")
        assign(name, GM.dredge)
        ls.abn.dredge[[i]] <- assign(name, GM.dredge)
}

# Dredge Biomass ####
p <- ncol(dt.rsp.bms)
ls.bms.dredge <- list()
for(i in 1:p) {
        #dt.exp <- dt.exp2[outlier.bms[[i]],]
        #dt.exp$y <- dt.rsp.bms[outlier.bms[[i]],i, with=F]
        dt.exp$y <- dt.rsp.bms[,i, with=F]
        clusterExport(clust, varlist=c("dt.exp", "con"))
        GM.dredge <- pdredge(ls.bms.Mglobal[[i]],  subset=smat, cluster=clust)
        #fixed=c("age_class", "samcam"),
        name <- paste("Ew_bms.dredge",i,names(dt.rsp.bms)[i], sep = ".")
        assign(name, GM.dredge)
        ls.bms.dredge[[i]] <- assign(name, GM.dredge)
}

# Dredge Biodiversity ####
p <- ncol(dt.rsp.bdv)-1
ls.bdv.dredge <- list()
for(i in 1:p) {
        dt.exp <- dt.exp2[outlier.bdv[[i]],]
        dt.exp$y <- dt.rsp.bdv[outlier.bdv[[i]],i, with=F]
        #dt.exp$y <- dt.rsp.bdv[,i, with=F]
        clusterExport(clust, varlist=c("dt.exp", "con"))
        GM.dredge <- pdredge(ls.bdv.Mglobal[[i]],  subset=smat, cluster=clust)
        #fixed=c("age_class", "samcam"),
        name <- paste("Ew_bv.dredge",i,names(dt.rsp.bdv)[i], sep = ".")
        assign(name, GM.dredge)
        ls.bdv.dredge[[i]] <- assign(name, GM.dredge)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Save #####
 save(list=c("ls.abn.dredge", "ls.abn.Mglobal", "ls.bms.dredge", "ls.bms.Mglobal", "ls.bdv.dredge", "ls.bdv.Mglobal"), 
      file="Analysis/DredgeRdata/F2_EW_PdredgeAll_NewCLimate.rda")


# candidate models, careful thinking:

# The only time models
formula(y ~ age_class*samcam)
formula(y ~ age_class + samcam)
formula(y ~ age_class)
formula(y ~ samcam)

# The only soil models
formula(y ~ scl.pH + scl.mc + scl.cn + scl.clay)
formula(y ~ scl.pH + scl.mc)
formula(y ~ scl.pH*scl.mc)
formula(y ~ I(scl.pH^2) + I(scl.mc^2) + I(scl.cn^2) + scl.clay) # unimodal relationship
formula(y ~ scl.pH + scl.mc + scl.cn + I(scl.pH^2) + I(scl.mc^2) + I(scl.cn^2) + scl.clay)

# The weather models
formula(y ~ scl.ats1 + scl.prec1 + scl.rad1 + scl.hum1)
formula(y ~ I(scl.ats1^2) + I(scl.prec1^2) + I(scl.rad1^2) + I(scl.hum1^2))
formula(y ~ I(scl.ats1 + scl.prec1 + scl.rad1 + scl.hum1 + scl.ats1^2) + I(scl.prec1^2) + I(scl.rad1^2) + I(scl.hum1^2))

# The global model

# The NULL model

# derivative models
# ats1 and mc are correlated!

# Candidate models for ancad

formula(y ~ age_class + samcam + scl.ats1 + I(scl.ats^2) + (1|field.ID) + offset(log(area)))
formula(y ~ age_class + samcam + scl.prec1 + I(scl.ats^2) + (1|field.ID) + offset(log(area)))
formula(y ~ age_class + samcam + scl.ats1 + I(scl.ats^2) + scl.prec1 + (1|field.ID) + offset(log(area)))
formula(y ~ age_class + samcam + scl.ats1 + scl.cn + (1|field.ID) + offset(log(area)))
formula(y ~ age_class + samcam + scl.ats1 + scl.pH + (1|field.ID) + offset(log(area)))
formula(y ~ age_class + samcam + scl.ats1 + scl.cn + scl.pH + (1|field.ID) + offset(log(area)))



f1 <- formula(y ~ age_class + samcam + I(scl.ats1^2) + scl.prec1 + (1|field.ID) + offset(log(area)))
f2 <- formula(y ~ age_class + samcam + I(scl.ats1^2) + scl.clay + scl.prec1 + (1|field.ID) + offset(log(area)))
f3 <- formula(y ~ age_class + samcam + I(scl.ats1^2) + scl.cn + scl.prec1 + (1|field.ID) + offset(log(area)))
f4 <- formula(y ~ age_class + samcam + I(scl.ats1^2) + scl.pH + scl.prec1 + (1|field.ID) + offset(log(area)))
f5 <- formula(y ~ age_class + samcam + I(scl.ats1^2) + scl.ats1 + scl.prec1 + (1|field.ID) + offset(log(area)))
f6 <- formula(y ~ age_class + samcam + I(scl.ats1^2) + scl.clay + scl.pH + scl.prec1 + (1|field.ID) + offset(log(area)))
f7 <- formula(y ~ age_class + samcam + I(scl.ats1^2) + scl.clay + scl.cn + scl.prec1 + (1|field.ID) + offset(log(area)))
f8 <- formula(y ~ age_class + samcam + I(scl.ats1^2) + scl.clay + scl.ats1 + scl.prec1 + (1|field.ID) + offset(log(area)))
f9 <- formula(y ~ age_class + samcam + I(scl.ats1^2) + scl.clay + (1|field.ID) + offset(log(area)))
f10 <- formula(y ~ age_class + samcam + I(scl.ats1^2) + scl.cn + scl.pH + scl.prec1 + (1|field.ID) + offset(log(area)))


# Candidate models for ancad.bm



