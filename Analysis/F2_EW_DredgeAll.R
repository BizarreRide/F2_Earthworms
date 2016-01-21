#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# F2_Eearthworms
# GLMM with Dredge
# Quentin Schorpp
# 29.12.2015
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setwd("C:/Users/Quentin/Documents/git_repositories/F2_Earthworms")

# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Data/GatherSource/F2_EW_MakeLikeFile.R")

# Additional Function:
is.correlated <- function(i, j, data, conf.level = .95, cutoff = .5, ...) {
        if(j >= i) return(NA)
        ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
        ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}
# Suppose we want to have set of models that exclude combinations of colinear
# variables, that are significantly (p < 0.05) correlated, with Pearson
# correlation coefficient larger than r = 0.5.

vCorrelated <- Vectorize(is.correlated, c("i", "j")) # Need vectorized function to use with 'outer'
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


dt.exp<- as.data.table(data[,c("age_class", "samcam", "field.ID", "location", "area")])
dt.exp <- cbind(dt.exp, std.var)
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



# Global Models ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Global Models Abundance ####
p <- ncol(dt.rsp.abn)
ls.abn.Mglobal <- list()

for(i in 1:p){
#         dt.exp2 <- dt.exp[outlier.abn[[i]],]
#         dt.exp2$y <- dt.rsp[outlier.abn[[i]],i, with=F]
        dt.exp$y <- dt.rsp.abn[,i, with=F]
        Mglobal <- glmer(y ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.ats1 + I(scl.ats1^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID),data=dt.exp,family=poisson, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
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
#         dt.exp2 <- dt.exp[outlier.bms[[i]],]
#         dt.exp2$y <- dt.rsp[outlier.bms[[i]],i, with=F]
        dt.exp$y <- dt.rsp.bms[,i, with=F]
        Mglobal <- lmer(log1p(y) ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.ats1 + I(scl.ats1^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID),data=dt.exp, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
        #Mglobal <- glmer(log1p(y) ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.ats1 + I(scl.ats1^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=dt.exp, family=gaussian(link="log")) #, control=glmerControl(optimizer="bobyqa"))
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
#         dt.exp2 <- dt.exp[outlier.bdv[[i]],]
#         dt.exp2$y <- dt.rsp[outlier.bdv[[i]],i, with=F]
        dt.exp$y <- dt.rsp.bdv[,i, with=F]
        Mglobal <- lmer(y ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.ats1 + I(scl.ats1^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID),data=dt.exp, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
        print(summary(Mglobal))
        print(overdisp_fun(Mglobal))
        name <- paste("Ew_Model",i,names(dt.rsp.bdv)[i], sep = ".")
        assign(name, Mglobal)
        ls.bdv.Mglobal[[i]] <- assign(name, Mglobal)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Dredge Models # demo(dredge.subset) ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
#         dt.exp2 <- dt.exp[outlier.abn[[i]],]
#         dt.exp2$y <- dt.rsp[outlier.abn[[i]],i, with=F]
        dt.exp$y <- dt.rsp.abn[,i, with=F]
        GM.dredge <- dredge(ls.abn.Mglobal[[i]],  subset = smat) #, m.lim=c(0,5), extra=c("R^2", "adjR^2"))
        #fixed=c("age_class", "samcam"),
        name <- paste("Ew_abn.dredge",i,names(dt.rsp.abn)[i], sep = ".")
        assign(name, GM.dredge)
        ls.abn.dredge[[i]] <- assign(name, GM.dredge)
}

# Dredge Biomass ####
p <- ncol(dt.rsp.bms)
ls.bms.dredge <- list()
for(i in 1:p) {
#         dt.exp2 <- dt.exp[outlier.bms[[i]],]
#         dt.exp2$y <- dt.rsp[outlier.bms[[i]],i, with=F]
        dt.exp$y <- dt.rsp.bms[,i, with=F]
        GM.dredge <- dredge(ls.bms.Mglobal[[i]],  subset=smat) #, m.lim=c(0,5), extra=c("R^2", "adjR^2"))
        #fixed=c("age_class", "samcam"),
        name <- paste("Ew_bms.dredge",i,names(dt.rsp.bms)[i], sep = ".")
        assign(name, GM.dredge)
        ls.bms.dredge[[i]] <- assign(name, GM.dredge)
}

# Dredge Biodiversity ####
p <- ncol(dt.rsp.bdv)-1
ls.bdv.dredge <- list()
for(i in 1:p) {
#         dt.exp2 <- dt.exp[outlier.bdv[[i]],]
#         dt.exp2$y <- dt.rsp[outlier.bdv[[i]],i, with=F]
        dt.exp$y <- dt.rsp.bdv[,i, with=F]
        GM.dredge <- dredge(ls.bdv.Mglobal[[i]],  subset=smat) #, m.lim=c(0,5)), extra=c("R^2", "adjR^2"))
        #fixed=c("age_class", "samcam"),
        name <- paste("Ew_bv.dredge",i,names(dt.rsp.bdv)[i], sep = ".")
        assign(name, GM.dredge)
        ls.bdv.dredge[[i]] <- assign(name, GM.dredge)
}



save(list=c("ls.abn.dredge", "ls.abn.Mglobal", "ls.bms.dredge", "ls.bms.Mglobal", "ls.bdv.dredge", "ls.bdv.Mglobal"), 
     file="C:/Users/Quentin/Documents/git_repositories/F2_Earthworms/Analysis/OutputTables/F2_EW_DredgeAll1.rda")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# QUESTIONS: 

# A.) Limit Number of estimated parameters to 5 (6,7,...) ?
#     Including or excluding m.lim=c(0.5), 
#     RUle of Thumb: 1 explanator for at least 10 (15) samples
#     soft limits related to statistical power and good statistical practice.
# http://stats.stackexchange.com/questions/12854/maximum-number-of-independent-variables-that-can-be-entered-into-a-multiple-regr
# http://stats.stackexchange.com/questions/10079/rules-of-thumb-for-minimum-sample-size-for-multiple-regression

# B.) Include age_class*samcam in all models, since it is the model we want to test?

# c.) Power Analysis (Software: G.Power 3, Power Primer)


# Special case for ancad and ancad.bm ####

dt.exp$y <- dt.rsp[,"ancad", with=F]

M.1 <- glmer(y ~ age_class + samcam + scl.ats1 + I(scl.ats1^2) + (1|field.ID),data=dt.exp,family=poisson, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
M.2 <- glmer(y ~ age_class + samcam + I(scl.ats1^2) + scl.prec1 + (1|field.ID),data=dt.exp,family=poisson, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
M.3 <- glmer(y ~ age_class + samcam + scl.ats1 + I(scl.ats1^2) + scl.prec1 + (1|field.ID),data=dt.exp,family=poisson, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
M.4 <- glmer(y ~ age_class + samcam + scl.ats1 + scl.cn + (1|field.ID),data=dt.exp,family=poisson, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
M.5 <- glmer(y ~ age_class + samcam + scl.ats1 + scl.pH + (1|field.ID),data=dt.exp,family=poisson, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
M.6 <- glmer(y ~ age_class + samcam + scl.ats1 + scl.cn + scl.pH + (1|field.ID),data=dt.exp,family=poisson, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
M.7 <- glmer(y ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.ats1 + I(scl.ats1^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID),data=dt.exp,family=poisson, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
M.8 <- glmer(y ~ 1 + (1|field.ID),data=dt.exp,family=poisson, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

Mselect.ancad <- model.sel(M.1,M.2,M.3,M.4,M.5, M.6, M.7, M.8)


dt.exp$y <- dt.rsp[,"ancad.bm", with=F]

M.1 <- lmer(log1p(y) ~ age_class + scl.ats1 + I(scl.ats1^2) + (1|field.ID),data=dt.exp, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
M.2 <- lmer(log1p(y) ~ age_class + I(scl.ats1^2) + (1|field.ID),data=dt.exp, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
M.3 <- lmer(log1p(y) ~ age_class + scl.ats1 + I(scl.ats1^2) + scl.prec1  + (1|field.ID),data=dt.exp, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
M.4 <- lmer(log1p(y) ~ age_class + scl.pH + (1|field.ID),data=dt.exp, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
M.5 <- lmer(log1p(y) ~ age_class + scl.ats1 + scl.pH + scl.prec1  + (1|field.ID),data=dt.exp, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
M.6 <- lmer(log1p(y) ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.ats1 + I(scl.ats1^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID),data=dt.exp, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
M.7 <- lmer(log1p(y) ~ 1 + (1|field.ID),data=dt.exp, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

Mselect.ancad.bm <- model.sel(M.1,M.2,M.3,M.4,M.5, M.6, M.7)

save(list=c("Mselect.ancad", "Mselect.ancad.bm"), 
     file="C:/Users/Quentin/Documents/git_repositories/F2_Earthworms/Analysis/OutputTables/F2_EW_Mselect.rda")



