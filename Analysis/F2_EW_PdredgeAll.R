#???????????????????????????
# F2_Eearthworms
# GLMM with Dredge
# Quentin Schorpp
# 29.12.2015
#???????????????????????????

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

# Global Models Abundance ####
p <- ncol(dt.rsp.abn)
ls.abn.Mglobal <- list()

for(i in 1:p){
        dt.exp$y <- dt.rsp.abn[,i, with=F]
        Mglobal <- glmer(y ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.ats1 + I(scl.ats1^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=dt.exp,family=poisson) #, control=glmerControl(optimizer="bobyqa"))
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
        dt.exp$y <- dt.rsp.bms[,i, with=F]
        Mglobal <- lmer(log1p(y) ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.ats1 + I(scl.ats1^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=dt.exp) #, control=glmerControl(optimizer="bobyqa"))
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
        dt.exp$y <- dt.rsp.bdv[,i, with=F]
        Mglobal <- lmer(y ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.ats1 + I(scl.ats1^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=dt.exp) #, control=glmerControl(optimizer="bobyqa"))
        print(summary(Mglobal))
        print(overdisp_fun(Mglobal))
        name <- paste("Ew_Model",i,names(dt.rsp.bdv)[i], sep = ".")
        assign(name, Mglobal)
        ls.bdv.Mglobal[[i]] <- assign(name, Mglobal)
}


# Dredge Models # demo(dredge.subset)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Generate Cluster
library(parallel)
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 4), type = clusterType))
clusterEvalQ(clust, c(library("lme4")))
clusterExport(clust, "dt.exp")


# Subsets of models excluded from dredge:
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


# i <- as.vector(smat == FALSE & !is.na(smat))
# 
# sexpr <-parse(text = paste("!(", paste("(",
#                                        nm[col(smat)[i]], " && ",
#                                        nm[row(smat)[i]], ")",
#                                        sep = "", collapse = " || "), ")"))
# sexpr <- parse(text = paste(sexpr, "&& dc(scl.ats1, I(scl.ats1^2)) && dc(scl.mc, I(scl.mc^2))"))



options(na.action = na.fail)

# Dredge Abundance
p <- ncol(dt.rsp.abn)
ls.abn.dredge <- list()
for(i in 1:p) {
        dt.exp$y <- dt.rsp.abn[,i, with=F]
        GM.dredge <- pdredge(ls.abn.Mglobal[[i]],  subset=smat, cluster=clust)
        #fixed=c("age_class", "samcam"),
        name <- paste("Ew_abn.dredge",i,names(dt.rsp.abn)[i], sep = ".")
        assign(name, GM.dredge)
        ls.abn.dredge[[i]] <- assign(name, GM.dredge)
}


# Dredge Biomass
p <- ncol(dt.rsp.bms)
ls.bms.dredge <- list()
for(i in 1:p) {
        dt.exp$y <- dt.rsp.bms[,i, with=F]
        GM.dredge <- pdredge(ls.bms.Mglobal[[i]],  subset=smat, cluster=clust)
        #fixed=c("age_class", "samcam"),
        name <- paste("Ew_bms.dredge",i,names(dt.rsp.abn)[i], sep = ".")
        assign(name, GM.dredge)
        ls.bms.dredge[[i]] <- assign(name, GM.dredge)
}


# Dredge Biodiversity
p <- ncol(dt.rsp.bdv)-1
ls.bdv.dredge <- list()
for(i in 1:p) {
        dt.exp$y <- dt.rsp.bdv[,i, with=F]
        GM.dredge <- pdredge(ls.bdv.Mglobal[[i]],  subset=smat, cluster=clust)
        #fixed=c("age_class", "samcam"),
        name <- paste("Ew_bv.dredge",i,names(dt.rsp.abn)[i], sep = ".")
        assign(name, GM.dredge)
        ls.bdv.dredge[[i]] <- assign(name, GM.dredge)
}

save(list=c("ls.abn.dredge", "ls.abn.Mglobal", "ls.bms.dredge", "ls.bms.Mglobal", "ls.bdv.dredge", "ls.bdv.Mglobal"), 
     file="C:/Users/Quentin/Documents/git_repositories/F2_Earthworms/Analysis/OutputTables/F2_EW_PdredgeAll.rda")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










