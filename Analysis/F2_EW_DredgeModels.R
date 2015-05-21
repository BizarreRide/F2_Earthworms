#################
# F2_Eearthworms
# Global Model Dredge
# Quentin Schorpp
# 19.05.2015
#################

##############################################################
## Create subset matrix for dredge() ####

# Suppose we want to have set of models that exclude combinations of colinear
# variables, that are significantly (p < 0.05) correlated, with Pearson
# correlation coefficient larger than r = 0.5.

std.var2 <- std.var[,c("scl.ats1","scl.cn", "scl.prec1", "scl.clay", "scl.mc", "scl.pH")]
std.var2 <- cbind(std.var2, "I(scl.mc^2)"=std.var$scl.mc^2,"I(scl.ats1^2)" = std.var$scl.ats1^2)#,age_class=data$age_class, age_class=data$samcam)

is.correlated <- function(i, j, data, conf.level = .95, cutoff = .5, ...) {
  if(j >= i) return(NA)
  ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}

# Need vectorized function to use with 'outer'
vCorrelated <- Vectorize(is.correlated, c("i", "j"))

# Create logical matrix
smat <- outer(1:8, 1:8, vCorrelated, data = std.var2)

nm <- colnames(std.var2[1:8])

dimnames(smat) <- list(nm, nm)

# Although the squared terms seem not to correlate between mc and ats, we won't fit them within the same model
smat[5,1] <- FALSE
smat[7,1]<- FALSE
smat[8,5]<- FALSE
smat[8,7]<- FALSE

smat
##############################################################



##############################################################
## Global Models for Abundance ####

# With glmer()
anc.glob1 <- glmer(anc ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))
anc.glob2 <- glmer(anc ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))

endad.glob1 <- glmer(endad ~ age_class*samcam + I(scl.mc^2) + scl.mc*scl.pH          + scl.pH*scl.cn + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))
endad.glob2 <- glmer(endad ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))

#endo.glob1 <- glmer(endo ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))
#endo.glob2 <- glmer(endo ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))

options(na.action = "na.fail")
endad.dredge1 <- dredge(endad.glob1)
endad.dredge2 <- dredge(endad.glob2)
anc.dredge1 <- dredge(anc.glob1)
anc.dredge2 <- dredge(anc.glob2)
options(na.action = "na.omit")

save(endad.dredge1, endad.dredge2, anc.dredge1, anc.dredge2, file="Analysis/F2_EW_glmerDredge.RData")
##############################################################

##############################################################
## Multimodel averaging and subset selection ####

head(anc.dredge1,10)
anc.avgmod1.2d4 <- model.avg(anc.dredge1, subset = delta < 4)
summary(anc.avgmod1.d4)

head(anc.dredge2,10)
anc.avgmod2.d4 <- model.avg(anc.dredge2, subset = delta < 4)
summary(anc.avgmod2.d4)

head(endad.dredge1,10)
endad.avgmod1.d4 <- model.avg(endad.dredge1, subset = delta < 4)
summary(endad.avgmod1.d4)

head(endad.dredge2,10)
endad.avgmod2.d4 <- model.avg(endad.dredge2, subset = delta < 4)
summary(endad.avgmod2.d4)
##############################################################

##############################################################
## With glmmadmb() ####

anc.glob1 <- glmmadmb(anc ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="poisson")
anc.glob2 <- glmmadmb(anc ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="poisson")

endad.glob1 <- glmmadmb(endad ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.ph*scl.cn + scl.prec1 + scl.clay + crop + (1|field.ID) + offset(log(area)),data=data,family="poisson")
endad.glob2 <- glmmadmb(endad ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="poisson")

#endo.glob1 <- glmmadmb(endo ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="poisson")
#endo.glob2 <- glmmadmb(endo ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="poisson")

anc.dredge1.2   <- dredge(anc.glob1)
anc.dredge2.2   <- dredge(anc.glob2)

endad.dredge1.2 <- dredge(endad.glob1)
endad.dredge2.2 <- dredge(endad.glob2)

save(endad.dredge1.2, endad.dredge2.2, anc.dredge1.2, anc.dredge2.2, file="Analysis/F2_EW_glmmadmbDredge.RData")
##############################################################



##############################################################
## Global Models for Biomass ####
# Biomass // 2nd Try; log normal distribution
# preferred over gamma since the variance seems not to be positively correlated with the mean.

# Anecic
anc.bm.glob1 <- lmer(log1p(anc.bm) ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc:scl.pH + + scl.ats1 + I(scl.ats1^2) + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data)

getAllTerms(anc.bm.glob3)

options(na.action = "na.fail")
anc.bm.dredge1 <- dredge(anc.bm.glob3, subset= smat)
options(na.action = "na.omit")


# Endogeic
endo.bm.glob1 <- lmer(log1p(endo.bm) ~ age_class*samcam + scl.ats1 + I(scl.ats1^2) + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data)

options(na.action = "na.fail")
endo.bm.dredge1 <- dredge(endo.bm.glob1, subset=smat)
options(na.action = "na.omit")


# Endogeic adult
endad.bm.glob1 <- lmer(log1p(endad.bm) ~ age_class*samcam + scl.ats1 + I(scl.ats1^2) + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data)

options(na.action = "na.fail")
endad.bm.dredge1 <- dredge(endad.bm.glob1, subset=smat)
options(na.action = "na.omit")

save(anc.bm.dredge1,endad.bm.dredge1,endo.bm.dredge1, file="Analysis/F2_EW_BiomassDredge2.RData")
##############################################################

##############################################################
## Multimodel averaging and subset selection ####

head(anc.bm.dredge1,10)
anc.bm.avgmod1.d4 <- model.avg(anc.bm.dredge1, subset = delta < 4)
summary(anc.bm.avgmod1.d4)

head(endo.bm.dredge1,10)
summary(endo.bm.dredge1)
endo.bm.avgmod1.d4 <- model.avg(endo.bm.dredge1, subset = delta < 4)
summary(endo.bm.avgmod1.d4)

head(endad.bm.dredge1,10)
endad.bm.avgmod1.d4 <- model.avg(endad.bm.dredge1, subset = delta < 4)
summary(endad.bm.avgmod1.d4)
##############################################################






# Biomass // First Try
# Anecic
anc.bm.glob1.gs <- glmmadmb(I(anc.bm+0.0001) ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data, family="gaussian")
anc.bm.glob1.gm <- glmmadmb(I(anc.bm+0.0001) ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="gamma")
anc.bm.glob2.gs <- glmmadmb(I(anc.bm+0.0001) ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="gaussian")
anc.bm.glob2.gm <- glmmadmb(I(anc.bm+0.0001) ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="gamma")

anc.bm.dredge1.gs <- dredge(anc.bm.glob1.gs)
anc.bm.dredge1.gm <- dredge(anc.bm.glob1.gm)
anc.bm.dredge2.gs <- dredge(anc.bm.glob2.gs)
anc.bm.dredge2.gm <- dredge(anc.bm.glob2.gm)

#Endogeic
endo.bm.glob1.gs <- glmmadmb(I(endo.bm+0.0001) ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data, family="gaussian")
endo.bm.glob1.gm <- glmmadmb(I(endo.bm+0.0001) ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="gamma")
endo.bm.glob2.gs <- glmmadmb(I(endo.bm+0.0001) ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="gaussian")
endo.bm.glob2.gm <- glmmadmb(I(endo.bm+0.0001) ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="gamma")

endo.bm.dredge1.gs <- dredge(endo.bm.glob1.gs)
endo.bm.dredge1.gm <- dredge(endo.bm.glob1.gm)
endo.bm.dredge2.gs <- dredge(endo.bm.glob2.gs)
endo.bm.dredge2.gm <- dredge(endo.bm.glob2.gm)

#Endogeic adult
endad.bm.glob1.gs <- glmmadmb(I(endad.bm+0.0001) ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data, family="gaussian")
endad.bm.glob1.gm <- glmmadmb(I(endad.bm+0.0001) ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="gamma")
endad.bm.glob2.gs <- glmmadmb(I(endad.bm+0.0001) ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="gaussian")
endad.bm.glob2.gm <- glmmadmb(I(endad.bm+0.0001) ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="gamma")

endad.bm.dredge1.gs <- dredge(endad.bm.glob1.gs)
endad.bm.dredge1.gm <- dredge(endad.bm.glob1.gm)
endad.bm.dredge2.gs <- dredge(endad.bm.glob2.gs)
endad.bm.dredge2.gm <- dredge(endad.bm.glob2.gm)

save(anc.bm.dredge1.gs,anc.bm.dredge1.gm,anc.bm.dredge2.gs,anc.bm.dredge2.gm, endad.bm.dredge1.gs,endad.bm.dredge1.gm,endad.bm.dredge2.gs,endad.bm.dredge2.gm,endo.bm.dredge1.gs,endo.bm.dredge1.gm,endo.bm.dredge2.gs,endo.bm.dredge2.gm, file="Analysis/F2_EW_BiomassDredge.RData")

head(anc.bm.dredge1.gs,10)
anc.bm.avgmod1.gsd4 <- model.avg(anc.bm.dredge1.gs, subset = delta < 4)
summary(anc.bm.avgmod1.gsd4)

head(anc.bm.dredge1.gm,30)
anc.bm.avgmod1.gmd4 <- model.avg(anc.bm.dredge1.gm, subset = delta < 4)
summary(anc.bm.avgmod2.gmd4)

head(anc.bm.dredge2.gs,10)
anc.bm.avgmod2.gsd4 <- model.avg(anc.bm.dredge2.gs, subset = delta < 4)
summary(anc.bm.avgmod1.gsd4)

head(anc.bm.dredge2.gm,30)
anc.bm.avgmod2.gmd4 <- model.avg(anc.bm.dredge2.gm, subset = delta < 4)
summary(anc.bm.avgmod2.gmd4)


# I gave up on this strategy, since the addition of 0.0001 is generally concerned as faulty!




