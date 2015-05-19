#################
# F2_Eearthworms
# Global Model Dredge
# Quentin Schorpp
# 19.05.2015
#################




##############################################################
# With glmer()


endad.glob1 <- glmer(endad ~ age_class*samcam + I(scl.mc^2) + scl.mc*scl.pH          + scl.pH*scl.cn + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))
endad.glob2 <- glmer(endad ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))


#endo.glob1 <- glmer(endo ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))
#endo.glob2 <- glmer(endo ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))


anc.glob1 <- glmer(anc ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))
anc.glob2 <- glmer(anc ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family=poisson, control=glmerControl(optimizer="bobyqa"))


options(na.action = "na.fail")
endad.dredge1 <- dredge(endad.glob1)
endad.dredge2 <- dredge(endad.glob2)
anc.dredge1 <- dredge(anc.glob1)
anc.dredge2 <- dredge(anc.glob2)
options(na.action = "na.omit")

save(endad.dredge1, endad.dredge2, anc.dredge1, anc.dredge2, file="Analysis/F2_EW_glmerDredge.RData")


summary(endad.dredge1.2)

head(endad.dredge1,10)
endad.avgmod1.d4 <- model.avg(endad.dredge1, subset = delta < 4)
summary(endad.avgmod1.d4)

head(endad.dredge2,10)
endad.avgmod2.d4 <- model.avg(endad.dredge2, subset = delta < 4)
summary(endad.avgmod2.d4)

head(anc.dredge1,10)
anc.avgmod1.2d4 <- model.avg(anc.dredge1, subset = delta < 4)
summary(anc.avgmod1.d4)

head(anc.dredge2,10)
anc.avgmod2.d4 <- model.avg(anc.dredge2, subset = delta < 4)
summary(anc.avgmod2.d4)
##############################################################



##############################################################
# With glmmadmb()


endad.glob1 <- glmmadmb(endad ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.ph*scl.cn + scl.prec1 + scl.clay + crop + (1|field.ID) + offset(log(area)),data=data,family="poisson")
endad.glob2 <- glmmadmb(endad ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="poisson")


#endo.glob1 <- glmmadmb(endo ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="poisson")
#endo.glob2 <- glmmadmb(endo ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="poisson")


anc.glob1 <- glmmadmb(anc ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="poisson")
anc.glob2 <- glmmadmb(anc ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="poisson")



endad.dredge1.1 <- dredge(endad.glob1)
endad.dredge2.2 <- dredge(endad.glob2)
anc.dredge1.2   <- dredge(anc.glob1)
anc.dredge2.2   <- dredge(anc.glob2)
##############################################################


# Biomass
# Anecic
anc.bm.glob1.gs <- glmmadmb(I(anc.bm+0.0001) ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data, family="gaussian")
anc.bm.glob1.gm <- glmmadmb(I(anc.bm+0.0001) ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="gamma")
anc.bm.glob2.gs <- glmmadmb(I(anc.bm+0.0001) ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="gaussian")
anc.bm.glob2.gm <- glmmadmb(I(anc.bm+0.0001) ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="gamma")

anc.bm.dredge1.gs <- dredge(anc.bm.glob1.gs)
anc.bm.dredge1.gm <- dredge(anc.bm.glob2.gm)
anc.bm.dredge2.gs <- dredge(anc.bm.glob1.gs)
anc.bm.dredge2.gm <- dredge(anc.bm.glob2.gm)

#Endogeic
endo.bm.glob1.gs <- glmmadmb(I(endo.bm+0.0001) ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data, family="gaussian")
endo.bm.glob1.gm <- glmmadmb(I(endo.bm+0.0001) ~ age_class*samcam + scl.mc + I(scl.mc^2) + scl.mc*scl.pH + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="gamma")
endo.bm.glob2.gs <- glmmadmb(I(endo.bm+0.0001) ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="gaussian")
endo.bm.glob2.gm <- glmmadmb(I(endo.bm+0.0001) ~ age_class*samcam + scl.ats1 + I(scl.ats1^2)             + scl.pH*scl.cn  + scl.prec1 + scl.clay + (1|field.ID) + offset(log(area)),data=data,family="gamma")

endo.bm.dredge1.gs <- dredge(endo.bm.glob1.gs)
endo.bm.dredge1.gm <- dredge(endo.bm.glob2.gm)
endo.bm.dredge2.gs <- dredge(endo.bm.glob1.gs)
endo.bm.dredge2.gm <- dredge(endo.bm.glob2.gm)

save(anc.bm.dredge1.gs,anc.bm.dredge1.gm,anc.bm.dredge2.gs,anc.bm.dredge2.gm, endo.bm.dredge1.gs,endo.bm.dredge1.gm,endo.bm.dredge2.gs,endo.bm.dredge2.gm, file="Analysis/F2_EW_BiomassDredge.RData")