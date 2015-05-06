---
title: "Anecic Abundance"
author:
  - name: Quentin Schorpp
    email: qurntin.schorpp@ti.bund.de
date: "2015-05-05"
output: 
      html_document:
          toc: true
          toc_depth: 2
          number_sections: true
          fig_caption: true
          theme: cerulean
          highlight: haddock
          
---



# Introduction  

I investigated how earthworm communities change in perennial energy farming systems. To assess the development in the long term, I chose an approach called Space-for-Time Substitution (Pickett, 1989). Instead of sampling sites over a long time period, I studied different aged sites. The sites were cropped either with the perennial bioenergy plant cup-plant (*Silphium perfoliatum*) or with annual Maize (*Zea mays*). The Main question of this investigation is wether older cup plant fields possess higher structural and functional biodiversity than younger ones and to compare cup-plant to maize. To distinguish the age-effect from environmental influences, several parameters regarding soil properties, climate and Management were collected.
Sampling took place 3-4 times in two consecutive years. In total 18 sites were surveyed, comprising 4 cup-plant age classes and maize. The sites were located at seven different locations, comprising fields from research and business. Each cup-plant age-class had three replicates, whereas maize plants were switched in the second year of sampling due to crop rotation and gained six replicates. To sample earthworms 4 pits with 0.25 m² were digged up to 20 cm in depth and both handsorting of soil and chemical extractioon with AITC took place. These 4 pits were randomly chosen at the field, hence they're nested in field.
9 earthworm species from 3 different ecological groups were identified and weighed.

Here I present the modeling of the abundance data of earthworms from the ecological group of anecic earthworms. 

# Materials and methods  

I use R with RStudio to perform all analyses. For creating reports i use knitr with RStudio.

## Modeling procedure:

* I start with a global model, with a full set of parameters. 
  
    + The choice of parameters took place according to preknowledge  about the biology of earthworms. 
    + Management parameters, are not included up to now. 

* Since Abundance is count data, i choose poisson distribution in the first place, and switch to negative binomial if overdipersion can be detected.

* As random effects i chose field. I could have chosen sampling campaign, too and analyse the model with crossed random effects structure, but I'm also interested in change over the time of sampling.

* I included an offset "area", since i want to present the data per m² not per 0.25 m².

* several global models were formulated if there are different combinations of parameters, meaningful for analysis, but not able to analyse at once, i.e. because of  collinearity.

* Then I use Multimodel Inference/Averaging to find the model with the lowest AIC and determine it as the best model.

* Model validation is performed by 

    + plotting residuals against all covariates in the model as well as against all covariates not in the model
    + plotting a histogram of the residual
    + checking for overdispersion

* after having found a model with good model fit, i present the data for the means of environmental covariates and plot barplots for the age classes.

* errorbars are calculated, using the inverse link function with confidence intervals.

* Post-Hoc Multicomparisons are performed using Tukey Test with the glht() function.

# Analysis






### Data Pre-Processing  
#### Standardize Explanatory Variables  




### Anecic Abundance
********************

#### Assess variability in random effects
<img src="./AnecicAbundance_files/figure-html/BoxPlots1.png" title="plot of chunk BoxPlots" alt="plot of chunk BoxPlots" style="display: block; margin: auto;" /><img src="./AnecicAbundance_files/figure-html/BoxPlots2.png" title="plot of chunk BoxPlots" alt="plot of chunk BoxPlots" style="display: block; margin: auto;" /><img src="./AnecicAbundance_files/figure-html/BoxPlots3.png" title="plot of chunk BoxPlots" alt="plot of chunk BoxPlots" style="display: block; margin: auto;" />
  

#### GLMM using glmmADMB


##### Model Formulation

Find the best model using function dredge(). Begin with three different global models.

```r
# Global Model

# str(data)

# Model Formulation
anc.glob1 <- glmmadmb(anc ~ age_class*samcam + scl.mc + scl.pH*scl.cn  +I(scl.hum1^2) + (1|field.ID) + offset(log(area)),data=data,family="poisson")
anc.glob2 <- glmmadmb(anc ~ age_class*samcam + I(scl.mc^2) + scl.pH*scl.cn  +I(scl.hum1^2) + (1|field.ID) + offset(log(area)),data=data,family="poisson")
anc.glob3 <- glmmadmb(anc ~ age_class*samcam + I(scl.ats1^2) + scl.pH*scl.cn  +I(scl.hum1^2) + (1|field.ID) + offset(log(area)),data=data,family="poisson")
# offsset is used due to the personal advice by T. Onkelinx:
# You better use an offset if you want to express the model in terms of m². Just add offset(log(0.25)) to the model. 

# summary(anc.glob3)

# Overdispersion
E1 <- resid(anc.glob3, type="pearson")
N <- nrow(data)
p <- length(coef(anc.glob3)) # The Plus "+1" used for determining the number of parameters(p) is due to the "k" of negative 
                              # binomial distribution
Dispersion <- sum(E1^2)/(N-p)
Dispersion
```

```
## [1] 1.37
```

```r
## Multimodel averaging ####
#anc.dredge1 <- dredge(anc.glob1)
#anc.dredge2 <- dredge(anc.glob2)
#anc.dredge3 <- dredge(anc.glob3)
#head(anc.dredge3,10)
#head(anc.dredge2,10)
#head(anc.dredge,10)
#####
```

Fit the best model

```
## 
## Call:
## glmmadmb(formula = anc ~ age_class * samcam + I(scl.ats1^2) + 
##     (1 | field.ID) + offset(log(area)), data = data, family = "poisson")
## 
## AIC: 736 
## 
## Coefficients:
##                     Estimate Std. Error z value Pr(>|z|)    
## (Intercept)            2.220      0.191   11.60  < 2e-16 ***
## age_class.L            1.985      0.435    4.57 0.000005 ***
## age_class.Q           -0.590      0.401   -1.47  0.14136    
## age_class.C           -0.249      0.394   -0.63  0.52845    
## age_class^4            0.660      0.353    1.87  0.06134 .  
## samcam2                0.679      0.197    3.45  0.00056 ***
## samcam3                0.225      0.181    1.24  0.21575    
## I(scl.ats1^2)         -0.423      0.100   -4.23 0.000023 ***
## age_class.L:samcam2   -0.316      0.543   -0.58  0.56046    
## age_class.Q:samcam2   -0.220      0.457   -0.48  0.63077    
## age_class.C:samcam2    0.872      0.362    2.41  0.01611 *  
## age_class^4:samcam2   -0.848      0.249   -3.41  0.00065 ***
## age_class.L:samcam3    0.155      0.510    0.30  0.76065    
## age_class.Q:samcam3   -0.210      0.441   -0.48  0.63428    
## age_class.C:samcam3    0.524      0.361    1.45  0.14605    
## age_class^4:samcam3   -0.500      0.257   -1.95  0.05168 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Number of observations: total=180, field.ID=18 
## Random effect variance(s):
## Group=field.ID
##             Variance StdDev
## (Intercept)   0.2499 0.4999
## 
## 
## Log-likelihood: -351
```

```
##                   Df Sum Sq Mean Sq F value  Pr(>F)    
## age_class          4    761     190   15.63 7.4e-11 ***
## samcam             2     67      34    2.75   0.067 .  
## I(scl.ats1^2)      1    332     332   27.27 5.3e-07 ***
## age_class:samcam   8     90      11    0.92   0.501    
## Residuals        164   1997      12                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

A former "best model" with humidity as significant term and age_class as ordered factor, found using dredge was
anc.best <- glmmadmb(anc ~ age_class + scl.dsamcam + I(scl.ats1^2) +I(scl.hum1^2) + (1|field.ID) + offset(log(area)) ,data=data,family="poisson")

##### Model Validation

```r
confint(anc.best)
```

```
##                       2.5 %   97.5 %
## (Intercept)          1.8445  2.59482
## age_class.L          1.1327  2.83670
## age_class.Q         -1.3753  0.19612
## age_class.C         -1.0214  0.52427
## age_class^4         -0.0314  1.35111
## samcam2              0.2932  1.06479
## samcam3             -0.1310  0.58003
## I(scl.ats1^2)       -0.6195 -0.22731
## age_class.L:samcam2 -1.3810  0.74849
## age_class.Q:samcam2 -1.1158  0.67632
## age_class.C:samcam2  0.1617  1.58191
## age_class^4:samcam2 -1.3354 -0.36079
## age_class.L:samcam3 -0.8437  1.15429
## age_class.Q:samcam3 -1.0733  0.65405
## age_class.C:samcam3 -0.1827  1.23159
## age_class^4:samcam3 -1.0029  0.00363
```

```r
coefplot2(anc.best)
```

<img src="./AnecicAbundance_files/figure-html/ModelValidation1_anc1.png" title="plot of chunk ModelValidation1_anc" alt="plot of chunk ModelValidation1_anc" style="display: block; margin: auto;" />

```r
## Model validation ####
E1 <- resid(anc.best, type="pearson")
F1 <- fitted(anc.best, type="response")

# Dataset with all variables IN the model
env1 <- cbind(E1, F1, data[,c("anc", "age_class", "samcam", "ats1")]) # should "ats1" be squared?

# covariates NOT in the model
env0 <- cbind(E1, F1,data[,c("mc", "pH", "cn", "clay", "ats2", "hum1","dgO")])

# plot residual versus all covariates in the model
par(mfrow=c(2,3),
    mar=c(3.8,4,1,3))
for (i in 1:length(env1)) {
  scatter.smooth(env1[,i],env1[,1], xlab=colnames(env1)[i])
  abline(h=0, lty=2, col="red")
}
```

<img src="./AnecicAbundance_files/figure-html/ModelValidation1_anc2.png" title="plot of chunk ModelValidation1_anc" alt="plot of chunk ModelValidation1_anc" style="display: block; margin: auto;" />

```r
par(mfrow=c(1,1))

# plot residual versus all covariates NOT in the model
par(mfrow=c(3,3),
    mar=c(3.8,4,1,3))
for (i in 1:length(env0)) {
  scatter.smooth(env0[,i],env0[,1], xlab=colnames(env0)[i])
  abline(h=0, lty=2, col="red")
}
```

<img src="./AnecicAbundance_files/figure-html/ModelValidation1_anc3.png" title="plot of chunk ModelValidation1_anc" alt="plot of chunk ModelValidation1_anc" style="display: block; margin: auto;" />

```r
par(mfrow=c(1,1))

# Overdispersion
N <- nrow(data)
p <- length(coef(anc.best)) # The Plus "+1" used for determining the number of parameters(p) is due to the "k" of negative 
                        # binomial distribution
Dispersion <- sum(E1^2)/(N-p)
Dispersion
```

```
## [1] 1.33
```

```r
#####

overdisp_fun(anc.best)
```

```
##     chisq     ratio       rdf         p 
## 217.80558   1.33623 163.00000   0.00267
```

```r
# Deviance
logLik(anc.best) 
```

```
## 'log Lik.' -351 (df=17)
```

```r
# Fraction of explained variation in the response variable
100*((anc.best$null.deviance-anc.best$deviance)/anc.best$null.deviance) # failed
```

```
## numeric(0)
```

```r
# Plots of Predicted Values
par(mfrow=c(2,2))
plot(data$age_class,predict(anc.best, type="response"))
plot(data$samcam,predict(anc.best, type="response"))
plot(data$ats1^2,predict(anc.best, type="response"))
plot(data$hum1^2,predict(anc.best, type="response"))
```

<img src="./AnecicAbundance_files/figure-html/ModelValidation1_anc4.png" title="plot of chunk ModelValidation1_anc" alt="plot of chunk ModelValidation1_anc" style="display: block; margin: auto;" />

```r
# Plots of fitted Values
par(mfrow=c(2,2))
plot(data$age_class,F1)
plot(data$samcam,F1)
plot(data$ats1^2,F1)

par(lo)
```

<img src="./AnecicAbundance_files/figure-html/ModelValidation1_anc5.png" title="plot of chunk ModelValidation1_anc" alt="plot of chunk ModelValidation1_anc" style="display: block; margin: auto;" />

Model validation 2

```r
# Model Validation 

E <- resid(anc.best, type="pearson")

par(mfrow=c(2,2))
par(mar=c(5,5,2,3))
plot(E~fitted(anc.best), ylab="Pearson Residuals", xlab="Fitted values") 
abline(h=0, lty=2)
hist(E, main = "", breaks = 20, cex.lab = 1.5, xlab = "Pearson Residuals")
plot(E ~ data$age_class, ylab="Pearson Residuals", xlab="Age Class")
abline(h=0, lty=2)
plot(E ~ data$ats1, ylab="Pearson Residuals", xlab="Age Class")
abline(h=0, lty=2)
```

<img src="./AnecicAbundance_files/figure-html/ModelValidation2_anc.png" title="plot of chunk ModelValidation2_anc" alt="plot of chunk ModelValidation2_anc" style="display: block; margin: auto;" />

##### Post-Hoc Multicomparisons

```r
data$ia.acl.smc <- interaction(data$age_class, data$samcam)
anc.tukey <- glmmadmb(anc ~ ia.acl.smc + I(scl.ats1^2) + (1|field.ID) + offset(log(area)) ,data=data,family="poisson")
summary(anc.tukey)
```

```
## 
## Call:
## glmmadmb(formula = anc ~ ia.acl.smc + I(scl.ats1^2) + (1 | field.ID) + 
##     offset(log(area)), data = data, family = "poisson")
## 
## AIC: 736 
## 
## Coefficients:
##                        Estimate Std. Error z value Pr(>|z|)    
## (Intercept)               0.807      0.529    1.53  0.12700    
## ia.acl.smcB_Sp_young.1    0.470      0.673    0.70  0.48450    
## ia.acl.smcC_Sp_int1.1     2.201      0.614    3.58  0.00034 ***
## ia.acl.smcD_Sp_int2.1     2.040      0.626    3.26  0.00111 ** 
## ia.acl.smcE_Sp_old.1      2.353      0.625    3.76  0.00017 ***
## ia.acl.smcA_Cm.2          0.385      0.802    0.48  0.63148    
## ia.acl.smcB_Sp_young.2    2.265      0.629    3.60  0.00032 ***
## ia.acl.smcC_Sp_int1.2     2.390      0.612    3.90 0.000095 ***
## ia.acl.smcD_Sp_int2.2     2.532      0.619    4.09 0.000044 ***
## ia.acl.smcE_Sp_old.2      2.889      0.619    4.67 0.000003 ***
## ia.acl.smcA_Cm.3         -0.211      0.750   -0.28  0.77806    
## ia.acl.smcB_Sp_young.3    1.272      0.634    2.01  0.04476 *  
## ia.acl.smcC_Sp_int1.3     2.180      0.614    3.55  0.00039 ***
## ia.acl.smcD_Sp_int2.3     2.277      0.613    3.72  0.00020 ***
## ia.acl.smcE_Sp_old.3      2.670      0.610    4.38 0.000012 ***
## I(scl.ats1^2)            -0.423      0.100   -4.23 0.000023 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Number of observations: total=180, field.ID=18 
## Random effect variance(s):
## Group=field.ID
##             Variance StdDev
## (Intercept)   0.2499 0.4999
## 
## 
## Log-likelihood: -351
```

```r
# Pairwise comparisons
anc.pairwise <- glht(anc.tukey, mcp(ia.acl.smc = "Tukey"))
anc.pw.ci <- confint(anc.pairwise)

names(anc.pw.ci)
```

```
##  [1] "model"       "linfct"      "rhs"         "coef"        "vcov"       
##  [6] "df"          "alternative" "type"        "focus"       "confint"
```

```r
anc.pw.sig <- which(anc.pw.ci$confint[,2]>0)
data.frame(names(anc.pw.sig))
```

```
##       names.anc.pw.sig.
## 1  C_Sp_int1.1 - A_Cm.1
## 2   E_Sp_old.1 - A_Cm.1
## 3 B_Sp_young.2 - A_Cm.1
## 4  C_Sp_int1.2 - A_Cm.1
## 5  D_Sp_int2.2 - A_Cm.1
## 6   E_Sp_old.2 - A_Cm.1
## 7  C_Sp_int1.3 - A_Cm.1
## 8  D_Sp_int2.3 - A_Cm.1
## 9   E_Sp_old.3 - A_Cm.1
```

```r
# Plot Errorbars
ggplot(anc.pw.ci, aes(y = lhs, x = exp(estimate), xmin = exp(lwr), xmax = exp(upr))) + 
  geom_errorbarh() + 
  geom_point() + 
  geom_vline(xintercept = 1) +
  mytheme
```

<img src="./AnecicAbundance_files/figure-html/PostHocMulticomparisons_anc1.png" title="plot of chunk PostHocMulticomparisons_anc" alt="plot of chunk PostHocMulticomparisons_anc" style="display: block; margin: auto;" />

```r
# plot confidence intervals

par(mar=c(2,15,2,2))
plot(anc.pw.ci) # only maize is significantly different from all SIlphie fields. Within SIlphie there are no differences
```

<img src="./AnecicAbundance_files/figure-html/PostHocMulticomparisons_anc2.png" title="plot of chunk PostHocMulticomparisons_anc" alt="plot of chunk PostHocMulticomparisons_anc" style="display: block; margin: auto;" />

##### Extract Predictions and plot predictions with error bars








