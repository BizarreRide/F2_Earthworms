#################
# F2_Eearthworms
# Load Data and Packages
# Quentin Schorpp
# 07.05.2015
#################


##############################################################
# Load Data

source("Data/GatherSource/Gather1.R") # Includes Correction for Watercontent
source("Data/GatherSource/Gather2.R")
source("Data/GatherSource/Gather3.R") # no adult endogeics separated
source("Data/GatherSource/CNRatio.R") # only for 1st order
##############################################################

##############################################################
# Load packages and functions

source("Analysis/F2_EW_RequiredPackages.R")
source("Data/GatherSource/Functions/coldiss.R")
source("Data/GatherSource/Functions/panelutils.R")
source("Data/GatherSource/Functions/HighStatLibV8.R")
source("Data/GatherSource/ggplotTheme.R")


lo <- par(mfrow=c(1,1))

overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

##############################################################
