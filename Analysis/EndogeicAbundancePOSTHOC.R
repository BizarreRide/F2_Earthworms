#################
# F2_Eearthworms
# Post-Hoc Multicomparisons for Endogeic Abundance
# Quentin Schorpp
# 18.05.2015
#################



##############################################################
# Post-Hoc Multicomparisons for age class x samcam, 
# interaction effect - reduced contrast matrix ####

# Model with the interaction term
data$ia.acl.smc <- interaction(data$age_class, data$samcam)
endad.tukey <- glmer(endad ~ ia.acl.smc + scl.prec1 + scl.mc*scl.pH + (1|field.ID) + offset(log(area)) ,data=data,family="poisson", control=glmerControl(optimizer="bobyqa"))
# summary(endad.tukey)

# Create contrast matrix
source("Analysis/F2_EW_ContrastMatrix.R")

# Pairwise comparisons (with interaction term)
endad.pairwise <- glht(endad.tukey, linfct=mcp(ia.acl.smc = cm1))
endad.pw.ci <- confint(endad.pairwise)
summary(endad.pairwise, test=adjusted(type="fdr"))

# Confidence intervals including 0
endad.pw.sig <- which(endad.pw.ci$confint[,2]>0)
data.frame(names(endad.pw.sig))

# Plot Errorbars
phfig1 <- ggplot(endad.pw.ci, aes(y = lhs, x = exp(estimate), xmin = exp(lwr), xmax = exp(upr))) + 
  geom_errorbarh() + 
  geom_point() + 
  geom_vline(xintercept = 1) +
  mytheme
phfig1

# plot confidence intervals
par(mar=c(2,15,2,2))
plot(endad.pw.ci) 
##############################################################



##############################################################
# Post-Hoc Multicomparisons for age class, 
# main effect ####

# Pairwise comparisons (without interaction term)
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
plot(endad.pw.ci) 
##############################################################



##############################################################
# Post-Hoc Multicomparisons for samcam, 
# main effect ####

# Pairwise comparisons (without interaction term)
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
plot(endad.pw.ci) 
##############################################################
