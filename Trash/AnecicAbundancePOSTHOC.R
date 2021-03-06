#################
# F2_Eearthworms
# Post-Hoc Multicomparisons for Anecic Abundance
# Quentin Schorpp
# 18.05.2015
#################



##############################################################
# Post-Hoc Multicomparisons for age class x samcam, 
# interaction effect - reduced contrast matrix ####

# Model with the interaction term
data$ia.acl.smc <- interaction(data$age_class, data$samcam)
anc.tukey <- glmer(anc ~ ia.acl.smc + I(scl.ats1^2) + scl.prec1 + (1|field.ID) + offset(log(area)) ,data=data,family="poisson", control=glmerControl(optimizer="bobyqa"))
# summary(anc.tukey)

# Create contrast matrix
source("Analysis/F2_EW_ContrastMatrix.R")

# Pairwise comparisons (with interaction term)
anc.pairwise <- glht(anc.tukey, linfct=mcp(ia.acl.smc = cm1))
anc.pw.ci <- confint(anc.pairwise)
summary(anc.pairwise, test=adjusted(type="fdr"))

# Confidence intervals including 0
anc.pw.sig <- which(anc.pw.ci$confint[,2]>0)
data.frame(names(anc.pw.sig))

# Plot Errorbars
phfig1 <- ggplot(anc.pw.ci, aes(y = lhs, x = exp(estimate), xmin = exp(lwr), xmax = exp(upr))) + 
            geom_errorbarh() + 
            geom_point() + 
            geom_vline(xintercept = 1) +
            mytheme
phfig1

# plot confidence intervals
par(mar=c(2,15,2,2))
  plot(anc.pw.ci) 
##############################################################



##############################################################
# Post-Hoc Multicomparisons for age class, 
# main effect ####

# Pairwise comparisons (without interaction term)
anc.pairwise <- glht(anc.best, mcp(age_class = "Tukey"))
anc.pw.ci <- confint(anc.pairwise)
summary(anc.pairwise, test=adjusted(type="none"))

# Confidence intervals including 0
anc.pw.sig <- which(anc.pw.ci$confint[,2]>0)
data.frame(names(anc.pw.sig))


# Plot Errorbars
ggplot(anc.pw.ci, aes(y = lhs, x = exp(estimate), xmin = exp(lwr), xmax = exp(upr))) + 
  geom_errorbarh() + 
  geom_point() + 
  geom_vline(xintercept = 1) +
  mytheme

# plot confidence intervals
par(mar=c(2,15,2,2))
plot(anc.pw.ci) 
##############################################################



##############################################################
# Post-Hoc Multicomparisons for samcam, 
# main effect ####

# Pairwise comparisons (without interaction term)
anc.pairwise <- glht(anc.best, mcp(samcam = "Tukey"))
anc.pw.ci <- confint(anc.pairwise)

# Confidence intervals including 0
anc.pw.sig <- which(anc.pw.ci$confint[,2]>0)
data.frame(names(anc.pw.sig))


# Plot Errorbars
ggplot(anc.pw.ci, aes(y = lhs, x = exp(estimate), xmin = exp(lwr), xmax = exp(upr))) + 
  geom_errorbarh() + 
  geom_point() + 
  geom_vline(xintercept = 1) +
  mytheme

# plot confidence intervals
par(mar=c(2,15,2,2))
plot(anc.pw.ci) 
##############################################################
