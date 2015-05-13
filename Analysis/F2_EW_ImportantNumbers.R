#################
# F2_Earthworms
# Calculation of some important numbers
# Quentin Schorpp
# 07.05.2015
#################


# How high was the fraction of A. longa on total anecic abundance in that single field (PO_2009), 
# where A. longa was found frequently?

sub <- subset(data, subset=field.ID==9, select=c("field.ID","LTR", "ALO"))
Alo.frac <- sum(sub$ALO)/(sum(sub$LTR)+sum(sub$ALO))


# How high was the fraction of adult anecic Individuals on total anecic abundance 
# in age_classes for each sampling campaign?

anc.mean <- with(data2, tapply(anc, list(samcam,age_class), mean))
ancad.mean <- with(data2, tapply(ancad, list(samcam,age_class), mean))

ancad.frac <- ancad.mean/anc.mean

ancad.frac-1
