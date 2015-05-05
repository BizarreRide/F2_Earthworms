#############
# Create Dummies for Excel Pivot tables
# Quentin Schorpp
# 05.05.2015
# F2_Earthworms
##############

dum <- read.delim("Data/F2_EW_Dummies.txt")

library(reshape2)

dum_melt <- melt(dum, id.vars=1:7, value.name="weights", variable="species")

write.csv(dum_melt, "Data/F2_EW_Dummies2.csv")

