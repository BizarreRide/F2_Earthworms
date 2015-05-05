#########
# Only adult Endogeics
# Quentin Schorpp
# 30.04.2015
##############

Seaparate abundance data of endogeic earthworms into adult and juveniles

EndoAdult <- read.delim("~/git_repositories/F2_Earthworms/Data/F2_EW_EndoAdult.txt")

data$endad <- rowSums(EndoAdult[,c("AROad", "ACAad", "ACHad", "ENDad","OCY","OLA")])
