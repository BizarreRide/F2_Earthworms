#########
# Only adult Endogeics
# Quentin Schorpp
# 30.04.2015
##############

EndoAdult <- read.delim("~/git_repositories/F2_Earthworms/Data/RW_Alle Aufnahmen.txt")

data$endad <- rowSums(EndoAdult[,c("AROad", "ACAad", "ACHad", "ENDad","OCY","OLA")])
