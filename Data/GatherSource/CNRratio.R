###############
# C and N Measurements
# data in %
# Calculation of CN ratio
# Claculation of average C and N content
###############

# Load data
CNpercent <- read.delim("Data/CNpercent.txt")

summary(CNpercent$Feld)
str(CNpercent)

# create field ID to maintain row order
CNpercent$field.ID <- rep(1:45,each=5)

# calculate C:N ratio
CNpercent$cn <- CNpercent$C/CNpercent$N

# build averages
cn <- aggregate(. ~ field.ID, CNpercent,mean)

cn <- cn[rep(seq_len(nrow(cn)), each=4),]

data$cn <- cn$cn
data$c <- cn$C
data$n <- cn$N
