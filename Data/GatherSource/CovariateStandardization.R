#################
# F2_Eearthworms
# Standardize covariates
# Quentin Schorpp
# 08.05.2015
#################




##############################################################

### Data Pre-Processing  ####

## Standardize all explanatory variables ####

# extract all numerical variables
std.var <- data[sapply(data,is.numeric)]
# exclude all response variables
std.var <- std.var[,-c(7:30,49:61)]


# standardize, each column by {x - mean(x)) / sd(x)}
std.var <- as.data.frame(scale(std.var, center=TRUE, scale=TRUE))
colnames(std.var) <- paste("scl", colnames(std.var), sep = ".") # rename columns

data <- cbind(data, std.var)
##############################################################



