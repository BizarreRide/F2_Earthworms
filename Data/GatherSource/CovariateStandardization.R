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

# Create variable d(delta)samcam as the time difference of each sampling date to the first date, that each field was sampled; first date is Zero for all fields. ####
data$dsamcam <- with(data, round(c(rep(0,60),date[61:120]-date[1:60],date[121:180]-date[1:60]),0))

# create covariate for offset, area = area of the pit where the soil-core was taken from
data$area <- rep(0.25,180)

rm(std.var)
