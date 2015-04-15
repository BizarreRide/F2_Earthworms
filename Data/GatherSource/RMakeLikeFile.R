#############
# Example R make.like file
# Quentin Schorpp
# Updated 15 April 2015
##############

# Set working directory
setwd("D:/Quentin_Schorpp/Arbeitsprozess/git_repositories/F2_Earthworms/Data/")

# Gather and clean up raw data files.
source("GatherSource/Gather1.R") # 1st Order data [180 rows]

source("GatherSource/Gather2.R") # 2nd Order data [45 rows]

source("GatherSource/Gather3.R") # 3rd Order data [15 rows]

# Merge cleaned data frames into data frame object CleanedData
source("GatherSource/MergeData.R")  # The MergeData.R file merges the data frames and saves the output data frame as CSV formatted file
                                    # In my Example it could probably create a RData file with all the different datasets?? However .RData files should be ignored by Github
