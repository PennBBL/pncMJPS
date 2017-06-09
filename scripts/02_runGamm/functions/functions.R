# Create a function which will run a gam model
# given an df and what to grep for
# It will return a list, this list will include
# a 2 x 2 table with N of signifianct interaction terms
# the ROI's that had the signifianct interaction and the 
runGamModel <- function(dataFrame, grepPattern, qualityIndex){
  # First find all fo the columns that we will loop through
  columnIndex <- grep(grepPattern, names(dataFrame))
  qualityCol <- grep(qualityIndex, names(dataFrame))
  # Now declare our statics
  ageCol <- which(names(dataFrame)=='ageAtScan1')
  sexCol <- which(names(dataFrame)=='sex')
  mjCol <- which(names(dataFrame)=='marcat')
  psCol <- which(names(dataFrame)=='goassessDxpmr7')
  # Declare som values for which will be returned from the gam function
  mjUserPS <- NULL
  mjUserTD <- NULL
  nonUserPS <- NULL
  nonUserTD <- NULL
  # Now loop through each region and run the model 
  for(i in columnIndex){
    # Get the roi name
    roiName <- names(dataFrame)[i]
    # First run the model 
    tmp <- gam(dataFrame[,i] ~ s(dataFrame[,ageCol]) + dataFrame[,sexCol] + dataFrame[,mjCol]*dataFrame[,psCol] + dataFrame[,qualityCol])
    # Now grap our interaction terms
    summaryVals <- summary(tmp)
    sigTerms <- summaryVals$p.table[7:11,4]
  }
}
