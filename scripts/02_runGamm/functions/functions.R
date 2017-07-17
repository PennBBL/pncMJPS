# Create a function which will run a gam model
# given an df and what to grep for
# It will return a list, this list will include
# a 2 x 2 table with N of signifianct interaction terms
# the ROI's that had the signifianct interaction and the 
runGamModel <- function(dataFrame, grepPattern, qualityIndex){
  # First find all fo the columns that we will loop through
  columnIndex <- grep(grepPattern, names(dataFrame))
  qualityCol <- grep(qualityIndex, names(dataFrame))[1]
  # Now declare our statics
  ageCol <- which(names(dataFrame)=='ageAtScan1')
  sexCol <- which(names(dataFrame)=='sex')
  mjCol <- which(names(dataFrame)=='marcat')
  psCol <- which(names(dataFrame)=='goassessDxpmr7')
  # Declare som values for which will be returned from the gam function
  output <- NULL
  # Now loop through each region and run the model 
  for(i in columnIndex){
    # Get the roi name
    roiName <- names(dataFrame)[i]
    # First run the model 
    tmp <- gam(dataFrame[,i] ~ s(dataFrame[,ageCol]) + dataFrame[,sexCol] + dataFrame[,mjCol]*dataFrame[,psCol] + dataFrame[,qualityCol])
    # Now grap our interaction terms
    summaryVals <- summary(tmp)
    sigTerms <- summaryVals$p.table[8:11,4]
    dim(sigTerms) <- c(1, 4)
    # Now export the variables
    allVals <- cbind(roiName, sigTerms)
    output <- rbind(output, allVals)
  }
  # Now fix the column names
  colnames(output) <- c('ROI', 'PS-Nonuser', 'PS-User', 'TD-Nonuser', 'TD-User')
  return(output)
}


# Return a 2x2 matrix given a p threshold (or perform fdr correction)
# for the number of signifianct interactions from th output of runGamModel
returnTable <- function(inputGamVals){
    # First FDR correct all of the columns
    newValues <- apply(inputGamVals[,2:5], 2, function(x) p.adjust(x, method='fdr'))
    
    # Now find the total number of signifianct regions
    binaryVals <-  matrix(0, nrow=dim(newValues)[1], ncol=dim(newValues)[2])
    
    # Now trun the values less then our 1 .05 value to 1 in the binary matrix
    binaryVals[which(newValues <=.05)] <- 1
    
    # Now prepare our output
    output <- matrix(NA, 2, 2)
    output[1,1] <- length(which(binaryVals[,1]==1))
    output[1,2] <- length(which(binaryVals[,2]==1))
    output[2,1] <- length(which(binaryVals[,3]==1))
    output[2,2] <- length(which(binaryVals[,4]==1))
    # Now return the output
    return(output)
}

runGamModel <- function(dataFrame, grepPattern, qualityIndex){
    # First find all fo the columns that we will loop through
    columnIndex <- grep(grepPattern, names(dataFrame))
    qualityCol <- grep(qualityIndex, names(dataFrame))[1]
    # Now declare our statics
    ageCol <- which(names(dataFrame)=='ageAtScan1')
    sexCol <- which(names(dataFrame)=='sex')
    mjCol <- which(names(dataFrame)=='marcat')
    psCol <- which(names(dataFrame)=='goassessDxpmr7')
    # Declare som values for which will be returned from the gam function
    output <- NULL
    # Now loop through each region and run the model
    for(i in columnIndex){
        # Get the roi name
        roiName <- names(dataFrame)[i]
        formValue <- as.formula(paste(roiName, '~ s(ageAtScan1) + sex + race + marcat*goassessDxpmr7 +', qualityIndex))
        # First run the model
        tmp <- gam(formValue, data=dataFrame)
        # Now grap our interaction terms
        summaryVals <- summary(tmp)
        sigTerms <- anova(tmp)$pTerms.pv[6]
        # Now combine the roi witht he p value
        sigTerms <- c(roiName, sigTerms)
        # Now export the variables
        output <- rbind(output, sigTerms)
    }
    # Now fix the column names
    return(output)
}
