runGamModel <- function(dataFrame, grepPattern, qualityIndex){
  # First find all fo the columns that we will loop through
  columnIndex <- grep(grepPattern, names(dataFrame))
  qualityCol <- grep(qualityIndex, names(dataFrame))[1]
  # Now declare our statics
  ageCol <- which(names(dataFrame)=='ageAtScan1')
  sexCol <- which(names(dataFrame)=='sex')
  mjCol <- which(names(dataFrame)=='dosage')
  psCol <- which(names(dataFrame)=='goassessDxpmr7')
  
  # Declare som values for which will be returned from the gam function
  output <- NULL
  # Now loop through each region and run the model 
  for(i in columnIndex){
    # Get the roi name
    roiName <- names(dataFrame)[i]
    # First run the model 
    tmp <- gam(dataFrame[,i] ~ s(dataFrame[,ageCol]) + dataFrame[,sexCol] + dataFrame[,mjCol]+dataFrame[,psCol] + dataFrame[,qualityCol])
    # Now grap our interaction terms
    summaryVals <- summary(tmp)
    sigTerms <- summaryVals$p.table[3,4]
    # Now export the variables
    allVals <- cbind(roiName, sigTerms, summaryVals$p.table[3,3])
    output <- rbind(output, allVals)
  }

  return(output)
}

runGamModelS <- function(dataFrame, grepPattern, qualityIndex){
  # First find all fo the columns that we will loop through
  columnIndex <- grep(grepPattern, names(dataFrame))
  qualityCol <- grep(qualityIndex, names(dataFrame))[1]
  # Now declare our statics
  ageCol <- which(names(dataFrame)=='ageAtScan1')
  sexCol <- which(names(dataFrame)=='sex')
  mjCol <- which(names(dataFrame)=='dosage')
  psCol <- which(names(dataFrame)=='goassessDxpmr7')
  
  # Declare som values for which will be returned from the gam function
  output <- NULL
  # Now loop through each region and run the model 
  for(i in columnIndex){
    # Get the roi name
    roiName <- names(dataFrame)[i]
    # First run the model 
    tmp <- gam(dataFrame[,i] ~ s(dataFrame[,ageCol]) + dataFrame[,sexCol] + s(dataFrame[,mjCol],k=4)+dataFrame[,psCol] + dataFrame[,qualityCol])
    # Now grap our interaction terms
    summaryVals <- summary(tmp)
    sigTerms <- summaryVals$p.table[3,4]
    # Now export the variables
    allVals <- cbind(roiName, sigTerms, summaryVals$p.table[3,3])
    output <- rbind(output, allVals)
  }

  return(output)
}

runGamModelG <- function(dataFrame, grepPattern, qualityIndex){
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
    tmp <- gam(dataFrame[,i]~s(dataFrame[,ageCol])+dataFrame[,sexCol]+dataFrame[,mjCol]+dataFrame[,psCol]+dataFrame[,qualityCol])
    # Now grap our interaction terms
    summaryVals <- summary(tmp)
    sigTerms <- summaryVals$p.table[3,4]
    # Now export the variables
    allVals <- cbind(roiName, sigTerms, summaryVals$p.table[3,3])
    output <- rbind(output, allVals)
  }

  return(output)
}


runGamModelI <- function(dataFrame, grepPattern, qualityIndex){
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
    tmp <- gam(dataFrame[,i]~dataFrame[,sexCol]+dataFrame[,ageCol]*dataFrame[,mjCol]+dataFrame[,psCol]+dataFrame[,qualityCol])
    tmp <- anova(tmp)
    # Now grap our interaction terms
    summaryVals <- summary(tmp)
    sigTerms <- tmp$pTerms.table[6,3]
    # Now export the variables
    allVals <- cbind(roiName, sigTerms)
    output <- rbind(output, allVals)
  }

  return(output)
}
