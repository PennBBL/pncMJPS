## Load library(s)
install_load('psych','ggplot2','pROC','ggrepel','caret','MatchIt','glmnet','foreach','useful','doParallel','utils')
source('figure6Functions.R')

## Now load the data
# Now start loading the data down here
# This should be run from the dta directory
mjData <- read.csv("../data/n9462_mj_ps_cnb_fortmm.csv")
mjData$dosage <- NA
mjData$dosage[which(mjData$marcat=='MJ Non-User')] <- 0
mjData$dosage[which(mjData$marcat=='MJ User' & mjData$mjpastyr=='')] <- 1
mjData$dosage[which(mjData$mjpastyr=="Less than once a month")] <- 2
mjData$dosage[which(mjData$mjpastyr=="About once a month")] <- 3
mjData$dosage[which(mjData$mjpastyr=="2-3 times a month")] <- 4
mjData$dosage[which(mjData$mjpastyr=="1-2 times a week")] <- 5
mjData$dosage[which(mjData$mjpastyr=="3-4 times a week")] <- 1
mjData$dosage[which(mjData$mjpastyr=="Everyday or nearly every day")] <- 1

# Now load imaging data
img.data <- read.csv('../data/imagingDataAll.csv')

# Now prepare all of the data
all.data <- merge(img.data, mjData)
all.data <- all.data[-which(all.data$dosage==1),]
all.data$usageBin <- 0
all.data$usageBin[all.data$dosage>1] <- 1
all.data <- all.data[-which(is.na(all.data$dosage)),]

# Now isolate genders
male.data <- all.data[which(all.data$sex==1),]
female.data <- all.data[which(all.data$sex==2),]

# Now create our age matched samples
# Starting with male
tmpDat <- male.data[c('bblid', 'scanid', 'usageBin', 'ageAtScan1', 'envSES', 'dti64Tsnr')]
tmpDat <- tmpDat[complete.cases(tmpDat),]
mod <- matchit(usageBin ~ ageAtScan1 + envSES, data=tmpDat, ratio=3, na.action=na.omit)
male.data.all <- male.data
male.data <- male.data[as.vector(mod$match.matrix),]
male.data <- rbind(male.data, male.data.all[which(male.data.all$usageBin==1),])
male.data.all.m <- male.data
propValueMale <- table(male.data$usageBin)[2]/sum(table(male.data$usageBin))
male.data.all.m$usageBinOrig <- male.data.all.m$usageBin
output <- male.data.all.m[,c('bblid', 'scanid')]
write.csv(output, "maleIDValues.csv", quote=F, row.names=F)

# Now do female
tmpDat <- female.data[c('bblid', 'scanid', 'usageBin', 'ageAtScan1', 'envSES', 'dti64Tsnr')]
tmpDat <- tmpDat[complete.cases(tmpDat),]
mod <- matchit(usageBin ~ ageAtScan1 + envSES, data=tmpDat, ratio=1, na.action=na.omit)
female.data.all <- female.data
female.data <- female.data[as.vector(mod$match.matrix),]
female.data <- rbind(female.data, female.data.all[which(female.data.all$usageBin==1),])
female.data.all.m <- female.data
propValueFemale <- table(female.data$usageBin)[2]/sum(table(female.data$usageBin))
female.data.all.m$usageBinOrig <- female.data.all.m$usageBin
output <- female.data.all.m[,c('bblid', 'scanid')]
write.csv(output, "femaleIDValues.csv", quote=F, row.names=F)

## Perform variable selection in each fold
male.data <- male.data.all.m[complete.cases(male.data.all.m[,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]),]
foldsToLoop <- createFolds(male.data$usageBin, 25)
cvPredVals <- rep(NA, length(male.data$usageBin))
selectVals <- list()
# Run variable selection in males
pb <- txtProgressBar(min=0, max=length(foldsToLoop), initial=0, style=3)
for(q in seq(1, length(foldsToLoop))){
    index <- foldsToLoop[[q]]
    #tmp <- runLassoforHiLo(x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), y=as.vector(male.data$usageBin[-index]), nofFolds=20, trainingIterations=10, nCor=8, alphaSequence=seq(0,1,.1))
    # Now get the selection index
    #colsToUse <- which(returnSelectionCol(tmp)>1)
    #selectVals[[q]] <- colsToUse
    inputX <- as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))])
    #inputX <- as.matrix(inputX[,colsToUse])
    outputX <- as.matrix(male.data[index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))])
    #outputX <- as.matrix(outputX[,colsToUse])
    # Now build a ridge regression model
    optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=inputX, alpha=0, parallel=F)
    lasModel1 <- glmnet(y=as.vector(male.data$usageBin[-index]), x=inputX, alpha=0, lambda=optLam$lambda.min)
    # Finally predict in the testing subjects
    cvPredVals[index] <- predict(lasModel1, newx=outputX, type='response')
    setTxtProgressBar(pb, q)
}
roc(male.data$usageBin ~ cvPredVals)
# Now build an ROC curve with these values
cvPredVals <- rep(NA, length(male.data$usageBin))
for(q in seq(1, length(foldsToLoop))){
    index <- foldsToLoop[[q]]
    
    # Now do the real labels
    optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, family="binomial", parallel=F)
    lasModel2 <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, lambda=optLam$lambda.min)
    
    cvPredVals[index] <- predict(lasModel2, newx=as.matrix(male.data[index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), type='response')
    
}
