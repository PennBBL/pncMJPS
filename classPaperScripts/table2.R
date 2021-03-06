## Load library(s)
install_load('psych','ggplot2','pROC','ggrepel','caret','randomForest','MatchIt','glmnet','useful')
source('figure5Functions.R')

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
mod <- matchit(usageBin ~ ageAtScan1 + envSES + dti64Tsnr, data=tmpDat, ratio=2, na.action=na.omit)
male.data.all <- male.data
male.data <- male.data[as.vector(mod$match.matrix),]
male.data <- rbind(male.data, male.data.all[which(male.data.all$usageBin==1),])
male.data.all.m <- male.data

# Now do female
tmpDat <- female.data[c('bblid', 'scanid', 'usageBin', 'ageAtScan1', 'envSES', 'dti64Tsnr')]
tmpDat <- tmpDat[complete.cases(tmpDat),]
mod <- matchit(usageBin ~ ageAtScan1 + envSES + dti64Tsnr, data=tmpDat, ratio=2, na.action=na.omit)
female.data.all <- female.data
female.data <- female.data[as.vector(mod$match.matrix),]
female.data <- rbind(female.data, female.data.all[which(female.data.all$usageBin==1),])
female.data.all.m <- female.data

# Now create a CV ridge reg model prediction stats
male.data <- male.data.all.m[complete.cases(male.data.all.m[,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]),]
foldsToLoop <- createFolds(male.data$usageBin, k=20)
cvPredVals <- rep(NA, length(male.data$usageBin))
for(q in seq(1, length(foldsToLoop))){
    index <- foldsToLoop[[q]]
    # Now do the real labels
    optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, family="binomial", parallel=F)
    lasModel2 <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, lambda=optLam$lambda.min)
    
    cvPredVals[index] <- predict(lasModel2, newx=as.matrix(male.data[index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), type='response')
    
}

# Now find the optimal cut off point
cutVal <- coords(roc(male.data$usageBin ~ cvPredVals), 'best')
# Now create the confusion matrix
male.data$usagePred <- 0
male.data$usagePred[cvPredVals<=cutVal[1]] <- 2
outTab <- table(male.data$usageBin, male.data$usagePred)
accVal <- sum(diag(outTab)) / sum(outTab)
write.csv(outTab, "cmNonVsUserMale.csv", quote=F, row.names=F)

## Now do females
female.data <- female.data.all.m[complete.cases(female.data.all.m[,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]),]
foldsToLoop <- createFolds(female.data$usageBin, 20)
cvPredVals <- rep(NA, length(female.data$usageBin))
for(q in seq(1, length(foldsToLoop))){
    index <- foldsToLoop[[q]]
    # Now do the real labels
    optLam <- cv.glmnet(y=as.vector(female.data$usageBin[-index]), x=as.matrix(female.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), alpha=0, family="binomial", parallel=F)
    lasModel2 <- glmnet(y=as.vector(female.data$usageBin[-index]), x=as.matrix(female.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), alpha=0, lambda=optLam$lambda.min)
    
    cvPredVals[index] <- predict(lasModel2, newx=as.matrix(female.data[index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), type='response')
    
}

# Now find the optimal cut off point
cutVal <- coords(roc(female.data$usageBin ~ cvPredVals), 'best')
# Now create the confusion matrix
female.data$usagePred <- 0
female.data$usagePred[cvPredVals<=cutVal[1]] <- 2
outTab <- table(female.data$usageBin, female.data$usagePred)
accVal <- sum(diag(outTab)) / sum(outTab)
write.csv(outTab, "cmNonVsUserFemale.csv", quote=F, row.names=F)

## Now do the user vs frequent analysis
mjData <- read.csv("../data/n9462_mj_ps_cnb_fortmm.csv")
mjData$dosage <- NA
mjData$dosage[which(mjData$marcat=='MJ Non-User')] <- 1
mjData$dosage[which(mjData$marcat=='MJ User' & mjData$mjpastyr=='')] <- 1
mjData$dosage[which(mjData$mjpastyr=="Less than once a month")] <- 2
mjData$dosage[which(mjData$mjpastyr=="About once a month")] <- 3
mjData$dosage[which(mjData$mjpastyr=="2-3 times a month")] <- 4
mjData$dosage[which(mjData$mjpastyr=="1-2 times a week")] <- 5
mjData$dosage[which(mjData$mjpastyr=="3-4 times a week")] <- 6
mjData$dosage[which(mjData$mjpastyr=="Everyday or nearly every day")] <- 7

# Now give us all the values
all.data <- merge(img.data, mjData)
all.data <- all.data[-which(all.data$dosage==1),]
all.data <- all.data[-which(is.na(all.data$dosage)),]
all.data$usageBin <- 0
all.data$usageBin[all.data$dosage>=6] <- 1
all.data$usageBinOrig <- all.data$usageBin

# Now prepare a sex specific values
male.data <- all.data[which(all.data$sex==1),]
male.data.all.m <- male.data
female.data <- all.data[which(all.data$sex==2),]
female.data.all.m <- female.data

# Now create a CV ridge reg model prediction stats
male.data <- male.data.all.m[complete.cases(male.data.all.m[,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]),]
foldsToLoop <- createFolds(male.data$usageBin, 20)
cvPredVals <- rep(NA, length(male.data$usageBin))
for(q in seq(1, length(foldsToLoop))){
    index <- foldsToLoop[[q]]
    # Now do the real labels
    optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, family="binomial", parallel=F)
    lasModel2 <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, lambda=optLam$lambda.min)
    
    cvPredVals[index] <- predict(lasModel2, newx=as.matrix(male.data[index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), type='response')
    
}

# Now find the optimal cut off point
cutVal <- coords(roc(male.data$usageBin ~ cvPredVals), 'best')
# Now create the confusion matrix
male.data$usagePred <- 0
male.data$usagePred[cvPredVals<=cutVal[1]] <- 2
outTab <- table(male.data$usageBin, male.data$usagePred)
accVal <- sum(diag(outTab)) / sum(outTab)
write.csv(outTab, "cmUserVsFreqMale.csv", quote=F, row.names=F)

## Now do females
female.data <- female.data.all.m[complete.cases(female.data.all.m[,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]),]
foldsToLoop <- createFolds(female.data$usageBin, 20)
cvPredVals <- rep(NA, length(female.data$usageBin))
for(q in seq(1, length(foldsToLoop))){
    index <- foldsToLoop[[q]]
    # Now do the real labels
    optLam <- cv.glmnet(y=as.vector(female.data$usageBin[-index]), x=as.matrix(female.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), alpha=0, family="binomial", parallel=F)
    lasModel2 <- glmnet(y=as.vector(female.data$usageBin[-index]), x=as.matrix(female.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), alpha=0, lambda=optLam$lambda.min)
    
    cvPredVals[index] <- predict(lasModel2, newx=as.matrix(female.data[index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), type='response')
    
}

# Now find the optimal cut off point
cutVal <- coords(roc(female.data$usageBin ~ cvPredVals), 'best')
# Now create the confusion matrix
female.data$usagePred <- 0
female.data$usagePred[cvPredVals<=cutVal[1]] <- 2
outTab <- table(female.data$usageBin, female.data$usagePred)
accVal <- sum(diag(outTab)) / sum(outTab)
write.csv(outTab, "cmUserVsFreqFemale.csv", quote=F, row.names=F)
