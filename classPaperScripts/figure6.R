## Load library(s)
install_load('psych','ggplot2','pROC','ggrepel','caret','MatchIt','glmnet','foreach','useful','doParallel','utils','reshape2')
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

# Now create our age matched samples
# Starting with male
tmpDat <- male.data[c('bblid', 'scanid', 'usageBin', 'ageAtScan1', 'envSES', 'dti64Tsnr')]
tmpDat <- tmpDat[complete.cases(tmpDat),]
mod <- matchit(usageBin ~ ageAtScan1 + envSES , data=tmpDat, ratio=3, na.action=na.omit)
male.data.all <- male.data
male.data <- male.data[as.vector(mod$match.matrix),]
male.data <- rbind(male.data, male.data.all[which(male.data.all$usageBin==1),])
male.data.all.m <- male.data
output <- male.data.all.m[,c('bblid', 'scanid')]
write.csv(output, "maleIDValues.csv", quote=F, row.names=F)

## Perform variable selection in each fold
male.data <- male.data.all.m[complete.cases(male.data.all.m[,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]),]
foldsToLoop <- createFolds(male.data$usageBin, 10)
# Now prepare an output matrix
variance.explore <- matrix(NA,48,length(foldsToLoop))
# And now prepare the values for t tests
tval.explore <- matrix(NA, 48, length(foldsToLoop))
seqVals <- grep('_dtitk_jhulabel', names(male.data))
# Run variable selection in males
for(q in seq(1, length(foldsToLoop))){
  index <- foldsToLoop[[q]]
  tmp <- runLassoforHiLo(x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), y=as.vector(male.data$usageBin)[-index], nofFolds=10, trainingIterations=100, nCor=8, alphaSequence=c(.5, 1))
  countVals <- returnSelectionCol(tmp)
  variance.explore[,q] <- countVals
  # Now produce a t statistic for these guys
  for(w in 1:length(seqVals)){
      z <- seqVals[w]
      tVal <-t.test(male.data[-index,z]~usageBin[-index], data=male.data)
      tVal <- tVal$statistic
      tval.explore[w,q] <- tVal
  }
  # Now explore variable selection using the selection count
}
rownames(variance.explore) <- names(countVals)
toPlot <- melt(variance.explore)
rownames(tval.explore) <- names(countVals)
toPlot2 <- melt(tval.explore)
# Now crate a violin plot
viol.plot.1 <- ggplot(toPlot, aes(x=Var1, y=value)) +
  geom_violin() +
  stat_summary(fun.y=mean, geom="point", shape=23) +
  stat_summary(fun.y=median, geom="point", size=2, color="red") +
  theme(axis.text.x=element_text(angle=90))
viol.plot.2 <- ggplot(toPlot2, aes(x=Var1, y=value)) +
  geom_violin() +
  stat_summary(fun.y=mean, geom="point", shape=23) +
  stat_summary(fun.y=median, geom="point", size=2, color="red") +
  theme(axis.text.x=element_text(angle=90))

# Now explore how much additional AUC each variable adds in a cv fashion
loopVals <- names(sort(rank(apply(variance.explore, 1, mean), ties.method='random'), decreasing=T))
outputAucVals <- NULL
for(i in 2:48){
    inputVals <- loopVals[1:i]
    cvPredVals <- rep(NA, length(male.data$usageBin))
    for(q in seq(1, length(foldsToLoop))){
        index <- foldsToLoop[[q]]
        # build a ridge model with the fake data
        optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,c(inputVals)]), alpha=0, family="binomial", parallel=F)
        lasModel1 <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,c(inputVals)]), alpha=0, lambda=optLam$lambda.min, family="binomial")
        # Now output the predicted values
        cvPredVals[index] <- predict(lasModel1, newx=as.matrix(male.data[index,c(inputVals)]), type='response')
    }
    tmpAUC <- pROC::auc(roc(male.data$usageBin~cvPredVals))
    print(tmpAUC)
    print(plot(roc(male.data$usageBin~cvPredVals)))
    outputAucValsRow <- c(i,tmpAUC)
    outputAucVals <- rbind(outputAucVals, outputAucValsRow)
    rm(cvPredVals)
}
