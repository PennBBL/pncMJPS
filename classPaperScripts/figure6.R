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
mod <- matchit(usageBin ~ ageAtScan1 + envSES + dti64Tsnr, data=tmpDat, ratio=3, na.action=na.omit)
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
mod <- matchit(usageBin ~ ageAtScan1 + envSES + dti64Tsnr, data=tmpDat, ratio=3, na.action=na.omit)
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
foldsToLoop <- createFolds(male.data$usageBin, 20)
# Now prepare an output matrix
variance.explore <- matrix(NA,48,length(foldsToLoop))
# Run variable selection in males
for(q in seq(1, length(foldsToLoop))){
  index <- foldsToLoop[[q]]
  tmp <- runLassoforHiLo(x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), y=as.vector(male.data$usageBin)[-index], nofFolds=20, trainingIterations=100, nCor=8, alphaSequence=c(0.05, .5, 1))
  countVals <- returnSelectionCol(tmp)
  variance.explore[,q] <- countVals
}
rownames(variance.explore) <- names(countVals)
toPlot <- melt(variance.explore)
# Now crate a violin plot
viol.plot.1 <- ggplot(toPlot, aes(x=Var1, y=value)) +
  geom_violin() +
  stat_summary(fun.y=mean, geom="point", shape=23) +
  stat_summary(fun.y=median, geom="point", size=2, color="red") +
  theme(axis.text.x=element_text(angle=90))

# Now get a t test between the groups
valsOut <- NULL
seqVals <- grep('_dtitk_jhulabel', names(male.data))
for(q in seqVals){
    tVal <-t.test(male.data[,q]~usageBin, data=male.data)
    rocVal <- roc(usageBin~male.data[,q], data=male.data)
    outputRow <- c(colnames(male.data)[q], as.numeric(rocVal$auc), tVal$statistic)
    valsOut <- rbind(valsOut, outputRow)
}
rownames(valsOut) <- NULL
# Now combine these metrics
allOut <- as.data.frame(valsOut)
countVals <- cbind(names(returnSelectionCol(tmp)), returnSelectionCol(tmp))
rownames(countVals) <- NULL
countVals <- as.data.frame(countVals)
allOut <- merge(allOut, countVals, by='V1')
allOut[,2:4] <- apply(allOut[,2:4], 2, function(x) as.numeric(as.character(x)))

# Now create a scatter plot for these values
outputPlot <- ggplot(allOut, aes(x=t, y=V2.y, ))
