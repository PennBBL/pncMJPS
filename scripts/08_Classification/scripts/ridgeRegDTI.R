#
#t1 <- read.csv('alffData.csv')
#t2 <- read.csv('cbfData.csv')
#t3 <- read.csv('ctData.csv')
#t4 <- read.csv('gmdData.csv')
#t5 <- read.csv('jhuFALabel.csv')
#t6 <- read.csv('jlfTRData.csv')
#t7 <- read.csv('rehoData.csv')
#t8 <- read.csv('volumeData.csv')

#allData <- merge(t1, t2, all=T)
#allData <- merge(allData, t3, all=T)
#allData <- merge(allData, t4, all=T)
#allData <- merge(allData, t5, all=T)
#allData <- merge(allData, t6, all=T)
#allData <- merge(allData, t7, all=T)
#allData <- merge(allData, t8, all=T)

## Load library(s)
install_load('psych', 'ggplot2', 'pROC', 'ggrepel', 'caret', 'randomForest', 'MatchIt', 'glmnet', 'doMC')
source('../functions/functions.R')

# Now start loading the data down here
# This should be run from the dta directory
mjData <- read.csv("../../../data/n9462_mj_ps_cnb_fortmm.csv")
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
img.data <- read.csv('../../../data/imagingDataAll.csv')
#img.data <- read.csv('../../../data/n1601_imagingDataDump_20180104.csv')

# Now give us all the values
all.data <- merge(img.data, mjData)
all.data <- all.data[-which(all.data$dosage==1),]

# Now prepare a sex specific values
male.data <- all.data[which(all.data$sex==1),]
#male.data <- male.data[which(male.data$dosage == 0 | male.data$dosage > 5),]
female.data <- all.data[which(all.data$sex==2),]

# Now add a binary matrix for use or no use
male.data$usageBin <- 0
male.data$usageBin[male.data$dosage>1] <- 1
male.data <- male.data[-which(is.na(male.data$dosage)),]

# Now match up our data
#tmpDat <- male.data[c('bblid', 'scanid', 'usageBin', 'ageAtGo1Scan', 'envSES')]
tmpDat <- male.data[c('bblid', 'scanid', 'usageBin', 'ageAtScan1', 'envSES')]
mod <- matchit(usageBin ~ ageAtScan1 + envSES, data=tmpDat, ratio=3, na.action=na.omit)
male.data.all <- male.data
male.data <- male.data[as.vector(mod$match.matrix),]
male.data <- rbind(male.data, male.data.all[which(male.data.all$usageBin==1),])
male.data.all.m <- male.data

# Now lets see how well we can build our model in a cross validated fashion
# This will be done within modality just to explore things
aucVals <- NULL
registerDoMC(4)
# tr
male.data <- male.data.all.m[complete.cases(male.data.all.m[,grep('dti_jlf_tr', names(male.data))]),]
#male.data <- male.data[,-grep('dti_jlf_tr_MeanTR', names(male.data))]
foldsToLoop <- createFolds(male.data$usageBin, table(male.data$usageBin)[2])
cvPredVals <- rep(NA, length(male.data$usageBin))
for(q in seq(1, length(foldsToLoop))){
    index <- foldsToLoop[[q]]
    #volMod <- buildStepROCModel(y=male.data$usageBin[-index], x=male.data[-index,grep('_jlf_vol_', names(male.data))], varAdd=6)
    #outModel <- as.formula(paste('usageBin~', volMod[dim(volMod)[1],2]))
    # Now build this model and
    # build a lasso model
    optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('dti_jlf_tr', names(male.data))]), alpha=0, family="binomial", parallel=T)
    print(optLam$lambda.min)
    lasModel <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('dti_jlf_tr', names(male.data))]), alpha=0, lambda=optLam$lambda.min)
    print(lasModel$beta)
    #tmpModel <- glm(outModel, data=male.data[-index,], family=binomial())
    cvPredVals[index] <- predict(lasModel, newx=as.matrix(male.data[index,grep('dti_jlf_tr_', names(male.data))]), type='response')
    
}
plot(roc(male.data$usageBin ~ cvPredVals))
pROC::auc(roc(male.data$usageBin ~ cvPredVals))
aucVals <- rbind(aucVals, c('tr', pROC::auc(roc(male.data$usageBin ~ cvPredVals))))

# Now prepare t values and roc values across these people
valsOut <- NULL
seqVals <- grep('dti_jlf_tr', names(male.data))
for(q in seqVals){
  tVal <-t.test(male.data[,q]~usageBin, data=male.data)
  rocVal <- roc(usageBin~male.data[,q], data=male.data)
  outputRow <- c(colnames(male.data)[q], as.numeric(rocVal$auc), tVal$statistic)
  valsOut <- rbind(valsOut, outputRow)
}

# Now create a confusion matrix at our "best" cut off value
cutVal <- coords(roc(male.data$usageBin ~ cvPredVals), 'best')
# Now grab our confusion matrix
male.data$usagePred <- 0
male.data$usagePred[cvPredVals<=cutVal[1]] <- 2
table(male.data$usageBin, male.data$usagePred)

# Now find our true log odds from a modality regressed data set
mod.reg.ds <- male.data.all.m[complete.cases(male.data.all.m[,grep('dti_jlf_tr', names(male.data))]),]
mod.reg.ds <- mod.reg.ds[,-grep('dti_jlf_tr_MeanTR', names(mod.reg.ds))]
mod.reg.ds <- regressWithinModality(mod.reg.ds, 'dti_jlf_tr')
colVals <- grep('dti_jlf_tr', names(mod.reg.ds))
mod.reg.ds <- mod.reg.ds[,c(1027, colVals)]
mod.reg.ds[,2:length(colVals)] <- as.matrix(scale(mod.reg.ds[,2:length(colVals)]))
# Now create the model
outMod <- glm(usageBin~., data=mod.reg.ds)
toWrite <- summary(outMod)
write.csv(toWrite$coefficients, 'maleTRCoefValuesWM.csv', quote=F, row.names=T)

# Now do the sme but exclude TBV and WM values
mod.reg.ds <- male.data.all.m[complete.cases(male.data.all.m[,grep('dti_jlf_tr', names(male.data))]),]
mod.reg.ds <- mod.reg.ds[,-grep('dti_jlf_tr_MeanTR', names(mod.reg.ds))]
mod.reg.ds <- mod.reg.ds[,-grep('_Lobe_WM', names(mod.reg.ds))]
mod.reg.ds <- regressWithinModality(mod.reg.ds, 'dti_jlf_tr')
colVals <- grep('dti_jlf_tr', names(mod.reg.ds))
mod.reg.ds <- mod.reg.ds[,c(991, colVals)]
mod.reg.ds[,2:length(colVals)] <- as.matrix(scale(mod.reg.ds[,2:length(colVals)]))
# Now create the model
outMod <- glm(usageBin~., data=mod.reg.ds)
toWrite <- summary(outMod)
write.csv(toWrite$coefficients, 'maleTRCoefValues.csv', quote=F, row.names=T)

# FA
male.data <- male.data.all.m[complete.cases(male.data.all.m[,grep('dti_dtitk_jhulabel_fa', names(male.data))]),]
foldsToLoop <- createFolds(male.data$usageBin, table(male.data$usageBin)[2])
cvPredVals <- rep(NA, length(male.data$usageBin))
for(q in seq(1, length(foldsToLoop))){
    index <- foldsToLoop[[q]]
    #volMod <- buildStepROCModel(y=male.data$usageBin[-index], x=male.data[-index,grep('_jlf_vol_', names(male.data))], varAdd=6)
    #outModel <- as.formula(paste('usageBin~', volMod[dim(volMod)[1],2]))
    # Now build this model and
    # build a lasso model
    optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('dti_dtitk_jhulabel_fa', names(male.data))]), alpha=0, family="binomial", parallel=T)
    print(optLam$lambda.min)
    lasModel <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('dti_dtitk_jhulabel_fa', names(male.data))]), alpha=0, lambda=optLam$lambda.min)
    #tmpModel <- glm(outModel, data=male.data[-index,], family=binomial())
    cvPredVals[index] <- predict(lasModel, newx=as.matrix(male.data[index,grep('dti_dtitk_jhulabel_fa', names(male.data))]), type='response')
    
}
plot(roc(male.data$usageBin ~ cvPredVals))
pROC::auc(roc(male.data$usageBin ~ cvPredVals))
aucVals <- rbind(aucVals, c('fa', pROC::auc(roc(male.data$usageBin ~ cvPredVals))))

# Now do the fa labels
seqVals <- grep('_dtitk_jhulabel', names(male.data))
for(q in seqVals){
    tVal <-t.test(male.data[,q]~usageBin, data=male.data)
    rocVal <- roc(usageBin~male.data[,q], data=male.data)
    outputRow <- c(colnames(male.data)[q], as.numeric(rocVal$auc), tVal$statistic)
    valsOut <- rbind(valsOut, outputRow)
}

# Now write the output
write.csv(valsOut, 'tValsandROCValsnonVsUser.csv', quote=F, row.names=F)

aucVals <- as.data.frame(aucVals)
aucVals$V2 <- as.numeric(as.character(aucVals$V2))
# Now plot the auc Values
aucPlot <- ggplot(aucVals, aes(x=V1, y=as.numeric(as.character(V2)))) +
 geom_col()
pdf('nonUserVsInFreqUser.pdf')
print(aucPlot)
dev.off()

# Now write the color maps and all of that good stuff
writeColorTableandKey(inputData=valsOut,inputColumn=2,outName='allValsA',minTmp=c(-1,0),maxTmp=c(.45,.8))
writeColorTableandKey(inputData=valsOut,inputColumn=3,outName='allValsT',minTmp=c(-3,0),maxTmp=c(0,3))

# Now do the same thing but for our in freq vs freq smokers
mjData <- read.csv("../../../data/n9462_mj_ps_cnb_fortmm.csv")
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

# Now prepare a sex specific values
male.data <- all.data[which(all.data$sex==2),]
#male.data <- male.data[which(male.data$dosage == 0 | male.data$dosage > 5),]
female.data <- all.data[which(all.data$sex==1),]

# Now add a binary matrix for use or no use
male.data$usageBin <- 0
male.data$usageBin[male.data$dosage>=6] <- 1
male.data.all.m <- male.data

# Now see about classification
aucVals <- NULL
registerDoMC(4)
# tr
male.data <- male.data.all.m[complete.cases(male.data.all.m[,grep('dti_jlf_tr', names(male.data))]),]
male.data <- male.data[,-grep('dti_jlf_tr_MeanTR', names(male.data))]
foldsToLoop <- createFolds(male.data$usageBin, table(male.data$usageBin)[2])
cvPredVals <- rep(NA, length(male.data$usageBin))
for(q in seq(1, length(foldsToLoop))){
    index <- foldsToLoop[[q]]
    #volMod <- buildStepROCModel(y=male.data$usageBin[-index], x=male.data[-index,grep('_jlf_vol_', names(male.data))], varAdd=6)
    #outModel <- as.formula(paste('usageBin~', volMod[dim(volMod)[1],2]))
    # Now build this model and
    # build a lasso model
    optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('dti_jlf_tr', names(male.data))]), alpha=0, family="binomial", parallel=T)
    print(optLam$lambda.min)
    lasModel <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('dti_jlf_tr', names(male.data))]), alpha=0, lambda=optLam$lambda.min)
    print(lasModel$beta)
    #tmpModel <- glm(outModel, data=male.data[-index,], family=binomial())
    cvPredVals[index] <- predict(lasModel, newx=as.matrix(male.data[index,grep('dti_jlf_tr_', names(male.data))]), type='response')
    
}
plot(roc(male.data$usageBin ~ cvPredVals))
pROC::auc(roc(male.data$usageBin ~ cvPredVals))
aucVals <- rbind(aucVals, c('tr', pROC::auc(roc(male.data$usageBin ~ cvPredVals))))

# Now prepare t values and roc values across these people
valsOut <- NULL
seqVals <- grep('dti_jlf_tr', names(male.data))
for(q in seqVals){
    tVal <-t.test(male.data[,q]~usageBin, data=male.data)
    rocVal <- roc(usageBin~male.data[,q], data=male.data)
    outputRow <- c(colnames(male.data)[q], as.numeric(rocVal$auc), tVal$statistic)
    valsOut <- rbind(valsOut, outputRow)
}

# FA
male.data <- male.data.all.m[complete.cases(male.data.all.m[,grep('dti_dtitk_jhulabel_fa', names(male.data))]),]
foldsToLoop <- createFolds(male.data$usageBin, table(male.data$usageBin)[2])
cvPredVals <- rep(NA, length(male.data$usageBin))
for(q in seq(1, length(foldsToLoop))){
    index <- foldsToLoop[[q]]
    #volMod <- buildStepROCModel(y=male.data$usageBin[-index], x=male.data[-index,grep('_jlf_vol_', names(male.data))], varAdd=6)
    #outModel <- as.formula(paste('usageBin~', volMod[dim(volMod)[1],2]))
    # Now build this model and
    # build a lasso model
    optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('dti_dtitk_jhulabel_fa', names(male.data))]), alpha=0, family="binomial", parallel=T)
    print(optLam$lambda.min)
    lasModel <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('dti_dtitk_jhulabel_fa', names(male.data))]), alpha=0, lambda=optLam$lambda.min)
    #tmpModel <- glm(outModel, data=male.data[-index,], family=binomial())
    cvPredVals[index] <- predict(lasModel, newx=as.matrix(male.data[index,grep('dti_dtitk_jhulabel_fa', names(male.data))]), type='response')
    
}
plot(roc(male.data$usageBin ~ cvPredVals))
pROC::auc(roc(male.data$usageBin ~ cvPredVals))
aucVals <- rbind(aucVals, c('fa', pROC::auc(roc(male.data$usageBin ~ cvPredVals))))

# Now do the fa labels
seqVals <- grep('_dtitk_jhulabel', names(male.data))
for(q in seqVals){
    tVal <-t.test(male.data[,q]~usageBin, data=male.data)
    rocVal <- roc(usageBin~male.data[,q], data=male.data)
    outputRow <- c(colnames(male.data)[q], as.numeric(rocVal$auc), tVal$statistic)
    valsOut <- rbind(valsOut, outputRow)
}

# Now write the output
write.csv(valsOut, 'tValsandROCValsuserVsFreqUser.csv', quote=F, row.names=F)

aucVals <- as.data.frame(aucVals)
aucVals$V2 <- as.numeric(as.character(aucVals$V2))
# Now plot the auc Values
aucPlot <- ggplot(aucVals, aes(x=V1, y=as.numeric(as.character(V2)))) +
geom_col()

pdf('userVsFreqUser.pdf')
print(aucPlot)
dev.off()

# Now write the color maps and all of that good stuff
writeColorTableandKey(inputData=valsOut,inputColumn=2,outName='allValsUvFA',minTmp=c(-1,0),maxTmp=c(.45,.8))
writeColorTableandKey(inputData=valsOut,inputColumn=3,outName='allValsUvFT',minTmp=c(-3,0),maxTmp=c(0,3))

