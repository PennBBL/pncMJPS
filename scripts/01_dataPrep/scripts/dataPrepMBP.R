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

# Now start loading the data down here
# This should be run from the dta directory
mjData <- read.csv("n9462_mj_ps_cnb_fortmm.csv")
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
img.data <- read.csv('./imagingDataAll.csv')

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

# Now match up our data
tmpDat <- male.data[c('bblid', 'scanid', 'usageBin', 'ageAtGo1Scan', 'envSES')]
mod <- matchit(usageBin ~ ageAtGo1Scan + envSES, data=tmpDat, ratio=3, na.action=na.omit)
male.data.all <- male.data
male.data <- male.data[as.vector(mod$match.matrix),]
male.data <- rbind(male.data, male.data.all[which(male.data.all$usageBin==1),])
male.data.all.m <- male.data


# Now run through every colume and give us a t value for every JLF region
tValsMaleOut <- NULL
colVals <- grep('_jlf_', names(male.data))
colVals <- append(colVals, grep('dtitk_jhulabel', names(male.data)))
for(z in colVals){
  tVal <-t.test(male.data[,z]~usageBin, data=male.data)
  outputRow <- c(colnames(male.data)[z], as.numeric(tVal$statistic), as.numeric(tVal$p.value))
  tValsMaleOut <- rbind(tValsMaleOut, outputRow)
}
tValsMaleOut <- tValsMaleOut[order(tValsMaleOut[,3]),]
tValsMaleSig <- tValsMaleOut[which(tValsMaleOut[,3]<.05),]
rownames(tValsMaleSig) <- NULL
tValsMaleSig <- as.data.frame(tValsMaleSig)

# Now do the same with ROC and find the corellation between roc and
rocValsMaleOut <- NULL
for(z in colVals){
  rocVal <- roc(usageBin~male.data[,z], data=male.data)
  outputRow <- c(colnames(male.data)[z], as.numeric(rocVal$auc))
  rocValsMaleOut <- rbind(rocValsMaleOut, outputRow)
}
rownames(rocValsMaleOut) <- NULL
rocValsMaleOut <- as.data.frame(rocValsMaleOut)

foobar <- merge(tValsMaleSig, rocValsMaleOut, by='V1')
foobar$V2.x <- as.numeric(as.character(foobar$V2.x))
foobar$V3 <- as.numeric(as.character(foobar$V3))
foobar$V2.y <- as.numeric(as.character(foobar$V2.y))
# Now seperate the positive and negative values
foobarPos <- foobar[which(foobar$V2.x>0),]
foobarNeg <- foobar[which(foobar$V2.x<0),]

# Now plot these values in a scatter plot
outPlot <- ggplot(foobarNeg, aes(x=V2.y, y=V2.x)) +
  geom_point() +
  geom_smooth(method=lm) +
  geom_label_repel(aes(label=V1,size=2)) +#,box.padding=unit(1,"lines"),point.padding=unit(1,"lines")) +
  xlab("AUC") +
  ylab("t value") +
  #geom_hline(yintercept = 0 , linetype=3) +
  #geom_vline(xintercept = 0 , linetype=3) +
  theme(legend.position="none") +
  geom_text(aes(x=-Inf, y=Inf, hjust=0, vjust=1, label=cor(foobarPos$V2.y,foobarPos$V2.x)))

# Now prepare the t values here
tValsMaleOut <- tValsMaleOut[order(tValsMaleOut[,3]),]

# Now prepare a bar graph with these values
plotData <- tValsMaleOut[which(tValsMaleOut[,3]<.05),]


# Now lets see how well we can build our model in a cross validated fashion
# This will be done within modality just to explore things
aucVals <- NULL
male.data <- male.data.all.m[complete.cases(male.data.all.m[,grep('_jlf_vol_', names(male.data))]),]
foldsToLoop <- createFolds(male.data$usageBin, 56)
cvPredVals <- rep(NA, length(male.data$usageBin))
registerDoMC(4)
for(q in seq(1, length(foldsToLoop))){
  index <- foldsToLoop[[q]]
  #volMod <- buildStepROCModel(y=male.data$usageBin[-index], x=male.data[-index,grep('_jlf_vol_', names(male.data))], varAdd=6)
  #outModel <- as.formula(paste('usageBin~', volMod[dim(volMod)[1],2]))
  # Now build this model and
  # build a lasso model
  optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('_jlf_vol_', names(male.data))]), alpha=.55, family="binomial", parallel=T)
  print(optLam$lambda.min)
  lasModel <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('_jlf_vol_', names(male.data))]), alpha=.55, lambda=optLam$lambda.min)
  print(lasModel$beta)
  #tmpModel <- glm(outModel, data=male.data[-index,], family=binomial())
  print(q)
  cvPredVals[index] <- predict(lasModel, newx=as.matrix(male.data[index,grep('_jlf_vol_', names(male.data))]), type='response')
  
}
plot(roc(male.data$usageBin ~ cvPredVals))
aucVals <- rbind(aucVals, c('vol', pROC::auc(roc(male.data$usageBin ~ cvPredVals))))
cvPredValsVol <- cvPredVals

# Now try cbf
male.data <- male.data.all.m[complete.cases(male.data.all.m[,grep('pcasl_jlf_cbf', names(male.data))]),]
foldsToLoop <- createFolds(male.data$usageBin, table(male.data$usageBin)[2])
cvPredVals <- rep(NA, length(male.data$usageBin))
for(q in seq(1, length(foldsToLoop))){
    index <- foldsToLoop[[q]]
    #volMod <- buildStepROCModel(y=male.data$usageBin[-index], x=male.data[-index,grep('_jlf_vol_', names(male.data))], varAdd=6)
    #outModel <- as.formula(paste('usageBin~', volMod[dim(volMod)[1],2]))
    # Now build this model and
    # build a lasso model
    optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('pcasl_jlf_cbf', names(male.data))]), alpha=.55, family="binomial", parallel=T)
    print(optLam$lambda.min)
    lasModel <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('pcasl_jlf_cbf', names(male.data))]), alpha=.55, lambda=optLam$lambda.min)
    #tmpModel <- glm(outModel, data=male.data[-index,], family=binomial())
    cvPredVals[index] <- predict(lasModel, newx=as.matrix(male.data[index,grep('pcasl_jlf_cbf_', names(male.data))]), type='response')
    
}
plot(roc(male.data$usageBin ~ cvPredVals))
pROC::auc(roc(male.data$usageBin ~ cvPredVals))
aucVals <- rbind(aucVals, c('cbf', pROC::auc(roc(male.data$usageBin ~ cvPredVals))))

# tr
male.data <- male.data.all.m[complete.cases(male.data.all.m[,grep('dti_jlf_tr', names(male.data))]),]
foldsToLoop <- createFolds(male.data$usageBin, table(male.data$usageBin)[2])
cvPredVals <- rep(NA, length(male.data$usageBin))
for(q in seq(1, length(foldsToLoop))){
    index <- foldsToLoop[[q]]
    #volMod <- buildStepROCModel(y=male.data$usageBin[-index], x=male.data[-index,grep('_jlf_vol_', names(male.data))], varAdd=6)
    #outModel <- as.formula(paste('usageBin~', volMod[dim(volMod)[1],2]))
    # Now build this model and
    # build a lasso model
    optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('dti_jlf_tr', names(male.data))]), alpha=.55, family="binomial", parallel=T)
    print(optLam$lambda.min)
    lasModel <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('dti_jlf_tr', names(male.data))]), alpha=.55, lambda=optLam$lambda.min)
    print(lasModel$beta)
    #tmpModel <- glm(outModel, data=male.data[-index,], family=binomial())
    cvPredVals[index] <- predict(lasModel, newx=as.matrix(male.data[index,grep('dti_jlf_tr_', names(male.data))]), type='response')
    
}
plot(roc(male.data$usageBin ~ cvPredVals))
pROC::auc(roc(male.data$usageBin ~ cvPredVals))
aucVals <- rbind(aucVals, c('tr', pROC::auc(roc(male.data$usageBin ~ cvPredVals))))

# gmd
male.data <- male.data.all.m[complete.cases(male.data.all.m[,grep('mprage_jlf_gmd', names(male.data))]),]
foldsToLoop <- createFolds(male.data$usageBin, table(male.data$usageBin)[2])
cvPredVals <- rep(NA, length(male.data$usageBin))
for(q in seq(1, length(foldsToLoop))){
    index <- foldsToLoop[[q]]
    #volMod <- buildStepROCModel(y=male.data$usageBin[-index], x=male.data[-index,grep('_jlf_vol_', names(male.data))], varAdd=6)
    #outModel <- as.formula(paste('usageBin~', volMod[dim(volMod)[1],2]))
    # Now build this model and
    # build a lasso model
    optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('mprage_jlf_gmd', names(male.data))]), alpha=.55, family="binomial", parallel=T)
    print(optLam$lambda.min)
    lasModel <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('mprage_jlf_gmd', names(male.data))]), alpha=.55, lambda=optLam$lambda.min)
    print(lasModel$beta)
    #tmpModel <- glm(outModel, data=male.data[-index,], family=binomial())
    cvPredVals[index] <- predict(lasModel, newx=as.matrix(male.data[index,grep('mprage_jlf_gmd_', names(male.data))]), type='response')
    
}
plot(roc(male.data$usageBin ~ cvPredVals))
pROC::auc(roc(male.data$usageBin ~ cvPredVals))
aucVals <- rbind(aucVals, c('gmd', pROC::auc(roc(male.data$usageBin ~ cvPredVals))))

# ct
male.data <- male.data.all.m[complete.cases(male.data.all.m[,grep('mprage_jlf_ct', names(male.data))]),]
foldsToLoop <- createFolds(male.data$usageBin, table(male.data$usageBin)[2])
cvPredVals <- rep(NA, length(male.data$usageBin))
for(q in seq(1, length(foldsToLoop))){
    index <- foldsToLoop[[q]]
    #volMod <- buildStepROCModel(y=male.data$usageBin[-index], x=male.data[-index,grep('_jlf_vol_', names(male.data))], varAdd=6)
    #outModel <- as.formula(paste('usageBin~', volMod[dim(volMod)[1],2]))
    # Now build this model and
    # build a lasso model
    optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('mprage_jlf_ct', names(male.data))]), alpha=.55, family="binomial", parallel=T)
    print(optLam$lambda.min)
    lasModel <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('mprage_jlf_ct', names(male.data))]), alpha=.55, lambda=optLam$lambda.min)
    print(lasModel$beta)
    #tmpModel <- glm(outModel, data=male.data[-index,], family=binomial())
    cvPredVals[index] <- predict(lasModel, newx=as.matrix(male.data[index,grep('mprage_jlf_ct_', names(male.data))]), type='response')
    
}
plot(roc(male.data$usageBin ~ cvPredVals))
pROC::auc(roc(male.data$usageBin ~ cvPredVals))
aucVals <- rbind(aucVals, c('ct', pROC::auc(roc(male.data$usageBin ~ cvPredVals))))

# reho
male.data <- male.data.all.m[complete.cases(male.data.all.m[,grep('rest_jlf_reho', names(male.data))]),]
foldsToLoop <- createFolds(male.data$usageBin, table(male.data$usageBin)[2])
cvPredVals <- rep(NA, length(male.data$usageBin))
for(q in seq(1, length(foldsToLoop))){
    index <- foldsToLoop[[q]]
    #volMod <- buildStepROCModel(y=male.data$usageBin[-index], x=male.data[-index,grep('_jlf_vol_', names(male.data))], varAdd=6)
    #outModel <- as.formula(paste('usageBin~', volMod[dim(volMod)[1],2]))
    # Now build this model and
    # build a lasso model
    optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('rest_jlf_reho', names(male.data))]), alpha=.55, family="binomial", parallel=T)
    print(optLam$lambda.min)
    lasModel <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('rest_jlf_reho', names(male.data))]), alpha=.55, lambda=optLam$lambda.min)
    #tmpModel <- glm(outModel, data=male.data[-index,], family=binomial())
    cvPredVals[index] <- predict(lasModel, newx=as.matrix(male.data[index,grep('rest_jlf_reho_', names(male.data))]), type='response')
    
}
plot(roc(male.data$usageBin ~ cvPredVals))
pROC::auc(roc(male.data$usageBin ~ cvPredVals))
aucVals <- rbind(aucVals, c('reho', pROC::auc(roc(male.data$usageBin ~ cvPredVals))))

# ALFF
male.data <- male.data.all.m[complete.cases(male.data.all.m[,grep('rest_jlf_alff', names(male.data))]),]
foldsToLoop <- createFolds(male.data$usageBin, table(male.data$usageBin)[2])
cvPredVals <- rep(NA, length(male.data$usageBin))
for(q in seq(1, length(foldsToLoop))){
    index <- foldsToLoop[[q]]
    #volMod <- buildStepROCModel(y=male.data$usageBin[-index], x=male.data[-index,grep('_jlf_vol_', names(male.data))], varAdd=6)
    #outModel <- as.formula(paste('usageBin~', volMod[dim(volMod)[1],2]))
    # Now build this model and
    # build a lasso model
    optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('rest_jlf_alff', names(male.data))]), alpha=.55, family="binomial", parallel=T)
    print(optLam$lambda.min)
    lasModel <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('rest_jlf_alff', names(male.data))]), alpha=.55, lambda=optLam$lambda.min)
    #tmpModel <- glm(outModel, data=male.data[-index,], family=binomial())
    cvPredVals[index] <- predict(lasModel, newx=as.matrix(male.data[index,grep('rest_jlf_alff_', names(male.data))]), type='response')
    
}
plot(roc(male.data$usageBin ~ cvPredVals))
pROC::auc(roc(male.data$usageBin ~ cvPredVals))
aucVals <- rbind(aucVals, c('alff', pROC::auc(roc(male.data$usageBin ~ cvPredVals))))

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
    optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('dti_dtitk_jhulabel_fa', names(male.data))]), alpha=.55, family="binomial", parallel=T)
    print(optLam$lambda.min)
    lasModel <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('dti_dtitk_jhulabel_fa', names(male.data))]), alpha=.55, lambda=optLam$lambda.min)
    #tmpModel <- glm(outModel, data=male.data[-index,], family=binomial())
    cvPredVals[index] <- predict(lasModel, newx=as.matrix(male.data[index,grep('dti_dtitk_jhulabel_fa', names(male.data))]), type='response')
    
}
plot(roc(male.data$usageBin ~ cvPredVals))
pROC::auc(roc(male.data$usageBin ~ cvPredVals))
aucVals <- rbind(aucVals, c('fa', pROC::auc(roc(male.data$usageBin ~ cvPredVals))))




aucVals <- as.data.frame(aucVals)
aucVals$V2 <- as.numeric(as.character(aucVals$V2))
# Now plot the auc Values
aucPlot <- ggplot(aucVals, aes(x=V1, y=as.numeric(as.character(V2)))) +
 geom_col()

pdf('test.pdf')
print(aucPlot)
dev.off()
