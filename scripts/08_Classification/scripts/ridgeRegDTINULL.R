## Load library(s)
install_load('psych', 'ggplot2', 'pROC', 'ggrepel', 'caret', 'randomForest', 'MatchIt', 'glmnet', 'foreach', 'useful', 'doParallel')
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
male.data <- all.data[which(all.data$sex==2),]
#male.data <- male.data[which(male.data$dosage == 0 | male.data$dosage > 5),]
female.data <- all.data[which(all.data$sex==1),]

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
propValue <- table(male.data$usageBin)[2]/sum(table(male.data$usageBin))
male.data.all.m$usageBinOrig <- male.data.all.m$usageBin

# Now create a random binary vector and see how well
# we can classify these two groups
# We are going to run this 100 times and then plot a heat ROC map
allAUC <- NULL
cl <- makeCluster(8)
registerDoParallel(cl)
allAUC <- foreach(z=seq(1,100), .combine=rbind) %dopar%{
  # Load library(s)
  install_load('glmnet', 'caret', 'pROC')
  # Create a random binary outcome
  male.data.all.m$usageBin <- rbinom(dim(male.data.all.m)[1], 1, propValue)
  # Now lets see how well we can build our model in a cross validated fashion
  # This will be done within modality just to explore things
  aucVals <- NULL
  # tr
  male.data <- male.data.all.m[complete.cases(male.data.all.m[,grep('dti_jlf_tr', names(male.data))]),]
  male.data <- male.data[,-grep('dti_jlf_tr_MeanTR', names(male.data))]
  foldsToLoop <- createFolds(male.data$usageBin, table(male.data$usageBin)[2])
  cvPredVals <- rep(NA, length(male.data$usageBin))
  cvPredValsReal <- rep(NA, length(male.data$usageBin))
  for(q in seq(1, length(foldsToLoop))){
      # First do the random folds
      index <- foldsToLoop[[q]]
      # Now build this model and
      # build a lasso model
      optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('dti_jlf_tr', names(male.data))]), alpha=0, family="binomial", parallel=F)
      lasModel <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('dti_jlf_tr', names(male.data))]), alpha=0, lambda=optLam$lambda.min)
      cvPredVals[index] <- predict(lasModel, newx=as.matrix(male.data[index,grep('dti_jlf_tr_', names(male.data))]), type='response')
      # Now do the real folds
      optLam <- cv.glmnet(y=as.vector(male.data$usageBinOrig[-index]), x=as.matrix(male.data[-index,grep('dti_jlf_tr', names(male.data))]), alpha=0, family="binomial", parallel=F)
      lasModel <- glmnet(y=as.vector(male.data$usageBinOrig[-index]), x=as.matrix(male.data[-index,grep('dti_jlf_tr', names(male.data))]), alpha=0, lambda=optLam$lambda.min)
      cvPredValsReal[index] <- predict(lasModel, newx=as.matrix(male.data[index,grep('dti_jlf_tr_', names(male.data))]), type='response')
  }
  # Now export the auc value to our bootstrapped AUC holder
  #outputValues <- rocdata(grp=binary.flip(male.data$usageBin), pred=cvPredVals)$roc
  #outputValues <- cbind(rocdata(grp=binary.flip(male.data$usageBin),pred=cvPredVals)$roc,rep(z,length(cvPredVals)))
  cm <- as.matrix(table(male.data$usageBin, male.data$usageBinOrig))
  diag <- diag(cm)
  n <- sum(cm)
  outputValues <- cbind(z, pROC::auc(roc(male.data$usageBin ~ cvPredVals)), pROC::auc(roc(male.data$usageBinOrig ~ cvPredValsReal)), sum(diag)/n)
  outputValues
}
colnames(allAUC) <- c('count', 'fake', 'real', 'accuracyOfFake')
# Kill the cluster
stopCluster(cl)

# Now calculate a t.test between the null and the real distribution
outStatMD <- t.test(allAUC[,2], allAUC[,3], alternative='less', paired=T)

# FA
allAUC <- NULL
cl <- makeCluster(8)
registerDoParallel(cl)
allAUC <- foreach(z=seq(1,100), .combine=rbind) %dopar%{
    # Load library(s)
    install_load('glmnet', 'caret', 'pROC')
    # Create a random binary outcome
    male.data.all.m$usageBin <- rbinom(dim(male.data.all.m)[1], 1, propValue)
    # Now lets see how well we can build our model in a cross validated fashion
    # This will be done within modality just to explore things
    aucVals <- NULL
    # tr
    male.data <- male.data.all.m[complete.cases(male.data.all.m[,grep('dti_dtitk_jhulabel_fa', names(male.data))]),]
    foldsToLoop <- createFolds(male.data$usageBin, table(male.data$usageBin)[2])
    cvPredVals <- rep(NA, length(male.data$usageBin))
    cvPredValsReal <- rep(NA, length(male.data$usageBin))
    for(q in seq(1, length(foldsToLoop))){
        # First do the random folds
        index <- foldsToLoop[[q]]
        # Now build this model and
        # build a lasso model
        optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('dti_dtitk_jhulabel_fa', names(male.data))]), alpha=0, family="binomial", parallel=F)
        lasModel <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('dti_dtitk_jhulabel_fa', names(male.data))]), alpha=0, lambda=optLam$lambda.min)
        cvPredVals[index] <- predict(lasModel, newx=as.matrix(male.data[index,grep('dti_dtitk_jhulabel_fa', names(male.data))]), type='response')
        # Now do the real folds
        optLam <- cv.glmnet(y=as.vector(male.data$usageBinOrig[-index]), x=as.matrix(male.data[-index,grep('dti_dtitk_jhulabel_fa', names(male.data))]), alpha=0, family="binomial", parallel=F)
        lasModel <- glmnet(y=as.vector(male.data$usageBinOrig[-index]), x=as.matrix(male.data[-index,grep('dti_dtitk_jhulabel_fa', names(male.data))]), alpha=0, lambda=optLam$lambda.min)
        cvPredValsReal[index] <- predict(lasModel, newx=as.matrix(male.data[index,grep('dti_dtitk_jhulabel_fa', names(male.data))]), type='response')
    }
    # Now export the auc value to our bootstrapped AUC holder
    #outputValues <- rocdata(grp=binary.flip(male.data$usageBin), pred=cvPredVals)$roc
    #outputValues <- cbind(rocdata(grp=binary.flip(male.data$usageBin),pred=cvPredVals)$roc,rep(z,length(cvPredVals)))
    cm <- as.matrix(table(male.data$usageBin, male.data$usageBinOrig))
    diag <- diag(cm)
    n <- sum(cm)
    outputValues <- cbind(z, pROC::auc(roc(male.data$usageBin ~ cvPredVals)), pROC::auc(roc(male.data$usageBinOrig ~ cvPredValsReal)), sum(diag)/n)
    outputValues
}
colnames(allAUC) <- c('count', 'fake', 'real', 'accuracyOfFake')

# Kill the cluster
stopCluster(cl)

# Now calculate a t.test between the null and the real distribution
outStatFA <- wilcox.test(x=allAUC[,2], y=allAUC[,3], paired=T, alternative='less')

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
male.data.all.m$usageBinOrig <- male.data.all.m$usageBin

# Now create a random binary vector and see how well
# we can classify these two groups
# We are going to run this 100 times and then plot a heat ROC map
allAUC <- NULL
cl <- makeCluster(8)
registerDoParallel(cl)
allAUC <- foreach(z=seq(1,1000), .combine=rbind, .errorhandling='remove') %dopar%{
    # Load library(s)
    install_load('glmnet', 'caret', 'pROC')
    # Create a random binary outcome
    male.data.all.m$usageBin <- rbinom(dim(male.data.all.m)[1], 1, propValue)
    # Now lets see how well we can build our model in a cross validated fashion
    # This will be done within modality just to explore things
    aucVals <- NULL
    # tr
    male.data <- male.data.all.m[complete.cases(male.data.all.m[,grep('dti_jlf_tr', names(male.data))]),]
    male.data <- male.data[,-grep('dti_jlf_tr_MeanTR', names(male.data))]
    foldsToLoop <- createFolds(male.data$usageBin, table(male.data$usageBin)[2])
    foldsToLoopReal <- createFolds(male.data$usageBinOrig, table(male.data$usageBinOrig)[2])
    cvPredVals <- rep(NA, length(male.data$usageBin))
    cvPredValsReal <- rep(NA, length(male.data$usageBin))
    for(q in seq(1, length(foldsToLoop))){
        # First do the random folds
        index <- foldsToLoop[[q]]
        index2 <- foldsToLoopReal[[q]]
        # Now build this model and
        # build a lasso model
        optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('dti_jlf_tr', names(male.data))]), alpha=0, family="binomial", parallel=F)
        lasModel <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('dti_jlf_tr', names(male.data))]), alpha=0, lambda=optLam$lambda.min)
        cvPredVals[index] <- predict(lasModel, newx=as.matrix(male.data[index,grep('dti_jlf_tr_', names(male.data))]), type='response')
        # Now do the real folds
        optLam <- cv.glmnet(y=as.vector(male.data$usageBinOrig[-index2]), x=as.matrix(male.data[-index2,grep('dti_jlf_tr', names(male.data))]), alpha=0, family="binomial", parallel=F)
        lasModel <- glmnet(y=as.vector(male.data$usageBinOrig[-index2]), x=as.matrix(male.data[-index2,grep('dti_jlf_tr', names(male.data))]), alpha=0, lambda=optLam$lambda.min)
        cvPredValsReal[index2] <- predict(lasModel, newx=as.matrix(male.data[index2,grep('dti_jlf_tr_', names(male.data))]), type='response')
    }
    # Now export the auc value to our bootstrapped AUC holder
    #outputValues <- rocdata(grp=binary.flip(male.data$usageBin), pred=cvPredVals)$roc
    #outputValues <- cbind(rocdata(grp=binary.flip(male.data$usageBin),pred=cvPredVals)$roc,rep(z,length(cvPredVals)))
    cm <- as.matrix(table(male.data$usageBin, male.data$usageBinOrig))
    diag <- diag(cm)
    n <- sum(cm)
    outputValues <- cbind(z, pROC::auc(roc(male.data$usageBin ~ cvPredVals)), pROC::auc(roc(male.data$usageBinOrig ~ cvPredValsReal)), sum(diag)/n)
    outputValues
}
colnames(allAUC) <- c('count', 'fake', 'real', 'accuracyOfFake')
# Kill the cluster
stopCluster(cl)

# Now calculate a t.test between the null and the real distribution
outStatMDFreq <- wilcox.test(allAUC[,2], allAUC[,3], alternative='less', paired=T)

# FA
allAUC <- NULL
cl <- makeCluster(8)
registerDoParallel(cl)
allAUC <- foreach(z=seq(1,1000), .combine=rbind, .errorhandling='remove') %dopar%{
    # Load library(s)
    install_load('glmnet', 'caret', 'pROC')
    # Create a random binary outcome
    male.data.all.m$usageBin <- rbinom(dim(male.data.all.m)[1], 1, propValue)
    # Now lets see how well we can build our model in a cross validated fashion
    # This will be done within modality just to explore things
    aucVals <- NULL
    # tr
    male.data <- male.data.all.m[complete.cases(male.data.all.m[,grep('dti_dtitk_jhulabel_fa', names(male.data))]),]
    foldsToLoop <- createFolds(male.data$usageBin, table(male.data$usageBin)[2])
    cvPredVals <- rep(NA, length(male.data$usageBin))
    cvPredValsReal <- rep(NA, length(male.data$usageBin))
    for(q in seq(1, length(foldsToLoop))){
        # First do the random folds
        index <- foldsToLoop[[q]]
        # Now build this model and
        # build a lasso model
        optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('dti_dtitk_jhulabel_fa', names(male.data))]), alpha=0, family="binomial", parallel=F)
        lasModel <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,grep('dti_dtitk_jhulabel_fa', names(male.data))]), alpha=0, lambda=optLam$lambda.min)
        cvPredVals[index] <- predict(lasModel, newx=as.matrix(male.data[index,grep('dti_dtitk_jhulabel_fa', names(male.data))]), type='response')
        # Now do the real folds
        optLam <- cv.glmnet(y=as.vector(male.data$usageBinOrig[-index]), x=as.matrix(male.data[-index,grep('dti_dtitk_jhulabel_fa', names(male.data))]), alpha=0, family="binomial", parallel=F)
        lasModel <- glmnet(y=as.vector(male.data$usageBinOrig[-index]), x=as.matrix(male.data[-index,grep('dti_dtitk_jhulabel_fa', names(male.data))]), alpha=0, lambda=optLam$lambda.min)
        cvPredValsReal[index] <- predict(lasModel, newx=as.matrix(male.data[index,grep('dti_dtitk_jhulabel_fa', names(male.data))]), type='response')
    }
    # Now export the auc value to our bootstrapped AUC holder
    #outputValues <- rocdata(grp=binary.flip(male.data$usageBin), pred=cvPredVals)$roc
    #outputValues <- cbind(rocdata(grp=binary.flip(male.data$usageBin),pred=cvPredVals)$roc,rep(z,length(cvPredVals)))
    cm <- as.matrix(table(male.data$usageBin, male.data$usageBinOrig))
    diag <- diag(cm)
    n <- sum(cm)
    outputValues <- cbind(z, pROC::auc(roc(male.data$usageBin ~ cvPredVals)), pROC::auc(roc(male.data$usageBinOrig ~ cvPredValsReal)), sum(diag)/n)
    outputValues
}
colnames(allAUC) <- c('count', 'fake', 'real', 'accuracyOfFake')

# Kill the cluster
stopCluster(cl)

# Now calculate a t.test between the null and the real distribution
outStatFAFreq <- wilcox.test(x=allAUC[,2], y=allAUC[,3], paired=T, alternative='less')
save.image(file='./fileName.RData')
