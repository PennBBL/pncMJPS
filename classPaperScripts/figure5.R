## Load library(s)
install_load('psych', 'ggplot2', 'pROC', 'ggrepel', 'caret', 'randomForest', 'MatchIt', 'glmnet', 'foreach', 'useful', 'doParallel')
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
tmpDat <- male.data[c('bblid', 'scanid', 'usageBin', 'ageAtScan1', 'envSES')]
mod <- matchit(usageBin ~ ageAtScan1 + envSES, data=tmpDat, ratio=3, na.action=na.omit)
male.data.all <- male.data
male.data <- male.data[as.vector(mod$match.matrix),]
male.data <- rbind(male.data, male.data.all[which(male.data.all$usageBin==1),])
male.data.all.m <- male.data
propValueMale <- table(male.data$usageBin)[2]/sum(table(male.data$usageBin))
male.data.all.m$usageBinOrig <- male.data.all.m$usageBin

# Now do female
tmpDat <- female.data[c('bblid', 'scanid', 'usageBin', 'ageAtScan1', 'envSES')]
mod <- matchit(usageBin ~ ageAtScan1 + envSES, data=tmpDat, ratio=3, na.action=na.omit)
female.data.all <- female.data
female.data <- female.data[as.vector(mod$match.matrix),]
female.data <- rbind(female.data, female.data.all[which(female.data.all$usageBin==1),])
female.data.all.m <- female.data
propValueFemale <- table(female.data$usageBin)[2]/sum(table(female.data$usageBin))
female.data.all.m$usageBinOrig <- female.data.all.m$usageBin

# Now lets make our bootstrapped male ROC curves
cl <- makeCluster(8)
registerDoParallel(cl)
allAUCM <- foreach(z=seq(1,300), .combine=rbind) %dopar%{
    # Load library(s)
    install_load('glmnet', 'caret', 'pROC', 'useful')
    # Create a random binary outcome
    male.data.all.m$usageBin <- rbinom(dim(male.data.all.m)[1], 1, propValueMale)
    # Now lets see how well we can build our model in a cross validated fashion
    male.data <- male.data.all.m[complete.cases(male.data.all.m[,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]),]
    foldsToLoop <- createFolds(male.data$usageBin, table(male.data$usageBin)[2])
    cvPredVals <- rep(NA, length(male.data$usageBin))
    cvPredValsReal <- rep(NA, length(male.data$usageBin))
    for(q in seq(1, length(foldsToLoop))){
        index <- foldsToLoop[[q]]
        # build a ridge model with the fake data
        optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, family="binomial", parallel=T)
        lasModel1 <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, lambda=optLam$lambda.min)
        
        # Now do the real labels
        optLam <- cv.glmnet(y=as.vector(male.data$usageBinOrig[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, family="binomial", parallel=T)
        lasModel2 <- glmnet(y=as.vector(male.data$usageBinOrig[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, lambda=optLam$lambda.min)
        
        cvPredValsReal[index] <- predict(lasModel2, newx=as.matrix(male.data[index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), type='response')
        cvPredVals[index] <- predict(lasModel1, newx=as.matrix(male.data[index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), type='response')
        
    }
    # Now export the auc value to our bootstrapped AUC holder
    outputValues <- cbind(rocdata(grp=binary.flip(male.data$usageBin),pred=cvPredVals)$roc,rep(z,length(cvPredVals)), rep('Fake', length(cvPredVals)))
    colnames(outputValues) <- c('x', 'y', 'Fold', 'Status')
    outputValues2 <- cbind(rocdata(grp=binary.flip(male.data$usageBinOrig),pred=cvPredValsReal)$roc,rep(z,length(cvPredVals)), rep('Real', length(cvPredVals)))
    colnames(outputValues2) <- c('x', 'y', 'Fold', 'Status')
    outputValues <- rbind(outputValues, outputValues2)
    outputValues
}

# Now prepare a plot
allAUCM$facPlot <- paste(allAUCM$Fold, allAUCM$Status)
p1 <- ggplot(allAUCM, aes(x = x, y = y, group=facPlot, col=Status)) +
  geom_line(alpha=1/10, size=3) +
  theme_bw() +
  scale_fill_manual(values=c("Real"="Blue", "Fake"="Red")) +
  geom_abline (intercept = 0, slope = 1) +
  scale_x_continuous("1-Specificity") +
  scale_y_continuous("Sensitivity") +
  theme(legend.position="none")

# Now do the female none vs usage
allAUCF <- foreach(z=seq(1,300), .combine=rbind) %dopar%{
    # Load library(s)
    install_load('glmnet', 'caret', 'pROC', 'useful')
    # Create a random binary outcome
    female.data.all.m$usageBin <- rbinom(dim(female.data.all.m)[1], 1, propValueMale)
    # Now lets see how well we can build our model in a cross validated fashion
    female.data <- female.data.all.m[complete.cases(female.data.all.m[,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]),]
    foldsToLoop <- createFolds(female.data$usageBin, table(female.data$usageBin)[2])
    cvPredVals <- rep(NA, length(female.data$usageBin))
    cvPredValsReal <- rep(NA, length(female.data$usageBin))
    for(q in seq(1, length(foldsToLoop))){
        index <- foldsToLoop[[q]]
        # build a ridge model with the fake data
        optLam <- cv.glmnet(y=as.vector(female.data$usageBin[-index]), x=as.matrix(female.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), alpha=0, family="binomial", parallel=T)
        lasModel1 <- glmnet(y=as.vector(female.data$usageBin[-index]), x=as.matrix(female.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), alpha=0, lambda=optLam$lambda.min)
        
        # Now do the real labels
        optLam <- cv.glmnet(y=as.vector(female.data$usageBinOrig[-index]), x=as.matrix(female.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), alpha=0, family="binomial", parallel=T)
        lasModel2 <- glmnet(y=as.vector(female.data$usageBinOrig[-index]), x=as.matrix(female.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), alpha=0, lambda=optLam$lambda.min)
        
        cvPredValsReal[index] <- predict(lasModel2, newx=as.matrix(female.data[index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), type='response')
        cvPredVals[index] <- predict(lasModel1, newx=as.matrix(female.data[index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), type='response')
        
    }
    # Now export the auc value to our bootstrapped AUC holder
    outputValues <- cbind(rocdata(grp=binary.flip(female.data$usageBin),pred=cvPredVals)$roc,rep(z,length(cvPredVals)), rep('Fake', length(cvPredVals)))
    colnames(outputValues) <- c('x', 'y', 'Fold', 'Status')
    outputValues2 <- cbind(rocdata(grp=binary.flip(female.data$usageBinOrig),pred=cvPredValsReal)$roc,rep(z,length(cvPredVals)), rep('Real', length(cvPredVals)))
    colnames(outputValues2) <- c('x', 'y', 'Fold', 'Status')
    outputValues <- rbind(outputValues, outputValues2)
    outputValues
}

# Now prepare a plot
allAUCF$facPlot <- paste(allAUCF$Fold, allAUCF$Status)
p2 <- ggplot(allAUCF, aes(x = x, y = y, group=facPlot, col=Status)) +
  geom_line(alpha=1/10, size=3) +
  theme_bw() +
  scale_fill_manual(values=c("Real"="Blue", "Fake"="Red")) +
  geom_abline (intercept = 0, slope = 1) +
  scale_x_continuous("1-Specificity") +
  scale_y_continuous("Sensitivity") +
  theme(legend.position="none")

## Now move on to the user vs frequent user!
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
propValueMale <- table(male.data$usageBin)[2]/sum(table(male.data$usageBin))
female.data <- all.data[which(all.data$sex==2),]
female.data.all.m <- female.data
propValueFemale <- table(female.data$usageBin)[2]/sum(table(female.data$usageBin))

# Now run the loop for the males
allAUCM <- foreach(z=seq(1,300), .combine=rbind, .errorhandling='remove') %dopar%{
    # Load library(s)
    install_load('glmnet', 'caret', 'pROC', 'useful')
    # Create a random binary outcome
    male.data.all.m$usageBin <- rbinom(dim(male.data.all.m)[1], 1, propValueMale)
    # Now lets see how well we can build our model in a cross validated fashion
    male.data <- male.data.all.m[complete.cases(male.data.all.m[,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]),]
    foldsToLoop <- createFolds(male.data$usageBin, table(male.data$usageBin)[2])
    cvPredVals <- rep(NA, length(male.data$usageBin))
    cvPredValsReal <- rep(NA, length(male.data$usageBin))
    for(q in seq(1, length(foldsToLoop))){
        index <- foldsToLoop[[q]]
        # build a ridge model with the fake data
        optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, family="binomial", parallel=T)
        lasModel1 <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, lambda=optLam$lambda.min)
        
        # Now do the real labels
        optLam <- cv.glmnet(y=as.vector(male.data$usageBinOrig[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, family="binomial", parallel=T)
        lasModel2 <- glmnet(y=as.vector(male.data$usageBinOrig[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, lambda=optLam$lambda.min)
        
        cvPredValsReal[index] <- predict(lasModel2, newx=as.matrix(male.data[index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), type='response')
        cvPredVals[index] <- predict(lasModel1, newx=as.matrix(male.data[index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), type='response')
        
    }
    # Now export the auc value to our bootstrapped AUC holder
    outputValues <- cbind(rocdata(grp=binary.flip(male.data$usageBin),pred=cvPredVals)$roc,rep(z,length(cvPredVals)), rep('Fake', length(cvPredVals)))
    colnames(outputValues) <- c('x', 'y', 'Fold', 'Status')
    outputValues2 <- cbind(rocdata(grp=binary.flip(male.data$usageBinOrig),pred=cvPredValsReal)$roc,rep(z,length(cvPredVals)), rep('Real', length(cvPredVals)))
    colnames(outputValues2) <- c('x', 'y', 'Fold', 'Status')
    outputValues <- rbind(outputValues, outputValues2)
    outputValues
}

# Now prepare a plot
allAUCM$facPlot <- paste(allAUCM$Fold, allAUCM$Status)
p3 <- ggplot(allAUCM, aes(x = x, y = y, group=facPlot, col=Status)) +
  geom_line(alpha=1/10, size=3) +
  theme_bw() +
  scale_fill_manual(values=c("Real"="Blue", "Fake"="Red")) +
  geom_abline (intercept = 0, slope = 1) +
  scale_x_continuous("1-Specificity") +
scale_y_continuous("Sensitivity") +
theme(legend.position="none")

# Now onto the females
allAUCF <- foreach(z=seq(1,300), .combine=rbind, .errorhandling='remove') %dopar%{
    # Load library(s)
    install_load('glmnet', 'caret', 'pROC', 'useful')
    # Create a random binary outcome
    female.data.all.m$usageBin <- rbinom(dim(female.data.all.m)[1], 1, propValueFemale)
    # Now lets see how well we can build our model in a cross validated fashion
    female.data <- female.data.all.m[complete.cases(female.data.all.m[,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]),]
    foldsToLoop <- createFolds(female.data$usageBin, table(female.data$usageBin)[2])
    cvPredVals <- rep(NA, length(female.data$usageBin))
    cvPredValsReal <- rep(NA, length(female.data$usageBin))
    for(q in seq(1, length(foldsToLoop))){
        index <- foldsToLoop[[q]]
        # build a ridge model with the fake data
        optLam <- cv.glmnet(y=as.vector(female.data$usageBin[-index]), x=as.matrix(female.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), alpha=0, family="binomial", parallel=T)
        lasModel1 <- glmnet(y=as.vector(female.data$usageBin[-index]), x=as.matrix(female.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), alpha=0, lambda=optLam$lambda.min)
        
        # Now do the real labels
        optLam <- cv.glmnet(y=as.vector(female.data$usageBinOrig[-index]), x=as.matrix(female.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), alpha=0, family="binomial", parallel=T)
        lasModel2 <- glmnet(y=as.vector(female.data$usageBinOrig[-index]), x=as.matrix(female.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), alpha=0, lambda=optLam$lambda.min)
        
        cvPredValsReal[index] <- predict(lasModel2, newx=as.matrix(female.data[index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), type='response')
        cvPredVals[index] <- predict(lasModel1, newx=as.matrix(female.data[index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), type='response')
        
    }
    # Now export the auc value to our bootstrapped AUC holder
    outputValues <- cbind(rocdata(grp=female.data$usageBin,pred=cvPredVals)$roc,rep(z,length(cvPredVals)), rep('Fake', length(cvPredVals)))
    colnames(outputValues) <- c('x', 'y', 'Fold', 'Status')
    outputValues2 <- cbind(rocdata(grp=female.data$usageBinOrig,pred=cvPredValsReal)$roc,rep(z,length(cvPredVals)), rep('Real', length(cvPredVals)))
    colnames(outputValues2) <- c('x', 'y', 'Fold', 'Status')
    outputValues <- rbind(outputValues, outputValues2)
    outputValues
}
# Kill the cluster
stopCluster(cl)

# Now prepare a plot
allAUCF$facPlot <- paste(allAUCF$Fold, allAUCF$Status)
p4 <- ggplot(allAUCF, aes(x = x, y = y, group=facPlot, col=Status)) +
  geom_line(alpha=1/10, size=3) +
  theme_bw() +
  scale_fill_manual(values=c("Real"="Blue", "Fake"="Red")) +
  geom_abline (intercept = 0, slope = 1) +
  scale_x_continuous("1-Specificity") +
  scale_y_continuous("Sensitivity") +
  theme(legend.position="none")

# Now combine all of the plots
pdf("testOut.pdf", height=20, width=20)
multiplot(p1, p2, p3, p4, cols=2)
dev.off()
