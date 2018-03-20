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

# Now lets make our bootstrapped male ROC curves
cl <- makeCluster(8)
registerDoParallel(cl)
allAUCM <- foreach(z=seq(1,10), .combine=rbind) %dopar%{
    # Load library(s)
    install_load('glmnet', 'caret', 'pROC', 'useful')
    # Create a random binary outcome
    male.data.all.m$usageBin <- rbinom(dim(male.data.all.m)[1], 1, propValueMale)
    # Now lets see how well we can build our model in a cross validated fashion
    male.data <- male.data.all.m[complete.cases(male.data.all.m[,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]),]
    foldsToLoop <- createFolds(male.data$usageBin, 25)
    cvPredVals <- rep(NA, length(male.data$usageBin))
    cvPredValsReal <- rep(NA, length(male.data$usageBin))
    for(q in seq(1, length(foldsToLoop))){
        index <- foldsToLoop[[q]]
        # build a ridge model with the fake data
        optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, family="binomial", parallel=F)
        lasModel1 <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, lambda=optLam$lambda.min)
        
        # Now do the real labels
        optLam <- cv.glmnet(y=as.vector(male.data$usageBinOrig[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, family="binomial", parallel=F)
        lasModel2 <- glmnet(y=as.vector(male.data$usageBinOrig[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, lambda=optLam$lambda.min)
        
        cvPredValsReal[index] <- predict(lasModel2, newx=as.matrix(male.data[index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), type='response')
        cvPredVals[index] <- predict(lasModel1, newx=as.matrix(male.data[index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), type='response')
        
    }
    # Now export the auc value to our bootstrapped AUC holder
    outputValues <- cbind(rocdata(grp=male.data$usageBin,pred=cvPredVals)$roc,rep(z,length(cvPredVals)), rep('Fake', length(cvPredVals)))
    colnames(outputValues) <- c('x', 'y', 'Fold', 'Status')
    outputValues2 <- cbind(rocdata(grp=male.data$usageBinOrig,pred=cvPredValsReal)$roc,rep(z,length(cvPredVals)), rep('Real', length(cvPredVals)))
    colnames(outputValues2) <- c('x', 'y', 'Fold', 'Status')
    outputValues3 <- rbind(c("AUC", pROC::auc(roc(male.data$usageBinOrig~cvPredValsReal)), "AUC", pROC::auc(roc(male.data$usageBin~cvPredVals))))
    colnames(outputValues3) <- c('x', 'y', 'Fold', 'Status')
    outputValues <- rbind(outputValues, outputValues2, outputValues3)
    outputValues
}
aucOutM <- allAUCM[grep("AUC", allAUCM[,1]),]
allAUCM <- allAUCM[-grep("AUC", allAUCM[,1]),]

# Now prepare a plot
allAUCM$facPlot <- paste(allAUCM$Fold, allAUCM$Status)
allAUCM$x <- as.numeric(as.character(allAUCM$x))
allAUCM$y <- as.numeric(as.character(allAUCM$y))
p1 <- ggplot(allAUCM, aes(x = x, y = y, group=facPlot, col=Status)) +
  geom_line(alpha=1/10, size=3) +
  theme_bw() +
  scale_fill_manual(values=c("Real"="Blue", "Fake"="Red")) +
  geom_abline (intercept = 0, slope = 1) +
  scale_x_continuous("1-Specificity") +
  scale_y_continuous("Sensitivity") +
  theme(legend.position="none")

# Now make a histogram for the AUC values
colnames(aucOutM) <- c("AUC", "Val", "AUC", "Val")
aucOutM[,2] <- round(as.numeric(as.character(aucOutM[,2])), digits=3)
out1 <- summary(aucOutM)[,2]
aucOutM[,4] <- round(as.numeric(as.character(aucOutM[,4])), digits=3)
# Now perform a t.test between the two AUC values
outDiffM <- t.test(aucOutM[,2] , aucOutM[,4], paired=T, alternative='greater')
# Now prepare the data for ggplot
tmp1 <- cbind(aucOutM[,1:2])
tmp1$Status <- "Real"
tmp2 <- cbind(aucOutM[,3:4])
tmp2$Status <- "Fake"
histData <- rbind(tmp1, tmp2)
h1 <- ggplot(histData) +
  geom_histogram(data=tmp1, aes(x=Val, fill='red'), alpha=.3) +
  geom_histogram(data=tmp2, aes(x=Val, fill='bue'), alpha=.3) +
  xlim(c(.35,.9)) +
  ylim(c(0,150)) +
  theme_bw() +
  theme(legend.position="none")

# Now do the female none vs usage
allAUCF <- foreach(z=seq(1,1000), .combine=rbind) %dopar%{
    # Load library(s)
    install_load('glmnet', 'caret', 'pROC', 'useful')
    # Create a random binary outcome
    female.data.all.m$usageBin <- rbinom(dim(female.data.all.m)[1], 1, propValueMale)
    # Now lets see how well we can build our model in a cross validated fashion
    female.data <- female.data.all.m[complete.cases(female.data.all.m[,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]),]
    foldsToLoop <- createFolds(female.data$usageBin, 20)
    cvPredVals <- rep(NA, length(female.data$usageBin))
    cvPredValsReal <- rep(NA, length(female.data$usageBin))
    for(q in seq(1, length(foldsToLoop))){
        index <- foldsToLoop[[q]]
        # build a ridge model with the fake data
        optLam <- cv.glmnet(y=as.vector(female.data$usageBin[-index]), x=as.matrix(female.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), alpha=0, family="binomial", parallel=T)
        lasModel1 <- glmnet(y=as.vector(female.data$usageBin[-index]), x=as.matrix(female.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), alpha=0, lambda=optLam$lambda.min)
        
        # Now do the real labels
        optLam <- cv.glmnet(y=as.vector(female.data$usageBinOrig[-index]), x=as.matrix(female.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), alpha=0, family="binomial", parallel=F)
        lasModel2 <- glmnet(y=as.vector(female.data$usageBinOrig[-index]), x=as.matrix(female.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), alpha=0, lambda=optLam$lambda.min)
        
        cvPredValsReal[index] <- predict(lasModel2, newx=as.matrix(female.data[index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), type='response')
        cvPredVals[index] <- predict(lasModel1, newx=as.matrix(female.data[index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), type='response')
        
    }
    roc(female.data$usageBinOrig ~ cvPredValsReal)
    # Now export the auc value to our bootstrapped AUC holder
    outputValues <- cbind(rocdata(grp=female.data$usageBin,pred=cvPredVals)$roc,rep(z,length(cvPredVals)), rep('Fake', length(cvPredVals)))
    colnames(outputValues) <- c('x', 'y', 'Fold', 'Status')
    outputValues2 <- cbind(rocdata(grp=female.data$usageBinOrig,pred=cvPredValsReal)$roc,rep(z,length(cvPredVals)), rep('Real', length(cvPredVals)))
    colnames(outputValues2) <- c('x', 'y', 'Fold', 'Status')
    outputValues3 <- rbind(c("AUC", pROC::auc(roc(female.data$usageBinOrig~cvPredValsReal)), "AUC", pROC::auc(roc(female.data$usageBin~cvPredVals))))
    colnames(outputValues3) <- c('x', 'y', 'Fold', 'Status')
    outputValues <- rbind(outputValues, outputValues2, outputValues3)
    outputValues
}
aucOutF <- allAUCF[grep("AUC", allAUCF[,1]),]
allAUCF <- allAUCF[-grep("AUC", allAUCF[,1]),]

# Now prepare a plot
allAUCF$facPlot <- paste(allAUCF$Fold, allAUCF$Status)
allAUCF$x <- as.numeric(as.character(allAUCF$x))
allAUCF$y <- as.numeric(as.character(allAUCF$y))
p2 <- ggplot(allAUCF, aes(x = x, y = y, group=facPlot, col=Status)) +
  geom_line(alpha=1/10, size=3) +
  theme_bw() +
  scale_fill_manual(values=c("Real"="Blue", "Fake"="Red")) +
  geom_abline (intercept = 0, slope = 1) +
  scale_x_continuous("1-Specificity") +
  scale_y_continuous("Sensitivity") +
  theme(legend.position="none")

# Now make a histogram for the AUC values
colnames(aucOutF) <- c("AUC", "Val", "AUC", "Val")
aucOutF[,2] <- round(as.numeric(as.character(aucOutF[,2])), digits=3)
out2 <- summary(aucOutF)[,2]
aucOutF[,4] <- round(as.numeric(as.character(aucOutF[,4])), digits=3)
# Now perform a t.test between the two AUC values
outDiffF <- t.test(aucOutF[,2] , aucOutF[,4], paired=T, alternative='greater')
# Now prepare the data for ggplot
tmp1 <- cbind(aucOutF[,1:2])
tmp1$Status <- "Real"
tmp2 <- cbind(aucOutF[,3:4])
tmp2$Status <- "Fake"
histData <- rbind(tmp1, tmp2)
h2 <- ggplot(histData) +
  geom_histogram(data=tmp1, aes(x=Val, fill='red'), alpha=.3) +
  geom_histogram(data=tmp2, aes(x=Val, fill='blue'), alpha=.3) +
  xlim(c(.35,.9)) +
  ylim(c(0,150)) +
  theme_bw() +
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
allAUCM <- foreach(z=seq(1,1000), .combine=rbind, .errorhandling='remove') %dopar%{
    # Load library(s)
    install_load('glmnet', 'caret', 'pROC', 'useful')
    # Create a random binary outcome
    male.data.all.m$usageBin <- rbinom(dim(male.data.all.m)[1], 1, propValueMale)
    # Now lets see how well we can build our model in a cross validated fashion
    male.data <- male.data.all.m[complete.cases(male.data.all.m[,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]),]
    foldsToLoop <- createFolds(male.data$usageBin, 20)
    cvPredVals <- rep(NA, length(male.data$usageBin))
    cvPredValsReal <- rep(NA, length(male.data$usageBin))
    for(q in seq(1, length(foldsToLoop))){
        index <- foldsToLoop[[q]]
        # build a ridge model with the fake data
        optLam <- cv.glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, family="binomial", parallel=F)
        lasModel1 <- glmnet(y=as.vector(male.data$usageBin[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, lambda=optLam$lambda.min)
        
        # Now do the real labels
        optLam <- cv.glmnet(y=as.vector(male.data$usageBinOrig[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, family="binomial", parallel=T)
        lasModel2 <- glmnet(y=as.vector(male.data$usageBinOrig[-index]), x=as.matrix(male.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), alpha=0, lambda=optLam$lambda.min)
        
        cvPredValsReal[index] <- predict(lasModel2, newx=as.matrix(male.data[index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), type='response')
        cvPredVals[index] <- predict(lasModel1, newx=as.matrix(male.data[index,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]), type='response')
        
    }
    # Now export the auc value to our bootstrapped AUC holder
    outputValues <- cbind(rocdata(grp=male.data$usageBin,pred=cvPredVals)$roc,rep(z,length(cvPredVals)), rep('Fake', length(cvPredVals)))
    colnames(outputValues) <- c('x', 'y', 'Fold', 'Status')
    outputValues2 <- cbind(rocdata(grp=male.data$usageBinOrig,pred=cvPredValsReal)$roc,rep(z,length(cvPredVals)), rep('Real', length(cvPredVals)))
    colnames(outputValues2) <- c('x', 'y', 'Fold', 'Status')
    outputValues3 <- rbind(c("AUC", pROC::auc(roc(male.data$usageBinOrig~cvPredValsReal)), "AUC", pROC::auc(roc(male.data$usageBin~cvPredVals))))
    colnames(outputValues3) <- c('x', 'y', 'Fold', 'Status')
    outputValues <- rbind(outputValues, outputValues2, outputValues3)
    outputValues
}
aucOutM <- allAUCM[grep("AUC", allAUCM[,1]),]
allAUCM <- allAUCM[-grep("AUC", allAUCM[,1]),]

# Now prepare a plot
allAUCM$facPlot <- paste(allAUCM$Fold, allAUCM$Status)
allAUCM$x <- as.numeric(as.character(allAUCM$x))
allAUCM$y <- as.numeric(as.character(allAUCM$y))

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

# Now make a histogram for the AUC values
colnames(aucOutM) <- c("AUC", "Val", "AUC", "Val")
aucOutM[,2] <- round(as.numeric(as.character(aucOutM[,2])), digits=3)
out3 <- summary(aucOutM)[,2]
aucOutM[,4] <- round(as.numeric(as.character(aucOutM[,4])), digits=3)
# Now perform a t.test between the two AUC values
outDiffM <- t.test(aucOutM[,2] , aucOutM[,4], paired=T, alternative='greater')
# Now prepare the data for ggplot
tmp1 <- cbind(aucOutM[,1:2])
tmp1$Status <- "Real"
tmp2 <- cbind(aucOutM[,3:4])
tmp2$Status <- "Fake"
histData <- rbind(tmp1, tmp2)
h3 <- ggplot(histData) +
  geom_histogram(data=tmp1, aes(x=Val, fill='red'), alpha=.3) +
  geom_histogram(data=tmp2, aes(x=Val, fill='blue'), alpha=.3) +
  xlim(c(.35,.9)) +
  ylim(c(0,150)) +
  theme_bw() +
  theme(legend.position="none")

# Now onto the females
allAUCF <- foreach(z=seq(1,1000), .combine=rbind, .errorhandling='remove') %dopar%{
    # Load library(s)
    install_load('glmnet', 'caret', 'pROC', 'useful')
    # Create a random binary outcome
    female.data.all.m$usageBin <- rbinom(dim(female.data.all.m)[1], 1, propValueFemale)
    # Now lets see how well we can build our model in a cross validated fashion
    female.data <- female.data.all.m[complete.cases(female.data.all.m[,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]),]
    foldsToLoop <- createFolds(female.data$usageBin, 10)
    cvPredVals <- rep(NA, length(female.data$usageBin))
    cvPredValsReal <- rep(NA, length(female.data$usageBin))
    for(q in seq(1, length(foldsToLoop))){
        index <- foldsToLoop[[q]]
        # build a ridge model with the fake data
        optLam <- cv.glmnet(y=as.vector(female.data$usageBin[-index]), x=as.matrix(female.data[-index,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]), alpha=0, family="binomial", parallel=F)
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
    outputValues3 <- rbind(c("AUC", pROC::auc(roc(female.data$usageBinOrig~cvPredValsReal)), "AUC", pROC::auc(roc(female.data$usageBin~cvPredVals))))
    colnames(outputValues3) <- c('x', 'y', 'Fold', 'Status')
    outputValues <- rbind(outputValues, outputValues2, outputValues3)
    outputValues
}
# Kill the cluster
stopCluster(cl)
aucOutF <- allAUCF[grep("AUC", allAUCF[,1]),]
allAUCF <- allAUCF[-grep("AUC", allAUCF[,1]),]

# Now prepare a plot
allAUCF$facPlot <- paste(allAUCF$Fold, allAUCF$Status)
allAUCF$x <- as.numeric(as.character(allAUCF$x))
allAUCF$y <- as.numeric(as.character(allAUCF$y))

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

# Now make a histogram for the AUC values
colnames(aucOutF) <- c("AUC", "Val", "AUC", "Val")
aucOutF[,2] <- round(as.numeric(as.character(aucOutF[,2])), digits=3)
out4 <- summary(aucOutF)[,2]
aucOutF[,4] <- round(as.numeric(as.character(aucOutF[,4])), digits=3)
# Now perform a t.test between the two AUC values
outDiffF <- t.test(aucOutF[,2] , aucOutF[,4], paired=T, alternative='greater')
# Now prepare the data for ggplot
tmp1 <- cbind(aucOutF[,1:2])
tmp1$Status <- "Real"
tmp2 <- cbind(aucOutF[,3:4])
tmp2$Status <- "Fake"
histData <- rbind(tmp1, tmp2)
h4 <- ggplot(histData) +
  geom_histogram(data=tmp1, aes(x=Val, fill='red'), alpha=.3) +
  geom_histogram(data=tmp2, aes(x=Val, fill='blue'), alpha=.3) +
  xlim(c(.35,.9)) +
  ylim(c(0,150)) +
  theme_bw() +
  theme(legend.position="none")

# Now combine all of the plots
png("figure5.png", height=20, width=20, units='in', res=300)
multiplot(p1, p2, p3, p4, cols=2)
dev.off()

# Now export the histograms
pdf("histImages.pdf")
h1
h2
h3
h4
dev.off()

# Now export auc values
toWrite <- rbind(out1, out2, out3, out4)
write.csv(toWrite, 'aucStats.csv', quote=F)
