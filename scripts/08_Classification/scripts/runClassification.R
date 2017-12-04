## Load library(s)
source("/home/arosen/adroseHelperScripts/R/afgrHelpFunc.R")
install_load('gbm', 'adabag', 'randomForest', 'caret', 'gbm', 'pROC', 'doMC', 'plyr')

## Load the data - this has to be age regressed data
vol.data <- read.csv('/data/joy/BBL/projects/pncMJPS/data/ageRegressedData/volumeData.csv')
ct.data <- read.csv('/data/joy/BBL/projects/pncMJPS/data/ageRegressedData/ctData.csv')
cc.data <- read.csv('/data/joy/BBL/projects/pncMJPS/data/ageRegressedData/ccData.csv')
gmd.data <- read.csv('/data/joy/BBL/projects/pncMJPS/data/ageRegressedData/gmdData.csv')
cbf.data <- read.csv('/data/joy/BBL/projects/pncMJPS/data/ageRegressedData/cbfData.csv')
reho.data <- read.csv('/data/joy/BBL/projects/pncMJPS/data/ageRegressedData/rehoData.csv')
alff.data <- read.csv('/data/joy/BBL/projects/pncMJPS/data/ageRegressedData/alffData.csv')
tr.data <- read.csv('/data/joy/BBL/projects/pncMJPS/data/ageRegressedData/jlfTRData.csv')
fa.data <- read.csv('/data/joy/BBL/projects/pncMJPS/data/ageRegressedData/jhuFALabel.csv')

# Now grab the MJ data 
mjData <- read.csv('/data/joy/BBL/projects/pncMJPS/data/n9462_mj_ps_cnb_fortmm.csv')
# Now create an ordinal variable for the MJ dosage
mjData$dosage <- NA
mjData$dosage[which(mjData$marcat=='MJ Non-User')] <- 0
mjData$dosage[which(mjData$marcat=='MJ User' & mjData$mjpastyr=='')] <- 1
mjData$dosage[which(mjData$mjpastyr=="Less than once a month")] <- 2
mjData$dosage[which(mjData$mjpastyr=="About once a month")] <- 3
mjData$dosage[which(mjData$mjpastyr=="2-3 times a month")] <- 4
mjData$dosage[which(mjData$mjpastyr=="1-2 times a week")] <- 5
mjData$dosage[which(mjData$mjpastyr=="3-4 times a week")] <- 6
mjData$dosage[which(mjData$mjpastyr=="Everyday or nearly every day")] <- 7
fakeData <- read.csv('/data/joy/BBL/projects/pncMJPS/data/fakesub_exclude.csv')

# Now merge the data and prepare it 
all.data <- merge(vol.data, ct.data)
all.data <- merge(all.data, cc.data)
all.data <- merge(all.data, gmd.data)
all.data <- merge(all.data, cbf.data)
all.data <- merge(all.data, reho.data)
all.data <- merge(all.data, alff.data)
all.data <- merge(all.data, tr.data)
all.data <- merge(all.data, fa.data)
all.data <- merge(all.data, mjData)

# Now apply exclusions
all.data <- all.data[-which(all.data$dosage==1),]
all.data <- all.data[complete.cases(all.data$dosage),]

# Now find our subjects with MJ data 
# Start with male data
all.data.male <- all.data[which(all.data$sex==1),]

# Now create a useage outcome ID
all.data.male$useageBin <- 0
all.data.male$useageBin[all.data.male$dosage>1] <- 1
outcome <- all.data.male$useageBin
dataAll <- all.data.male[,grep('_jlf_', names(all.data.male))]
dataAll <- cbind(outcome, dataAll)
rm(outcome)
dataAll$outcome <- as.factor(dataAll$outcome)
# Now create a train test index
index <- createFolds(dataAll$outcome, k=3, returnTrain=T, list=T)
tree1 <- train(y=dataAll$outcome[index[[1]]], x=dataAll[index[[1]],-1], method='rpart', tuneLength=20, metric="ROC", trControl = fitControl)
treeTest <- rpart(outcome ~ ., data=dataAll[index[[1]],], method='anova', control=rpart.control(minsplit=20,minbucket=15))
# Now train the hyperparameters in the training dataset
fitControl <- trainControl(method='repeatedcv', number=3, repeats=3, classProbs=T, summaryFunction=twoClassSummary)
#gbmGrid <- expand.grid(interaction.depth = (1:5) * 2, n.trees = (1:10)*25, shrinkage = c(.1, .001), n.minobsinnode=c(10, 15, 18, 20))
set.seed(16)
#gbmFit1 <- train(outcome~., data=dataAll[index,], distribution='bernoulli', method='gbm', trControl=fitControl, tuneGrid=gbmGrid)
#dataAll$outcome <- factor(dataAll$outcome)
#dataAll$outcome <- revalue(dataAll$outcome, c('0'='NotUse', '1'='User'))
#registerDoMC(cores = 3)
#gbmFit1 <- train(outcome~., data=dataAll, distribution='bernoulli', method='gbm', trControl=fitControl, tuneGrid=gbmGrid,metric="Kappa")
outPredVals <- rep(NA, dim(dataAll)[1])
for(i in 1:3){
  # First train a model
  tmp <- gbm(outcome ~ ., distribution='bernoulli', n.trees=75, n.minobsinnode=18, shrinkage=.01, data=dataAll[index[[i]],], cv.folds=3, interaction.depth=2)
  # Now predict in the left out
  outPredVals[-index[[i]]] <- predict(tmp, n.trees=50, newdata=dataAll[-index[[i]],], type='response')
}
dataAll$outcome <- as.character(dataAll$outcome)
dataAll$outcome <- revalue(dataAll$outcome, c('NotUse'='0', 'User'='1'))
tmp <- gbm(outcome ~ ., distribution = 'bernoulli', n.trees=100, n.minobsinnode=18, shrinkage=0.001, data=dataAll[index[[1]],],cv.folds=10, interaction.depth=3)

# Now I need to get the relative imporantce values
relImp <- relative.influence(tmp, n.trees=100)
toCheck <- which(relImp!=0)
relImp <- relImp[toCheck]
relImp <- relImp[order(relImp)]
# Now go through each of these and plot the difference 
pdf('varImp.pdf')
for(i in names(relImp)){
  sumVal <- summarySE(data=dataAll, groupvars='outcome', measurevar=i)
  # Now make the ggplot value
  bg1 <- ggplot(sumVal, aes(x=outcome, y=sumVal[,3], group=outcome)) +
    geom_bar(stat='identity',position=position_dodge(), width=.5) +
    labs(title=i)+ 
    geom_errorbar(aes(ymin=as.numeric(as.character(sumVal[,3]))-se, 
      ymax=as.numeric(as.character(sumVal[,3]))+se), 
      width = .1, position=position_dodge(.9))
  print(bg1)
}
plot(roc(dataAll$outcome[index] ~ predict(tmp, n.trees=100, type='response'))) 
plot(roc(dataAll$outcome[-index[[1]]] ~ predict(tmp, n.trees=81, type='response', newdata=dataAll[-index[[1]],]))) 
dev.off()

# Now do this all in never vs frequent
dataAll <- dataAll[which(all.data.male$dosage == 0 | all.data.male$dosage > 5),]
index <- createFolds(dataAll$outcome, k=5, returnTrain=T, list=T)[[1]]

dataAll$outcome <- as.character(dataAll$outcome)
fitControl <- trainControl(method='repeatedcv', number=3, repeats=3, classProbs=T, summaryFunction=twoClassSummary)
gbmGrid <- expand.grid(interaction.depth = (1:5) * 2, n.trees = (1:10)*25, shrinkage = c(.1, .001), n.minobsinnode=c(10, 15, 18, 20))
set.seed(16)
dataAll$outcome <- revalue(dataAll$outcome, c('NotUse'='0', 'User'='1'))
gbmFit1 <- train(outcome~., data=dataAll[index,], method='gbm', trControl=fitControl, tuneGrid=gbmGrid)

tmp <- gbm(outcome ~ ., distribution = 'bernoulli', n.trees=100, n.minobsinnode=18, shrinkage=0.001, data=dataAll[index,], cv.folds=10)
treeVal <- gbm.perf(tmp)
plot(roc(dataAll$outcome[index] ~ predict(tmp, n.trees=65, type='response')))
auc(roc(dataAll$outcome[index] ~ predict(tmp, n.trees=65, type='response'))) 
plot(roc(dataAll$outcome[-index] ~ predict(tmp, n.trees=65, type='response', newdata=dataAll[-index,]))) 
auc(roc(dataAll$outcome[-index] ~ predict(tmp, n.trees=65, type='response', newdata=dataAll[-index,]))) 

# Now I need to get the relative imporantce values
relImp <- relative.influence(tmp, n.trees=65)
toCheck <- which(relImp!=0)
relImp <- relImp[toCheck]
relImp <- relImp[order(relImp)]
# Now go through each of these and plot the difference 
pdf('varImpF.pdf')
for(i in names(relImp)){
  sumVal <- summarySE(data=dataAll, groupvars='outcome', measurevar=i)
  # Now make the ggplot value
  bg1 <- ggplot(sumVal, aes(x=outcome, y=sumVal[,3], group=outcome)) +
    geom_bar(stat='identity',position=position_dodge(), width=.5) +
    labs(title=i)+ 
    geom_errorbar(aes(ymin=as.numeric(as.character(sumVal[,3]))-se, 
      ymax=as.numeric(as.character(sumVal[,3]))+se), 
      width = .1, position=position_dodge(.9))
  print(bg1)
}
plot(roc(dataAll$outcome[index] ~ predict(tmp, n.trees=100, type='response'))) 
plot(roc(dataAll$outcome[-index] ~ predict(tmp, n.trees=100, type='response', newdata=dataAll[-index,]))) 
dev.off()
