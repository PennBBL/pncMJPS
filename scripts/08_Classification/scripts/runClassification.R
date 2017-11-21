## Load library(s)
source("/home/arosen/adroseHelperScripts/R/afgrHelpFunc.R")
install_load('gbm', 'adabag', 'randomForest', 'caret', 'gbm', 'pROC')

## Load the data - this has to be age regressed data
vol.data <- read.csv('/data/joy/BBL/projects/pncMJPS/data/ageRegressedData/volumeData.csv')
ct.data <- read.csv('/data/joy/BBL/projects/pncMJPS/data/ageRegressedData/ctData.csv')
cc.data <- read.csv('/data/joy/BBL/projects/pncMJPS/data/ageRegressedData/ccData.csv')
gmd.data <- read.csv('/data/joy/BBL/projects/pncMJPS/data/ageRegressedData/gmdData.csv')
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
index <- createFolds(dataAll$outcome, k=5, returnTrain=T, list=T)[[1]]

f# Now train the hyperparameters in the training dataset
fitControl <- trainControl(method='cv', number=5)
gbmGrid <- expand.grid(interaction.depth=c(1, 2, 3),
  n.trees=seq(1,100,5), 
  shrinkage=0.001, 
  n.minobsinnode=20)
set.seed(16)
gbmFit1 <- train(outcome~., data=dataAll[index,], distribution='bernoulli', method='gbm', trControl=fitControl, tuneGrid=gbmGrid)

dataAll$outcome <- as.character(dataAll$outcome)
tmp <- gbm(outcome ~ ., distribution = 'bernoulli', n.trees=100, n.minobsinnode=20, shrinkage=0.001, data=dataAll[index,],cv.folds=5, interaction.depth=1)
plot(roc(dataAll$outcome[index] ~ predict(tmp, n.trees=4, type='response'))) 
auc(roc(dataAll$outcome[index] ~ predict(tmp, n.trees=4, type='response'))) 
plot(roc(dataAll$outcome[-index] ~ predict(tmp, n.trees=4, type='response', newdata=dataAll[-index,])))
auc(roc(dataAll$outcome[-index] ~ predict(tmp, n.trees=20, type='response', newdata=dataAll[-index,])))

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
plot(roc(dataAll$outcome[-index] ~ predict(tmp, n.trees=100, type='response', newdata=dataAll[-index,]))) 
dev.off()

# Now do this all in never vs frequent
dataAll <- dataAll[which(all.data.male$dosage == 0 | all.data.male$dosage > 4),]
index <- createFolds(dataAll$outcome, k=5, returnTrain=T, list=T)[[1]]

dataAll$outcome <- as.character(dataAll$outcome)
fitControl <- trainControl(method='repeatedcv', number=5, repeats=5)
gbmGrid <- expand.grid(interaction.depth=c(1, 2, 3),
  n.trees=(1:30)*10, 
  shrinkage=0.001, 
  n.minobsinnode=18)
set.seed(16)
#gbmFit1 <- train(outcome~., data=dataAll[index,], method='gbm', trControl=fitControl, tuneGrid=gbmGrid)

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
