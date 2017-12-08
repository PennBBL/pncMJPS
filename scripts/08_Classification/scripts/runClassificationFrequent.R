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
all.data.male <- all.data.male[which(all.data.male$dosage == 0 | all.data.male$dosage > 5),]
all.data.male$useageBin <- 0
all.data.male$useageBin[all.data.male$dosage>1] <- 1
outcome <- all.data.male$useageBin
dataAll <- all.data.male[,grep('_jlf_', names(all.data.male))]
dataAll <- cbind(outcome, dataAll)
rm(outcome)
dataAll$outcome <- as.factor(dataAll$outcome)
dataAll$outcome <- revalue(dataAll$outcome, c('0'='NotUse', '1'='User'))
# Now create a train test index
index <- createFolds(dataAll$outcome, k=3, returnTrain=T, list=T)
fitControl <- trainControl(method = "repeatedcv",
                           repeats = 5,
                           classProbs = TRUE,
                           summaryFunction = twoClassSummary, 
		           search="random")
gridVal <- data.frame(mtry=c(1:7))
tree1 <- train(y=dataAll$outcome[index[[1]]], x=dataAll[index[[1]],-1],method='rf',metric="ROC",trControl=fitControl,verbose=T,ntree=1000, tuneGrid=gridVal)

roc(dataAll$outcome[index[[1]]] ~ predict(tree1, type='prob')[,1])
roc(dataAll$outcome[-index[[1]]] ~ predict(tree1, type='prob', newdata=dataAll[-index[[1]],])[,1])

# Now try gbm 
gbmGrid <- expand.grid(interaction.depth = (1:5), n.trees = (1:100), shrinkage = c(.1, .001), n.minobsinnode=c(10, 15, 18, 20))
gbmFit1 <- train(outcome~., data=dataAll[index[[1]],], distribution='bernoulli', method='gbm', trControl=fitControl, tuneGrid=gbmGrid)

roc(dataAll$outcome[index[[1]]] ~ predict(gbmFit1, type='prob')[,1])
roc(dataAll$outcome[-index[[1]]] ~ predict(gbmFit1, type='prob', newdata=dataAll[-index[[1]],])[,1])
