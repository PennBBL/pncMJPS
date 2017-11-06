## Load library(s)
source("/home/arosen/adroseHelperScripts/R/afgrHelpFunc.R")
install_load('gbm', 'adabag', 'randomForest')

## Load the data - this has to be age regressed data
vol.data <- read.csv('/data/joy/BBL/projects/pncMJPS/data/ageRegressedData/volumeData.csv')
ct.data <- read.csv('/data/joy/BBL/projects/pncMJPS/data/ageRegressedData/ctData.csv')
cc.data <- read.csv('/data/joy/BBL/projects/pncMJPS/data/ageRegressedData/ccData.csv')
gmd.data <- read.csv('/data/joy/BBL/projects/pncMJPS/data/ageRegressedData/gmdData.csv')
reho.data <- read.csv('/data/joy/BBL/projects/pncMJPS/data/ageRegressedData/rehoData.csv')
alff.data <- read.csv('/data/joy/BBL/projects/pncMJPS/data/ageRegressedData/alffData.csv')
tr.data <- read.csv('/data/joy/BBL/projects/pncMJPS/data/ageRegressedData/alffData.csv')

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

boosting.cv(outcome ~ ., data=dataAll, control=rpart.control(maxdepth=5, minsplit=15))

boosting(outcome ~ ., data=dataAll, control=rpart.control(maxdepth=5, minsplit=15))
