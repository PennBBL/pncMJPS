## Load library(s)
source('/home/arosen/adroseHelperScripts/R/afgrHelpFunc.R')

## Load data
mjData <- read.csv('/data/jux/BBL/projects/pncMJPS/data/n9462_mj_ps_cnb_fortmm.csv')
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
mjData$mjLabel <- NA
mjData$mjLabel[which(mjData$marcat=="MJ Non-User")] <- "NonUser"
mjData$mjLabel[which(mjData$marcat=="MJ User")] <- "User"
mjData$mjLabel[which(mjData$marcat=="MJ Frequent User")] <- "FreqUser"
mjData$mjBinLabel <- "NonUser"
mjData$mjBinLabel[which(mjData$mjLabel=="User")] <- "User"
mjData$mjBinLabel[which(mjData$mjLabel=="FreqUser")] <- "User"
mjData$mjLabel <- as.factor(mjData$mjLabel)

fakeData <- read.csv('/data/jux/BBL/projects/pncMJPS/data/fakesub_exclude.csv')
# Now create a new label which collapses ps op and td into a new label. This will be a factor w/ two levels
psData <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/clinical/n1601_diagnosis_dxpmr_20170509.csv')
psData <- psData[,-grep('goassessDxpmr4', names(psData))]
psData$pathLabel <- NA
psData$pathLabel[which(psData$goassessDxpmr7=='TD' | psData$goassessDxpmr7=='OP')] <- 'TDOP'
psData$pathLabel[which(psData$goassessDxpmr7=='PS')] <- 'PS'
psData$pathLabel <- as.factor(psData$pathLabel)

demoData <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_demographics_go1_20161212.csv')
psData <- merge(psData, demoData, by=c('bblid', 'scanid'))
psData <- merge(psData, mjData)

all.data <- read.csv('/data/jux/BBL/projects/pncMJPS/scripts/07_MJEffects/scripts/n1601_imagingDataDump_2018-04-04.csv')
all.data <- merge(all.data, psData)

# Now add the clinical bifactor scores
fac.data <- read.csv('/data/joy/BBL/studies/pnc/n9498_dataFreeze/clinical/n9498_goassess_itemwise_bifactor_scores_age_regressed_20170131.csv')
all.data <- merge(all.data, fac.data)

# Now also add our fs ct values
fs.values <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_freesurferCt_20180213.csv')
all.data <- merge(all.data, fs.values)
all.data$fsAvgCT <- (all.data$LThickness + all.data$RThickness)/2

# Now remove all subjects that have no marcat variable, are older than 14, do not pass structural QA, and endorsed fake drugs
all.data <- all.data[-which(all.data$marcat==''),]
all.data <- all.data[which((all.data$ageAtScan1/12)>=14),]
all.data <- all.data[-which(all.data$averageManualRating==0),]
all.data <- all.data[-which(all.data$bblid %in% fakeData$bblid),]
all.data <- all.data[-which(all.data$dosage==1),]

## Now save our RDS file
saveRDS(object=all.data, file="mjAnovaData.RDS")
