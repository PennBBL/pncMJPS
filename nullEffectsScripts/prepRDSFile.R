## Load library(s)
source('/home/arosen/adroseHelperScripts/R/afgrHelpFunc.R')
install_load('MASS','MatchIt')
## Load data
mjData <- read.csv('/data/jux/BBL/projects/pncMJPS/data/n9462_mj_ps_cnb_fortmm.csv')
mjData <- read.csv('/data/jux/BBL/projects/pncMJPS/data/n9498_go1_foradon_061518.csv')
# Now create an ordinal variable for the MJ dosage
mjData$dosage <- NA
mjData$dosage[which(mjData$marcat=='MJ Non-User')] <- 0
mjData$dosage[which(mjData$marcat=='MJ Occ User' & mjData$mjpastyr=='')] <- 1
mjData$dosage[which(mjData$mjpastyr=="Less than once a month")] <- 2
mjData$dosage[which(mjData$mjpastyr=="About once a month")] <- 3
mjData$dosage[which(mjData$mjpastyr=="2-3 times a month")] <- 4
mjData$dosage[which(mjData$mjpastyr=="1-2 times a week")] <- 5
mjData$dosage[which(mjData$mjpastyr=="3-4 times a week")] <- 6
mjData$dosage[which(mjData$mjpastyr=="Everyday or nearly every day")] <- 7
mjData$mjLabel <- NA
mjData$mjLabel[which(mjData$marcat=="MJ Non-User")] <- "NonUser"
mjData$mjLabel[which(mjData$marcat=="MJ Occ User")] <- "User"
mjData$mjLabel[which(mjData$marcat=="MJ Freq User")] <- "FreqUser"
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
all.data <- merge(all.data, psData, by=c('bblid', 'scanid'), suffixes=c(".all", ""))

# Now add the clinical factor scores
fac.data <- read.csv('/data/joy/BBL/studies/pnc/n9498_dataFreeze/clinical/n9498_goassess_itemwise_bifactor_scores_age_regressed_20170131.csv')
all.data <- merge(all.data, fac.data, by='bblid', suffixes=c(".all", ""))

# Now also add our fs ct values
fs.values <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_freesurferCt_20180213.csv')
all.data <- merge(all.data, fs.values)
all.data$fsAvgCT <- (all.data$LThickness + all.data$RThickness)/2

# Now do FS volume
fs.values <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_freesurferAsegVol_20161220.csv')
all.data <- merge(all.data, fs.values)

# Now remove all subjects that have no marcat variable, are older than 14, do not pass structural QA, and endorsed fake drugs
all.data <- all.data[which((all.data$ageAtScan1/12)>=14),]
all.data <- all.data[-which(all.data$marcat==''),]
#all.data <- all.data[-which(all.data$averageManualRating==0),]
all.data <- all.data[which(all.data$t1Exclude==0),]
all.data <- all.data[-which(all.data$bblid %in% fakeData$bblid),]

## Now save our RDS file
saveRDS(object=all.data, file="mjAnovaDataWithDosage1.RDS")
all.data <- all.data[-which(all.data$dosage==1),]
saveRDS(object=all.data, file="mjAnovaData.RDS")

## Now prepare table 1 our demographics and what not
percentFemale <- table(all.data$sex, all.data$marcat)[2,]/table(all.data$marcat)
ageVals <- summarySE(data=all.data, measurevar='ageAtScan1', groupvars=c('marcat'))
# Now do our group comparisons
t.val.oc.vs.non <- t.test(ageAtScan1 ~ marcat, data=all.data[-which(all.data$marcat=="MJ Freq User"),])
t.val.oc.vs.freq <- t.test(ageAtScan1 ~ marcat, data=all.data[-which(all.data$marcat=="MJ Non-User"),])
t.val.non.vs.freq <- t.test(ageAtScan1 ~ marcat, data=all.data[-which(all.data$marcat=="MJ Occ User"),])
aov.val <- aov(ageAtScan1 ~ marcat, data=all.data)
## Now prepare an output row with these data
ageRow <- t(c(ageVals[2,3:4]/12, ageVals[3,3:4]/12, ageVals[1,3:4]/12, t.val.oc.vs.non$statistic, t.val.oc.vs.non$p.value, t.val.non.vs.freq$statistic, t.val.non.vs.freq$p.value, t.val.oc.vs.freq$statistic, t.val.oc.vs.freq$p.value))

## Now do our proportions and export the n val row
prop.val.oc.vs.non <- prop.test(table(all.data$marcat, all.data$sex)[c(4,3),])
prop.val.oc.vs.freq <- prop.test(table(all.data$marcat, all.data$sex)[c(2,4),])
prop.val.freq.vs.non <- prop.test(table(all.data$marcat, all.data$sex)[c(2,3),])
prop.val.all <- prop.test(table(all.data$marcat, all.data$sex)[c(2,3,4),])
n.row <- c(ageVals[2,2], percentFemale[3], ageVals[3,2], percentFemale[4], ageVals[1,2], percentFemale[2], prop.val.oc.vs.non$statistic, prop.val.oc.vs.non$p.value, prop.val.oc.vs.freq$statistic, prop.val.oc.vs.freq$p.value, prop.val.freq.vs.non$statistic, prop.val.freq.vs.non$p.value)

## Now report race
out.table <- table(all.data$race2, all.data$marcat)
## Now test for race proportion differences
all.data$raceTest <- 1
all.data$raceTest[which(all.data$race2!=1)] <- 2
race.prop.test <- chisq.test(x=all.data$raceTest, y=all.data$marcat)

## Now do the wrat values
wratVals <- summarySE(data=all.data, groupvars='marcat', measurevar='wrat4crstd', na.rm=T)
t.val.oc.vs.non <- t.test(wrat4crstd ~ marcat, data=all.data[-which(all.data$marcat=="MJ Freq User"),])
t.val.oc.vs.freq <- t.test(wrat4crstd ~ marcat, data=all.data[-which(all.data$marcat=="MJ Non-User"),])
t.val.non.vs.freq <- t.test(wrat4crstd ~ marcat, data=all.data[-which(all.data$marcat=="MJ Occ User"),])
aov.val <- aov(wrat4crstd ~ marcat, data=all.data)
wrat.row <- t(c(wratVals[2,3:4]/12, wratVals[3,3:4]/12, wratVals[1,3:4]/12, t.val.oc.vs.non$statistic, t.val.oc.vs.non$p.value, t.val.non.vs.freq$statistic, t.val.non.vs.freq$p.value, t.val.oc.vs.freq$statistic, t.val.oc.vs.freq$p.value))

## Now onto age at first use
firstUse <- summarySE(data=all.data, measurevar='mj_firstuse', groupvars='marcat', na.rm=T)
t.val.oc.vs.freq <- t.test(mj_firstuse ~ marcat, data=all.data[-which(all.data$marcat=="MJ Non-User"),])

## Now lets do dosage and what not
dosagePercent <- table(all.data$dosage, all.data$marcat)

## Now report alcohol values
alc.table <- table(all.data$marcat, all.data$cnb_substance_alc_0501)
# now test our alcohol prop differences
all.data$alcTest <- 1
all.data$alcTest[which(all.data$cnb_substance_alc_0501<5)] <- 3
all.data$alcTest[which(all.data$cnb_substance_alc_0501>=5)] <- 2
alc.prop.test <- chisq.test(x=all.data$alcTest, y=all.data$marcat)
alc.prop.test.occ.vs.non <- chisq.test(x=all.data$alcTest[-which(all.data$marcat=='MJ Freq User')], y=all.data$marcat[-which(all.data$marcat=='MJ Freq User')])
alc.prop.test.freq.vs.non <- chisq.test(x=all.data$alcTest[-which(all.data$marcat=='MJ Occ User')], y=all.data$marcat[-which(all.data$marcat=='MJ Occ User')])
alc.prop.test.freq.vs.occ <- chisq.test(x=all.data$alcTest[-which(all.data$marcat=='MJ Non-User')], y=all.data$marcat[-which(all.data$marcat=='MJ Non-User')])

## Now onto our pathology factor scores
oapVals <- summarySE(data=all.data, groupvars='marcat', measurevar='overall_psychopathology_ar_4fact', na.rm=T)
t.val.oc.vs.non <- t.test(overall_psychopathology_ar_4fact ~ marcat, data=all.data[-which(all.data$marcat=="MJ Freq User"),])
t.val.oc.vs.freq <- t.test(overall_psychopathology_ar_4fact ~ marcat, data=all.data[-which(all.data$marcat=="MJ Non-User"),])
t.val.non.vs.freq <- t.test(overall_psychopathology_ar_4fact ~ marcat, data=all.data[-which(all.data$marcat=="MJ Occ User"),])
aov.val <- aov(overall_psychopathology_ar_4fact ~ marcat, data=all.data)

mapVals <- summarySE(data=all.data, groupvars='marcat', measurevar='mood_ar_4factor', na.rm=T)
t.val.oc.vs.non <- t.test(mood_ar_4factor ~ marcat, data=all.data[-which(all.data$marcat=="MJ Freq User"),])
t.val.oc.vs.freq <- t.test(mood_ar_4factor ~ marcat, data=all.data[-which(all.data$marcat=="MJ Non-User"),])
t.val.non.vs.freq <- t.test(mood_ar_4factor ~ marcat, data=all.data[-which(all.data$marcat=="MJ Occ User"),])
aov.val <- aov(mood_ar_4factor ~ marcat, data=all.data)

expVals <- summarySE(data=all.data, groupvars='marcat', measurevar='externalizing_ar_4factor', na.rm=T)
t.val.oc.vs.non <- t.test(externalizing_ar_4factor ~ marcat, data=all.data[-which(all.data$marcat=="MJ Freq User"),])
t.val.oc.vs.freq <- t.test(externalizing_ar_4factor ~ marcat, data=all.data[-which(all.data$marcat=="MJ Non-User"),])
t.val.non.vs.freq <- t.test(externalizing_ar_4factor ~ marcat, data=all.data[-which(all.data$marcat=="MJ Occ User"),])
aov.val <- aov(externalizing_ar_4factor ~ marcat, data=all.data)

pspVals <- summarySE(data=all.data, groupvars='marcat', measurevar='psychosis_ar_4factor', na.rm=T)
t.val.oc.vs.non <- t.test(psychosis_ar_4factor ~ marcat, data=all.data[-which(all.data$marcat=="MJ Freq User"),])
t.val.oc.vs.freq <- t.test(psychosis_ar_4factor ~ marcat, data=all.data[-which(all.data$marcat=="MJ Non-User"),])
t.val.non.vs.freq <- t.test(psychosis_ar_4factor ~ marcat, data=all.data[-which(all.data$marcat=="MJ Occ User"),])
aov.val <- aov(psychosis_ar_4factor ~ marcat, data=all.data)


## Now return an age, and sex matched sample
all.data$usageBin <- 0
all.data$usageBin[all.data$dosage>0] <- 1
tmpDat <- all.data[c('bblid', 'scanid', 'usageBin', 'ageAtScan1', 'sex')]
mod <- matchit(usageBin ~ ageAtScan1 + sex, data=tmpDat, ratio=1, na.action=na.omit)
all.data.out <- all.data[as.vector(mod$match.matrix),]
all.data.out <- rbind(all.data.out, all.data[which(all.data$usageBin==1),])
saveRDS(object=all.data.out, file="mjAnovaDataMatch1.RDS")
all.data$usageBin <- 0
all.data$usageBin[all.data$dosage>0] <- 1
tmpDat <- all.data[c('bblid', 'scanid', 'usageBin', 'ageAtScan1', 'sex')]
mod <- matchit(usageBin ~ ageAtScan1 + sex, data=tmpDat, ratio=2, na.action=na.omit)
all.data.out <- all.data[as.vector(mod$match.matrix),]
all.data.out <- rbind(all.data.out, all.data[which(all.data$usageBin==1),])
saveRDS(object=all.data.out, file="mjAnovaDataMatch2.RDS")
all.data$usageBin <- 0
all.data$usageBin[all.data$dosage>0] <- 1
tmpDat <- all.data[c('bblid', 'scanid', 'usageBin', 'ageAtScan1', 'sex')]
mod <- matchit(usageBin ~ ageAtScan1 + sex, data=tmpDat, ratio=3, na.action=na.omit)
all.data.out <- all.data[as.vector(mod$match.matrix),]
all.data.out <- rbind(all.data.out, all.data[which(all.data$usageBin==1),])
saveRDS(object=all.data.out, file="mjAnovaDataMatch3.RDS")
