source('/home/arosen/adroseHelperScripts/R/afgrHelpFunc.R')
install_load('mgcv', 'voxel')

# Load data
source('/data/joy/BBL/projects/pncMJPS/scripts/01_dataPrep/scripts/dataPrepNoFake.R')
cogData <- read.csv('/data/joy/BBL/projects/pncMJPS/data/PS_Cannabis/Marijuana.csv')
mjData <- read.csv('/data/joy/BBL/projects/pncMJPS/data/n9462_mj_ps_cnb_fortmm.csv')
goassessData <- read.csv('/data/joy/BBL/studies/pnc/n9498_dataFreeze/clinical/n9498_diagnosis_dxpmr_20161014.csv')
# Now create an ordinal variable for the MJ dosage
mjData$mjBinary <- mjData$marcat 
mjData$mjBinary[mjData$mjBinary=='MJ Frequent User'] <- 'MJ User'
mjData$dosage <- NA
mjData$dosage[which(mjData$marcat=='MJ Non-User')] <- 0
mjData$dosage[which(mjData$marcat=='MJ User' & mjData$mjpastyr=='')] <- 1
mjData$dosage[which(mjData$mjpastyr=="Less than once a month")] <- 2
mjData$dosage[which(mjData$mjpastyr=="About once a month")] <- 3
mjData$dosage[which(mjData$mjpastyr=="2-3 times a month")] <- 4
mjData$dosage[which(mjData$mjpastyr=="1-2 times a week")] <- 5
mjData$dosage[which(mjData$mjpastyr=="3-4 times a week")] <- 6

# Now get the imaging sample
bblidIndex <- strucData$bblid

# Now create our data frame
allData <- merge(mjData, cogData, by=intersect(names(mjData), names(cogData)))
allData <- merge(allData, goassessData, by=intersect(names(allData), names(goassessData)))

allData$mjBinary <- as.character(allData$mjBinary)
allData$mjBinary <- as.factor(allData$mjBinary)
allData$goassessDxpmr6 <- ordered(allData$goassessDxpmr6)

# Now rm the subjects that have not used any within the last year
allDataAll <- allData
allData <- allData[which(allData$dosage>1),]
allData <- allData[allData$bblid %in% bblidIndex,]
allData$race <- as.factor(allData$race)

# Now produce non linear effects across sex
allData <- allData[which(allData$goassessDxpmr6=='TD' | allData$goassessDxpmr6=='PS'),]
variableValues <- c(seq(21,33), seq(148,159))
pdf('cogDataAcrossSex.pdf')
for(i in variableValues){
  # Grab our names and make our formula
  nameValue <- names(allData)[i]
  formulaValue <- as.formula(paste(nameValue, '~ s(ageAtGo1Cnb) + s(dosage, k=4, by=goassessDxpmr6) + race + goassessDxpmr6 + s(dosage,k=4) + Sex22'))
  print(formulaValue)
  mod1 <- gam(formulaValue, data=allData)
  if(summary(mod1)$s.table[3,4] < .05){
    print(summary(mod1))
    tmp <- plotGAM(gamFit=mod1, smooth.cov='dosage', groupCovs=c('goassessDxpmr6'))
    print(tmp)
  }
}
dev.off()

# Now run through our variables that we want to check for non linear effects
allDataFreeze <- allData
allData <- allDataFreeze[which(allDataFreeze$Sex22=='Male'),]
allData <- allData[which(allData$goassessDxpmr6=='TD' | allData$goassessDxpmr6=='PS'),]

pdf('sigNonLinearDosageEffectsMaleOnlyLastYear.pdf')
for(i in variableValues){
  # Grab our names and make our formula
  nameValue <- names(allData)[i]
  formulaValue <- as.formula(paste(nameValue, '~ s(ageAtGo1Cnb) + s(dosage, k=4, by=goassessDxpmr6) + race + goassessDxpmr6 + s(dosage,k=4)'))
  print(formulaValue)
  mod1 <- gam(formulaValue, data=allData)
  if(summary(mod1)$s.table[3,4] < .05){
    print(summary(mod1))
    tmp <- plotGAM(gamFit=mod1, smooth.cov='dosage', groupCovs=c('goassessDxpmr6'))
    print(tmp)
  }
}
dev.off()

# Now do female
allData <- allDataFreeze[which(allDataFreeze$Sex22=='Female'),]
allData <- allData[which(allData$goassessDxpmr6=='TD' | allData$goassessDxpmr6=='PS'),]
pdf('sigNonLinearDosageEffectsFemaleOnlyLastYear.pdf')
for(i in variableValues){
  # Grab our names and make our formula
  nameValue <- names(allData)[i]
  formulaValue <- as.formula(paste(nameValue, '~ s(ageAtGo1Cnb) + s(dosage, k=4, by=goassessDxpmr6) + race + goassessDxpmr6 + s(dosage,k=4)'))
  print(formulaValue)
  mod1 <- gam(formulaValue, data=allData)
  if(summary(mod1)$s.table[2,4] < .05){
    print(summary(mod1))
    tmp <- plotGAM(gamFit=mod1, smooth.cov='dosage', groupCovs=c('goassessDxpmr6'))
    print(tmp)
  }
}
dev.off()

# Now look at PS in isolation
allData <- allDataFreeze[which(allDataFreeze$Sex22=='Male'),]
allData <- allData[which(allData$goassessDxpmr6=='PS'),]
pdf('sigNonLinearDosageEffectsMaleIsolatedPSOnlyLastYear.pdf')
for(i in variableValues){
  # Grab our names and make our formula
  nameValue <- names(allData)[i]
  formulaValue <- as.formula(paste(nameValue, '~ s(ageAtGo1Cnb) + s(dosage, k=4)'))
  print(formulaValue)
  mod1 <- gam(formulaValue, data=allData)
  if(summary(mod1)$s.table[2,4] < .05){
    tmp <- plotGAM(gamFit=mod1, smooth.cov='dosage')
    print(tmp)
  }
}
dev.off()

allData <- allDataFreeze[which(allDataFreeze$Sex22=='Female'),]
allData <- allData[which(allData$goassessDxpmr6=='PS'),]
pdf('sigNonLinearDosageEffectsFemaleIsolatedPSOnlyLastYear.pdf')
for(i in variableValues){
  # Grab our names and make our formula
  nameValue <- names(allData)[i]
  formulaValue <- as.formula(paste(nameValue, '~ s(ageAtGo1Cnb) + s(dosage, k=4)'))
  print(formulaValue)
  mod1 <- gam(formulaValue, data=allData)
  if(summary(mod1)$s.table[2,4] < .05){
    tmp <- plotGAM(gamFit=mod1, smooth.cov='dosage')
    print(tmp)
  }
}
dev.off()
