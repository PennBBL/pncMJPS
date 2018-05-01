## AFGR March 2018

## This script is going to be used to grab t stats from our marcat groups 
## isolated to the PS groups

## First load library(s)
source('/data/joy/BBL/projects/pncMJPS/scripts/01_dataPrep/scripts/dataPrepNoFake.R')
source('/home/arosen/adroseHelperScripts/R/afgrHelpFunc.R')

## First lets explore our demographics for all subbjects
ageVals <- summarySE(data=psOut, measurevar='ageAtScan1', groupvars='dosage')
raceVals <- table(psOut$race2, psOut$dosage)
sexVals <- table(psOut$sex, psOut$dosage)

## Now write these guys
write.csv(ageVals, "ageValues.csv", quote=F, row.names=F)
write.csv(raceVals, "raceValues.csv", quote=F, row.names=F)
write.csv(sexVals, "sexValues.csv", quote=F, row.names=F)

## Now make a age2 and age3 value
psOut$age <- scale(psOut$ageAtScan1)
psOut$age2 <- scale(psOut$ageAtScan1)^2
psOut$age3 <- scale(psOut$ageAtScan1)^3

## Now I am going to loop through every value and grab a t stat for group differences 
## Between our non users and occasional users
loopVals <- names(psOut)[grep("jlf", names(psOut))]
loopVals <- append(loopVals, names(psOut)[grep("dtitk_jhulabel_", names(psOut))])
loopVals <- append(loopVals, names(psOut)[grep("Nmf", names(psOut))])
psOut <- psOut[-which(psOut$marcat=="MJ Frequent User"),]
outputVals <- matrix(NA, length(loopVals), 3)
rowCheck <- 1
for(i in loopVals){
  formTemp <- as.formula(paste(i, "~age + age2 + age3 + sex + marcat + race2"))
  tmpMod <- lm(formTemp, data=psOut)
  outputRow <- c(i,as.numeric(summary(tmpMod)$coefficients['marcatMJ User',3]), as.numeric(summary(tmpMod)$coefficients['marcatMJ User',4]))
  outputVals[rowCheck,] <- outputRow
  rowCheck <- rowCheck+1
}

# Now perform FWE correction
addCol <- rep(NA, length(loopVals), 3)
grepVals <- c('vol', 'gmd', 'ct', 'cortcon', 'cbf', 'tr', 'jlf_fa', 'alff', 'reho', 'jhulabel', 'Ct', 'Ravens')
for(gV in grepVals){
  pVals <- as.numeric(outputVals[grep(gV, outputVals[,1]),3])
  addCol[grep(gV, outputVals[,1])] <- p.adjust(pVals, method='fdr')  
}
outputVals <- cbind(outputVals,addCol)

# Now do frequent vs non
source('/data/joy/BBL/projects/pncMJPS/scripts/01_dataPrep/scripts/dataPrepNoFake.R')
psOut$age <- scale(psOut$ageAtScan1)
psOut$age2 <- scale(psOut$ageAtScan1)^2
psOut$age3 <- scale(psOut$ageAtScan1)^3
psOut <- psOut[-which(psOut$marcat=="MJ User"),]
outputVals2 <- matrix(NA, length(loopVals), 3)
rowCheck <- 1
for(i in loopVals){
  formTemp <- as.formula(paste(i, "~age + age2 + age3 + sex + marcat"))
  tmpMod <- lm(formTemp, data=psOut)
  outputRow <- c(i,as.numeric(summary(tmpMod)$coefficients['marcatMJ Non-User',3]), as.numeric(summary(tmpMod)$coefficients['marcatMJ Non-User',4]))
  outputVals2[rowCheck,] <- outputRow
  rowCheck <- rowCheck+1
}

# Now perform FWE correction
addCol <- rep(NA, length(loopVals), 3)
grepVals <- c('vol', 'gmd', 'ct', 'cortcon', 'cbf', 'tr', 'jlf_fa', 'alff', 'reho', 'jhulabel', 'Ct', 'Ravens')
for(gV in grepVals){
  pVals <- as.numeric(outputVals2[grep(gV, outputVals[,1]),3])
  addCol[grep(gV, outputVals[,1])] <- p.adjust(pVals, method='fdr')  
}
outputVals2 <- cbind(outputVals2,addCol)

## Now do this analysis in the regions that are deemed as "important" by hi lo standards
cog.names <- c("F1_Exec_Comp_Cog_AccuracySummedES.csv")#,
#"F1_Slow_SpeedSummedES.csv",
#"F1_Social_Cognition_EfficiencySummedES.csv",
#"F2_Complex_Reasoning_EfficiencySummedES.csv",
#"F2_Memory_SpeedSummedES.csv",
#"F2_Social_Cog_AccuracySummedES.csv",
#"F3_Fast_SpeedSummedES.csv",
#"F3_Memory_AccuracySummedES.csv",
#"F3_Memory_EfficiencySummedES.csv",
#"F4_Executive_EfficiencySummedES.csv")
modal.names <- c('mprage_jlf_vol_', 'mprage_jlf_gmd_', 'mprage_jlf_ct_', 'mprage_jlf_cortcon_', 'pcasl_jlf_cbf_', 'dti_jlf_tr_', 'rest_jlf_alff_', 'rest_jlf_reho_')
# Now we want to loop through each of the cognitive metrics and read the csv, and find the ROI's that come out
for(x in c(.7,.8,.9,1,1.1,1.2)){
for(c in cog.names){
  tmpVals <- read.csv(c)
  tmpVals$sexAvg <- rowSums(tmpVals[,4:5])/2
  importReg <- tmpVals[which(tmpVals$sexAvg>x),]
  importReg$ROI_readable <- as.character(importReg$ROI_readable)
  print(dim(importReg))
  # Now go through each of our imaging modalities, and run FDR correction in the 
  # regions that are deemed as "important"
  outCol1 <- rep(NA, length(loopVals))
  outCol2 <- rep(NA, length(loopVals))
  for(m in modal.names){
    nameVals <- gsub(x=importReg$ROI_readable, pattern='mprage_jlf_vol_', replacement=m)
    nameVals <- append(nameVals, paste(m, "R_Pallidum", sep=''))
    index <- which(outputVals[,1] %in% nameVals)
    # Now apply FWE correction w/in those specific regions
    #hist(as.numeric(outputVals[index,3]))
    outCol1[index] <- p.adjust(as.numeric(outputVals[index,3]), method='fdr')
    #hist( p.adjust(as.numeric(outputVals[index,3]), method='fdr'))
    outCol2[index] <- p.adjust(as.numeric(outputVals2[index,3]), method='fdr')
  }
  outputVals <- cbind(outputVals, outCol1)
  outputVals2 <- cbind(outputVals2, outCol2)
}
}
## Now we need to write the output
toWrite1 <- outputVals[which(outputVals[,3]<.05),]
write.csv(toWrite1, "nonVsUser.csv", quote=F, row.names=F)
toWrite2 <- outputVals2[which(outputVals2[,3]<.05),]
write.csv(toWrite2, "nonVsFreq.csv", quote=F, row.names=F)
