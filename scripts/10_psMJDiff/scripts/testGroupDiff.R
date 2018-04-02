## AFGR March 2018

## This script is going to be used to grab t stats from our marcat groups 
## isolated to the PS groups

## First load libbrary(s)
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
outputVals <- matrix(NA, 1026, 3)
rowCheck <- 1
for(i in loopVals){
  formTemp <- as.formula(paste(i, "~age + age2 + age3 + sex + race2 + marcat"))
  tmpMod <- lm(formTemp, data=psOut)
  outputRow <- c(i,as.numeric(summary(tmpMod)$coefficients[7,3]), as.numeric(summary(tmpMod)$coefficients[7,4]))
  outputVals[rowCheck,] <- outputRow
  rowCheck <- rowCheck+1
}

# Now perform FWE correction
addCol <- rep(NA, 1026, 3)
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
outputVals2 <- matrix(NA, 1026, 3)
rowCheck <- 1
for(i in loopVals){
  formTemp <- as.formula(paste(i, "~age + age2 + age3 + sex + race2 + marcat"))
  tmpMod <- lm(formTemp, data=psOut)
  outputRow <- c(i,as.numeric(summary(tmpMod)$coefficients[7,3]), as.numeric(summary(tmpMod)$coefficients[7,4]))
  outputVals2[rowCheck,] <- outputRow
  rowCheck <- rowCheck+1
}

# Now perform FWE correction
addCol <- rep(NA, 1026, 3)
grepVals <- c('vol', 'gmd', 'ct', 'cortcon', 'cbf', 'tr', 'jlf_fa', 'alff', 'reho', 'jhulabel', 'Ct', 'Ravens')
for(gV in grepVals){
  pVals <- as.numeric(outputVals2[grep(gV, outputVals[,1]),3])
  addCol[grep(gV, outputVals[,1])] <- p.adjust(pVals, method='fdr')  
}
outputVals2 <- cbind(outputVals2,addCol)

## Now we need to write the output
toWrite1 <- outputVals[which(outputVals[,3]<.05),]
write.csv(toWrite1, "nonVsUser.csv", quote=F, row.names=F)
toWrite2 <- outputVals2[which(outputVals2[,3]<.05),]
write.csv(toWrite2, "nonVsFreq.csv", quote=F, row.names=F)


## Now explore the same analyses in a matched sample
source('/home/arosen/adroseHelperScripts/R/afgrHelpFunc.R')

