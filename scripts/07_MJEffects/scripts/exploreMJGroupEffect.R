# Load all data and library(s)
cogValues <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/cnb/n1601_cnb_factor_scores_tymoore_20151006.csv')
source('/data/joy/BBL/projects/pncMJPS/scripts/01_dataPrep/scripts/dataPrepNoFake.R')
source('/data/joy/BBL/projects/pncMJPS/scripts/07_MJEffects/functions/functions.R')
install_load('mgcv', 'voxel', 'ggplot2')

# Merge our data values 
strucData <- merge(strucData, cogValues, by='bblid')
cbfData <- merge(cbfData, cogValues, by='bblid')
restData <- merge(restData, cogValues, by='bblid')

# Now run the gam's for the structural data 
# Collapse all fo the Users into MJ User bin for marcat 
StrucDataFreeze <- strucData
# First run all subjs who have used within the last year
strucData <- strucData[-which(strucData$dosage==1),]
strucData <- strucData[-which(strucData$marcat==levels(strucData$marcat)[1]),]
# Now find the number of significant interactions for each group
sigVol <- runGamModelG(strucData, 'mprage_jlf_vol', 'averageManualRating')
sigVolN <- sigVol[which(sigVol[,2]<.05),]
sigCT <- runGamModelG(strucData, 'mprage_jlf_ct', 'averageManualRating')
sigCTN <- sigCT[which(sigCT[,2]<.05),]
sigGMD <- runGamModelG(strucData, 'mprage_jlf_gmd', 'averageManualRating')
sigGMDN <- sigGMD[which(sigGMD[,2]<.05),]
sigCC <- runGamModelG(strucData, 'mprage_jlf_cortcon', 'averageManualRating')
sigCCN <- sigCC[which(sigCC[,2]<.05),]

# Now find the corellation between these regions and our cognitive variables
outAll <- NULL
for(c in sigCTN[,1]){
  outC <- NULL
  for(corVal in 486:511){
    corStat <- cor(strucData[c], strucData[,corVal], use='complete')
    print(corStat)
    outC <- cbind(outC, corStat)
  }
  outAll <- rbind(outAll, outC)
}
colnames(outAll) <- names(strucData)[486:511]
outAllCT <- outAll

# Now run CBF
CbfDataFreeze <- cbfData
cbfData <- cbfData[-which(cbfData$dosage==1),]
cbfData <- cbfData[-which(cbfData$marcat==levels(cbfData$marcat)[1]),]
sigCBF <- runGamModelG(cbfData, 'pcasl_jlf_cbf', 'pcaslTSNR')
sigCBFN <- sigCBF[which(sigCBF[,2]<.05),]


# Now do rest
RestDataFreeze <- restData
restData <- restData[-which(restData$dosage==1),]
restData <- restData[-which(restData$marcat==levels(restData$marcat)[1]),]
sigReho <- runGamModelG(restData, 'rest_jlf_reho', 'restRelMeanRMSMotion')
sigRehoN <- sigReho[which(sigReho[,2]<.05),]
sigAlff <- runGamModelG(restData, 'rest_jlf_alff', 'restRelMeanRMSMotion')
sigAlffN <- sigAlff[which(sigAlff[,2]<.05),]

# Now do the cor values
outAll <- NULL
for(c in sigRehoN[,1]){
  outC <- NULL
  for(corVal in 279:304){
    corStat <- cor(restData[c], restData[,corVal], use='complete')
    print(corStat)
    outC <- cbind(outC, corStat)
  }
  outAll <- rbind(outAll, outC)
}
colnames(outAll) <- names(restData)[279:304]
outAllREHO <- outAll

# Now do alff
# Now do the cor values
outAll <- NULL
for(c in sigAlffN[,1]){
  outC <- NULL
  for(corVal in 279:304){
    corStat <- cor(restData[c], restData[,corVal], use='complete')
    print(corStat)
    outC <- cbind(outC, corStat)
  }
  outAll <- rbind(outAll, outC)
}
colnames(outAll) <- names(restData)[279:304]
outAllALFF <- outAll

# Now do TR
sigTR <- runGamModelG(dtiData, 'dti_jlf_tr_', 'dti64Tsnr')
sigTRN <- sigTR[which(sigTR[,2]<.05),]

# Now do FA
sigFA <- runGamModelG(faData, 'dti_dtitk_jhutract_fa', 'dti64Tsnr')
sigFA <- rbind(sigFA, runGamModelG(faData, 'dti_jlf_fa_', 'dti64Tsnr'))
sigFAN <- sigFA[which(sigFA[,2]<.05),] 


# Now write the output
output <- rbind(sigVolN, sigCTN, sigGMDN, sigCBFN, sigRehoN, sigAlffN, sigTRN, sigFAN)
write.csv(output, 'nomSigValsGroup.csv', quote=F, row.names=F)

output <- rbind(outAllCT, outAllCBF, outAllREHO, outAllALFF)
write.csv(output, 'corValsFromNomSig.csv', quote=F, row.names=T)
