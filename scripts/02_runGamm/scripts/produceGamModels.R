# AFGR June 2017
# This script is going to be used to run the GAM models for the PS-MJ study
# The steps include:
#	1.) Running GAM's for each modality for each ROI
#		a. The model will look like:"brainData ~ s(age) + PS*MJ + sex + dataQuality"
# 	2.) Reporting significant intgeractions from PS & MJ

# Load all data and library(s)
source('/data/joy/BBL/projects/pncMJPS/scripts/01_dataPrep/scripts/dataPrepNoFake.R')
source('/data/joy/BBL/projects/pncMJPS/scripts/02_runGamm/functions/functions.R')
install_load('mgcv', 'voxel', 'ggplot2')

# Now run the gam's for the structural data 
# Collapse all fo the Users into MJ User bin for marcat 
strucData$marcat[strucData$marcat=='MJ Frequent User'] <- 'MJ User'
#strucData <- strucData[-which(strucData$marcat=='MJ User'),]
strucData <- strucData[-which(strucData$marcat==levels(strucData$marcat)[1]),]
strucData <- strucData[-which(strucData$goassessDxpmr7==levels(strucData$goassessDxpmr7)[1]),]
# Now find the number of significant interactions for each group
sigVol <- runMainEffect(strucData, 'mprage_jlf_vol', 'averageManualRating')
sigVolN <- length(which(p.adjust(sigVol[,2], method='fdr')<.05))
plotMainEffects(pValueInfo=sigVol, dataFrame=strucData, pdfName='volMainMarcat.pdf', QC='averageManualRating')
sigCT <- runMainEffect(strucData, 'mprage_jlf_ct', 'averageManualRating')
sigCTN <- length(which(p.adjust(sigCT[,2], method='fdr')<.05))
plotMainEffects(pValueInfo=sigCT, dataFrame=strucData, pdfName='ctMainMarcat.pdf', QC='averageManualRating')
sigGMD <- runMainEffect(strucData, 'mprage_jlf_gmd', 'averageManualRating')
sigGMDN <- length(which(p.adjust(sigGMD[,2], method='fdr')<.05))
plotMainEffects(pValueInfo=sigGMD, dataFrame=strucData, pdfName='gmdMainMarcat.pdf', QC='averageManualRating')

# Now do the sttah parc
stathParc$marcat[stathParc$marcat=='MJ Frequent User'] <- 'MJ User'
stathParc <- stathParc[-which(stathParc$marcat==levels(stathParc$marcat)[1]),]
stathParc <- stathParc[-which(stathParc$goassessDxpmr7==levels(stathParc$goassessDxpmr7)[1]),]
stathParc <- merge(strucData, stathParc, by=c('bblid', 'scanid'))
colnames(stathParc) <- gsub(x=colnames(stathParc), pattern='.y', replacement='')
stathGmd <- runMainEffect(stathParc, 'NZMean_', 'averageManualRating')
stathGMDN <- length(which(p.adjust(stathGmd[,2], method='fdr')<.05))
stathGMDN <-  which(stathGmd[,2]<.05)
plotMainEffects(pValueInfo=stathGmd, dataFrame=stathParc, pdfName='gmdStathMainMarcat.pdf', QC='averageManualRating')
# Now produce the nominally significant ROI's
tmpVals <- unname(gsub(x=stathGmd[stathGMDN,1], pattern="NZMean_", replacement=''))


# Now do stath CT 
stathCT$marcat[stathCT$marcat=='MJ Frequent User'] <- 'MJ User'
stathCT <- stathCT[-which(stathCT$marcat==levels(stathCT$marcat)[1]),]
stathCT <- stathCT[-which(stathCT$goassessDxpmr7==levels(stathCT$goassessDxpmr7)[1]),]
stathCT <- merge(strucData, stathCT, by=intersect(names(stathCT),names(strucData)))
colnames(stathCT) <- gsub(x=colnames(stathCT), pattern='.y', replacement='')
sigStathCT <- runMainEffect(stathCT, 'NZMean_', 'averageManualRating')

# Now run CBF
cbfData$marcat[cbfData$marcat=="MJ Frequent User"] <- 'MJ User'
#cbfData <- cbfData[-which(cbfData$marcat=='MJ User'),]
cbfData <- cbfData[-which(cbfData$marcat==levels(cbfData$marcat)[1]),]
cbfData <- cbfData[-which(cbfData$goassessDxpmr7==levels(cbfData$goassessDxpmr7)[1]),]
#cbfData$marcat[cbfData$marcat=='MJ Frequent User'] <- 'MJ User'
sigCBF <- runMainEffect(cbfData, 'pcasl_jlf_cbf', 'pcaslTSNR')
sigCBFN <- length(which(p.adjust(sigCBF[,2], method='fdr')<.05))
plotMainEffects(pValueInfo=sigCBF, dataFrame=cbfData, pdfName='cbfMainMarcat.pdf', QC='pcaslTSNR')

# Now run the rest data
# These don't run becvause we do not have any PS-Users
restData <- restData[-which(restData$marcat==levels(restData$marcat)[1]),]
#restData <- restData[-which(restData$goassessDxpmr7==levels(restData$goassessDxpmr7)[1]),]
restData$marcat[restData$marcat=='MJ Frequent User'] <- 'MJ User'
sigReho <- runMainEffect(restData, 'rest_jlf_reho', 'restRelMeanRMSMotion')
sigRehoN <- length(which(p.adjust(sigReho[,2], method='fdr')<.05))
sigAlff <- runMainEffect(restData, 'rest_jlf_alff', 'restRelMeanRMSMotion')
sigAlffN <- length(which(p.adjust(sigAlff[,2], method='fdr')<.05))

# Now run the DTI data
sigDTI <- runMainEffect(dtiData, 'dti_jlf_tr', 'dti64Tsnr')
sigDTIN <- length(which(p.adjust(sigDTI[,2], method='fdr')<.05))

# Now assess the non linear impact of MJ dosgae on the PS group
valsToLoop <- grep('mprage_jlf_gmd', names(strucData))
strucDataNon <- strucData[which(strucData$goassessDxpmr7=='PS'),]
strucDataNon <- strucDataNon[complete.cases(strucDataNon$dosage),]
pdf('nonLinearRelationshipInGMD.pdf')
for(i in valsToLoop){
  formulaVal <- as.formula(paste(names(strucDataNon)[i], '~ s(ageAtScan1) + s(dosage, k=7) + sex + averageManualRating'))
  mod1 <- gam(formulaVal, data=strucDataNon)
  if(summary(mod1)$s.table[2,4] < .05){
    print(summary(mod1))
    tmp <- plotGAM(gamFit=mod1, smooth.cov='dosage')
    print(tmp)
  }
}
dev.off()

valsToLoop <- grep('mprage_jlf_ct', names(strucData))
strucDataNon <- strucData[which(strucData$goassessDxpmr7=='PS'),]
strucDataNon <- strucDataNon[complete.cases(strucDataNon$dosage),]
pdf('nonLinearRelationshipInCT.pdf')
for(i in valsToLoop){
  formulaVal <- as.formula(paste(names(strucDataNon)[i], '~ s(ageAtScan1) + s(dosage, k=7) + sex + averageManualRating'))
  mod1 <- gam(formulaVal, data=strucDataNon)
  if(summary(mod1)$s.table[2,4] < .05){
    print(summary(mod1))
    tmp <- plotGAM(gamFit=mod1, smooth.cov='dosage')
    print(tmp)
  }
}
dev.off()
