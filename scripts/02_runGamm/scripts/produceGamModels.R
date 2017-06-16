# AFGR June 2017
# This script is going to be used to run the GAM models for the PS-MJ study
# The steps include:
#	1.) Running GAM's for each modality for each ROI
#		a. The model will look like:"brainData ~ s(age) + PS*MJ + sex + dataQuality"
# 	2.) Reporting significant intgeractions from PS & MJ

# Load all data and library(s)
source('/data/joy/BBL/projects/pncMJPS/scripts/01_dataPrep/scripts/dataPrepNoFake.R')
source('/data/joy/BBL/projects/pncMJPS/scripts/02_runGamm/functions/functions.R')
install_load('mgcv')

# Now run the gam's for the structural data 
# Collapse all fo the Users into MJ User bin for marcat 
strucData$marcat[strucData$marcat=='MJ Frequent User'] <- 'MJ User'
strucData <- strucData[-which(strucData$marcat==levels(strucData$marcat)[1]),]
strucData <- strucData[-which(strucData$goassessDxpmr7==levels(strucData$goassessDxpmr7)[1]),]
# Now find the number of significant interactions for each group
sigVol <- runGamModel(strucData, 'mprage_jlf_vol', 'averageManualRating')
sigVolN <- length(which(p.adjust(sigVol[,2], method='fdr')<.05))
sigCT <- runGamModel(strucData, 'mprage_jlf_ct', 'averageManualRating')
sigCTN <- length(which(p.adjust(sigCT[,2], method='fdr')<.05))
sigGMD <- runGamModel(strucData, 'mprage_jlf_gmd', 'averageManualRating')
sigGMDN <- length(which(p.adjust(sigGMD[,2], method='fdr')<.05))

# Now run CBF
cbfData$marcat[cbfData$marcat=='MJ Frequent User'] <- 'MJ User'
sigCBF <- runGamModel(cbfData, 'pcasl_jlf_cbf', 'pcaslTSNR')
sigCBFN <- length(which(p.adjust(sigCBF[,2], method='fdr')<.05))

# Now run the rest data
# These don't run becvause we do not have any PS-Users
restData$marcat[restData$marcat=='MJ Frequent User'] <- 'MJ User'
sigReho <- runGamModel(restData, 'rest_jlf_reho', 'restRelMeanRMSMotion')
sigRehoN <- length(which(p.adjust(sigReho[,2], method='fdr')<.05))
sigAlff <- runGamModel(restData, 'rest_jlf_alff', 'restRelMeanRMSMotion')
sigAlffN <- length(which(p.adjust(sigAlff[,2], method='fdr')<.05))

# Now run the DTI data
sigDTI <- runGamModel(dtiData, 'dti_jlf_tr', 'dti64Tsnr')
sigDTIN <- length(which(p.adjust(sigDTI[,2], method='fdr')<.05))
