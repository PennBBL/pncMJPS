# AFGR June 2017
# This script will be used to prepare all of the image and demographic data for the MJ-PS project
# The various steps - roughly - for this script include:
#	1.) Loading all required CSV's
#	2.) combining the csv's 
#	3.) Applying any exlcusion criteria


## Load any functions / library(s)
source('/home/arosen/adroseHelperScripts/R/afgrHelpFunc.R')

## Load the required data
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
psData <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/clinical/n1601_diagnosis_dxpmr_20170509.csv')
psData[which(psData$goassessDxpmr7=='OP'),] <- "TD"
demoData <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_demographics_go1_20161212.csv')
psData <- merge(psData, demoData, by=c('bblid', 'scanid'))
volData <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionVol_20170412.csv')
gmdData <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAtroposIntersectionGMD_20170410.csv')
ctData <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionCT_20170331.csv')
ccData <- read.csv('/data/joy/BBL/studies/pnc/n2416_dataFreeze/neuroimaging/t1struct/n2416_jlfAntsCTIntersectionCortCon_20170814.csv')
t1QAData <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_t1QaData_20170306.csv')
cbfData <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/asl/n1601_jlfAntsCTIntersectionPcaslValues_20170403.csv')
cbfQAData <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/asl/n1601_PcaslQaData_20170403.csv')
trData <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/dti/n1601_jlfTRValues_20170411.csv')
trData2 <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/dti/n1601_jlfWmLobesTRValues_20170405.csv')
trData <- merge(trData, trData2)
rm(trData2)
dtiQAData <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/dti/n1601_dti_qa_20170301.csv')
alffData <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/rest/n1601_jlfALFFValues_20170714.csv')
rehoData <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/rest/n1601_jlfReHoValues_20170714.csv')
restQAData <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/rest/n1601_RestQAData_20170714.csv')
faData <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/dti/n1601_JHUTractFA_20170321.csv')
faData2 <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/dti/n1601_jlfWmLobesFAValues_20170405.csv')
faData <- merge(faData, faData2)
rm(faData2)

## Now apply all restrictions 
#Start with struc
strucData <- merge(volData, gmdData, by=c('bblid', 'scanid'))
strucData <- merge(strucData, ctData, by=c('bblid', 'scanid'))
strucData <- merge(strucData, t1QAData, by=c('bblid', 'scanid'))
strucData <- merge(strucData, ccData)
strucData <- strucData[which(strucData$averageManualRating!=0),]
strucData <- merge(strucData, mjData, by='bblid')
strucData <- merge(strucData, psData, by=c('bblid', 'scanid'))
bblidIndex <- strucData$bblid
bblidIndex <- bblidIndex[which(bblidIndex%in%fakeData$bblid=='FALSE')]
strucData <- strucData[-which(strucData$goassessDxpmr7==levels(strucData$goassessDxpmr7)[1]),]

# Now do cbfData
cbfData <- merge(cbfQAData, cbfData, by=c('bblid', 'scanid'))
cbfData <- cbfData[which(cbfData$pcaslExclude==0 & cbfData$bblid %in% bblidIndex),]
cbfData <- merge(cbfData, psData, by=c('bblid', 'scanid'))
cbfData <- merge(cbfData, mjData, by='bblid')
cbfData <- cbfData[-which(cbfData$goassessDxpmr7==levels(cbfData$goassessDxpmr7)[1]),]

# Now do dti data
dtiData <- merge(dtiQAData, trData, by=c('bblid', 'scanid'))
dtiData <- dtiData[which(dtiData$dti64Exclude!=1),]
dtiData <- merge(dtiData, mjData, by='bblid')
dtiData <- merge(dtiData, psData, by=c('bblid', 'scanid'))

# Now do FA data
faData <- merge(dtiQAData, faData)
faData <- faData[which(faData$dti64Exclude!=1),]
faData <- merge(faData, mjData)
faData <- merge(faData, psData)

# Now do the rest data
restData <- merge(alffData, rehoData, by=c('bblid', 'scanid'))
restData <- merge(restQAData, restData, by=c('bblid', 'scanid'))
restData <- restData[which(restData$restExclude==0 & restData$bblid %in% bblidIndex),]
restData <- merge(restData, mjData, by='bblid')
restData <- merge(restData, psData, by=c('bblid', 'scanid'))

# Now rm all variables we won't need
rm(mjData, volData, gmdData, ctData, t1QAData, cbfQAData, trData, dtiQAData, alffData, rehoData, restQAData)
