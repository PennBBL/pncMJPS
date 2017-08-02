# Load all data and library(s)
source('/data/joy/BBL/projects/pncMJPS/scripts/01_dataPrep/scripts/dataPrepNoFake.R')
source('/data/joy/BBL/projects/pncMJPS/scripts/02_runGamm/functions/functions.R')

# Now run the gam's for the structural data 
# Collapse all fo the Users into MJ User bin for marcat 
strucData$marcat[strucData$marcat=='MJ Frequent User'] <- 'MJ User'
#strucData <- strucData[-which(strucData$marcat=='MJ User'),]
strucData <- strucData[-which(strucData$marcat==levels(strucData$marcat)[1]),]
strucData <- strucData[-which(strucData$goassessDxpmr7==levels(strucData$goassessDxpmr7)[1]),]

# Prepare the paths 
pathVals <- paste("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/voxelwiseMaps_antsCt/", strucData$scanid, "_CorticalThicknessNormalizedToTemplate2mm.nii.gz", sep='')
pathVals2 <- paste("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/voxelwiseMaps_gmd/", strucData$scanid, "_atropos3class_prob02SubjToTemp2mm.nii.gz", sep='')


bblid <- strucData$bblid
age <- strucData$ageAtScan1
sex <- strucData$sex
qaVal <- strucData$averageManualRating
mjDosage <- strucData$dosage
psTD <- strucData$goassessDxpmr7
race <- strucData$race2
race[which(race == 3)] <- 2
marBin <- strucData$marcat
quality <- strucData$averageManualRating


output <- as.data.frame(cbind(bblid, pathVals, pathVals2, age, sex, mjDosage, psTD, race, marBin, quality))
output$age <- as.numeric(as.character(output$age))
output$pathVals <- as.character(output$pathVals)
output$pathVals2 <- as.character(output$pathVals2)
output$marBin <- ordered(output$marBin)
output$mjDosage <- ordered(output$mjDosage)
output$psTD <- ordered(output$psTD)
output$quality <- as.factor(output$quality)
output$inclusion <- 1
output$inclusionFrequent <- 0
output$inclusionFrequent[which(output$mjDosage == 0 | output$mjDosage == 5 | output$mjDosage == 6)] <- 1 


saveRDS(output, "ctMJPSTDInteraction.rds")


# Now do the cbf RDS file 
cbfData$marcat[cbfData$marcat=='MJ Frequent User'] <- 'MJ User'
cbfData <- cbfData[-which(cbfData$marcat==levels(cbfData$marcat)[1]),]
cbfData <- cbfData[-which(cbfData$goassessDxpmr7==levels(cbfData$goassessDxpmr7)[1]),]
pathVals3 <- paste("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/asl/voxelwiseMaps_cbf/", cbfData$scanid, "_asl_quant_ssT1Std.nii.gz", sep='')

bblid <- cbfData$bblid
age <- cbfData$ageAtScan1
sex <- cbfData$sex
qaVal <- cbfData$averageManualRating
mjDosage <- cbfData$dosage
psTD <- cbfData$goassessDxpmr7
race <- cbfData$race2
race[which(race == 3)] <- 2
marBin <- cbfData$marcat
quality <- cbfData$pcaslRelMeanRMSMotion

output <- as.data.frame(cbind(bblid, pathVals3, age, sex, mjDosage, psTD, race, marBin, quality))
output$age <- as.numeric(as.character(output$age))
output$pathVals3 <- as.character(output$pathVals3)
output$marBin <- ordered(output$marBin)
output$mjDosage <- ordered(output$mjDosage)
output$psTD <- ordered(output$psTD)
output$quality <- as.numeric(output$quality)
output$inclusion <- 1
output$inclusionFrequent <- 0
output$inclusionFrequent[which(output$mjDosage == 0 | output$mjDosage == 5 | output$mjDosage == 6)] <- 1 

saveRDS(output, "cbfMJPSTDInteraction.rds")

