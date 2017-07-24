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


output <- as.data.frame(cbind(bblid, pathVals, age, sex, mjDosage, psTD, race, marBin, quality))
output$age <- as.numeric(as.character(output$age))
output$marBin <- ordered(output$marBin)
output$mjDosage <- ordered(output$mjDosage)
output$psTD <- ordered(output$psTD)
output$quality <- as.factor(output$quality)
output$inclusion <- 1

saveRDS(output, "ctMJPSTDInteraction.rds")
