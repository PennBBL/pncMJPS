## This script will prepare all of the cognitition and imaging data
## It will return one massive RDS file
## Begin by loading any library(s)
install_load('psych','readstata13')

## Now load all of the data
all_cog <- read.dta13("./cannabis_psychosis_cnb_foradon.dta")
all_cog <- all_cog[,-grep('mar', names(all_cog))]
all_mar <- read.csv('./marcatSource.csv')
all_cog <- merge(all_cog, all_mar)#,all=T)
all_img <- read.csv("./n1601_imagingDataDump_2018-09-20.csv")
## Apply the age restriction to to the img data
all_img <- all_img[which(all_img$ageAtScan1/12 > 13.9),]
all_img <- all_img[which(all_img$bblid %in% all_cog$bblid),]
all_fake <- read.csv("./fakesub_exclude.csv")

## Now begin by mergin all of the data
all_data <- merge(all_cog, all_img, all=T,by='bblid',suffixes = c("",".y"))

## Now apply the fake endorsement restriction
all_data <- all_data[-which(all_data$bblid %in% all_fake$bblid),]

## Now check against the old img data
old_img <-  readRDS("./mjAnovaData.RDS")
mjData <- read.csv("./n9498_go1_foradon_061518.csv")
mjData$dosage <- NA
mjData$dosage[which(mjData$marcat=='MJ Non-User')] <- 0
mjData$dosage[which(mjData$marcat=='MJ Occ User' & mjData$mjpastyr=='')] <- 1
mjData$dosage[which(mjData$mjpastyr=="Less than once a month")] <- 2
mjData$dosage[which(mjData$mjpastyr=="About once a month")] <- 3
mjData$dosage[which(mjData$mjpastyr=="2-3 times a month")] <- 4
mjData$dosage[which(mjData$mjpastyr=="1-2 times a week")] <- 5
mjData$dosage[which(mjData$mjpastyr=="3-4 times a week")] <- 6
mjData$dosage[which(mjData$mjpastyr=="Everyday or nearly every day")] <- 7
old_mj <- mjData
rm(mjData)
old_mj <- old_mj[which(old_mj$ageatcnb1/12 > 13.9),]
flaggedSubjs <-  which(old_img$bblid %in% all_data$bblid == 'FALSE')
flaggedSubjs2 <- which(!is.na(all_data$scanid) & all_data$bblid %in% old_img$bblid == 'FALSE')

## Remove any dosage == 1 subjects
all_data <- all_data[-which(all_data$bblid %in% old_mj$bblid[old_mj$dosage==1]),]

## Now check for the same flagged subjs
all_data_img <- all_data[!is.na(all_data$scanid),]

## Now find which subjects from the struc paper are not in the
## new cog dataset
discrep_subjs <- all_data_img[which(!all_data_img$bblid %in% old_img$bblid),]
discrep_subjs <- discrep_subjs[-which(discrep_subjs$bblid %in% old_mj$bblid[which(old_mj$dosage==1)]),]

## These subjects are all dosage == 1 variables
## So lets remove them
#all_data <- all_data[-which(all_data$bblid %in% discrep_subjs$bblid),]

## Now ensure our factors are factors!
all_data$marcat <- factor(all_data$marcat)
all_data$race2 <- factor(all_data$race2)
all_data$race2.y <- factor(all_data$race2.y)
all_data$psBinary <- 'NO'
all_data$psBinary[which(all_data$goassessdxpmr6=='PS')] <- 'YES'
all_data$mjbinary <- 'NO'
all_data$mjbinary[which(all_data$marcat!=1)] <- 'YES'

# DO the same for the img data
all_data_img$marcat <- factor(all_data_img$marcat)
all_data_img$race2 <- factor(all_data_img$race2)
all_data_img$race2.y <- factor(all_data_img$race2.y)
all_data_img$psBinary <- 'NO'
all_data_img$psBinary[which(all_data_img$goassessdxpmr6=='PS')] <- 'YES'
all_data_img$mjbinary <- 'NO'
all_data_img$mjbinary[which(all_data_img$marcat!=1)] <- 'YES'

## Now write the full dataset
saveRDS(all_data, file="mjPSCogImg.RDS")
saveRDS(all_data_img, file="mjPSCogImgIsol.RDS")
