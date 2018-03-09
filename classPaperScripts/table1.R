# Library(s)
install_load('MatchIt')

# Now start loading the data down here
# This should be run from the dta directory
mjData <- read.csv("../data/n9462_mj_ps_cnb_fortmm.csv")
mjData$dosage <- NA
mjData$dosage[which(mjData$marcat=='MJ Non-User')] <- 0
mjData$dosage[which(mjData$marcat=='MJ User' & mjData$mjpastyr=='')] <- 1
mjData$dosage[which(mjData$mjpastyr=="Less than once a month")] <- 2
mjData$dosage[which(mjData$mjpastyr=="About once a month")] <- 3
mjData$dosage[which(mjData$mjpastyr=="2-3 times a month")] <- 4
mjData$dosage[which(mjData$mjpastyr=="1-2 times a week")] <- 5
mjData$dosage[which(mjData$mjpastyr=="3-4 times a week")] <- 1
mjData$dosage[which(mjData$mjpastyr=="Everyday or nearly every day")] <- 1

# Now load imaging data
img.data <- read.csv('../data/imagingDataAll.csv')

# Now prepare all of the data
all.data <- merge(img.data, mjData)
all.data <- all.data[-which(all.data$dosage==1),]
all.data$usageBin <- 0
all.data$usageBin[all.data$dosage>1] <- 1
all.data <- all.data[-which(is.na(all.data$dosage)),]

# Now isolate genders
male.data <- all.data[which(all.data$sex==1),]
female.data <- all.data[which(all.data$sex==2),]

# Now create our age matched samples
# Starting with male
tmpDat <- male.data[c('bblid', 'scanid', 'usageBin', 'ageAtScan1', 'envSES')]
mod <- matchit(usageBin ~ ageAtScan1 + envSES, data=tmpDat, ratio=3, na.action=na.omit)
male.data.all <- male.data
male.data <- male.data[as.vector(mod$match.matrix),]
male.data <- rbind(male.data, male.data.all[which(male.data.all$usageBin==1),])
male.data <- male.data[complete.cases(male.data[,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]),]

# Now do female
tmpDat <- female.data[c('bblid', 'scanid', 'usageBin', 'ageAtScan1', 'envSES')]
mod <- matchit(usageBin ~ ageAtScan1 + envSES, data=tmpDat, ratio=3, na.action=na.omit)
female.data.all <- female.data
female.data <- female.data[as.vector(mod$match.matrix),]
female.data <- rbind(female.data, female.data.all[which(female.data.all$usageBin==1),])
female.data <- female.data[complete.cases(female.data[,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]),]

# Before I wirte this though I have to get the deomgraphics from the freq users
demo.vals.male <- summarySE(data=male.data, measurevar="ageAtScan1", groupvars="marcat")
demo.vals.female <- summarySE(data=female.data, measurevar="ageAtScan1", groupvars="marcat")

## Now do the frequent users
mjData <- read.csv("../data/n9462_mj_ps_cnb_fortmm.csv")
mjData$dosage <- NA
mjData$dosage[which(mjData$marcat=='MJ Non-User')] <- 1
mjData$dosage[which(mjData$marcat=='MJ User' & mjData$mjpastyr=='')] <- 1
mjData$dosage[which(mjData$mjpastyr=="Less than once a month")] <- 2
mjData$dosage[which(mjData$mjpastyr=="About once a month")] <- 3
mjData$dosage[which(mjData$mjpastyr=="2-3 times a month")] <- 4
mjData$dosage[which(mjData$mjpastyr=="1-2 times a week")] <- 5
mjData$dosage[which(mjData$mjpastyr=="3-4 times a week")] <- 6
mjData$dosage[which(mjData$mjpastyr=="Everyday or nearly every day")] <- 7

# Now give us all the values
all.data <- merge(img.data, mjData)
all.data <- all.data[-which(all.data$dosage==1),]
all.data <- all.data[-which(is.na(all.data$dosage)),]
all.data$usageBin <- 0
all.data$usageBin[all.data$dosage>=6] <- 1
all.data$usageBinOrig <- all.data$usageBin

# Now prepare a sex specific values
male.data <- all.data[which(all.data$sex==1),]
male.data <- male.data[complete.cases(male.data[,c(grep('dti_dtitk_jhulabel_fa', names(male.data)))]),]
female.data <- all.data[which(all.data$sex==2),]
female.data <- female.data[complete.cases(female.data[,c(grep('dti_dtitk_jhulabel_fa', names(female.data)))]),]
sex <- rep('Male', 3)
demo.vals.male <- rbind(demo.vals.male, summarySE(data=male.data, measurevar="ageAtScan1", groupvars="marcat")[c(1),])
demo.vals.male <- cbind(sex, demo.vals.male)
sex <- rep('Female', 3)
demo.vals.female <- rbind(demo.vals.female, summarySE(data=female.data, measurevar="ageAtScan1", groupvars="marcat")[c(1),])
demo.vals.female <- cbind(sex, demo.vals.female)
demo.vals <- rbind(demo.vals.male, demo.vals.female)

# Now write the table
write.csv(demo.vals, "demographicVals.csv", quote=F, row.names=F)
