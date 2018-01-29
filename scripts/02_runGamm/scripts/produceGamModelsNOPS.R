## Load library(s)
install_load('psych', 'ggplot2', 'pROC', 'ggrepel', 'caret', 'randomForest', 'MatchIt', 'glmnet', 'doMC', 'voxel', 'mgcv')

# Now start loading the data down here
# This should be run from the dta directory
mjData <- read.csv("n9462_mj_ps_cnb_fortmm.csv")
mjData$dosage <- NA
mjData$dosage[which(mjData$marcat=='MJ Non-User')] <- 0
mjData$dosage[which(mjData$marcat=='MJ User' & mjData$mjpastyr=='')] <- 1
mjData$dosage[which(mjData$mjpastyr=="Less than once a month")] <- 2
mjData$dosage[which(mjData$mjpastyr=="About once a month")] <- 3
mjData$dosage[which(mjData$mjpastyr=="2-3 times a month")] <- 4
mjData$dosage[which(mjData$mjpastyr=="1-2 times a week")] <- 5
mjData$dosage[which(mjData$mjpastyr=="3-4 times a week")] <- 6
mjData$dosage[which(mjData$mjpastyr=="Everyday or nearly every day")] <- 7

# Now load imaging data
img.data <- read.csv('./imagingDataAll.csv')
#img.data <- read.csv('./n1601_imagingDataDump_20180104.csv')

# Now give us all the values
all.data <- merge(img.data, mjData)
all.data <- all.data[-which(all.data$dosage==1),]
all.data <- all.data[-which(is.na(all.data$dosage)),]

# Now prepare a sex specific values
male.data <- all.data[which(all.data$sex==1),]
#male.data <- male.data[which(male.data$dosage == 0 | male.data$dosage > 5),]
female.data <- all.data[which(all.data$sex==2),]

# Now add a binary matrix for use or no use
male.data$usageBin <- 0
male.data$usageBin[male.data$dosage>1] <- 1

# Now create a matched cohort here
tmpDat <- male.data[c('bblid', 'scanid', 'usageBin', 'ageAtScan1', 'envSES')]
mod <- matchit(usageBin ~ ageAtScan1 + envSES, data=tmpDat, ratio=3, na.action=na.omit)
male.data.all <- male.data
male.data <- male.data[as.vector(mod$match.matrix),]
male.data <- rbind(male.data, male.data.all[which(male.data.all$usageBin==1),])
male.data.all.m <- male.data

# Now loop through every region and produce a gam for it
colVals <- grep('_jlf_', names(male.data))
colVals <- append(colVals, grep('_dtitk_', names(male.data)))
pdf('nonLinearTrends.pdf')
for(z in colVals){
    tVal <-t.test(male.data[,z]~usageBin, data=male.data)
    outputRow <- c(colnames(male.data)[z], as.numeric(tVal$statistic), as.numeric(tVal$p.value))
    # First make a gam model
    formVal <- as.formula(paste(colnames(male.data)[z], "~s(ageAtScan1)+sex+s(dosage,k=3)"))
    mod <- gam(formVal, data=male.data)
    sumVal <- summary(mod)
    if(sumVal$s.table[2,4] < .1){
      print(plotGAM(mod, smooth.cov='dosage'))
    }
}
dev.off()

# Now do the same but create bar graphs
pdf('nonLinearTrendsBarGraph.pdf')
male.data.ss <- male.data[which(male.data$marcat!=""),]
male.data.ss$marcat <- factor(male.data.ss$marcat, levels(male.data.ss$marcat)[c(3, 4, 2)])
male.data.ss[,colVals] <- as.matrix(scale(male.data.ss[,colVals]))
for(z in colVals){
    tVal <-t.test(male.data[,z]~usageBin, data=male.data)
    outputRow <- c(colnames(male.data)[z], as.numeric(tVal$statistic), as.numeric(tVal$p.value))
    # First make a gam model
    formVal <- as.formula(paste(colnames(male.data)[z], "~s(ageAtScan1)+sex+s(dosage,k=3)"))
    mod <- gam(formVal, data=male.data.ss)
    sumVal <- summary(mod)
    if(sumVal$s.table[2,4] < .1){
        tmpVals <- summarySE(data=male.data.ss, groupvars='marcat', measurevar=colnames(male.data)[z], na.rm=T)
        val <-colnames(male.data)[z]
        outPlot <- ggplot(tmpVals, aes(y=tmpVals[,3], x=marcat)) +
        geom_bar(stat='identity', position=position_dodge(), width=.5) +
          labs(x='marcat', y=val) +
          geom_errorbar(aes(ymin=tmpVals[,3]-se, ymax=tmpVals[,3]+se),
            width = .1, position=position_dodge(.9)) +
          theme_bw()
        print(outPlot)
    }
}
dev.off()
