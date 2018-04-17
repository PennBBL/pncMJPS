# AFGR Aprin 2018
# This script is gfoing to be used to assess variance across our three marcat labels.
# It is going to go from a very coarse, whole brain, lobular, and then specific JLF labels.
# I am going to run an anova for our marcat factor for each of these metrics.

## Load library(s)
source('/home/arosen/adroseHelperScripts/R/afgrHelpFunc.R')
install_load('psych', 'pwr', 'ggplot2')

## Load data
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
mjData$mjLabel <- NA
mjData$mjLabel[which(mjData$marcat=="MJ Non-User")] <- "NonUser"
mjData$mjLabel[which(mjData$marcat=="MJ User")] <- "User"
mjData$mjLabel[which(mjData$marcat=="MJ Frequent User")] <- "FreqUser"
mjData$mjBinLabel <- "NonUser"
mjData$mjBinLabel[which(mjData$mjLabel=="User")] <- "User"
mjData$mjBinLabel[which(mjData$mjLabel=="FreqUser")] <- "User"
mjData$mjLabel <- as.factor(mjData$mjLabel)

fakeData <- read.csv('/data/joy/BBL/projects/pncMJPS/data/fakesub_exclude.csv')
# Now create a new label which collapses ps op and td into a new label. This will be a factor w/ two levels
psData <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/clinical/n1601_diagnosis_dxpmr_20170509.csv')
psData <- psData[,-grep('goassessDxpmr4', names(psData))]
psData$pathLabel <- NA
psData$pathLabel[which(psData$goassessDxpmr7=='TD' | psData$goassessDxpmr7=='OP')] <- 'TDOP'
psData$pathLabel[which(psData$goassessDxpmr7=='PS')] <- 'PS'
psData$pathLabel <- as.factor(psData$pathLabel)

demoData <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_demographics_go1_20161212.csv')
psData <- merge(psData, demoData, by=c('bblid', 'scanid'))
psData <- merge(psData, mjData)

all.data <- read.csv('/data/joy/BBL/projects/pncMJPS/scripts/07_MJEffects/scripts/n1601_imagingDataDump_2018-04-04.csv')
all.data <- merge(all.data, psData)
# Now also add our fs ct values
fs.values <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_freesurferCt_20180213.csv')
all.data <- merge(all.data, fs.values)
all.data$fsAvgCT <- (all.data$LThickness + all.data$RThickness)/2

# Now clean the marcat variable
all.data <- all.data[-which(all.data$marcat==''),]

# Now first run through some power analyses
# First do the power test for the summary metrics
pwr.anova.test(k=3, n=41, f=.01, sig.level=.05)

## Now go through the summary metrics and check for any 3x1 differences
summaryMetrics <- names(all.data)[c(1540,471,1550:1552)]
outputVals <- NULL
for(s in summaryMetrics){
  # build a lm in all of the data
  tmpMod <- lm(all.data[,s] ~ ageAtScan1 + sex + averageManualRating + marcat + goassessDxpmr7, data=all.data)
  aovMod <- aov(tmpMod)
  tmpVals <- matrix(unlist(summary(aovMod)), ncol=5, nrow=6)[4,4:5]
  tmpRow <- c(s, tmpVals)
  outputVals <- rbind(outputVals, tmpRow)
}

# Now do the lobular values
summaryMetrics <- names(all.data)[c(1553:1600)]
for(s in summaryMetrics){
  # build a lm in all of the data
  tmpMod <- lm(all.data[,s] ~ ageAtScan1 + sex + averageManualRating + marcat + goassessDxpmr7, data=all.data)
  aovMod <- aov(tmpMod)
  tmpVals <- matrix(unlist(summary(aovMod)), ncol=5, nrow=6)[4,4:5]
  tmpRow <- c(s, tmpVals)
  outputVals <- rbind(outputVals, tmpRow)
}

# Now do individual ROI's
summaryMetrics <- names(all.data)[c(107:245,255:352,353:470,472:569)]
for(s in summaryMetrics){
  # build a lm in all of the data
  tmpMod <- lm(all.data[,s] ~ ageAtScan1 + sex + averageManualRating + marcat + goassessDxpmr7, data=all.data)
  aovMod <- aov(tmpMod)
  tmpVals <- matrix(unlist(summary(aovMod)), ncol=5, nrow=6)[4,4:5]
  tmpRow <- c(s, tmpVals)
  outputVals <- rbind(outputVals, tmpRow)
}

# Now do the FS values down here
summaryMetrics <- names(all.data)[c(1668:1735,1740,1741)]
for(s in summaryMetrics){
  # build a lm in all of the data
  tmpMod <- lm(all.data[,s] ~ ageAtScan1 + sex + averageManualRating + marcat + goassessDxpmr7, data=all.data)
  aovMod <- aov(tmpMod)
  tmpVals <- matrix(unlist(summary(aovMod)), ncol=5, nrow=6)[4,4:5]
  tmpRow <- c(s, tmpVals)
  outputVals <- rbind(outputVals, tmpRow)
}

# Now add in our multiple comparisions and what not
pValFDR <- rep(NA, dim(outputVals)[1])
pValFDR[6:17] <- p.adjust(as.numeric(outputVals[6:17,3]), method='fdr')
pValFDR[18:29] <- p.adjust(as.numeric(outputVals[18:29,3]), method='fdr')
pValFDR[30:41] <- p.adjust(as.numeric(outputVals[30:41,3]), method='fdr')
pValFDR[42:53] <- p.adjust(as.numeric(outputVals[42:53,3]), method='fdr')
pValFDR[54:192] <- p.adjust(as.numeric(outputVals[54:192,3]), method='fdr')
pValFDR[193:290] <- p.adjust(as.numeric(outputVals[193:290,3]), method='fdr')
pValFDR[291:408] <- p.adjust(as.numeric(outputVals[291:408,3]), method='fdr')
pValFDR[409:506] <- p.adjust(as.numeric(outputVals[409:506,3]), method='fdr')
pValFDR[507:574] <- p.adjust(as.numeric(outputVals[507:574,3]), method='fdr')
outputVals <- cbind(outputVals, pValFDR)

# Now write the outputVals
colnames(outputVals) <- c('ROI', 'F Stat', 'p Value', 'Q Value')
rownames(outputVals) <- NULL
write.csv(outputVals, "threeByOneAnovaMarcat.csv", quote=F, row.names=F)

## Now I need to go about plotting the nominally significant differences
roiNames <- outputVals[which(as.numeric(outputVals[,3])<.05),1]
all.data$marcat <- factor(all.data$marcat, levels=c('MJ Non-User', 'MJ User', 'MJ Frequent User'))
# Now go through each roi and plot the differences after controlling for all of our covariates
# We are also going to do the 3 group difference tests so I am going to have to prepare a matrix 
# for all of these t values and associated p values
groupDiffMatrix <- matrix(NA, nrow=length(roiNames), ncol=10)
colnames(groupDiffMatrix) <- c('roi', 'N-Mt', 'N-Mp', 'N-Ft', 'N-Fp', 'M-Fp', 'M-Ft', 'N-Mpfdr', 'N-Fpfdr', 'M-Fpfdr')
index <- 1
pdf('anovaGroupDiffExplore.pdf')
for(n in roiNames){
  # First create the regressed values
  tmpData <- all.data[complete.cases(all.data[,n]),]
  #Now create the regressed values
  tmpModel <- as.formula(paste(n, "~ageAtScan1+sex+race2+averageManualRating+goassessDxpmr7"))
  tmpMod <- lm(tmpModel, data=tmpData)
  tmpData$tmpVals <- residuals(tmpMod)
  foo <- summarySE(data=tmpData, groupvars='marcat', measurevar='tmpVals')
  foo$marcat <- factor(foo$marcat, levels=c('MJ Non-User', 'MJ User', 'MJ Frequent User'))
  # Now fill out our groupDiffMatrix row
  test1 <- t.test(tmpVals~marcat, data=tmpData[-which(tmpData$marcat=='MJ Frequent User'),])
  test2 <- t.test(tmpVals~marcat, data=tmpData[-which(tmpData$marcat=='MJ User'),])
  test3 <- t.test(tmpVals~marcat, data=tmpData[-which(tmpData$marcat=='MJ Non-User'),])
  fdrValues <- p.adjust(c(test1$p.value, test2$p.value, test3$p.value), method='fdr')
  outRow <- c(n,test1$statistic,test1$p.value,test2$statistic,test2$p.value,test3$statistic,test3$p.value, fdrValues)
  groupDiffMatrix[index,] <- outRow
  # Now create a plot for these values
  tmpPlot <- ggplot(foo, aes(x=marcat, y=tmpVals)) + 
                           geom_bar(stat="identity", position=position_dodge(), size=.1) + 
                           geom_errorbar(aes(ymin=tmpVals-se, ymax=tmpVals+se), 
                           width = .2, position=position_dodge(.9)) + 
                           ylab(n)
			   
  print(tmpPlot)
  index <- index+1
}
dev.off()
write.csv(groupDiffMatrix, "~/groupDiffMat.csv", quote=F, row.names=F)
