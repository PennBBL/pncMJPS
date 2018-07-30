# AFGR April 2018
# This script is gfoing to be used to assess variance across our three marcat labels.
# It is going to go from a very coarse, whole brain, lobular, and then specific JLF labels.
# I am going to run an anova for our marcat factor for each of these metrics.

## Load library(s)
#source('/home/arosen/adroseHelperScripts/R/afgrHelpFunc.R')
install_load('psych', 'pwr', 'ggplot2', 'caret', 'mgcv')

## Load data
all.data <- readRDS('mjAnovaData.RDS')
all.data$marcat <- factor(all.data$marcat, levels=c("MJ Non-User", "MJ Occ User", "MJ Freq User"))
alc.data <- read.csv('alcData.csv')
all.data <- merge(all.data, alc.data)
all.data <- all.data[complete.cases(all.data$alchoholZScore),]

## Now create an alcohol usage factor
cut.vals <- quantile(all.data$alchoholZScore, na.rm=T, probs=c(.33, .67))
all.data$alcFactor <- NA
all.data$alcFactor[which(all.data$alchoholZScore<=cut.vals[1])] <- "LowAlcGroup"
all.data$alcFactor[which(all.data$alchoholZScore>cut.vals[1] & all.data$alchoholZScore<=cut.vals[2])] <- "MidAlcGroup"
all.data$alcFactor[which(all.data$alchoholZScore>cut.vals[2])] <- "HiAlcGroup"
all.data$alcFactor <- factor(all.data$alcFactor, levels=c("LowAlcGroup", "MidAlcGroup", "HiAlcGroup"))

# Now create a scaled alcohol variable
all.data$alcScale <- log(range01(all.data$alchoholZScore)+1)

## Now go through the summary metrics and check for any 3x1 differences
summaryMetrics <- names(all.data)[c(1540,471,1550:1552)]
outputVals <- NULL
for(s in summaryMetrics){
    # build a lm in all of the data
    tmpMod <- gam(all.data[,s] ~ s(ageAtScan1) + sex + averageManualRating + marcat + race2 + overall_psychopathology_ar_4factor + alcScale, data=all.data)
    aovMod <- anova.gam(tmpMod)
    tmpVals <- aovMod$pTerms.table['marcat',c('F', 'p-value')]
  tmpRow <- c(s, tmpVals)
  outputVals <- rbind(outputVals, tmpRow)
}

# Now do the lobular values
summaryMetrics <- names(all.data)[c(1553:1600)]
for(s in summaryMetrics){
  # build a lm in all of the data
  tmpMod <- gam(all.data[,s] ~ s(ageAtScan1) + sex + averageManualRating + marcat + factor(race2) + overall_psychopathology_ar_4factor + alcScale, data=all.data)
  aovMod <- anova.gam(tmpMod)
  tmpVals <- aovMod$pTerms.table['marcat',c('F', 'p-value')]
  tmpRow <- c(s, tmpVals)
  outputVals <- rbind(outputVals, tmpRow)
}

# Now do individual ROI's
summaryMetrics <- names(all.data)[c(107:245,255:352,353:470,472:569)]
for(s in summaryMetrics){
  # build a lm in all of the data
  tmpMod <- gam(all.data[,s] ~ s(ageAtScan1) + sex + averageManualRating + marcat + factor(race2) + overall_psychopathology_ar_4factor + alcScale, data=all.data)
  aovMod <- anova.gam(tmpMod)
  tmpVals <- aovMod$pTerms.table['marcat',c('F', 'p-value')]
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
#pValFDR[507:574] <- p.adjust(as.numeric(outputVals[507:574,3]), method='fdr')
outputVals <- cbind(outputVals, pValFDR)

## Now run FDR w/in specfic subsets for each of our modalities
##Create the index to use for our lobe values
lobeIndex <- c(1,2,3,4,5,6,11,12)
volIndex <- c(3,4,5,6,8,9,15,16,21:26,30,31,34,35,44,45,46,47,48,49,60,61,62,63,64,65,66,67,70,71,76,77,84,85,86,87,88,89,100,101,110,111,124,125)
gmdIndex <- c(1:4,6,7,10,11,12,13,21,22,25,26,35,36,37,38,45,46,51,52,53,54,55,56,61,62,67,68,75,76,77,78,79,80,91,92,101,102,115,116)
ctIndex <- gmdIndex[-c(1:10)]
ctIndex <- ctIndex-20
pValFDRLS <- rep(NA, dim(outputVals)[1])
pValFDRLS[c(6:17)[lobeIndex]] <- p.adjust(as.numeric(outputVals[6:17,3])[lobeIndex], method='fdr')
pValFDRLS[c(18:29)[lobeIndex]] <- p.adjust(as.numeric(outputVals[18:29,3])[lobeIndex], method='fdr')
pValFDRLS[c(30:41)[lobeIndex]] <- p.adjust(as.numeric(outputVals[30:41,3])[lobeIndex], method='fdr')
pValFDRLS[c(42:53)[lobeIndex]] <- p.adjust(as.numeric(outputVals[42:53,3])[lobeIndex], method='fdr')
pValFDRLS[c(54:192)[volIndex]] <- p.adjust(as.numeric(outputVals[54:192,3][volIndex]), method='fdr')
pValFDRLS[c(193:290)[ctIndex]] <- p.adjust(as.numeric(outputVals[193:290,3][ctIndex]), method='fdr')
pValFDRLS[c(291:408)[gmdIndex]] <- p.adjust(as.numeric(outputVals[291:408,3][gmdIndex]), method='fdr')
pValFDRLS[c(409:506)[ctIndex]] <- p.adjust(as.numeric(outputVals[409:506,3][ctIndex]), method='fdr')
#pValFDRLS[c(507:574)] <- p.adjust(as.numeric(outputVals[507:574,3]), method='fdr')
outputVals <- cbind(outputVals, pValFDRLS)

# Now write the outputVals
#colnames(outputVals) <- c('ROI', 'F Stat', 'p Value', 'Q Value', 'Q Value Manual Selection')
rownames(outputVals) <- NULL
write.csv(outputVals, "threeByOneAnovaMarcatAlc.csv", quote=F, row.names=F)

## Now I need to go about plotting the nominally significant differences
roiNames <- as.character(outputVals[,1])#[which(outputVals[,3]<.05),1])
all.data$marcat <- factor(all.data$marcat, levels=c('MJ Non-User', 'MJ Occ User', 'MJ Freq User'))
# Now go through each roi and plot the differences after controlling for all of our covariates
# We are also going to do the 3 group difference tests so I am going to have to prepare a matrix 
# for all of these t values and associated p values
groupDiffMatrix <- matrix(NA, nrow=length(roiNames), ncol=15)
colnames(groupDiffMatrix) <- c('roi', 'kruskall-wallis-stat', 'kruskall-wallis-p.val','N-Mt', 'N-Mp', 'N-Ft', 'N-Fp', 'M-Fp', 'M-Ft', 'N-Mpfdr', 'N-Fpfdr', 'M-Fpfdr', 'N-Md', 'N-Fd', 'M-Fd')
index <- 1
pdf('anovaGroupDiffExplore.pdf')
for(n in roiNames){
  # First create the regressed values
  tmpData <- all.data[complete.cases(all.data[,c(n, "alchoholZScore")]),]
  tmpData$race2 <- factor(tmpData$race2)
  #Now create the regressed values
  tmpModel <- as.formula(paste(n, "~s(ageAtScan1)+sex+race2+averageManualRating+overall_psychopathology_ar_4factor+alcScale"))
  tmpMod <- gam(tmpModel, data=tmpData)
  tmpData$tmpVals <- as.numeric(scale(residuals(tmpMod)))
  foo <- summarySE(data=tmpData, groupvars='marcat', measurevar='tmpVals')
  foo$marcat <- factor(foo$marcat, levels=c('MJ Non-User', 'MJ Occ User', 'MJ Freq User'))
  # Now fill out our groupDiffMatrix row
  test1 <- t.test(tmpVals~marcat, data=tmpData[-which(tmpData$marcat=='MJ Freq User'),])
  test2 <- t.test(tmpVals~marcat, data=tmpData[-which(tmpData$marcat=='MJ Occ User'),])
  test3 <- t.test(tmpVals~marcat, data=tmpData[-which(tmpData$marcat=='MJ Non-User'),])
  fdrValues <- p.adjust(c(test1$p.value, test2$p.value, test3$p.value), method='fdr')
  # Now calculate cohen d's
  occ.vs.non.d <- effsize::cohen.d(tmpVals~marcat, data=tmpData[-which(tmpData$marcat=='MJ Freq User'),])$estimate
  non.vs.freq.d <- effsize::cohen.d(tmpVals~marcat, data=tmpData[-which(tmpData$marcat=='MJ Occ User'),])$estimate
  freq.vs.occ.d <- effsize::cohen.d(tmpVals~marcat, data=tmpData[-which(tmpData$marcat=='MJ Non-User'),])$estimate
  # Now do a nonparametric diff in ranks test
  krusk.val <- kruskal.test(tmpVals ~marcat, data = tmpData)


  outRow <- c(n,krusk.val$statistic, krusk.val$p.value, test1$statistic,test1$p.value,test2$statistic,test2$p.value,test3$statistic,test3$p.value, fdrValues, occ.vs.non.d, non.vs.freq.d, freq.vs.occ.d)
  
  groupDiffMatrix[index,] <- outRow
  # Now create a plot for these values
  tmpPlot <- ggplot(foo, aes(x=marcat, y=tmpVals)) + 
                           geom_bar(stat="identity", position=position_dodge(), size=.1) + 
                           geom_errorbar(aes(ymin=tmpVals-se, ymax=tmpVals+se), 
                           width = .2, position=position_dodge(.9)) +
			   #geom_jitter(data=tmpData, aes(x=marcat, y=tmpVals),alpha=.1, position=position_jitter(0.2)) +
                           ylab(n)
			   
               #print(tmpPlot)
  tmpPlot <- ggplot(foo, aes(x=marcat, y=tmpVals)) + 
                           geom_bar(stat="identity", position=position_dodge(), size=.1) + 
                           geom_errorbar(aes(ymin=tmpVals-se, ymax=tmpVals+se), 
                           width = .2, position=position_dodge(.9)) +
			   geom_jitter(data=tmpData, aes(x=marcat, y=tmpVals),alpha=.1, position=position_jitter(0.2)) +
                           ylab(n)
			   
               #print(tmpPlot)
  index <- index+1
}
dev.off()


## Now look fdr correct all of the t values for the frequent vs non group comparisons
pValFDR <- rep(NA, dim(outputVals)[1])
pValFDR[6:17] <- p.adjust(as.numeric(groupDiffMatrix[6:17,3]), method='fdr')
pValFDR[18:29] <- p.adjust(as.numeric(groupDiffMatrix[18:29,3]), method='fdr')
pValFDR[30:41] <- p.adjust(as.numeric(groupDiffMatrix[30:41,3]), method='fdr')
pValFDR[42:53] <- p.adjust(as.numeric(groupDiffMatrix[42:53,3]), method='fdr')
pValFDR[54:192] <- p.adjust(as.numeric(groupDiffMatrix[54:192,3]), method='fdr')
pValFDR[193:290] <- p.adjust(as.numeric(groupDiffMatrix[193:290,3]), method='fdr')
pValFDR[291:408] <- p.adjust(as.numeric(groupDiffMatrix[291:408,3]), method='fdr')
pValFDR[409:506] <- p.adjust(as.numeric(groupDiffMatrix[409:506,3]), method='fdr')
pValFDRKW <- pValFDR
groupDiffMatrix <- cbind(groupDiffMatrix, pValFDRKW)

## Now apply fdr correction to the kruskall wallis p values
# Now add in our multiple comparisions and what not
pValFDR <- rep(NA, dim(outputVals)[1])
pValFDR[6:17] <- p.adjust(as.numeric(groupDiffMatrix[6:17,7]), method='fdr')
pValFDR[18:29] <- p.adjust(as.numeric(groupDiffMatrix[18:29,7]), method='fdr')
pValFDR[30:41] <- p.adjust(as.numeric(groupDiffMatrix[30:41,7]), method='fdr')
pValFDR[42:53] <- p.adjust(as.numeric(groupDiffMatrix[42:53,7]), method='fdr')
pValFDR[54:192] <- p.adjust(as.numeric(groupDiffMatrix[54:192,7]), method='fdr')
pValFDR[193:290] <- p.adjust(as.numeric(groupDiffMatrix[193:290,7]), method='fdr')
pValFDR[291:408] <- p.adjust(as.numeric(groupDiffMatrix[291:408,7]), method='fdr')
pValFDR[409:506] <- p.adjust(as.numeric(groupDiffMatrix[409:506,7]), method='fdr')
pValFDRNonVsFrequentT <- pValFDR
groupDiffMatrix <- cbind(groupDiffMatrix, pValFDRNonVsFrequentT)
groupDiffFull <- cbind(outputVals, groupDiffMatrix)

## Now write out all of these values
write.csv(groupDiffFull, "groupDiffMatAlc.csv", quote=F, row.names=F)

# Now write out the regions with nominally significant F stats
write.csv(groupDiffFull[which(groupDiffFull[,3]<.05),], "nominallySigValsAlc.csv", quote=F, row.names=F)


## Now test the interactions
summaryMetrics <- names(all.data)[c(1540,471,1550:1552)]
outputVals <- NULL
for(s in summaryMetrics){
    # build a lm in all of the data
    tmpMod <- lm(all.data[,s] ~ poly(ageAtScan1,3)*marcat + sex + averageManualRating + factor(race2) + overall_psychopathology_ar_4factor + alcScale, data=all.data)
    aovMod <- anova(tmpMod)
    tmpVals <- aovMod['poly(ageAtScan1, 3):marcat',c("F value","Pr(>F)")]
    tmpRow <- c(s, tmpVals)
    outputVals <- rbind(outputVals, tmpRow)
}

# Now do the lobular values
summaryMetrics <- names(all.data)[c(1553:1600)]
for(s in summaryMetrics){
    # build a lm in all of the data
    tmpMod <- lm(all.data[,s] ~ poly(ageAtScan1,3)*marcat + sex + averageManualRating + factor(race2) + overall_psychopathology_ar_4factor + alcScale, data=all.data)
    aovMod <- anova(tmpMod)
    tmpVals <- aovMod['poly(ageAtScan1, 3):marcat',c("F value","Pr(>F)")]
    tmpRow <- c(s, tmpVals)
    outputVals <- rbind(outputVals, tmpRow)
}

# Now do individual ROI's
summaryMetrics <- names(all.data)[c(107:245,255:352,353:470,472:569)]
for(s in summaryMetrics){
    # build a lm in all of the data
    tmpMod <- lm(all.data[,s] ~ poly(ageAtScan1,3)*marcat + sex + averageManualRating + factor(race2) + overall_psychopathology_ar_4factor + alcScale, data=all.data)
    aovMod <- anova(tmpMod)
    tmpVals <- aovMod['poly(ageAtScan1, 3):marcat',c("F value","Pr(>F)")]
    tmpRow <- c(s, tmpVals)
    outputVals <- rbind(outputVals, tmpRow)
}
write.csv(outputVals, "interactionFandPStatsAlc.csv", quote=F, row.names=F)
