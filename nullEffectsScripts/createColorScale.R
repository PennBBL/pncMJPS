# AFGR September 2018

## This script will be used to prepare f maps
## for the marcat variable in the PNC sample

## Load library(s)
#source('/home/arosen/adroseHelperScripts/R/afgrHelpFunc.R')
install_load('psych', 'pwr', 'ggplot2', 'caret', 'mgcv', 'reshape', 'lmerTest')
source("../../../jlfVisualizer/scripts/Rfunction/makeITKSnapColorTable.R")

## Load data
all.data <- readRDS('mjAnovaData.RDS')
all.data$marcat <- factor(all.data$marcat, levels=c("MJ Non-User", "MJ Occ User", "MJ Freq User"))

## Now do individual ROI's
# volume
summaryMetrics <- names(all.data)[c(107:245)]
outputVals <- NULL
for(s in summaryMetrics){
    # build a lm in all of the data
    tmpMod <- gam(all.data[,s] ~ s(ageAtScan1) + sex + averageManualRating + marcat + factor(race2) + overall_psychopathology_ar_4factor, data=all.data)
    aovMod <- anova.gam(tmpMod)
    # Now get our group t stats
    all.data$tmpVals <- residuals(gam(all.data[,s] ~ s(ageAtScan1) + sex + averageManualRating + factor(race2) + overall_psychopathology_ar_4factor, data=all.data, na.action=na.exclude))
    test1 <- t.test(tmpVals~marcat, data=all.data[-which(all.data$marcat=='MJ Freq User'),])
    test2 <- t.test(tmpVals~marcat, data=all.data[-which(all.data$marcat=='MJ Occ User'),])
    test3 <- t.test(tmpVals~marcat, data=all.data[-which(all.data$marcat=='MJ Non-User'),])


    tmpVals <- aovMod$pTerms.table['marcat',c('F', 'p-value')]
    tmpVals[2] <- 1-tmpVals[2]
    tmpRow <- c(s,tmpVals,test1$statistic,1-test1$p.value,test2$statistic,1-test2$p.value,test3$statistic,1-test3$p.value)
    outputVals <- rbind(outputVals, tmpRow)
    
}
# Now write the output files
writeColorTableandKey(inputColumn=2, inputData=outputVals, maxTmp=c(0,8), minTmp=c(-.8, 0), outName="volumef")
writeColorTableandKey(inputColumn=3, inputData=outputVals, maxTmp=c(0,1), minTmp=c(1,1.5), outName="volumep")
writeColorTableandKey(inputColumn=4, inputData=outputVals, maxTmp=c(0,4), minTmp=c(-4, 0), outName="volumetno")
writeColorTableandKey(inputColumn=5, inputData=outputVals, maxTmp=c(0,1), minTmp=c(1,1.5), outName="volumepno")
writeColorTableandKey(inputColumn=6, inputData=outputVals, maxTmp=c(0,4), minTmp=c(-4, 0), outName="volumetfo")
writeColorTableandKey(inputColumn=7, inputData=outputVals, maxTmp=c(0,1), minTmp=c(1,1.5), outName="volumepfo")
writeColorTableandKey(inputColumn=8, inputData=outputVals, maxTmp=c(0,4), minTmp=c(-4, 0), outName="volumetof")
writeColorTableandKey(inputColumn=9, inputData=outputVals, maxTmp=c(0,1), minTmp=c(1,1.5), outName="volumepof")



# CT ,255:352
summaryMetrics <- names(all.data)[c(255:352)]
outputVals <- NULL
for(s in summaryMetrics){
    # build a lm in all of the data
    tmpMod <- gam(all.data[,s] ~ s(ageAtScan1) + sex + averageManualRating + marcat + factor(race2) + overall_psychopathology_ar_4factor, data=all.data)
    aovMod <- anova.gam(tmpMod)
    # Now get our group t stats
    all.data$tmpVals <- residuals(gam(all.data[,s] ~ s(ageAtScan1) + sex + averageManualRating + factor(race2) + overall_psychopathology_ar_4factor, data=all.data, na.action=na.exclude))
    test1 <- t.test(tmpVals~marcat, data=all.data[-which(all.data$marcat=='MJ Freq User'),])
    test2 <- t.test(tmpVals~marcat, data=all.data[-which(all.data$marcat=='MJ Occ User'),])
    test3 <- t.test(tmpVals~marcat, data=all.data[-which(all.data$marcat=='MJ Non-User'),])
    
    
    tmpVals <- aovMod$pTerms.table['marcat',c('F', 'p-value')]
    tmpVals[2] <- 1-tmpVals[2]
    tmpRow <- c(s,tmpVals,test1$statistic,1-test1$p.value,test2$statistic,1-test2$p.value,test3$statistic,1-test3$p.value)
    outputVals <- rbind(outputVals, tmpRow)
    
}
# Now write the output files
writeColorTableandKey(inputColumn=2, inputData=outputVals, maxTmp=c(0,8), minTmp=c(-.8, 0), outName="ctf")
writeColorTableandKey(inputColumn=3, inputData=outputVals, maxTmp=c(0,1), minTmp=c(1,1.5), outName="ctp")
writeColorTableandKey(inputColumn=4, inputData=outputVals, maxTmp=c(0,4), minTmp=c(-4, 0), outName="cttno")
writeColorTableandKey(inputColumn=5, inputData=outputVals, maxTmp=c(0,1), minTmp=c(1,1.5), outName="ctpno")
writeColorTableandKey(inputColumn=6, inputData=outputVals, maxTmp=c(0,4), minTmp=c(-4, 0), outName="cttfo")
writeColorTableandKey(inputColumn=7, inputData=outputVals, maxTmp=c(0,1), minTmp=c(1,1.5), outName="ctpfo")
writeColorTableandKey(inputColumn=8, inputData=outputVals, maxTmp=c(0,4), minTmp=c(-4, 0), outName="cttof")
writeColorTableandKey(inputColumn=9, inputData=outputVals, maxTmp=c(0,1), minTmp=c(1,1.5), outName="ctpof")


# GMD ,353:470
summaryMetrics <- names(all.data)[c(353:470)]
outputVals <- NULL
for(s in summaryMetrics){
    # build a lm in all of the data
    tmpMod <- gam(all.data[,s] ~ s(ageAtScan1) + sex + averageManualRating + marcat + factor(race2) + overall_psychopathology_ar_4factor, data=all.data)
    aovMod <- anova.gam(tmpMod)
    # Now get our group t stats
    all.data$tmpVals <- residuals(gam(all.data[,s] ~ s(ageAtScan1) + sex + averageManualRating + factor(race2) + overall_psychopathology_ar_4factor, data=all.data, na.action=na.exclude))
    test1 <- t.test(tmpVals~marcat, data=all.data[-which(all.data$marcat=='MJ Freq User'),])
    test2 <- t.test(tmpVals~marcat, data=all.data[-which(all.data$marcat=='MJ Occ User'),])
    test3 <- t.test(tmpVals~marcat, data=all.data[-which(all.data$marcat=='MJ Non-User'),])
    
    
    tmpVals <- aovMod$pTerms.table['marcat',c('F', 'p-value')]
    tmpVals[2] <- 1-tmpVals[2]
    tmpRow <- c(s,tmpVals,test1$statistic,1-test1$p.value,test2$statistic,1-test2$p.value,test3$statistic,1-test3$p.value)
    outputVals <- rbind(outputVals, tmpRow)
    
}
# Now write the output files
writeColorTableandKey(inputColumn=2, inputData=outputVals, maxTmp=c(0,8), minTmp=c(-.8, 0), outName="gmdf")
writeColorTableandKey(inputColumn=3, inputData=outputVals, maxTmp=c(0,1), minTmp=c(1,1.5), outName="gmdp")
writeColorTableandKey(inputColumn=4, inputData=outputVals, maxTmp=c(0,4), minTmp=c(-4, 0), outName="gmdtno")
writeColorTableandKey(inputColumn=5, inputData=outputVals, maxTmp=c(0,1), minTmp=c(1,1.5), outName="gmdpno")
writeColorTableandKey(inputColumn=6, inputData=outputVals, maxTmp=c(0,4), minTmp=c(-4, 0), outName="gmdtfo")
writeColorTableandKey(inputColumn=7, inputData=outputVals, maxTmp=c(0,1), minTmp=c(1,1.5), outName="gmdpfo")
writeColorTableandKey(inputColumn=8, inputData=outputVals, maxTmp=c(0,4), minTmp=c(-4, 0), outName="gmdtof")
writeColorTableandKey(inputColumn=9, inputData=outputVals, maxTmp=c(0,1), minTmp=c(1,1.5), outName="gmdpof")

# Here is the scp call
# scp *-* arosen@chead:/data/jux/BBL/projects/pncMJPS/data/colorAndKey/
