## Load library(s)
#source('/home/arosen/adroseHelperScripts/R/afgrHelpFunc.R')
install_load('psych','ggplot2','caret','equivalence', 'mgcv','TOSTER','foreach','doParallel','plyr')

## Now load the data
all.data <- readRDS('mjAnovaData.RDS')
all.data$marcat <- factor(all.data$marcat, levels=c("MJ Non-User", "MJ Occ User", "MJ Freq User"))

# Delcare our variables of interest
roiNames <- colnames(all.data)[c(1540,471,1550:1552,1553:1600,107:245,255:352,353:470,472:569)]

# Declare our output matrix
groupDiffMatrix <- matrix(NA, nrow=length(roiNames), ncol=15)
colnames(groupDiffMatrix) <- c('roi', 'kruskall-wallis-stat', 'kruskall-wallis-p.val','N-Mt', 'N-Mp', 'N-Ft', 'N-Fp', 'M-Fp', 'M-Ft', 'N-Mpfdr', 'N-Fpfdr', 'M-Fpfdr', 'nonVsOccD', 'nonVsFreqD', 'occVsFreqD')
index <- 1
# Now run through every ROI and calculate our effect sizes
for(n in roiNames){
    # First create the regressed values
    tmpData <- all.data[complete.cases(all.data[,n]),]
    #Now create the regressed values
    tmpModel <- as.formula(paste(n, "~s(ageAtScan1)+sex+factor(race2)+averageManualRating+overall_psychopathology_ar_4factor"))
    tmpMod <- gam(tmpModel, data=tmpData)
    tmpData$tmpVals <- as.numeric(scale(residuals(tmpMod)))
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
    # now export em
    outRow <- c(n,krusk.val$statistic, krusk.val$p.value, test1$statistic,test1$p.value,test2$statistic,test2$p.value,test3$statistic,test3$p.value, fdrValues, occ.vs.non.d, non.vs.freq.d, freq.vs.occ.d)
    
    groupDiffMatrix[index,] <- outRow
    index <- index+1
}

## Now isolate only the values we want to plot for cortical thickness
toPlot <- groupDiffMatrix[c(grep("mprage_jlf_ct_", groupDiffMatrix[,1]),grep("mprage_jlfLobe_ct_", groupDiffMatrix[,1])),c(1,13,14,15)]
toPlot <- as.data.frame(toPlot)
toPlot[,2:4] <- apply(toPlot[,2:4], 2, function(x) as.numeric(as.character(x)))
# Now I need to attach lobe distinctions for all of the cortical ROI's
lobarValues <- c(100,101,102,103,104,105,106,107,108,109,104,105,108,109,100,101,104,105,104,105,122,123,104,105,108,109,122,123,108,109,104,105,100,101,104,105,104,105,108,109,104,105,106,107,104,105,104,105,122,123,108,109,108,109,104,105,104,105,100,101,106,107,100,101,102,103,106,107,106,107,104,105,122,123,104,105,122,123,104,105,104,105,104,105,106,107,108,109,106,107,122,123,122,123,104,105,122,123)
toPlot$lobeValue <- NA
for(i in names(table(lobarValues))){
    toPlot[seq(2,99)[which(lobarValues==i)],'lobeValue'] <- i
}
toPlot$lobeFactor <- factor(toPlot$lobeValue, levels=c(104,105,102,103,100,101,106,107,122,123,108,109))
toPlot$lobeFactor <- revalue(toPlot$lobeFactor, c("104"="Frontal", "105"="Frontal", "102"="Insular", "103"="Insular", "100"="Limbic","101"="Limbic", "106"="Parietal", "107"="Parietal", "122"="Temporal", "123"="Temporal", "108"="Occipital", "109"="Occipital"))
toPlot$hemisphere <- "Right"
toPlot$hemisphere[grep("_L_", toPlot$roi)] <- "Left"
## Now make our plots!
toPlot[,2:4] <- apply(toPlot[,2:4], 2, function(x) as.numeric(as.character(x)))
out.plot.o.n <- ggplot(toPlot[complete.cases(toPlot$lobeValue),], aes(x=lobeFactor, y=nonVsOccD, grouplobeFactor)) +
  geom_violin() +
  geom_point(position= position_jitter(), aes(shape=hemisphere)) +
  coord_cartesian(ylim=c(-.6, .6)) +
  facet_grid(.~lobeFactor, scales="free", space="free_x") +
  ylab("Effect Size (d)") +
xlab("Lobe") +
theme_bw() +
theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(), legend.position="none") +
ggtitle("Non-User - Occasional")

out.plot.f.n <- ggplot(toPlot[complete.cases(toPlot$lobeValue),], aes(x=lobeFactor, y=occVsFreqD, grouplobeFactor)) +
geom_violin() +
geom_point(position= position_jitter(), aes(shape=hemisphere)) +
coord_cartesian(ylim=c(-.6, .6)) +
facet_grid(.~lobeFactor, scales="free", space="free_x") +
ylab("") +
xlab("Lobe") +
theme_bw() +
theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(), axis.text.y = element_blank(),axis.ticks.y=element_blank()) +
ggtitle("Non-User - Frequent")
pdf("cohenDPlot.pdf", width=18, height=12)
multiplot(out.plot.o.n, out.plot.f.n,cols=2)
dev.off()
