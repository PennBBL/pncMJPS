# AFGR June 2017
# This script will be used to produce bar graphs for the GMD and CT ROI's
# that were returned w/ signifianct interactions

# Load any library(s)
source('/home/arosen/adroseHelperScripts/R/afgrHelpFunc.R')
source('/data/joy/BBL/projects/pncMJPS/scripts/02_runGamm/functions/functions.R')
#install_load('ggplot2', 'grid', 'gridExtra', 'scales', 'mgcv')
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)
library(mgcv)

# Declare any functions
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

# Now load the data
source('/data/joy/BBL/projects/pncMJPS/scripts/01_dataPrep/scripts/dataPrepNoFake.R')
# Now limit to only frequent users
strucData <- strucData[-which(strucData$dosage==2 |strucData$dosage==3 | strucData$dosage==4),]
cbfData <- cbfData[-which(cbfData$dosage==2 |cbfData$dosage==3 | cbfData$dosage==4),]
restData <- restData[-which(restData$dosage==2 |restData$dosage==3 | restData$dosage==4),]
stathParc <- stathParc[-which(stathParc$dosage==2 |stathParc$dosage==3 | stathParc$dosage==4),]
stathCT <- stathCT[-which(stathCT$dosage==2 | stathCT$dosage==3 | stathCT$dosage==4),]

# Now produce the bar graphs for all nominally sig GMD ROI's

# Now loop through each one and produce a bar graph for it
strucData$marcat[strucData$marcat=='MJ Frequent User'] <- 'MJ User'
strucData <- strucData[-which(strucData$marcat==levels(strucData$marcat)[1]),]
sigGMD <- runGamModel(strucData, 'mprage_jlf_gmd', 'averageManualRating')
sigGMD <- sigGMD[which(sigGMD[,2]<.05),]
#pdf('gmdInteractionsFreq.pdf')
#for(i in varsOfInterest){
#  mainTitle <- names(strucData)[i]
#  formulaValue <- as.formula(paste(mainTitle, '~ageAtScan1+ageAtScan1^2+sex'))
#  strucData[,i] <- scale(residuals(lm(formulaValue, data=strucData)))
#  foo <- summarySE(strucData, measurevar=mainTitle, groupvars=c('marcat','goassessDxpmr7') , na.rm=T)
#  barPlotToPrint <- ggplot(foo, aes(x=factor(marcat), y=foo[,4], fill=goassessDxpmr7)) + 
#                           geom_bar(stat="identity", position=position_dodge(), size=.1) + 
#                           geom_errorbar(aes(ymin=foo[,4]-se, ymax=foo[,4]+se), 
#                           width = .2, position=position_dodge(.9)) + 
#                           ggtitle(mainTitle) +
#                           ylab('Mean GMD Value')
#  print(barPlotToPrint) 
#}
#dev.off()

# Now onto volume
sigVOL <- runGamModel(strucData, 'mprage_jlf_vol', 'averageManualRating')
sigVOL <- sigVOL[which(sigVOL[,2] <.05),]
pdf('volInteractions.pdf')
for(i in sigVOL[,1]){
  mainTitle <- i
  formulaValue <- as.formula(paste(mainTitle, '~s(ageAtScan1)+sex'))
  strucData[,i] <- scale(residuals(gam(formulaValue, data=strucData)))
  foo <- summarySE(strucData, measurevar=mainTitle, groupvars=c('marcat','goassessDxpmr7') , na.rm=T)
  barPlotToPrint <- ggplot(foo, aes(x=factor(marcat), y=foo[,4], fill=goassessDxpmr7)) + 
                           geom_bar(stat="identity", position=position_dodge(), size=.1) + 
                           geom_errorbar(aes(ymin=foo[,4]-se, ymax=foo[,4]+se), 
                           width = .2, position=position_dodge(.9)) + 
                           ggtitle(mainTitle) +
                           ylab('Mean CT Value')#+ 
                           #scale_y_continuous(limits=c(4,5.2),oob=rescale_none)
  print(barPlotToPrint) 
}
dev.off()

# Now do this for where we see signifianct CT interactions
# This means we will have to find out which CT regions have signifianct interactions
sigCT <- runGamModel(strucData, 'mprage_jlf_ct', 'averageManualRating')
sigCT <- sigCT[which(sigCT[,2]<.05),]
pdf('ctInteractions.pdf')
for(i in sigCT[,1]){
  mainTitle <- i
  formulaValue <- as.formula(paste(mainTitle, '~s(ageAtScan1)+sex'))
  strucData[,i] <- scale(residuals(gam(formulaValue, data=strucData)))
  foo <- summarySE(strucData, measurevar=mainTitle, groupvars=c('marcat','goassessDxpmr7') , na.rm=T)
  barPlotToPrint <- ggplot(foo, aes(x=factor(marcat), y=foo[,4], fill=goassessDxpmr7)) + 
                           geom_bar(stat="identity", position=position_dodge(), size=.1) + 
                           geom_errorbar(aes(ymin=foo[,4]-se, ymax=foo[,4]+se), 
                           width = .2, position=position_dodge(.9)) + 
                           ggtitle(mainTitle) +
                           ylab('Mean CT Value')#+ 
                           #scale_y_continuous(limits=c(4,5.2),oob=rescale_none)
  print(barPlotToPrint) 
}
dev.off()

# Now do CBF
cbfData$marcat[cbfData$marcat=='MJ Frequent User'] <- 'MJ User'
cbfData <- cbfData[-which(cbfData$marcat==levels(cbfData$marcat)[1]),]
sigCBF <- runGamModel(cbfData, 'pcasl_jlf_cbf', 'pcaslTSNR')
sigCBF <- sigCBF[which(sigCBF[,2]<.05),]
pdf('cbfInteraction.pdf')
for(i in sigCBF[,1]){
  mainTitle <- i
  formulaValue <- as.formula(paste(mainTitle, '~s(ageAtScan1)+sex'))
  index <- which(!is.na(cbfData[,i]))
  cbfData[index,i] <- scale(residuals(gam(formulaValue, data=cbfData)))
  foo <- summarySE(cbfData, measurevar=mainTitle, groupvars=c('marcat','goassessDxpmr7') , na.rm=T)
  barPlotToPrint <- ggplot(foo, aes(x=factor(marcat), y=foo[,4], fill=goassessDxpmr7)) + 
                           geom_bar(stat="identity", position=position_dodge(), size=.1) + 
                           geom_errorbar(aes(ymin=foo[,4]-se, ymax=foo[,4]+se), 
                           width = .2, position=position_dodge(.9)) + 
                           ggtitle(mainTitle) +
                           ylab('Mean CBF Value')#+ 
                           #scale_y_continuous(limits=c(4,5.2),oob=rescale_none)
  print(barPlotToPrint) 
}
dev.off()

# Now onto rest 
restData$marcat[restData$marcat=='MJ Frequent User'] <- 'MJ User'
restData <- restData[-which(restData$marcat==levels(restData$marcat)[1]),]
sigReho <- runGamModel(restData, 'rest_jlf_reho', 'restRelMeanRMSMotion')
sigReho <- sigReho[which(sigReho[,2]<.05),]
pdf('rehoInteractions.pdf')
for(i in sigReho[1]){
  mainTitle <- i
  formulaValue <- as.formula(paste(mainTitle, '~s(ageAtScan1)+sex'))
  index <- which(!is.na(restData[,i]))
  restData[index,i] <- scale(residuals(gam(formulaValue, data=restData)))
  foo <- summarySE(restData, measurevar=mainTitle, groupvars=c('marcat','goassessDxpmr7') , na.rm=T)
  barPlotToPrint <- ggplot(foo, aes(x=factor(marcat), y=foo[,4], fill=goassessDxpmr7)) + 
                           geom_bar(stat="identity", position=position_dodge(), size=.1) + 
                           geom_errorbar(aes(ymin=foo[,4]-se, ymax=foo[,4]+se), 
                           width = .2, position=position_dodge(.9)) + 
                           ggtitle(mainTitle) +
                           ylab('Mean CT Value')#+ 
                           #scale_y_continuous(limits=c(4,5.2),oob=rescale_none)
  print(barPlotToPrint) 
}
dev.off()

sigAlff <- runGamModel(restData, 'rest_jlf_alff', 'restRelMeanRMSMotion')
sigAlff <- sigAlff[which(sigAlff[,2]<.05),]
pdf('alffInteractions.pdf')
for(i in sigAlff[1]){
  mainTitle <- i
  formulaValue <- as.formula(paste(mainTitle, '~s(ageAtScan1)+sex'))
  index <- which(!is.na(restData[,i]))
  restData[index,i] <- scale(residuals(gam(formulaValue, data=restData)))
  foo <- summarySE(restData, measurevar=mainTitle, groupvars=c('marcat','goassessDxpmr7') , na.rm=T)
  barPlotToPrint <- ggplot(foo, aes(x=factor(marcat), y=foo[,4], fill=goassessDxpmr7)) + 
                           geom_bar(stat="identity", position=position_dodge(), size=.1) + 
                           geom_errorbar(aes(ymin=foo[,4]-se, ymax=foo[,4]+se), 
                           width = .2, position=position_dodge(.9)) + 
                           ggtitle(mainTitle) +
                           ylab('Mean ALFF Value')#+ 
                           #scale_y_continuous(limits=c(4,5.2),oob=rescale_none)
  print(barPlotToPrint) 
}
dev.off()

# Now produce our figures for stathis regions
# First find the regions with nominally sig interactions
# Now do the sttah parc
tmp <- stathParc
#stathParc <- stathParc[-which(stathParc$marcat=='MJ User'),]
stathParc$marcat[stathParc$marcat=='MJ Frequent User'] <- 'MJ User'
stathParc <- stathParc[-which(stathParc$marcat==levels(stathParc$marcat)[1]),]
stathParc <- stathParc[-which(stathParc$goassessDxpmr7==levels(stathParc$goassessDxpmr7)[1]),]
#stathParc <- merge(strucData, stathParc, by=c('bblid', 'scanid'))
colnames(stathParc) <- gsub(x=colnames(stathParc), pattern='.y', replacement='')
stathGmd <- runGamModel(stathParc, 'NZMean_', 'averageManualRating')
stathGmd <- stathGmd[which(as.numeric(stathGmd[,2]) <.05),]
pdf('nominallySigStahGmd.pdf')
for(i in stathGmd[,1]){
  colIndex <- grep(i, names(stathParc))[1]
  colName <- i
  formulaValue <- as.formula(paste(colName, '~s(ageAtScan1)+sex'))
  stathParc[,colIndex] <- scale(residuals(gam(formulaValue, data=stathParc)))
  foo <- summarySE(stathParc, measurevar=colName, groupvars=c('marcat','goassessDxpmr7') , na.rm=T)
  barPlotToPrint <- ggplot(foo, aes(x=factor(marcat), y=foo[,4], fill=goassessDxpmr7)) + 
                           geom_bar(stat="identity", position=position_dodge(), size=.1) + 
                           geom_errorbar(aes(ymin=foo[,4]-se, ymax=foo[,4]+se), 
                           width = .2, position=position_dodge(.9)) + 
                           ggtitle(colName) +
                           ylab('Mean GMD Value')#+ 
                           #scale_y_continuous(limits=c(4,5.2),oob=rescale_none)
  print(barPlotToPrint)   
}
dev.off()

# Now do stath CT 
stathCT$marcat[stathCT$marcat=='MJ Frequent User'] <- 'MJ User'
stathCT <- stathCT[-which(stathCT$marcat==levels(stathCT$marcat)[1]),]
stathCT <- stathCT[-which(stathCT$goassessDxpmr7==levels(stathCT$goassessDxpmr7)[1]),]
stathCT <- merge(strucData, stathCT, by=intersect(names(stathCT),names(strucData)))
colnames(stathCT) <- gsub(x=colnames(stathCT), pattern='.y', replacement='')
sigStathCT <- runGamModel(stathCT, 'NZMean_', 'averageManualRating')
sigStathCT <- sigStathCT[which(as.numeric(sigStathCT[,2]) < .05),]
pdf('nominallySigStahCT.pdf')
for(i in sigStathCT[,1]){
  colName <- i
  colIndex <- grep(colName, names(stathCT))[1]
  formulaValue <- as.formula(paste(colName, '~s(ageAtScan1)+sex'))
  stathCT[,colIndex] <- scale(residuals(gam(formulaValue, data=stathCT)))
  foo <- summarySE(stathCT, measurevar=colName, groupvars=c('marcat','goassessDxpmr7') , na.rm=T)
  barPlotToPrint <- ggplot(foo, aes(x=factor(marcat), y=foo[,4], fill=goassessDxpmr7)) + 
                           geom_bar(stat="identity", position=position_dodge(), size=.1) + 
                           geom_errorbar(aes(ymin=foo[,4]-se, ymax=foo[,4]+se), 
                           width = .2, position=position_dodge(.9)) + 
                           ggtitle(colName) +
                           ylab('Mean CT Value')#+ 
                           #scale_y_continuous(limits=c(4,5.2),oob=rescale_none)
  print(barPlotToPrint)   
}
dev.off()
