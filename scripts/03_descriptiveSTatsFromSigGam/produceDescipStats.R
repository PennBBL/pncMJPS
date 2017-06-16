# AFGR June 2017
# This script will be used to produce bar graphs for the GMD and CT ROI's
# that were returned w/ signifianct interactions

# Load any library(s)
source('/home/arosen/adroseHelperScripts/R/afgrHelpFunc.R')
source('/data/joy/BBL/projects/pncMJPS/scripts/02_runGamm/functions/functions.R')
install_load('ggplot2', 'grid', 'gridExtra', 'scales', 'mgcv')

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

# Now produce the bar graphs for all GMD ROI's
varsOfInterest <- grep('mprage_jlf_gmd', names(strucData))
# Now loop through each one and produce a bar graph for it
strucData$marcat[strucData$marcat=='MJ Frequent User'] <- 'MJ User'
strucData <- strucData[-which(strucData$marcat==levels(strucData$marcat)[1]),]
strucData <- strucData[-which(strucData$goassessDxpmr7==levels(strucData$goassessDxpmr7)[1]),]
pdf('gmdInteractions.pdf')
for(i in varsOfInterest){
  mainTitle <- names(strucData)[i]
  formulaValue <- as.formula(paste(mainTitle, '~ageAtScan1+ageAtScan1^2+sex'))
  strucData[,i] <- scale(residuals(lm(formulaValue, data=strucData)))
  foo <- summarySE(strucData, measurevar=mainTitle, groupvars=c('marcat','goassessDxpmr7') , na.rm=T)
  barPlotToPrint <- ggplot(foo, aes(x=factor(marcat), y=foo[,4], fill=goassessDxpmr7)) + 
                           geom_bar(stat="identity", position=position_dodge(), size=.1) + 
                           geom_errorbar(aes(ymin=foo[,4]-se, ymax=foo[,4]+se), 
                           width = .2, position=position_dodge(.9)) + 
                           ggtitle(mainTitle) +
                           ylab('Mean GMD Value')
  print(barPlotToPrint) 
}
dev.off()

# Now do this for where we see signifianct CT interactions
# This means we will have to find out which CT regions have signifianct interactions
sigCT <- runGamModel(strucData, 'mprage_jlf_ct', 'averageManualRating')
sigCT <- sigCT[which(p.adjust(sigCT[,2], method='fdr')<.05),]
strucData <- strucData[-which(strucData$marcat==levels(strucData$marcat)[1]),]
strucData <- strucData[-which(strucData$goassessDxpmr7==levels(strucData$goassessDxpmr7)[1]),]
pdf('ctInteractions.pdf')
for(i in sigCT[,1]){
  mainTitle <- i
  formulaValue <- as.formula(paste(mainTitle, '~ageAtScan1+ageAtScan1^2+sex'))
  strucData[,i] <- scale(residuals(lm(formulaValue, data=strucData)))
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
cbfData <- cbfData[-which(cbfData$goassessDxpmr7==levels(cbfData$goassessDxpmr7)[1]),]
sigCBF <- runGamModel(cbfData, 'pcasl_jlf_cbf', 'pcaslTSNR')
sigCBFN <- length(which(p.adjust(sigCBF[,2], method='fdr')<.05))
for(i in sigCBF[,1]){
  mainTitle <- i
  formulaValue <- as.formula(paste(mainTitle, '~ageAtScan1+ageAtScan1^2+sex'))
  cbfData[,i] <- scale(residuals(lm(formulaValue, data=cbfData)))
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
