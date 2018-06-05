# AFGR June 2018 

## This script will be used to test for equivalence in our marcat groups
## There are three group comparisons we have to make:
##	1. User - Non
##	2. Freq - Non
##	3. Freq - User
## Differences will be calculated by taking the bootstrapped differences in means
## This will give us a confidence interval as well..

## Load library(s)
source('/home/arosen/adroseHelperScripts/R/afgrHelpFunc.R')
install_load('psych','ggplot2','caret','equivalence')

## Now load the data
mjData <- read.csv('/data/jux/BBL/projects/pncMJPS/data/n9462_mj_ps_cnb_fortmm.csv')
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

fakeData <- read.csv('/data/jux/BBL/projects/pncMJPS/data/fakesub_exclude.csv')
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

all.data <- read.csv('/data/jux/BBL/projects/pncMJPS/scripts/07_MJEffects/scripts/n1601_imagingDataDump_2018-04-04.csv')
all.data <- merge(all.data, psData)

# Now add the clinical bifactor scores
fac.data <- read.csv('/data/joy/BBL/studies/pnc/n9498_dataFreeze/clinical/n9498_goassess_itemwise_bifactor_scores_age_regressed_20170131.csv')
all.data <- merge(all.data, fac.data)

# Now also add our fs ct values
fs.values <- read.csv('/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_freesurferCt_20180213.csv')
all.data <- merge(all.data, fs.values)
all.data$fsAvgCT <- (all.data$LThickness + all.data$RThickness)/2

# Now clean the marcat variable
all.data <- all.data[-which(all.data$marcat==''),]
all.data <- all.data[which((all.data$ageAtScan1/12)>=14),]
 
### Begin bootstrapping down here
## First create our resamples
tmp.folds <- createResample(y=all.data$marcat, 10000)

## Now regress out age and sex so we can compare our groups
orig <- all.data
base.model <- paste("ageAtScan1 + sex + factor(race2)")
vars.of.interest <- c(107:245, 255:352, 353:470,471,1540,1550:1552)
for(v in vars.of.interest){
  name.val <- names(all.data)[v]
  tmp.formula <- as.formula(paste(name.val, "~", base.model))
  tmp.col <- rep(NA, 1504)
  tmp.mod <- lm(tmp.formula, all.data)
  index <- names(residuals(tmp.mod))
  all.data[index,name.val] <- NA
  all.data[index,name.val] <- scale(residuals(tmp.mod))
}

## Now go through each value and find our bootstrapped mean differences
output.differences <- NULL
for(q in 1:length(tmp.folds)){
  tmpData <- all.data[tmp.folds[[q]],]
  output.vals.tmp <- NULL
  for(v in vars.of.interest){
    name.val <- names(all.data)[v]
    mean.values <- summarySE(data=tmpData, measurevar=name.val, groupvars='marcat', na.rm=T)
    # Now prepare the output data
    output.row <- c(name.val, t(mean.values[,name.val]))
    output.vals.tmp <- rbind(output.vals.tmp, output.row)
  }
  output.differences <- rbind(output.differences, output.vals.tmp)
  print(q)
}
colnames(output.differences) <- c('ROI', 'MJ Frequent', 'Non-User', 'MJ User')
rownames(output.differences) <- NULL
output.differences <- as.data.frame(output.differences)

## Now find our mean differences
output.differences$user.minus.non <- as.numeric(as.character(output.differences[,4])) - as.numeric(as.character(output.differences[,3]))
output.differences$freq.minus.non <- as.numeric(as.character(output.differences[,2])) - as.numeric(as.character(output.differences[,3]))
output.differences$freq.minus.user <- as.numeric(as.character(output.differences[,2])) - as.numeric(as.character(output.differences[,4]))

## Now write this data
write.csv(output.differences, "~/mjBSVals.csv", quote=F, row.names=F)

## NOw get our means and confidence intervals four our differences
mean.vals.1 <- summarySE(data=output.differences, groupvars='ROI', measurevar='user.minus.non')
mean.vals.2 <- summarySE(data=output.differences, groupvars='ROI', measurevar='freq.minus.non')
mean.vals.3 <- summarySE(data=output.differences, groupvars='ROI', measurevar='freq.minus.user')

## Now plot our values
pdf('user.minus.non.pdf')
for(r in mean.vals.1$ROI){
  ## First grab all of our values
  tmp.dat <- output.differences[which(output.differences$ROI==r),]
  mean.value <- mean.vals.1[which(mean.vals.1$ROI==r),'user.minus.non']
  ci.value <- mean.vals.1[which(mean.vals.1$ROI==r),'ci']
  ## Now plot our histogram for these values
  out.plot <- ggplot(tmp.dat) +
    geom_histogram(aes(user.minus.non), color='red', alpha=.3) +
    geom_vline(xintercept=mean.value, yintercept=0) +
    geom_vline(xintercept=mean.value+ci.value, yintercept=0) +
    geom_vline(xintercept=mean.value-ci.value, yintercept=0) +
    geom_vline(xintercept=-.3, linetype="dotted") +
    geom_vline(xintercept=.3, linetype="dotted") +
    ggtitle(r) +
    coord_cartesian(xlim=c(-.7, .7))
  print(out.plot)
}
dev.off()
pdf('freq.minus.non.pdf')
for(r in mean.vals.1$ROI){
  ## First grab all of our values
  tmp.dat <- output.differences[which(output.differences$ROI==r),]
  mean.value <- mean.vals.2[which(mean.vals.1$ROI==r),'freq.minus.non']
  ci.value <- mean.vals.2[which(mean.vals.1$ROI==r),'ci']
  ## Now plot our histogram for these values
  out.plot <- ggplot(tmp.dat) +
    geom_histogram(aes(freq.minus.non), color='red', alpha=.3) +
    geom_vline(xintercept=mean.value, yintercept=0) +
    geom_vline(xintercept=mean.value+ci.value, yintercept=0) +
    geom_vline(xintercept=mean.value-ci.value, yintercept=0) +
    geom_vline(xintercept=-.3, linetype="dotted") +
    geom_vline(xintercept=.3, linetype="dotted") +
    ggtitle(r) +
    coord_cartesian(xlim=c(-.7, .7))
  print(out.plot)
}
dev.off()
pdf('freq.minus.user.pdf')
for(r in mean.vals.1$ROI){
  ## First grab all of our values
  tmp.dat <- output.differences[which(output.differences$ROI==r),]
  mean.value <- mean.vals.3[which(mean.vals.1$ROI==r),'freq.minus.user']
  ci.value <- mean.vals.3[which(mean.vals.1$ROI==r),'ci']
  ## Now plot our histogram for these values
  out.plot <- ggplot(tmp.dat) +
    geom_histogram(aes(freq.minus.user), color='red', alpha=.3) +
    geom_vline(xintercept=mean.value, yintercept=0) +
    geom_vline(xintercept=mean.value+ci.value, yintercept=0) +
    geom_vline(xintercept=mean.value-ci.value, yintercept=0) +
    geom_vline(xintercept=-.3, linetype="dotted") +
    geom_vline(xintercept=.3, linetype="dotted") +
    ggtitle(r) +
    coord_cartesian(xlim=c(-.7, .7))
  print(out.plot)
}
dev.off()
