## Load library(s)
install_load('ggplot2', 'grid', 'gridExtra', 'MatchIt', 'mgcv', 'reshape2')

## Create a function
swr = function(string, nwrap=30) {
    paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)

## Load the data
mjData <- read.csv("../data/n9462_mj_ps_cnb_fortmm.csv")
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
img.data <- read.csv('../data/imagingDataAll.csv')
fa.data <- read.csv('../data/meanFAVals.csv')
img.data <- merge(img.data, fa.data)

# Now load the new names
newNames <- read.csv('../data/labelNames.csv')

# Now age regress the fa.data
img.data$age <- img.data$ageAtScan1
img.data$age2 <- img.data$age^2
img.data$age3 <- img.data$ageAtScan1^3
mod <- lm(dti_jlf_fa_MeanFA ~ age + age2 + age3, data=img.data)
img.data$dti_jlf_fa_MeanFAar <- residuals(mod)
img.data$mean.fa <- apply(img.data[,572:619], 1, function(x) mean(x, na.rm=T))

# Now prepare all of the data
all.data <- merge(img.data, mjData)
all.data <- all.data[-which(all.data$dosage==1),]
all.data.freq <- all.data[which(all.data$dosage>=6),]
all.data$usageBin <- 0
all.data$usageBin[all.data$dosage>1] <- 1
all.data <- all.data[-which(is.na(all.data$dosage)),]
all.data.freq$usageBin <- 1

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
male.data <- rbind(male.data, all.data.freq[which(all.data.freq$sex==1),])
male.data$dti_jlf_fa_MeanFAar <- as.matrix(scale(male.data$dti_jlf_fa_MeanFAar))[1:352]
male.data[,579:612] <- scale(male.data[,579:612])

# Now do female
tmpDat <- female.data[c('bblid', 'scanid', 'usageBin', 'ageAtScan1', 'envSES')]
mod <- matchit(usageBin ~ ageAtScan1 + envSES, data=tmpDat, ratio=3, na.action=na.omit)
female.data.all <- female.data
female.data <- female.data[as.vector(mod$match.matrix),]
female.data <- rbind(female.data, female.data.all[which(female.data.all$usageBin==1),])
female.data <- rbind(female.data, all.data.freq[which(all.data.freq$sex==2),])
female.data$dti_jlf_fa_MeanFAar <- as.matrix(scale(female.data$dti_jlf_fa_MeanFAar))[1:255]
female.data[,579:612] <- scale(female.data[,579:612])

## Now prepare our vaues for FA in the labels
outVals <- NULL
names(male.data)[572:619] <- as.character(newNames$Unclassified)
allVals <- names(male.data)[572:619]
vals <- c(4, 35, 36, 39, 21)
inVals <- allVals[vals]
for(q in inVals){
    tmpVal <- summarySE(male.data, measurevar=q, groupvars='marcat', na.rm=T)
    colnames(tmpVal)[3] <- 'mean'
    q <- gsub(pattern="dti_dtitk_jhulabel_fa_", replacement="", x=q)
    tmpVal$ROI <- q
    outVals <- rbind(outVals, tmpVal)
}
outVals$Gender <- 'Male'
male.fa.vals <- outVals
outVals <- NULL
names(female.data)[572:619] <- as.character(newNames$Unclassified)
for(q in inVals){
    tmpVal <- summarySE(female.data, measurevar=q, groupvars='marcat', na.rm=T)
    colnames(tmpVal)[3] <- 'mean'
    q <- gsub(pattern="dti_dtitk_jhulabel_fa_", replacement="", x=q)
    tmpVal$ROI <- q
    outVals <- rbind(outVals, tmpVal)
}
outVals$Gender <- 'Female'
female.fa.vals <- outVals
fa.vals <- rbind(male.fa.vals, female.fa.vals)
fa.vals$marcat <- factor(fa.vals$marcat, levels=c("MJ Non-User", "MJ User", "MJ Frequent User"))
fa.vals$Marcat <- fa.vals$marcat
fa.vals$ROI <- swr(fa.vals$ROI)

## Now plot these values
outPlot <- ggplot(fa.vals, aes(x=Marcat, y=mean, fill=Marcat)) +
  geom_bar(stat='identity', position=position_dodge(), size=.1) +
  labs(title="", x="", y="Mean FA (z-score)") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(.9)) +
  theme_bw() +
  facet_grid(Gender~ROI, scales="free") +
  theme(legend.position="right",
  axis.text.x = element_blank(),
  axis.text.y = element_text(size=16, face="bold"),
  axis.ticks.x=element_blank(),
  axis.title=element_text(size=16,face="bold", angle=180),
  strip.text.y = element_text(size = 16, angle = 270, face="bold"),
  strip.text.x = element_text(size = 16, angle = 0, face="bold"),
  panel.spacing = unit(2, "lines")) + scale_fill_grey()

# Now print the output
png('figure2-marcatMeanValues.png', height=12, width=20, units='in', res=300)
outPlot
dev.off()
q()

# Now make a gam w/ the dosage value
modFemale <- gam(mean.fa ~ s(dosage,k=2) + s(age), data=female.data)
plotGAM(gamFit=modFemale, smooth.cov='dosage')
modMale <- gam(mean.fa ~ s(dosage,k=2) + s(age), data=male.data)
plotGAM(gamFit=modMale, smooth.cov='dosage')
