## Load library(s)
install_load('readstata13','ggplot2')

## Load the data
data <- read.dta13('./cannabis_psychosis_cnb_foradon.dta')
## Now load the imaging data
img.data <- read.csv("~/Documents/dataRel/n2416_imagingDataDump_2018-04-22.csv")
img.data <- img.data[which(img.data$tpvalue==1),]
## Make a ps vs not ps group
data$psBinary <- 'NO'
data$psBinary[which(data$goassessdxpmr6=='PS')] <- 'YES'
## Now find all of our cognition means
cog.names <- names(data)[c(65,66,67,71,72,73,74)]

## Now go through a loop and prepare all of the values
## the first loop will be looking at an interaction between ps and mj
out.data.one <- NULL
for(i in cog.names){
    tmp.dat <- summarySE(groupvars=c('maruse','psBinary'), data=data, measurevar=i, na.rm=T)
    colnames(tmp.dat)[4] <- c('mean')
    tmp.dat$cog <- i
    out.data.one <- rbind(out.data.one, tmp.dat)
}

## Now clean the data
out.data.one <- out.data.one[which(out.data.one$maruse=='frequser' | out.data.one$maruse=='nonuser'),]
#out.data.one <- out.data.one[-which(out.data.one$goassessdxpmr6==''),]

## Now make the first graph comparing frequent users PS vs !PS
out.data.one$cog <- factor(out.data.one$cog, levels=cog.names)
out.data.one$maruse <- factor(out.data.one$maruse, levels=c('nonuser','frequser'))
out.plot.one <- ggplot(out.data.one, aes(x=cog, y=mean, group=maruse, color=maruse)) +
  geom_point(position=position_dodge(.3)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(.3)) +
  facet_grid(.~psBinary) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## Now do the same analysis in the isolated img cohort
data.img.reduce <- data[which(data$bblid %in% img.data$bblid),]

out.data.one <- NULL
for(i in cog.names){
    tmp.dat <- summarySE(groupvars=c('maruse','psBinary'), data=data.img.reduce, measurevar=i, na.rm=T)
    colnames(tmp.dat)[4] <- c('mean')
    tmp.dat$cog <- i
    out.data.one <- rbind(out.data.one, tmp.dat)
}


## Now clean the data
out.data.one <- out.data.one[which(out.data.one$maruse=='frequser' | out.data.one$maruse=='nonuser'),]
#out.data.one <- out.data.one[-which(out.data.one$goassessdxpmr6==''),]

## Now make the first graph comparing frequent users PS vs !PS
out.data.one$cog <- factor(out.data.one$cog, levels=cog.names)
out.data.one$maruse <- factor(out.data.one$maruse, levels=c('nonuser','frequser'))
out.plot.two <- ggplot(out.data.one, aes(x=cog, y=mean, group=maruse, color=maruse)) +
geom_point(position=position_dodge(.3)) +
geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(.3)) +
facet_grid(.~psBinary) +
theme(axis.text.x = element_text(angle = 90, hjust = 1))

## Now do the same analysis w/ age bins
## the first loop will be looking at an interaction between ps and mj
out.data.one <- NULL
for(i in cog.names){
    tmp.dat <- summarySE(groupvars=c('maruse','psBinary','age_cat'), data=data, measurevar=i, na.rm=T)
    colnames(tmp.dat)[5] <- c('mean')
    tmp.dat$cog <- i
    out.data.one <- rbind(out.data.one, tmp.dat)
}

out.data.one <- out.data.one[which(out.data.one$maruse=='frequser' | out.data.one$maruse=='nonuser'),]

## Now plot em
## Now make the first graph comparing frequent users PS vs !PS
out.data.one$cog <- factor(out.data.one$cog, levels=cog.names)
out.data.one$maruse <- factor(out.data.one$maruse, levels=c('nonuser','frequser'))
out.plot.three <- ggplot(out.data.one, aes(x=cog, y=mean, group=maruse, color=maruse)) +
geom_point(position=position_dodge(.3)) +
geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(.3)) +
facet_grid(age_cat~psBinary) +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
coord_cartesian(ylim=c(-.7,.6))


## the first loop will be looking at an interaction between ps and mj
out.data.one <- NULL
for(i in cog.names){
    tmp.dat <- summarySE(groupvars=c('maruse','psBinary','age_cat'), data=data.img.reduce, measurevar=i, na.rm=T)
    colnames(tmp.dat)[5] <- c('mean')
    tmp.dat$cog <- i
    out.data.one <- rbind(out.data.one, tmp.dat)
}

out.data.one <- out.data.one[which(out.data.one$maruse=='frequser' | out.data.one$maruse=='nonuser'),]

## Now plot em
## Now make the first graph comparing frequent users PS vs !PS
out.data.one$cog <- factor(out.data.one$cog, levels=cog.names)
out.data.one$maruse <- factor(out.data.one$maruse, levels=c('nonuser','frequser'))
out.plot.four <- ggplot(out.data.one, aes(x=cog, y=mean, group=maruse, color=maruse)) +
geom_point(position=position_dodge(.3)) +
geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(.3)) +
facet_grid(age_cat~psBinary) +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
coord_cartesian(ylim=c(-.7,.6))

## Now plot all of these
pdf('freqIsolCogVals.pdf', height=14, width=14)
out.plot.one
out.plot.two
out.plot.three
out.plot.four
dev.off()

### Now do all of the same plots, but put collapse marcat into binary
data$maruse[which(data$maruse=="user")] <- "frequser"

## Now go through a loop and prepare all of the values
## the first loop will be looking at an interaction between ps and mj
out.data.one <- NULL
for(i in cog.names){
    tmp.dat <- summarySE(groupvars=c('maruse','psBinary'), data=data, measurevar=i, na.rm=T)
    colnames(tmp.dat)[4] <- c('mean')
    tmp.dat$cog <- i
    out.data.one <- rbind(out.data.one, tmp.dat)
}

## Now clean the data
out.data.one <- out.data.one[which(out.data.one$maruse=='frequser' | out.data.one$maruse=='nonuser'),]
#out.data.one <- out.data.one[-which(out.data.one$goassessdxpmr6==''),]

## Now make the first graph comparing frequent users PS vs !PS
out.data.one$cog <- factor(out.data.one$cog, levels=cog.names)
out.data.one$maruse <- factor(out.data.one$maruse, levels=c('nonuser','frequser'))
out.plot.one <- ggplot(out.data.one, aes(x=cog, y=mean, group=maruse, color=maruse)) +
geom_point(position=position_dodge(.3)) +
geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(.3)) +
facet_grid(.~psBinary) +
theme(axis.text.x = element_text(angle = 90, hjust = 1))

## Now do the same analysis in the isolated img cohort
data.img.reduce <- data[which(data$bblid %in% img.data$bblid),]

out.data.one <- NULL
for(i in cog.names){
    tmp.dat <- summarySE(groupvars=c('maruse','psBinary'), data=data.img.reduce, measurevar=i, na.rm=T)
    colnames(tmp.dat)[4] <- c('mean')
    tmp.dat$cog <- i
    out.data.one <- rbind(out.data.one, tmp.dat)
}


## Now clean the data
out.data.one <- out.data.one[which(out.data.one$maruse=='frequser' | out.data.one$maruse=='nonuser'),]
#out.data.one <- out.data.one[-which(out.data.one$goassessdxpmr6==''),]

## Now make the first graph comparing frequent users PS vs !PS
out.data.one$cog <- factor(out.data.one$cog, levels=cog.names)
out.data.one$maruse <- factor(out.data.one$maruse, levels=c('nonuser','frequser'))
out.plot.two <- ggplot(out.data.one, aes(x=cog, y=mean, group=maruse, color=maruse)) +
geom_point(position=position_dodge(.3)) +
geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(.3)) +
facet_grid(.~psBinary) +
theme(axis.text.x = element_text(angle = 90, hjust = 1))

## Now do the same analysis w/ age bins
## the first loop will be looking at an interaction between ps and mj
out.data.one <- NULL
for(i in cog.names){
    tmp.dat <- summarySE(groupvars=c('maruse','psBinary','age_cat'), data=data, measurevar=i, na.rm=T)
    colnames(tmp.dat)[5] <- c('mean')
    tmp.dat$cog <- i
    out.data.one <- rbind(out.data.one, tmp.dat)
}

out.data.one <- out.data.one[which(out.data.one$maruse=='frequser' | out.data.one$maruse=='nonuser'),]

## Now plot em
## Now make the first graph comparing frequent users PS vs !PS
out.data.one$cog <- factor(out.data.one$cog, levels=cog.names)
out.data.one$maruse <- factor(out.data.one$maruse, levels=c('nonuser','frequser'))
out.plot.three <- ggplot(out.data.one, aes(x=cog, y=mean, group=maruse, color=maruse)) +
geom_point(position=position_dodge(.3)) +
geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(.3)) +
facet_grid(age_cat~psBinary) +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
coord_cartesian(ylim=c(-.7,.6))


## the first loop will be looking at an interaction between ps and mj
out.data.one <- NULL
for(i in cog.names){
    tmp.dat <- summarySE(groupvars=c('maruse','psBinary','age_cat'), data=data.img.reduce, measurevar=i, na.rm=T)
    colnames(tmp.dat)[5] <- c('mean')
    tmp.dat$cog <- i
    out.data.one <- rbind(out.data.one, tmp.dat)
}

out.data.one <- out.data.one[which(out.data.one$maruse=='frequser' | out.data.one$maruse=='nonuser'),]

## Now plot em
## Now make the first graph comparing frequent users PS vs !PS
out.data.one$cog <- factor(out.data.one$cog, levels=cog.names)
out.data.one$maruse <- factor(out.data.one$maruse, levels=c('nonuser','frequser'))
out.plot.four <- ggplot(out.data.one, aes(x=cog, y=mean, group=maruse, color=maruse)) +
geom_point(position=position_dodge(.3)) +
geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(.3)) +
facet_grid(age_cat~psBinary) +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
coord_cartesian(ylim=c(-.7,.6))

## Now plot all of these
pdf('freqCollapseCogVals.pdf', height=14, width=14)
out.plot.one
out.plot.two
out.plot.three
out.plot.four
dev.off()
