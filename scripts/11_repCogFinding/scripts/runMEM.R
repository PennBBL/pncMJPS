## Load library(s)
install_load('readstata13','ggplot2','lme4','lmerTest','reshape','psych','visreg')
source('../functions/functions.R')
## Load the data
data <- read.dta13('./cannabis_psychosis_cnb_foradon.dta')
## Now load the imaging data
img.data <- read.csv("~/Documents/dataRel/n2416_imagingDataDump_2018-04-22.csv")
img.data <- img.data[which(img.data$tpvalue==1),]
## Make a ps vs not ps group
data$psBinary <- 'NO'
data$psBinary[which(data$goassessdxpmr6=='PS')] <- 'YES'
data$mjbinary <- 'NO'
data$mjbinary[which(data$marcat!=1)] <- 'YES'

## Now melt the data frame, so it is ready for a MEM
char.vec <- c("bblid","envses","sex","race2","ageatcnb1","psBinary","mjbinary","psychosis_ar_4factor","marcat")
cog.names <- names(data)[c(71,72,73,74)]
x <- data[,c(char.vec,cog.names)]
xCog <- melt(x, id.vars=char.vec)
xCog$marcat <- factor(xCog$marcat)
## Now run the model
mod.cog <- lmerTest::lmer(value ~ ageatcnb1+sex+envses+mjbinary*psBinary*variable+(1|bblid),data=xCog,na.action=na.exclude)
mod.cog2 <- lmerTest::lmer(value ~ ageatcnb1+sex+envses+mjbinary*psychosis_ar_4factor*variable+(1|bblid),data=xCog,na.action=na.exclude)
mod.cog3 <- lmerTest::lmer(value ~ ageatcnb1+sex+envses+marcat*psychosis_ar_4factor*variable+(1|bblid),data=xCog,na.action=na.exclude)

## Now do the fine for accuracy
cog.names <- names(data)[c(65,66,67)]
x <- data[,c(char.vec,cog.names)]
xCog <- melt(x, id.vars=char.vec)
xCog$marcat <- factor(xCog$marcat)
## Now run the model
mod.cog <- lmerTest::lmer(value ~ ageatcnb1+sex+envses+mjbinary*psBinary*variable+(1|bblid),data=xCog,na.action=na.exclude)
mod.cog2 <- lmerTest::lmer(value ~ ageatcnb1+sex+envses+mjbinary*psychosis_ar_4factor*variable+(1|bblid),data=xCog,na.action=na.exclude)
mod.cog3 <- lmerTest::lmer(value ~ ageatcnb1+sex+envses+marcat*psychosis_ar_4factor*variable+(1|bblid),data=xCog,na.action=na.exclude)

## Now plot it
## Now do continous plots here
out.data.one <- NULL
data.tmp <- data
for(i in names(data)[c(65,66,67)]){
    data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + ageatcnb1 + envses")), data=data.tmp, na.action=na.exclude))
    data.tmp[,i] <- scale(data.tmp[,i])
}
tmp.data <- data.tmp[,c(char.vec,cog.names)]
xVol <- melt(tmp.data, id.vars=char.vec)
xVol <- xVol[-which(is.na(xVol$marcat)),]
out.plot.one <- ggplot(xVol, aes(x=psychosis_ar_4factor, y=value, group=marcat, color=factor(marcat))) +
geom_point() +
facet_grid(.~variable) +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
geom_smooth(method='lm')
pdf('marcatByPsychFactorInteractionCog.pdf', width=30, height=10)
out.plot.one
dev.off()

## Now do speed
cog.names <- names(data)[c(68,69,70)]
x <- data[,c(char.vec,cog.names)]
xCog <- melt(x, id.vars=char.vec)
xCog$marcat <- factor(xCog$marcat)
## Now run the model
mod.cog <- lmerTest::lmer(value ~ ageatcnb1+sex+envses+mjbinary*psBinary*variable+(1|bblid),data=xCog,na.action=na.exclude)
mod.cog2 <- lmerTest::lmer(value ~ ageatcnb1+sex+envses+mjbinary*psychosis_ar_4factor*variable+(1|bblid),data=xCog,na.action=na.exclude)
mod.cog3 <- lmerTest::lmer(value ~ ageatcnb1+sex+envses+marcat*psychosis_ar_4factor*variable+(1|bblid),data=xCog,na.action=na.exclude)


## Prep the imaging data
## Load the data
data <- read.dta13('./cannabis_psychosis_cnb_foradon.dta')
## Now load the imaging data
img.data <- read.csv("~/Documents/dataRel/n2416_imagingDataDump_2018-09-20.csv")
img.data <- img.data[which(img.data$tpvalue==1),]
## Make a ps vs not ps group
data$psBinary <- 'NO'
data$psBinary[which(data$goassessdxpmr6=='PS')] <- 'YES'
data$mjbinary <- 'NO'
data$mjbinary[which(data$marcat!=1)] <- 'YES'
img.data <- img.data[which(img.data$tpvalue==1),]
img.data <- merge(data, img.data, by='bblid',suffixes = c("",".y"))
img.data <- img.data[-which(is.na(img.data$marcat)),]
## Now run the imaging analyses
char.vec <- c("bblid","envses","sex","race2","scanageMonths","psBinary","mjbinary","psychosis_ar_4factor","marcat")
global.val <- c("mprage_jlf_vol_TBWM","mprage_jlf_gmd_MeanGMD","dti_jlf_tr_MeanWholeBrainTR")
x <- img.data[,c(char.vec,global.val)]
## Now scale the values w/in modality
x[,10:dim(x)[2]] <- scale(x[,10:dim(x)[2]])
xCog <- melt(x, id.vars=char.vec)
xCog$marcat <- factor(xCog$marcat)
## Now run the model
mod.glo <- nlme::lme(value~scanageMonths+sex+envses+mjbinary*psBinary*variable,random=~1|bblid,data=xCog,na.action=na.exclude)
mod.glo2 <- nlme::lme(value~scanageMonths+sex+envses+mjbinary*psychosis_ar_4factor*variable,random=~1|bblid,data=xCog,na.action=na.exclude)
mod.glo3 <- nlme::lme(value~scanageMonths+sex+envses+marcat*psychosis_ar_4factor*variable,random=~1|bblid,data=xCog,na.action=na.exclude)

## Do volume here
tmp <- names(img.data)[c(grep("mprage_jlfHiLoLobe_vol", names(img.data)))]
xVol <- img.data[,c(char.vec, tmp)]
xVol[,10:18] <- scale(xVol[,10:18])
xVol <- melt(xVol, id.vars=char.vec)
xVol <- xVol[complete.cases(xVol),]
mod.vol3 <- nlme::lme(value~scanageMonths+sex+envses+marcat*psychosis_ar_4factor*variable,random=~1|bblid,data=xVol,na.action=na.exclude)

## Now do all of the individual volume HiLo lobes
for(lobeVal in c(1:9)){
    colVals <- names(img.data)[349:487]
    dim(colVals) <- c(139,1)
    tmp <- colVals[which(apply(colVals, 1, findLobe)==lobeVal),]
    xVol <- img.data[,c(char.vec, tmp)]
    xVol[,10:dim(xVol)[2]] <- scale(xVol[,10:dim(xVol)[2]])
    if(lobeVal!=8){
        xVol <- averageLeftAndRightVol(xVol)
    }
    xVol <- melt(xVol, id.vars=char.vec)
    mod.vol.roi <-nlme::lme(value~scanageMonths+sex+envses+marcat*psychosis_ar_4factor*variable,random=~1|bblid,data=xVol,na.action=na.exclude)
    assign(paste0("mod.vol.roi_",lobeVal),mod.vol.roi)
}

## Now onto gmd
tmp <- names(img.data)[c(grep("mprage_jlfHiLoLobe_gmd", names(img.data)))]
xGmd <- img.data[,c(char.vec, tmp)]
xGmd[,10:17] <- scale(xGmd[,10:17])
xGmd <- melt(xGmd, id.vars=char.vec)
mod.gmd3 <- nlme::lme(value~scanageMonths+sex+envses+marcat*psychosis_ar_4factor*variable,random=~1|bblid,data=xGmd,na.action=na.exclude)


## Now check the DTI lobes
tmp <- names(img.data)[c(grep("dti_jlfHiLoLobe_tr", names(img.data)))][-9]
xTr <- img.data[,c(char.vec, tmp)]
xTr[,10:17] <- scale(xTr[,10:17])
xTr <- melt(xTr, id.vars=char.vec)
mod.tr3 <- nlme::lme(value~scanageMonths+sex+envses+marcat*psychosis_ar_4factor*variable,random=~1|bblid,data=xTr,na.action=na.exclude)
for(lobeVal in c(1:9)){
    colVals <- names(img.data)[956:1085]
    dim(colVals) <- c(130,1)
    tmp <- colVals[which(apply(colVals, 1, findLobe)==lobeVal),]
    xVol <- img.data[,c(char.vec, tmp)]
    xVol[,10:dim(xVol)[2]] <- scale(xVol[,10:dim(xVol)[2]])
    if(lobeVal!=8){
        #xVol <- averageLeftAndRightVol(xVol)
    }
    xVol <- melt(xVol, id.vars=char.vec)
    mod.tr.roi <-nlme::lme(value~scanageMonths+sex+envses+marcat*psychosis_ar_4factor*variable,random=~1|bblid,data=xVol,na.action=na.exclude)
    assign(paste0("mod.tr.roi_",lobeVal),mod.tr.roi)
    ## Now write the anove values
    outputFileName <- paste("anovaTRLobeValue", lobeVal, ".csv", sep='')
    to.write <- anova(mod.tr.roi)
    write.csv(to.write, outputFileName, quote=F)
}

## Now do FA
tmp <- names(img.data)[c(grep("dti_dtitk_jhutract_fa", names(img.data)))]
xFA <- img.data[,c(char.vec, tmp)]
#xFA <- averageLeftAndRight1(xFA)
xFA[,10:27] <- scale(xFA[,10:27])
xFA <- melt(xFA, id.vars=char.vec)
mod.fa3 <- nlme::lme(value~scanageMonths+sex+envses+marcat*psychosis_ar_4factor*variable,random=~1|bblid,data=xFA,na.action=na.exclude)

## Looks like we have some significant ROI level psy*MJ interactions
## Now we are going to have to plot these effects
## So let the plotting begin down here


## Lets begin here with the volume lobes
out.data.one <- NULL
for(i in names(img.data)[c(grep("mprage_jlfHiLoLobe_vol", names(img.data)))]){
    img.data.tmp <- img.data
    img.data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + scanageMonths")), data=img.data.tmp, na.action=na.exclude))
    img.data.tmp[,i] <- scale(img.data.tmp[,i])
    tmp.dat <- summarySE(groupvars=c('maruse','psBinary'), data=img.data.tmp, measurevar=i, na.rm=T)
    colnames(tmp.dat)[4] <- c('mean')
    tmp.dat$lobe <- i
    out.data.one <- rbind(out.data.one, tmp.dat)
}
out.data.one$maruse <- factor(out.data.one$maruse, levels=c('nonuser','user','frequser'))
out.plot.one <- ggplot(out.data.one, aes(x=lobe, y=mean, group=maruse, color=maruse)) +
geom_point(position=position_dodge(.3)) +
geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(.3)) +
facet_grid(psBinary~maruse) +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
coord_cartesian(ylim=c(-.7,.6))

## Now do continous plots here
out.data.one <- NULL
img.data.tmp <- img.data
for(i in names(img.data)[c(grep("mprage_jlfHiLoLobe_vol", names(img.data)))]){
    img.data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + scanageMonths")), data=img.data.tmp, na.action=na.exclude))
    img.data.tmp[,i] <- scale(img.data.tmp[,i])
}
tmp.data <- img.data.tmp[,c(char.vec, names(img.data)[c(grep("mprage_jlfHiLoLobe_vol", names(img.data)))])]
xVol <- melt(tmp.data, id.vars=char.vec)
out.plot.one <- ggplot(xVol, aes(x=psychosis_ar_4factor, y=value, group=marcat, color=factor(marcat))) +
geom_point() +
facet_grid(.~variable) +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
coord_cartesian(ylim=c(-.7,.6)) +
geom_smooth(method='lm')

## Now do mj binary
out.data.one <- NULL
for(i in names(img.data)[c(grep("mprage_jlfHiLoLobe_vol", names(img.data)))]){
    img.data.tmp <- img.data
    img.data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + scanageMonths")), data=img.data.tmp, na.action=na.exclude))
    img.data.tmp[,i] <- scale(img.data.tmp[,i])
    tmp.dat <- summarySE(groupvars=c('mjbinary','psBinary'), data=img.data.tmp, measurevar=i, na.rm=T)
    colnames(tmp.dat)[4] <- c('mean')
    tmp.dat$lobe <- i
    out.data.one <- rbind(out.data.one, tmp.dat)
}
out.plot.one <- ggplot(out.data.one, aes(x=lobe, y=mean, group=mjbinary, color=mjbinary)) +
geom_point(position=position_dodge(.3)) +
geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(.3)) +
facet_grid(psBinary~mjbinary) +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
coord_cartesian(ylim=c(-.7,.6))

## Now do TBV
i <- 'mprage_jlf_vol_TBV'
img.data.tmp <- img.data
img.data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + scanageMonths")), data=img.data.tmp, na.action=na.exclude))
img.data.tmp[,i] <- scale(img.data.tmp[,i])
out.plot.one <- ggplot(img.data.tmp, aes(x=psychosis_ar_4factor, y=mprage_jlf_vol_TBV, group=marcat, color=factor(marcat))) +
geom_point(position=position_dodge(.3)) +
geom_smooth(method='lm')

## Now do TR
out.data.one <- NULL
for(i in names(img.data)[c(grep("dti_jlfHiLoLobe_tr", names(img.data)))]){
    img.data.tmp <- img.data
    img.data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + scanageMonths")), data=img.data.tmp, na.action=na.exclude))
    img.data.tmp[,i] <- scale(img.data.tmp[,i])
    tmp.dat <- summarySE(groupvars=c('maruse','psBinary'), data=img.data.tmp, measurevar=i, na.rm=T)
    colnames(tmp.dat)[4] <- c('mean')
    tmp.dat$lobe <- i
    out.data.one <- rbind(out.data.one, tmp.dat)
}
out.data.one$maruse <- factor(out.data.one$maruse, levels=c('nonuser','user','frequser'))
out.plot.two <- ggplot(out.data.one, aes(x=lobe, y=mean, group=maruse, color=maruse)) +
geom_point(position=position_dodge(.3)) +
geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(.3)) +
facet_grid(psBinary~maruse) +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
coord_cartesian(ylim=c(-.7,.6))

out.data.one <- NULL
for(i in names(img.data)[c(grep("dti_jlfHiLoLobe_tr", names(img.data)))]){
    img.data.tmp <- img.data
    img.data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + scanageMonths")), data=img.data.tmp, na.action=na.exclude))
    img.data.tmp[,i] <- scale(img.data.tmp[,i])
    tmp.dat <- summarySE(groupvars=c('mjbinary','psBinary'), data=img.data.tmp, measurevar=i, na.rm=T)
    colnames(tmp.dat)[4] <- c('mean')
    tmp.dat$lobe <- i
    out.data.one <- rbind(out.data.one, tmp.dat)
}
out.plot.one <- ggplot(out.data.one, aes(x=lobe, y=mean, group=mjbinary, color=mjbinary)) +
geom_point(position=position_dodge(.3)) +
geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(.3)) +
facet_grid(psBinary~mjbinary) +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
coord_cartesian(ylim=c(-.7,.6))

i <- 'dti_jlf_tr_MeanWholeBrainTR'
img.data.tmp <- img.data
img.data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + scanageMonths + dti64Tsnr + envses")), data=img.data.tmp, na.action=na.exclude))
img.data.tmp[,i] <- scale(img.data.tmp[,i])
out.plot.one <- ggplot(img.data.tmp, aes(x=psychosis_ar_4factor, y=dti_jlf_tr_MeanWholeBrainTR, group=marcat, color=factor(marcat))) +
geom_point(position=position_dodge(.3)) +
geom_smooth(method='lm')

## Now do continous plots here
out.data.one <- NULL
img.data.tmp <- img.data
for(i in names(img.data)[c(grep("dti_jlfHiLoLobe_tr", names(img.data)))]){
    img.data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + scanageMonths")), data=img.data.tmp, na.action=na.exclude))
    img.data.tmp[,i] <- scale(img.data.tmp[,i])
}
tmp.data <- img.data.tmp[,c(char.vec, names(img.data)[c(grep("dti_jlfHiLoLobe_tr", names(img.data)))])]
xVol <- melt(tmp.data, id.vars=char.vec)
out.plot.one <- ggplot(xVol, aes(x=psychosis_ar_4factor, y=value, group=marcat, color=factor(marcat))) +
geom_point() +
facet_grid(.~variable) +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
coord_cartesian(ylim=c(-.7,.6)) +
geom_smooth(method='lm')
pdf('lobularMarcatVsPsych.pdf', width=30, height=10)
out.plot.one
dev.off()
