## Load library(s)
install_load('readstata13','ggplot2','lme4','lmerTest','reshape','psych','visreg')
source('../functions/functions.R')
## Load the data
data <- read.dta13('./cannabis_psychosis_cnb_foradon.dta')
## Now load the imaging data
img.data <- read.csv("~/Documents/dataRel/n1601_imagingDataDump_2018-09-20.csv")
img.data <- img.data[which(img.data$tpvalue==1),]
## Make a ps vs not ps group
data$psBinary <- 'NO'
data$psBinary[which(data$goassessdxpmr6=='PS')] <- 'YES'
data$mjbinary <- 'NO'
data$mjbinary[which(data$marcat!=1)] <- 'YES'
data <- data[-which(data$mjpastyr==1),]

## Now load the RDS file as our data
data <-readRDS("../../01_dataPrep/scripts/mjPSCogImg.RDS")

## Now melt the data frame, so it is ready for a MEM
char.vec <- c("bblid","envses","sex","race2","ageatcnb1","psBinary","mjbinary","psychosis_ar_4factor","marcat")
cog.names <- c("f1_social_cognition_efficiency","f2_complex_reasoning_efficiency","f3_memory_efficiency","f4_executive_efficiency")
x <- data[,c(char.vec,cog.names)]
xCog <- melt(x, id.vars=char.vec)
xCog$marcat <- factor(xCog$marcat)
## Now run the model
mod.cog <- lmerTest::lmer(value ~ ageatcnb1+sex+envses+mjbinary*psBinary*variable+(1|bblid),data=xCog,na.action=na.exclude)
mod.cog2 <- lmerTest::lmer(value ~ ageatcnb1+sex+envses+mjbinary*psychosis_ar_4factor*variable+(1|bblid),data=xCog,na.action=na.exclude)
mod.cog3 <- lmerTest::lmer(value ~ ageatcnb1+sex+envses+marcat*psychosis_ar_4factor*variable+(1|bblid),data=xCog,na.action=na.exclude)

## Now do the fine for accuracy
cog.names <- c("f1_exec_comp_cog_accuracy","f2_social_cog_accuracy","f3_memory_accuracy")
x <- data[,c(char.vec,cog.names)]
xCog <- melt(x, id.vars=char.vec)
xCog$marcat <- factor(xCog$marcat)
## Now run the model
mod.cog <- lmerTest::lmer(value ~ ageatcnb1+sex+envses+mjbinary*psBinary*variable+(1|bblid),data=xCog,na.action=na.exclude)
mod.cog2 <- lmerTest::lmer(value ~ ageatcnb1+sex+envses+mjbinary*psychosis_ar_4factor*variable+(1|bblid),data=xCog,na.action=na.exclude)
mod.cog3 <- lmerTest::lmer(value ~ ageatcnb1+sex+envses+(marcat+psychosis_ar_4factor+variable)^3+(1|bblid),data=xCog,na.action=na.exclude)
## Now plot it
## Now do continous plots here
out.data.one <- NULL
data.tmp <- data
for(i in c("f1_exec_comp_cog_accuracy","f2_social_cog_accuracy","f3_memory_accuracy")){
    data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + ageatcnb1")), data=data.tmp, na.action=na.exclude))
    data.tmp[,i] <- scale(data.tmp[,i])
}
tmp.data <- data.tmp[,c(char.vec,cog.names)]
xVol <- melt(tmp.data, id.vars=char.vec)
xVol <- xVol[-which(is.na(xVol$marcat)),]
out.plot.one <- ggplot(xVol, aes(x=psychosis_ar_4factor, y=value, group=marcat, color=factor(marcat))) +
#geom_point() +
facet_grid(.~variable) +
theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
geom_smooth(method='lm') +

pdf('marcatByPsychFactorInteractionCog.pdf', width=30, height=10)
out.plot.one
dev.off()

## It looks like this effect is being driven by F1 exec comp cog
## so I will focus in on that effect
data.tmp <- data
data.tmp$None <- residuals(lm(as.formula(paste("f1_exec_comp_cog_accuracy"," ~ sex + ageatcnb1")), data=data.tmp, na.action=na.exclude))
data.tmp$onlySES <- residuals(lm(as.formula(paste("f1_exec_comp_cog_accuracy"," ~ sex + ageatcnb1+envses")), data=data.tmp, na.action=na.exclude))
data.tmp$onlyRace <- residuals(lm(as.formula(paste("f1_exec_comp_cog_accuracy"," ~ sex + ageatcnb1+factor(race2)")), data=data.tmp, na.action=na.exclude))
data.tmp$Both <- residuals(lm(as.formula(paste("f1_exec_comp_cog_accuracy"," ~ sex + ageatcnb1+envses+factor(race2)")), data=data.tmp, na.action=na.exclude))
data.tmp$marcat <- factor(data.tmp$marcat,levels=c('NU','OU','FU'))
data.tmp <- data.tmp[-which(is.na(data.tmp$marcat)),]
data.tmp <- data.tmp[,c(char.vec,'None','onlySES','onlyRace','Both')]
data.tmp <- melt(data.tmp, id.vars=c(char.vec))
out.plot.one.cog <- ggplot(data.tmp, aes(x=psychosis_ar_4factor,y=value, group=marcat, color=marcat)) +
#    geom_point() +
    geom_smooth(method='lm') +
    xlab("Psychosis Loading") +
    ylab("Executive Function") +
    scale_colour_manual(name = "marcat",values=c("NU"="Red","OU"="Green","FU"="Blue")) +
    theme_bw() +
theme(legend.position="bottom",text = element_text(size=28,face='bold')) +
facet_grid(.~variable)
png("cogInteraction.png", width=32, height=12, units='in', res=300)
out.plot.one.cog
dev.off()
## Now do speed
cog.names <- c("f1_slow_speed","f2_memory_speed","f3_fast_speed")
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
img.data <- read.csv("~/Documents/dataRel/n1601_imagingDataDump_2018-09-20.csv")
old_img <-  readRDS("./mjAnovaData.RDS")
#img.data <- img.data[which(img.data$tpvalue==1),]
## Make a ps vs not ps group
data$psBinary <- 'NO'
data$psBinary[which(data$goassessdxpmr6=='PS')] <- 'YES'
data$mjbinary <- 'NO'
data$mjbinary[which(data$marcat!=1)] <- 'YES'
img.data <- img.data[img.data$bblid %in% old_img$bblid,]
img.data <- merge(data, img.data, by='bblid',suffixes = c("",".y"))
img.data <- img.data[-which(is.na(img.data$marcat)),]

## Now load the RDS
img.data <-readRDS("../../01_dataPrep/scripts/mjPSCogImgIsol.RDS")
#img.data <- img.data[-which(is.na(img.data$marcat)),]
## Now run the imaging analyses
char.vec <- c("bblid","envses","sex","race2","ageAtScan1","psBinary","mjbinary","psychosis_ar_4factor","marcat")
global.val <- c("mprage_jlf_vol_TBGM","dti_jlf_tr_MeanWholeBrainTR","pcasl_jlf_cbf_MeanGMCBF","mprage_jlf_gmd_MeanGMD","rest_jlf_alff_MeanALFF","rest_jlf_reho_MeanReho")
x <- img.data[,c(char.vec,global.val)]
## Now scale the values w/in modality
x[,10:dim(x)[2]] <- scale(x[,10:dim(x)[2]])
xCog <- melt(x, id.vars=char.vec)
xCog$marcat <- factor(xCog$marcat)
## Now run the model
mod.glo <- nlme::lme(value~ageAtScan1+sex+envses+mjbinary*psBinary*variable,random=~1|bblid,data=xCog,na.action=na.exclude)
mod.glo2 <- nlme::lme(value~ageAtScan1+sex+envses+mjbinary*psychosis_ar_4factor*variable,random=~1|bblid,data=xCog,na.action=na.exclude)
mod.glo3 <- nlme::lme(value~ageAtScan1+sex+envses+marcat*psychosis_ar_4factor*variable,random=~1|bblid,data=xCog,na.action=na.exclude)

## Do volume here
tmp <- names(img.data)[c(grep("mprage_jlfHiLoLobe_vol", names(img.data)))]
xVol <- img.data[,c(char.vec, tmp)]
xVol[,10:18] <- scale(xVol[,10:18])
xVol <- melt(xVol, id.vars=char.vec)
xVol <- xVol[complete.cases(xVol),]
mod.vol3 <- nlme::lme(value~ageAtScan1+sex+envses+marcat*psychosis_ar_4factor*variable,random=~1|bblid,data=xVol,na.action=na.exclude)

## Now do all of the individual volume HiLo lobes
pdf('regionEffects.pdf', width=60, height=12)
for(lobeVal in c(1:9)){
    colVals <- names(img.data)[357:495]
    dim(colVals) <- c(139,1)
    tmp <- colVals[which(apply(colVals, 1, findLobe)==lobeVal),]
    xVol <- img.data[,c(char.vec, tmp)]
    xVol[,10:dim(xVol)[2]] <- scale(xVol[,10:dim(xVol)[2]])
    if(lobeVal!=8){
        #xVol <- averageLeftAndRightVol(xVol)
    }
    xVol <- melt(xVol, id.vars=char.vec)
    mod.vol.roi <-nlme::lme(value~ageAtScan1+sex+envses+marcat*psychosis_ar_4factor*variable,random=~1|bblid,data=xVol,na.action=na.exclude)
    assign(paste0("mod.vol.roi_",lobeVal),mod.vol.roi)
    ## Now write the anove values
    outputFileName <- paste("anovaVolLobeValue", lobeVal, ".csv", sep='')
    to.write <- anova(mod.vol.roi)
    #write.csv(to.write, outputFileName, quote=F)
    ## Now if we have a significnat marcat interaction lets plot the ROI's
    ## FInd the minimum p value between our variables of interest
    min.p.val <- min(anova(mod.vol.roi)[c('marcat:psychosis_ar_4factor','marcat:psychosis_ar_4factor:variable'),'p-value'])
    if(min.p.val<.05){
        write.csv(to.write, outputFileName, quote=F)
        ## Now plot these ROI's after removing
        dim(colVals) <- c(139,1)
        tmp <- colVals[which(apply(colVals, 1, findLobe)==lobeVal),]
        xVol <- img.data[,c(char.vec, tmp)]
        xVol[,10:dim(xVol)[2]] <- scale(xVol[,10:dim(xVol)[2]])
        if(lobeVal!=8){
            #xVol <- averageLeftAndRightVol(xVol)
        }
        ## Now remove age effects
        xVol[,10:dim(xVol)[2]] <- apply(xVol[,10:dim(xVol)[2]],2,function(x) residuals(lm(x~ageAtScan1+sex,data=xVol,na.action=na.exclude)))
        xVol <- melt(xVol, id.vars=char.vec)
        xVol$marcat <- factor(xVol$marcat,levels=c('NU','OU','FU'))
        xVol <- xVol[-which(is.na(xVol$marcat)),]
        out.plot.one <- ggplot(xVol, aes(x=psychosis_ar_4factor, y=value, group=marcat, color=factor(marcat))) +
        #geom_point() +
        facet_grid(.~variable) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        geom_smooth(method='lm')
        print(out.plot.one)
    }
}
dev.off()
## Now onto gmd
tmp <- names(img.data)[c(grep("mprage_jlfHiLoLobe_gmd", names(img.data)))]
xGmd <- img.data[,c(char.vec, tmp)]
xGmd[,10:17] <- scale(xGmd[,10:17])
xGmd <- melt(xGmd, id.vars=char.vec)
mod.gmd3 <- nlme::lme(value~ageAtScan1+sex+envses+marcat*psychosis_ar_4factor*variable,random=~1|bblid,data=xGmd,na.action=na.exclude)

## Now check the DTI lobes
tmp <- names(img.data)[c(grep("dti_jlfHiLoLobe_tr", names(img.data)))]
xTr <- img.data[,c(char.vec, tmp)]
xTr[,10:17] <- scale(xTr[,10:17])
xTr <- melt(xTr, id.vars=char.vec)
mod.tr3 <- nlme::lme(value~ageAtScan1+sex+envses+marcat*psychosis_ar_4factor*variable,random=~1|bblid,data=xTr,na.action=na.exclude)
pdf('regionEffectsTR.pdf', width=60, height=12)
for(lobeVal in c(1:9)){
    colVals <- names(img.data)[965:1095]
    dim(colVals) <- c(131,1)
    tmp <- colVals[which(apply(colVals, 1, findLobe)==lobeVal),]
    xVol <- img.data[,c(char.vec, tmp)]
    xVol[,10:dim(xVol)[2]] <- scale(xVol[,10:dim(xVol)[2]])
    if(lobeVal!=8){
        #xVol <- averageLeftAndRightVol(xVol)
    }
    xVol <- melt(xVol, id.vars=char.vec)
    mod.tr.roi <-nlme::lme(value~ageAtScan1+sex+envses+marcat*psychosis_ar_4factor*variable,random=~1|bblid,data=xVol,na.action=na.exclude)
    assign(paste0("mod.tr.roi_",lobeVal),mod.tr.roi)
    ## Now write the anove values
    outputFileName <- paste("anovaTRLobeValue", lobeVal, ".csv", sep='')
    to.write <- anova(mod.tr.roi)
    print(paste(lobeVal, anova(mod.tr.roi)[c('marcat:psychosis_ar_4factor','marcat:psychosis_ar_4factor:variable'),]))
    #write.csv(to.write, outputFileName, quote=F)
    ## Now if we have a significnat marcat interaction lets plot the ROI's
    ## FInd the minimum p value between our variables of interest
    min.p.val <- min(anova(mod.tr.roi)[c('marcat:psychosis_ar_4factor','marcat:psychosis_ar_4factor:variable'),'p-value'])
    if(min.p.val<.005){
        write.csv(to.write, outputFileName, quote=F)
        ## Now plot these ROI's after removing
        dim(colVals) <- c(131,1)
        tmp <- colVals[which(apply(colVals, 1, findLobe)==lobeVal),]
        xVol <- img.data[,c(char.vec, tmp)]
        xVol[,10:dim(xVol)[2]] <- scale(xVol[,10:dim(xVol)[2]])
        if(lobeVal!=8){
            #xVol <- averageLeftAndRightVol(xVol)
        }
        ## Now remove age effects
        xVol[,10:dim(xVol)[2]] <- apply(xVol[,10:dim(xVol)[2]],2,function(x) residuals(lm(x~ageAtScan1+sex,data=xVol,na.action=na.exclude)))
        xVol <- melt(xVol, id.vars=char.vec)
        xVol$marcat <- factor(xVol$marcat,levels=c('NU','OU','FU'))
        xVol <- xVol[-which(is.na(xVol$marcat)),]
        out.plot.one <- ggplot(xVol, aes(x=psychosis_ar_4factor, y=value, group=marcat, color=factor(marcat))) +
        #geom_point() +
        facet_grid(.~variable) +
        #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        geom_smooth(method='lm')
        print(out.plot.one)
    }
}
dev.off()

## Now do reho
tmp <- names(img.data)[c(grep("rest_jlfHiLoLobe_reho", names(img.data)))]
xReho <- img.data[,c(char.vec, tmp)]
xReho[,10:17] <- scale(xReho[,10:17])
xReho <- melt(xReho, id.vars=char.vec)
mod.reho3 <- nlme::lme(value~ageAtScan1+sex+envses+marcat*psychosis_ar_4factor*variable,random=~1|bblid,data=xReho,na.action=na.exclude)

## Now do alff
tmp <- names(img.data)[c(grep("rest_jlfHiLoLobe_alff", names(img.data)))]
xAlff <- img.data[,c(char.vec, tmp)]
xAlff[,10:17] <- scale(xAlff[,10:17])
xAlff <- melt(xAlff, id.vars=char.vec)
mod.alff3 <- nlme::lme(value~ageAtScan1+sex+envses+marcat*psychosis_ar_4factor*variable,random=~1|bblid,data=xAlff,na.action=na.exclude)

## Looks like we have some significant ROI level psy*MJ interactions
## Now we are going to have to plot these effects
## So let the plotting begin down here
## Lets begin here with the volume lobes
out.data.one <- NULL
for(i in names(img.data)[c(grep("mprage_jlfHiLoLobe_vol", names(img.data)))]){
    img.data.tmp <- img.data
    img.data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + ageAtScan1 + envses")), data=img.data.tmp, na.action=na.exclude))
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
    img.data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + ageAtScan1")), data=img.data.tmp, na.action=na.exclude))
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
    img.data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + ageAtScan1")), data=img.data.tmp, na.action=na.exclude))
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
anova(lm(mprage_jlf_vol_TBV ~ ageatcnb1+sex+envses+(marcat+psychosis_ar_4factor)^2, data=data))
i <- 'mprage_jlf_vol_TBV'
img.data.tmp <- img.data
img.data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + ageAtScan1")), data=img.data.tmp, na.action=na.exclude))
img.data.tmp[,i] <- scale(img.data.tmp[,i])
img.data.tmp$marcat <- factor(img.data.tmp$marcat,levels=c('NU','OU','FU'))
img.data.tmp <- img.data.tmp[-which(is.na(img.data.tmp$marcat)),]
out.plot.one <- ggplot(img.data.tmp, aes(x=psychosis_ar_4factor, y=mprage_jlf_vol_TBV, group=marcat, color=factor(marcat))) +
scale_colour_manual(name = "marcat",values=c("NU"="Red","OU"="Green","FU"="Blue")) +
geom_smooth(method='lm') +
xlab("Psychosis Loading") +
ylab("Total Brain Volume") +
theme(legend.position="none",text = element_text(size=28,face='bold')) +
annotate("text", -Inf, Inf, label = "F = 6.0", hjust = 0, vjust = 25.5, size=12) +
annotate("text", -Inf, Inf, label = "P < .001", hjust = 0, vjust = 27.5, size=12)
## Now plot TBV and cognition
out.plot.two <- ggplot(img.data.tmp, aes(y=F1_Exec_Comp_Cog_Accuracy, x=mprage_jlf_vol_TBV)) +
geom_smooth(method='lm') +
geom_point() +
xlab("Executive Function") +
xlab("Total Brain Volume") +
theme(legend.position="none",text = element_text(size=28,face='bold')) +
annotate("text", -Inf, Inf, label = "r = 0.40", hjust = 0, vjust = 25.5, size=12) +
annotate("text", -Inf, Inf, label = "P < .001", hjust = 0, vjust = 27.5, size=12)
png("imageDataTBV.png", width=24, height=12, units='in', res=300)
multiplot(out.plot.one, out.plot.two, cols=2)
dev.off()

## Now the strongest volume interaction, the temporal lobe
anova(lm(mprage_jlfHiLoLobe_vol_Temporal ~ ageatcnb1+sex+envses+(marcat+psychosis_ar_4factor)^2, data=data))
i <- 'mprage_jlfHiLoLobe_vol_Temporal'
img.data.tmp <- img.data
img.data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + ageAtScan1")), data=img.data.tmp, na.action=na.exclude))
img.data.tmp[,i] <- scale(img.data.tmp[,i])
img.data.tmp$marcat <- factor(img.data.tmp$marcat,levels=c('NU','OU','FU'))
img.data.tmp <- img.data.tmp[-which(is.na(img.data.tmp$marcat)),]
out.plot.one <- ggplot(img.data.tmp, aes(x=psychosis_ar_4factor, y=mprage_jlfHiLoLobe_vol_Temporal, group=marcat, color=factor(marcat))) +
scale_colour_manual(name = "marcat",values=c("NU"="Red","OU"="Green","FU"="Blue")) +
geom_smooth(method='lm') +
xlab("Psychosis Loading") +
ylab("Temporal Lobe Volume") +
theme(legend.position="none",text = element_text(size=28,face='bold')) +
annotate("text", -Inf, Inf, label = "F = 8.4", hjust = 0, vjust = 12.5, size=12) +
annotate("text", -Inf, Inf, label = "P < .001", hjust = 0, vjust = 14.5, size=12)
# Now plot it
png("lobarEffectVol.png")
out.plot.one
dev.off()

## Now do TR
out.data.one <- NULL
for(i in names(img.data)[c(grep("dti_jlfHiLoLobe_tr", names(img.data)))]){
    img.data.tmp <- img.data
    img.data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + ageAtScan1")), data=img.data.tmp, na.action=na.exclude))
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
    img.data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + ageAtScan1")), data=img.data.tmp, na.action=na.exclude))
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
img.data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + ageAtScan1 + dti64Tsnr + envses")), data=img.data.tmp, na.action=na.exclude))
img.data.tmp[,i] <- scale(img.data.tmp[,i])
img.data.tmp$marcat <- factor(img.data.tmp$marcat,levels=c('NU','OU','FU'))
img.data.tmp <- img.data.tmp[-which(is.na(img.data.tmp$marcat)),]
anova(lm(dti_jlf_tr_MeanWholeBrainTR ~ ageatcnb1+sex+envses+(marcat+psychosis_ar_4factor)^2, data=data))
out.plot.one <- ggplot(img.data.tmp, aes(x=psychosis_ar_4factor, y=dti_jlf_tr_MeanWholeBrainTR, group=marcat, color=factor(marcat))) +
scale_colour_manual(name = "marcat",values=c("NU"="Red","OU"="Green","FU"="Blue")) +
geom_smooth(method='lm') +
xlab("Psychosis Loading") +
ylab("Whole Brain Mean Diffusivity") +
theme(legend.position="none",text = element_text(size=28,face='bold')) +
annotate("text", -Inf, Inf, label = "F = 4.7", hjust = 0, vjust = 25.5, size=12) +
annotate("text", -Inf, Inf, label = "P < .01", hjust = 0, vjust = 27.5, size=12)
## Now plot TBV and cognition
out.plot.two <- ggplot(img.data.tmp, aes(y=F1_Exec_Comp_Cog_Accuracy, x=dti_jlf_tr_MeanWholeBrainTR)) +
geom_smooth(method='lm') +
geom_point() +
xlab("Executive Function") +
xlab("Whole Brain Mean Diffusivity") +
theme(legend.position="none",text = element_text(size=28,face='bold')) +
annotate("text", -Inf, Inf, label = "r = -0.2", hjust = 0, vjust = 25.5, size=12) +
annotate("text", -Inf, Inf, label = "P < .001", hjust = 0, vjust = 27.5, size=12)
png("imageDataMD.png", width=24, height=12, units='in', res=300)
multiplot(out.plot.one, out.plot.two, cols=2)
dev.off()


## Now do the strongest lobar relationship
anova(lm(dti_jlfHiLoLobe_tr_Occipital ~ ageatcnb1+sex+envses+(marcat+psychosis_ar_4factor)^2, data=data))
i <- 'dti_jlfHiLoLobe_tr_Occipital'
img.data.tmp <- img.data
img.data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + ageAtScan1 + dti64Tsnr + envses")), data=img.data.tmp, na.action=na.exclude))
img.data.tmp[,i] <- scale(img.data.tmp[,i])
img.data.tmp$marcat <- factor(img.data.tmp$marcat,levels=c('NU','OU','FU'))
img.data.tmp <- img.data.tmp[-which(is.na(img.data.tmp$marcat)),]
anova(lm(dti_jlf_tr_MeanWholeBrainTR ~ ageatcnb1+sex+envses+(marcat+psychosis_ar_4factor)^2, data=data))
out.plot.one <- ggplot(img.data.tmp, aes(x=psychosis_ar_4factor, y=dti_jlfHiLoLobe_tr_Occipital, group=marcat, color=factor(marcat))) +
scale_colour_manual(name = "marcat",values=c("NU"="Red","OU"="Green","FU"="Blue")) +
geom_smooth(method='lm') +
xlab("Psychosis Loading") +
ylab("Occipital Lobe Mean Diffusivity") +
theme(legend.position="none",text = element_text(size=28,face='bold')) +
annotate("text", -Inf, Inf, label = "F = 5.3", hjust = 0, vjust = 12.5, size=12) +
annotate("text", -Inf, Inf, label = "P < .005", hjust = 0, vjust = 14.5, size=12)

png("lobarEffectMD.png")
out.plot.one
dev.off()

## Now do continous plots here
out.data.one <- NULL
img.data.tmp <- img.data
for(i in names(img.data)[c(grep("dti_jlfHiLoLobe_tr", names(img.data)))]){
    img.data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + ageAtScan1")), data=img.data.tmp, na.action=na.exclude))
    img.data.tmp[,i] <- scale(img.data.tmp[,i])
}
tmp.data <- img.data.tmp[,c(char.vec, names(img.data)[c(grep("dti_jlfHiLoLobe_tr", names(img.data)))])]
xVol <- melt(tmp.data, id.vars=char.vec)
out.plot.one <- ggplot(xVol, aes(x=psychosis_ar_4factor, y=value, group=marcat, color=factor(marcat))) +
#geom_point() +
facet_grid(.~variable) +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
coord_cartesian(ylim=c(-.7,.6)) +
geom_smooth(method='lm')
pdf('lobularMarcatVsPsych.pdf', width=30, height=10)
out.plot.one
dev.off()


## Now plot regional effects here
sig.lobe.index.vol <- c(2,3,5,6,7,9)
sig.lobe.index.tr <- c(1,2,5,7,8,9)
## Now go through all of the rois and identify those that belong to each lobe and
## plot the conitnous interaction between marcat and psychosis for each of these

