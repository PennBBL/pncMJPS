## Load library(s)
install_load('readstata13','ggplot2','lme4','lmerTest','reshape','psych','visreg')
source('../functions/functions.R')

## Load data
img.data <-readRDS("../../01_dataPrep/scripts/mjPSCogImgIsol.RDS")
char.vec <- c("bblid","envses","sex","race2","ageAtScan1","psBinary","mjbinary","psychosis_ar_4factor","marcat","averageManualRating")
global.val <- c("mprage_jlf_vol_TBV","dti_jlf_tr_MeanWholeBrainTR","pcasl_jlf_cbf_MeanGMCBF","mprage_jlf_gmd_MeanGMD","rest_jlf_alff_MeanALFF","rest_jlf_reho_MeanReho")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
## Run the mem
x <- img.data[,c(char.vec,global.val)]
x <- x[complete.cases(x$mprage_jlf_vol_TBV),]
## Now scale the values w/in modality
x[,11:dim(x)[2]] <- scale(x[,11:dim(x)[2]])
xCog <- melt(x, id.vars=char.vec)
xCog$marcat <- factor(xCog$marcat)
## Now run the model
mod.glo3 <- nlme::lme(value~ageAtScan1+sex+envses+race2+averageManualRating+marcat*psychosis_ar_4factor*variable,random=~1|bblid,data=xCog,na.action=na.exclude)

anova(lm(mprage_jlf_vol_TBV~ageAtScan1+sex+envses+averageManualRating+marcat*psychosis_ar_4factor, data=x, na.action=na.exclude))
anova(lm(dti_jlf_tr_MeanWholeBrainTR~ageAtScan1+sex+envses+averageManualRating+marcat*psychosis_ar_4factor, data=x, na.action=na.exclude))

## Now prepare table 1
n.val <- summarySE(data=x, measurevar='ageAtScan1', groupvars='marcat')[c('2','3','1'),'N']
mean.age <- round(summarySE(data=x, measurevar='ageAtScan1', groupvars='marcat')[c('2','3','1'),'ageAtScan1']/12,2)
mean.age.f.val <- anova(lm(ageAtScan1~marcat, data=x))
mean.age.sd <- round(summarySE(data=x, measurevar='ageAtScan1', groupvars='marcat')[c('2','3','1'),'sd']/12,2)
sex.perc.male <- round(summarySE(data=x[which(x$sex=='1'),], measurevar='ageAtScan1', groupvars='marcat')[c('2','3','1'),'N']/summarySE(data=x, measurevar='ageAtScan1', groupvars='marcat')[c('2','3','1'),'N'],2)
race.row.NU <- round((table(x$marcat, x$race2) / rowSums(table(x$marcat, x$race2)))[2,],2)
race.row.OU <- round((table(x$marcat, x$race2) / rowSums(table(x$marcat, x$race2)))[3,],2)
race.row.FU <- round((table(x$marcat, x$race2) / rowSums(table(x$marcat, x$race2)))[1,],2)
## Now load the wrat scores
wrat.scores <- read.csv('../../11_repCogFinding/scripts/n9498_cnb_wrat_scores_20161215.csv')
xWrat <- merge(x, wrat.scores)
wrat.vals <- round(summarySE(data=xWrat, measurevar='wrat4CrStd', groupvars='marcat',na.rm=T)[c('2','3','1'),'wrat4CrStd'],2)
wrat.sd <- round(summarySE(data=xWrat, measurevar='wrat4CrStd', groupvars='marcat',na.rm=T)[c('2','3','1'),'sd'],2)
## Now do the psychosis factor
psy.fac <- round(summarySE(data=x, measurevar='psychosis_ar_4factor', groupvars='marcat',na.rm=T)[c('2','3','1'),'psychosis_ar_4factor'],2)
psy.fac.sd <- round(summarySE(data=x, measurevar='psychosis_ar_4factor', groupvars='marcat',na.rm=T)[c('2','3','1'),'sd'],2)

## Now write the table
out.mat <- matrix(NA, ncol=3, nrow=8)
out.mat[1,] <- n.val
out.mat[2,] <- paste(mean.age, '(', mean.age.sd, ')', sep='')
out.mat[3,] <- sex.perc.male
out.mat[4:6,1] <- race.row.NU
out.mat[4:6,2] <- race.row.OU
out.mat[4:6,3] <- race.row.FU
out.mat[7,] <- paste(wrat.vals, '(', wrat.sd, ')', sep='')
out.mat[8,] <- paste(psy.fac, '(', psy.fac.sd, ')', sep='')
colnames(out.mat) <- c("Non-User", "Occasional User", "Frequent User")
rownames(out.mat) <- c("N","Age", "% Male", "Caucasian", "African-American", "Asian/Native American/Other", "Wrat Scores", "Psychosis Factor Score")
write.csv(out.mat, "VolumeTable.csv", quote=F)

## Now produce the global effects scatter plot
x$mprage_jlf_vol_TBV <- residuals(lm(mprage_jlf_vol_TBV~ageAtScan1+sex+envses+averageManualRating, data=x, na.action=na.exclude))

out.plot.one <- ggplot(x, aes(x=psychosis_ar_4factor, y=mprage_jlf_vol_TBV, group=marcat, color=factor(marcat))) +
theme(axis.text.x = element_text(angle = 0, hjust = 1),legend.position="none", text=element_text(size=28)) +
scale_colour_manual(values=cbbPalette) +
geom_smooth(method='lm')+#,alpha=0) +
#geom_point() +
xlab("PSYCHOSIS FACTOR SCORE") +
ylab("TBV") +
coord_cartesian(ylim=c(-1,1))
png('TBV.png')
print(out.plot.one)
dev.off()

## Do all of the volume figure production down here
tmp <- names(img.data)[c(grep("mprage_jlfHiLoLobe_vol", names(img.data)))]
xVol <- img.data[,c(char.vec, tmp)]
xVol[,10:18] <- scale(xVol[,10:18])
xVol <- melt(xVol, id.vars=char.vec)
xVol <- xVol[complete.cases(xVol),]
mod.vol3 <- nlme::lme(value~ageAtScan1+sex+envses+averageManualRating+marcat*psychosis_ar_4factor*variable,random=~1|bblid,data=xVol,na.action=na.exclude)
## Now plot the lobe effects
xVol <- img.data[,c(char.vec, tmp)]
## Now create one consolodated frontal lobe
xVol$mprage_jlfHiLoLobe_vol_Frontal <- xVol$mprage_jlfHiLoLobe_vol_Frontal_Orbital + xVol$mprage_jlfHiLoLobe_vol_Frontal_Dorsal
## ANd now remove the other frontal nonsense
xVol <- xVol[,-c(13,14)]
xVol[,11:18] <- scale(xVol[,11:18])
## Now find which lobes correct
plotLobes <- names(which(p.adjust(apply(xVol[,11:dim(xVol)[2]],2,function(x) anova(lm(x~ageAtScan1+sex+envses+race2+averageManualRating+marcat*psychosis_ar_4factor,data=xVol,na.action=na.exclude))['marcat:psychosis_ar_4factor','Pr(>F)']), method='fdr')<.05))
## Now remove covariates of non interest
xVol[,11:18] <- apply(xVol[,11:18], 2, function(x) residuals(lm(x~ageAtScan1+sex+envses+race2+averageManualRating, data=xVol, na.action=na.exclude)))
data.tmp <- melt(xVol, id.vars=char.vec)
## Now make the plots
yLab <- c("Limbic","Temporal","Parietal","Occipital","White Matter", "Frontal")
for(z in 1:length(plotLobes)){
    i <- plotLobes[z]
    xVol <- data.tmp[which(data.tmp$variable==i),]
    out.plot.one <- ggplot(xVol, aes(x=psychosis_ar_4factor, y=value, group=marcat, color=factor(marcat))) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1),legend.position="none", text=element_text(size=32)) +
    scale_colour_manual(values=cbbPalette) +
    geom_smooth(method='lm')+#,alpha=0) +
    #geom_point() +
    xlab("") +
    ylab(yLab[z]) +
    coord_cartesian(ylim=c(-1,1))
    out.plot.name <- paste(i, ".png", sep='')
    png(out.plot.name)
    print(out.plot.one)
    dev.off()
}


## Now here I will produce the regional effect sizes for each significant lobe
source('../../../../../jlfVisualizer/scripts/Rfunction/makeITKSnapColorTable.R')
## Now prepare an output table with all of the f stats
## We will be comparing the occasional users vs the non users here
## so remove the frequent users


## Now get all of the F stats
out.mat.no.rm <- NULL
out.mat.rm.FU <- NULL
out.mat.rm.OU <- NULL
out.mat.rm.NU <- NULL
for(roi in names(img.data)[358:496]){
    ## Run the model
    tmp.form <- as.formula(paste(roi, "~ageAtScan1+sex+envses+race2+averageManualRating+marcat*psychosis_ar_4factor", sep=''))
    mod.out.all <- anova(lm(tmp.form, data=img.data, na.action=na.exclude))['marcat:psychosis_ar_4factor',c('F value','Pr(>F)')]
    mod.out <- anova(lm(tmp.form, data=img.data[-which(img.data$marcat=='FU'),], na.action=na.exclude))['marcat:psychosis_ar_4factor',c('F value','Pr(>F)')]
    mod.out.two <- anova(lm(tmp.form, data=img.data[-which(img.data$marcat=='OU'),], na.action=na.exclude))['marcat:psychosis_ar_4factor',c('F value','Pr(>F)')]
    mod.out.three <- anova(lm(tmp.form, data=img.data[-which(img.data$marcat=='NU'),], na.action=na.exclude))['marcat:psychosis_ar_4factor',c('F value','Pr(>F)')]
    ## Now grab the lobe value
    lobeVal <- findLobe(roi)
    
    ## Now prepare the data
    out.row.all <- cbind(roi,lobeVal,mod.out.all)
    out.row.one <- cbind(roi,lobeVal,mod.out)
    out.row.two <- cbind(roi, lobeVal, mod.out.two)
    out.row.three <- cbind(roi, lobeVal, mod.out.three)
    
    ## Now make the output row
    out.mat.no.rm <- rbind(out.mat.no.rm, out.row.all)
    out.mat.rm.FU <- rbind(out.mat.rm.FU, out.row.one)
    out.mat.rm.OU <- rbind(out.mat.rm.OU, out.row.two)
    out.mat.rm.NU <- rbind(out.mat.rm.NU, out.row.three)
}

## Now write the color map with all ROI's
#out.mat.rm.FU <- cbind(out.mat.rm.FU, rep(NA, dim(out.mat.rm.FU)[1]))
#out.mat.rm.OU <- cbind(out.mat.rm.OU, rep(NA, dim(out.mat.rm.FU)[1]))
#out.mat.rm.NU <- cbind(out.mat.rm.NU, rep(NA, dim(out.mat.rm.FU)[1]))


writeColorTableandKey(inputColumn=3, inputData=out.mat.rm.FU, maxTmp=c(0,17), minTmp=c(-.8, 0), outName="volOUvNUf")
writeColorTableandKey(inputColumn=3, inputData=out.mat.rm.OU, maxTmp=c(0,17), minTmp=c(-.8, 0), outName="volFUvNUf")
writeColorTableandKey(inputColumn=3, inputData=out.mat.rm.NU, maxTmp=c(0,17), minTmp=c(-.8, 0), outName="volFUvOUf")

## Now copy these over to chead
system("scp vol*f*csv arosen@chead:/data/jux/BBL/projects/pncMJPS/data/colorAndKey/")
system("scp vol*f*txt arosen@chead:/data/jux/BBL/projects/pncMJPS/data/colorAndKey/")


## Now combine all of these into one table
out.mat.no.rm <- as.data.frame(out.mat.no.rm)
out.mat.rm.FU <- as.data.frame(out.mat.rm.FU)
out.mat.rm.OU <- as.data.frame(out.mat.rm.OU)
out.mat.rm.NU <- as.data.frame(out.mat.rm.NU)
out.sup.table <- merge(out.mat.no.rm, out.mat.rm.FU, by=c('roi','lobeVal'),suffixes = c(".NoRm",".NoFU"))
out.sup.table <- merge(out.sup.table, out.mat.rm.OU, by=c('roi','lobeVal'))
colnames(out.sup.table)[7:8] <- paste(names(out.sup.table)[7:8], '.NoOU', sep='')
out.sup.table <- merge(out.sup.table, out.mat.rm.NU, by=c('roi','lobeVal'))
colnames(out.sup.table)[9:10] <- paste(names(out.sup.table)[9:10], '.NoNU', sep='')

## Now save the table with any nom sig values
index <- unique(unlist(apply(out.sup.table[,c(4,6,8,10)], 2, function(x) which(x<.005))))
out.sup.table <- out.sup.table[index,]
write.csv(out.sup.table, "allNomSigSupValsVol.csv", quote=F, row.names=F)
