## Load library(s)
install_load('readstata13','ggplot2','lme4','lmerTest','reshape','psych','visreg')
source('../functions/functions.R')

## Load data
img.data <-readRDS("../../01_dataPrep/scripts/mjPSCogImgIsol.RDS")
char.vec <- c("bblid","envses","sex","race2","ageAtScan1","psBinary","mjbinary","psychosis_ar_4factor","marcat","dti64Tsnr")
global.val <- c("mprage_jlf_vol_TBV","dti_jlf_tr_MeanWholeBrainTR","pcasl_jlf_cbf_MeanGMCBF","mprage_jlf_gmd_MeanGMD","rest_jlf_alff_MeanALFF","rest_jlf_reho_MeanReho")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
## Run the mem
x <- img.data[,c(char.vec,global.val)]
## Now scale the values w/in modality
x[,11:dim(x)[2]] <- scale(x[,11:dim(x)[2]])
xCog <- melt(x, id.vars=char.vec)
xCog$marcat <- factor(xCog$marcat)
## Now run the model
mod.glo3 <- nlme::lme(value~ageAtScan1+sex+envses+race2+dti64Tsnr+marcat*psychosis_ar_4factor*variable,random=~1|bblid,data=xCog,na.action=na.exclude)

## Now produce the global effects scatter plot
x$dti_jlf_tr_MeanWholeBrainTR <- residuals(lm(dti_jlf_tr_MeanWholeBrainTR~ageAtScan1+sex+envses+dti64Tsnr, data=x, na.action=na.exclude))


out.plot.one <- ggplot(x, aes(x=psychosis_ar_4factor, y=dti_jlf_tr_MeanWholeBrainTR, group=marcat, color=factor(marcat))) +
theme(axis.text.x = element_text(angle = 0, hjust = 1),legend.position="none", text=element_text(size=32)) +
scale_colour_manual(values=cbbPalette) +
geom_smooth(method='lm') +
xlab("") +
ylab("Whole Brain MD") +
coord_cartesian(ylim=c(-1,1))
png('TBMD.png')
print(out.plot.one)
dev.off()

## Do all of the MD figure production down here
tmp <- names(img.data)[c(grep("dti_jlfHiLoLobe_tr", names(img.data)))]
xVol <- img.data[,c(char.vec, tmp)]
xVol[,11:19] <- scale(xVol[,11:19])
xVol <- melt(xVol, id.vars=char.vec)
xVol <- xVol[complete.cases(xVol),]
mod.vol3 <- nlme::lme(value~ageAtScan1+sex+envses+dti64Tsnr+marcat*psychosis_ar_4factor*variable,random=~1|bblid,data=xVol,na.action=na.exclude)
## Now plot the lobe effects
xVol <- img.data[,c(char.vec, tmp)]
xVol[,11:19] <- scale(xVol[,11:19])
plotLobes <- names(which(apply(xVol[,11:dim(xVol)[2]],2,function(x) anova(lm(x~ageAtScan1+sex+envses+dti64Tsnr+race2+marcat*psychosis_ar_4factor,data=xVol,na.action=na.exclude))['marcat:psychosis_ar_4factor','Pr(>F)'])<.05))
plotLobes <- names(xVol)[c(11,12,15,17,18,19)]
yLab <- c("Basal Ganglia","Limbic","Temporal","Occipital","Cerebellum", "White Matter")
## Remove confounds
xVol[,11:19] <- apply(xVol[,11:19], 2, function(x) residuals(lm(x~ageAtScan1+sex+envses+race2+dti64Tsnr, data=xVol, na.action=na.exclude)))
data.tmp <- melt(xVol, id.vars=char.vec)
for(z in 1:length(plotLobes)){
    i <- plotLobes[z]
    xVol <- data.tmp[which(data.tmp$variable==i),]
    out.plot.one <- ggplot(xVol, aes(x=psychosis_ar_4factor, y=value, group=marcat, color=factor(marcat))) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1),legend.position="none", text=element_text(size=32)) +
    scale_colour_manual(values=cbbPalette) +
    geom_smooth(method='lm') +
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
img.data.freeze <- img.data
img.data <- img.data[-which(img.data$marcat=='OU'),]

## Now get all of the F stats
out.mat <- NULL
for(roi in names(img.data)[965:1095]){
    ## Run the model
    tmp.form <- as.formula(paste(roi, "~ageAtScan1+sex+envses+race2+averageManualRating+marcat*psychosis_ar_4factor", sep=''))
    mod.out <- anova(lm(tmp.form, data=img.data, na.action=na.exclude))['marcat:psychosis_ar_4factor',c('F value','Pr(>F)')]
    ## Now grab the lobe value
    lobeVal <- findLobe(roi)
    out.row <- cbind(roi,lobeVal,mod.out)
    out.mat <- rbind(out.mat, out.row)
}
out.mat <- cbind(out.mat, rep(NA, dim(out.mat)[1]))
## Now go through the sig lobes and find the ROIs we want to keep
for(lobeVal in c(1,3:4,6:9)){
    ## Grab the lobes we want to work with
    vals.of.int <- out.mat[which(out.mat[,2]==lobeVal),]
    index <- which(out.mat[,2]==lobeVal)
    out.mat[index,5] <- p.adjust(vals.of.int[,4], method='fdr')
    
}
## Now write the color map with all ROI's
to.write <- out.mat[which(out.mat[,2]==1|out.mat[,2]==2|out.mat[,2]==3|out.mat[,2]==5|out.mat[,2]==6|out.mat[,2]==7|out.mat[,2]==8|out.mat[,2]==9),]
writeColorTableandKey(inputColumn=3, inputData=to.write, maxTmp=c(0,17), minTmp=c(-.8, 0), outName="trf")

## Now do the FDR sig vals
to.write.sig <- to.write[which(to.write[,5]<.05),]
writeColorTableandKey(inputColumn=3, inputData=to.write.sig, maxTmp=c(0,17), minTmp=c(-.8, 0), outName="trfsig")

## Now copy these over to chead
system("scp trf* arosen@chead:/data/jux/BBL/projects/pncMJPS/data/colorAndKey/")
