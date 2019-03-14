## Load library(s)
install_load('readstata13','ggplot2','lme4','lmerTest','reshape','psych','visreg')

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
global.val <- c("mprage_jlf_vol_TBV","mprage_jlf_gmd_MeanGMD","dti_jlf_tr_MeanWholeBrainTR")
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
mod.vol3 <- nlme::lme(value~scanageMonths+sex+envses+marcat*psychosis_ar_4factor*variable,random=~1|bblid,data=xVol,na.action=na.exclude)

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

## Now do FA
tmp <- names(img.data)[c(grep("dti_dtitk_jhutract_fa", names(img.data)))]
xFA <- img.data[,c(char.vec, tmp)]
#xFA <- averageLeftAndRight1(xFA)
xFA[,10:27] <- scale(xFA[,10:27])
xFA <- melt(xFA, id.vars=char.vec)
mod.fa3 <- nlme::lme(value~scanageMonths+sex+envses+marcat*psychosis_ar_4factor*variable,random=~1|bblid,data=xFA,na.action=na.exclude)
