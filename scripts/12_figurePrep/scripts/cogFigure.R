## Load library(s)
install_load('readstata13','ggplot2','lme4','lmerTest','reshape','psych','visreg')
source('../functions/functions.R')
## Load the data
data <-readRDS("../../01_dataPrep/scripts/mjPSCogImg.RDS")

## Now declare some static variables
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
char.vec <- c("bblid","envses","sex","race2","ageatcnb1","psBinary","mjbinary","psychosis_ar_4factor","marcat","mj_firstuse")
cog.names <- c("f1_exec_comp_cog_accuracy","f2_social_cog_accuracy","f3_memory_accuracy","overall_accuracy")

## Now create the data set
x <- data[,c(char.vec,cog.names)]
xCog <- melt(x, id.vars=char.vec)
xCog$marcat <- factor(xCog$marcat)

## Now run the model
mod.cog3 <- lmerTest::lmer(value ~ ageatcnb1+sex+envses+factor(race2)+(marcat+psychosis_ar_4factor+variable)^3+(1|bblid),data=xCog,na.action=na.exclude)

## Now run individual models
p.adjust(
c(anova(lm(f1_exec_comp_cog_accuracy~ageatcnb1+sex+envses+factor(race2)+marcat*psychosis_ar_4factor, data=x))['marcat:psychosis_ar_4factor','Pr(>F)'],
anova(lm(f2_social_cog_accuracy~ageatcnb1+sex+envses+factor(race2)+marcat*psychosis_ar_4factor, data=x))['marcat:psychosis_ar_4factor','Pr(>F)'],
anova(lm(f3_memory_accuracy~ageatcnb1+sex+envses+factor(race2)+marcat*psychosis_ar_4factor, data=x))['marcat:psychosis_ar_4factor','Pr(>F)'],
anova(lm(overall_accuracy~ageatcnb1+sex+envses+factor(race2)+marcat*psychosis_ar_4factor, data=x))['marcat:psychosis_ar_4factor','Pr(>F)']), method='fdr')

## Now explore a GAF interaction
gaf.data <- read.csv('../../11_repCogFinding/scripts/ERS_7-20-2018.csv')
x.gaf <- merge(x, gaf.data, by='bblid')
char.vec2 <- c("bblid","envses","sex.x","race2","ageatcnb1","psBinary","mjbinary","psychosis_ar_4factor","marcat")
cog.names <- c("gaf001")
## Now look for an interaction with gaf
x2 <- x.gaf[,c(char.vec2,cog.names)]
lm.check <- lm(gaf001 ~ ageatcnb1+sex.x+envses+factor(race2)+marcat*psychosis_ar_4factor, data=x2)

## Now prepare table 1
n.val <- summarySE(data=x, measurevar='ageatcnb1', groupvars='marcat')[c('2','3','1'),'N']
mean.age <- round(summarySE(data=x, measurevar='ageatcnb1', groupvars='marcat')[c('2','3','1'),'ageatcnb1']/12,2)
mean.age.sd <- round(summarySE(data=x, measurevar='ageatcnb1', groupvars='marcat')[c('2','3','1'),'sd']/12,2)
sex.perc.male <- round(summarySE(data=x[which(x$sex=='1'),], measurevar='ageatcnb1', groupvars='marcat')[c('2','3','1'),'N']/summarySE(data=x, measurevar='ageatcnb1', groupvars='marcat')[c('2','3','1'),'N'],2)
sex.perc.male.n <- summarySE(data=x[which(x$sex=='1'),], measurevar='ageatcnb1', groupvars='marcat')[c('2','3','1'),'N']
sex.perc.female.n <- summarySE(data=x[which(x$sex=='2'),], measurevar='ageatcnb1', groupvars='marcat')[c('2','3','1'),'N']

race.row.NU <- round((table(x$marcat, x$race2) / rowSums(table(x$marcat, x$race2)))[2,],2)
race.row.NU.n <-  table(x$marcat, x$race2)[2,]
race.row.OU <- round((table(x$marcat, x$race2) / rowSums(table(x$marcat, x$race2)))[3,],2)
race.row.OU.n <- table(x$marcat, x$race2)[3,]
race.row.FU <- round((table(x$marcat, x$race2) / rowSums(table(x$marcat, x$race2)))[1,],2)
race.row.FU.n <- table(x$marcat, x$race2)[1,]
## Now load the wrat scores
wrat.scores <- read.csv('../../11_repCogFinding/scripts/n9498_cnb_wrat_scores_20161215.csv')
xWrat <- merge(x, wrat.scores)
wrat.vals <- round(summarySE(data=xWrat, measurevar='wrat4CrStd', groupvars='marcat',na.rm=T)[c('2','3','1'),'wrat4CrStd'],2)
wrat.sd <- round(summarySE(data=xWrat, measurevar='wrat4CrStd', groupvars='marcat',na.rm=T)[c('2','3','1'),'sd'],2)
## Now do age at first use
age.first <- round(summarySE(data=x, measurevar='mj_firstuse', groupvars='marcat', na.rm=T)[c('2','3','1'),'mj_firstuse'],2)
age.first.sd <- round(summarySE(data=x, measurevar='mj_firstuse', groupvars='marcat', na.rm=T)[c(2,3,1),'sd'],2)
## Now do the psychosis factor
psy.fac <- round(summarySE(data=x, measurevar='psychosis_ar_4factor', groupvars='marcat',na.rm=T)[c('2','3','1'),'psychosis_ar_4factor'],2)
psy.fac.sd <- round(summarySE(data=x, measurevar='psychosis_ar_4factor', groupvars='marcat',na.rm=T)[c('2','3','1'),'sd'],2)

## Now obtain p values for all of these variables
age.p.val <- round(anova(lm(ageatcnb1 ~ marcat, data=x))['marcat', 'Pr(>F)'],8)
sex.p.val <- round(as.numeric(chisq.test(table(x$sex, x$marcat))['p.value']),8)
race.p.val <- round(as.numeric(chisq.test(table(x$race2, x$marcat))['p.value']),5)
wrat.p.val <- round(anova(lm(wrat4CrStd ~ marcat, data=xWrat))['marcat', 'Pr(>F)'],5)
age.first.p.val <- round(anova(lm(mj_firstuse~marcat, data=x))['marcat','Pr(>F)'],8)
psy.p.val <- round(anova(lm(psychosis_ar_4factor ~ marcat, data=xWrat))['marcat', 'Pr(>F)'],5)

## Now write the table
out.mat <- matrix(NA, ncol=4, nrow=9)
out.mat[1,] <- c(n.val, 'NA')
out.mat[2,] <- c(paste(mean.age, '(', mean.age.sd, ')', sep=''), age.p.val)
out.mat[3,1:3] <- paste(sex.perc.male.n, "/", sex.perc.female.n, " (", sex.perc.male, ")", sep='')
out.mat[3,4] <- as.numeric(sex.p.val)
out.mat[4:6,1] <- paste(race.row.NU.n, " (",race.row.NU,")", sep='')
out.mat[4:6,2] <- paste(race.row.OU.n, " (",race.row.OU,")", sep='')
out.mat[4:6,3] <- paste(race.row.FU.n, " (",race.row.FU,")", sep='')
out.mat[4:6,4] <- as.numeric(race.p.val)
out.mat[7,] <- c(paste(wrat.vals, '(', wrat.sd, ')', sep=''),wrat.p.val)
out.mat[8,] <- c(paste(age.first, '(', age.first.sd, ')', sep=''),age.first.p.val)
out.mat[9,] <- c(paste(psy.fac, '(', psy.fac.sd, ')', sep=''),psy.p.val)
colnames(out.mat) <- c("Non-User", "Occasional User", "Frequent User", "p")
rownames(out.mat) <- c("N","Age", "% Male", "Caucasian", "African-American", "Asian/Native American/Other", "Wrat Scores", "Age first use","Psychosis Factor Score")
write.csv(out.mat, "CognitiveTable.csv", quote=F)


## Now regress out covariates of non interest
out.data.one <- NULL
data.tmp <- data
for(i in c("f1_exec_comp_cog_accuracy","f2_social_cog_accuracy","f3_memory_accuracy","overall_accuracy")){
    data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + ageatcnb1 + envses + race2")), data=data.tmp, na.action=na.exclude))
    data.tmp[,i] <- scale(data.tmp[,i])
}


## Now create our plots
char.vec <- c("bblid","envses","sex","race2","ageatcnb1","psBinary","mjbinary","psychosis_ar_4factor","marcat")
cog.names <- c("f1_exec_comp_cog_accuracy","f2_social_cog_accuracy","f3_memory_accuracy","overall_accuracy")
tmp.data <- data.tmp[,c(char.vec,cog.names)]
xVol <- melt(tmp.data, id.vars=char.vec)
xVol <- xVol[-which(is.na(xVol$marcat)),]
out.plot.one <- ggplot(xVol[which(xVol$variable=='f1_exec_comp_cog_accuracy'),], aes(x=psychosis_ar_4factor, y=value, group=marcat, color=factor(marcat))) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1),legend.position="none", text=element_text(size=32)) +
    scale_colour_manual(values=cbbPalette) +
    geom_smooth(method='lm') +
    xlab("") +
    ylab("") +
    coord_cartesian(ylim=c(-.6,.6))

out.plot.two <- ggplot(xVol[which(xVol$variable=='f2_social_cog_accuracy'),], aes(x=psychosis_ar_4factor, y=value, group=marcat, color=factor(marcat))) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1),legend.position="none", text=element_text(size=32)) +
    scale_colour_manual(values=cbbPalette) +
    geom_smooth(method='lm') +
    xlab("") +
    ylab("") +
    coord_cartesian(ylim=c(-.6,.6))

out.plot.three <- ggplot(xVol[which(xVol$variable=='f3_memory_accuracy'),], aes(x=psychosis_ar_4factor, y=value, group=marcat, color=factor(marcat))) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1),legend.position="none", text=element_text(size=32)) +
    scale_colour_manual(values=cbbPalette) +
    geom_smooth(method='lm') +
    xlab("") +
    ylab("") +
    coord_cartesian(ylim=c(-.6,.6))

out.plot.four <- ggplot(xVol[which(xVol$variable=='overall_accuracy'),], aes(x=psychosis_ar_4factor, y=value, group=marcat, color=factor(marcat))) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1),legend.position="none", text=element_text(size=32)) +
    scale_colour_manual(values=cbbPalette) +
    geom_smooth(method='lm') +
    xlab("") +
    ylab("") +
    coord_cartesian(ylim=c(-.6,.6))


## Now write these figures
png('f1Exec.png')
out.plot.one
dev.off()
png('f2Social.png')
out.plot.two
dev.off()
png('f3Memory.png')
out.plot.three
dev.off()
png('f4Overall.png')
out.plot.four
dev.off()

