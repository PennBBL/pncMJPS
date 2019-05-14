## Load library(s)
install_load('readstata13','ggplot2','lme4','lmerTest','reshape','psych','visreg')
source('../functions/functions.R')
## Load the data
data <-readRDS("../../01_dataPrep/scripts/mjPSCogImg.RDS")

## Now declare some static variables
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
char.vec <- c("bblid","envses","sex","race2","ageatcnb1","psBinary","mjbinary","psychosis_ar_4factor","marcat")
cog.names <- c("f1_exec_comp_cog_accuracy","f2_social_cog_accuracy","f3_memory_accuracy","overall_accuracy")

## Now create the data set
x <- data[,c(char.vec,cog.names)]
xCog <- melt(x, id.vars=char.vec)
xCog$marcat <- factor(xCog$marcat)

## Now run the model
mod.cog3 <- lmerTest::lmer(value ~ ageatcnb1+sex+envses+factor(race2)+(marcat+psychosis_ar_4factor+variable)^3+(1|bblid),data=xCog,na.action=na.exclude)


## Now regress out covariates of non interest
out.data.one <- NULL
data.tmp <- data
for(i in c("f1_exec_comp_cog_accuracy","f2_social_cog_accuracy","f3_memory_accuracy","overall_accuracy")){
    data.tmp[,i] <- residuals(lm(as.formula(paste(i," ~ sex + ageatcnb1 + envses + race2")), data=data.tmp, na.action=na.exclude))
    data.tmp[,i] <- scale(data.tmp[,i])
}


## Now create our plots
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

