## Load library(s)
install_load('readstata13','ggplot2','lme4','lmerTest','reshape','psych','visreg')
source('../functions/functions.R')
## Load the data
data <- read.dta13('./cannabis_psychosis_cnb_foradon.dta')
## Now load the imaging data
img.data <- read.csv("~/Documents/dataRel/n1601_imagingDataDump_2018-09-20.csv")
img.data <- img.data[which(img.data$tpvalue==1),]
## Now load the RDS
data <-readRDS("../../01_dataPrep/scripts/mjPSCogImg.RDS")
img.data <-readRDS("../../01_dataPrep/scripts/mjPSCogImgIsol.RDS")


## Now create an isolated data set
vars.of.int <- c("f1_exec_comp_cog_accuracy","overall_accuracy","psychosis_ar_4factor","mprage_jlfHiLoLobe_vol_Basal_Ganglia","mprage_jlfHiLoLobe_vol_Limbic","mprage_jlfHiLoLobe_vol_Frontal_Orbital","mprage_jlfHiLoLobe_vol_Frontal_Dorsal","mprage_jlfHiLoLobe_vol_Temporal","mprage_jlfHiLoLobe_vol_Parietal","mprage_jlfHiLoLobe_vol_Occipital","mprage_jlfHiLoLobe_vol_Cerebellum","mprage_jlfHiLoLobe_vol_White_Matter","dti_jlfHiLoLobe_tr_Basal_Ganglia","dti_jlfHiLoLobe_tr_Limbic","dti_jlfHiLoLobe_tr_Frontal_Orbital","dti_jlfHiLoLobe_tr_Frontal_Dorsal","dti_jlfHiLoLobe_tr_Temporal","dti_jlfHiLoLobe_tr_Parietal","dti_jlfHiLoLobe_tr_Occipital","dti_jlfHiLoLobe_tr_Cerebellum","dti_jlfHiLoLobe_tr_White_Matter","envses","ageatcnb1","ageAtScan1","sex")

## Now create the data frame of interest
data.of.int <- data[,vars.of.int]
data.of.int <- apply(data.of.int[,1:21], 2, function(x) scale(residuals(lm(x~envses+ageatcnb1+sex, data=data.of.int, na.action=na.exclude)))[,])

## Now create the heat matrices
tmp.data <- data.of.int
to.rm <- c('.x', '_CorrTraits', 'mprage_jlfHiLoLobe_','mprage_jlfHiLoLobe_','mprage_jlfHiLoLobe_','pcasl_jlfHiLoLobe_','dti_jlfHiLoLobe_','rest_jlfHiLoLobe_','rest_jlfHiLoLobe_','dti_dtitk_jhutract_','_Lobe','_Efficiency')
for(i in to.rm){colnames(tmp.data) <- gsub(x=colnames(tmp.data), pattern=i, replacement='',perl=FALSE)}
colnames(tmp.data)[1:3] <- c('Executive Complex','Overall Accuracy','Psychosis Factor')
tmp.vals <- colnames(tmp.data)
cor.vals <- cor(tmp.data, use='pairwise')
diag(cor.vals) <- NA
sig.index <- corr.test(tmp.data, adjust='none')$p
cor.vals.sig <- cor.vals
cor.vals.sig[sig.index>.05] <- NA
cor.vals <- melt(cor.vals)
cor.vals.sig <- melt(cor.vals.sig)
cor.vals$X1 <- factor(cor.vals$X1, levels=tmp.vals)
cor.vals$X2 <- factor(cor.vals$X2, levels=tmp.vals)
cor.vals <- cor.vals[which(cor.vals$X2=="Executive Complex" | cor.vals$X2=="Overall Accuracy" | cor.vals$X2=="Psychosis Factor"),]
cor.vals.sig <- cor.vals.sig[which(cor.vals.sig$X2=="Executive Complex" | cor.vals.sig$X2=="Overall Accuracy" | cor.vals.sig$X2=="Psychosis Factor"),]
cor.vals.sig$X1 <- factor(cor.vals.sig$X1, levels=tmp.vals)
cor.vals.sig$X2 <- factor(cor.vals.sig$X2, levels=tmp.vals)
tmp <- ggplot(data=cor.vals, aes(x=X1, y=X2,fill=value)) +
geom_tile() +
theme(text=element_text(size=20), axis.text.x = element_text(angle = 45, hjust = 1, face="bold"),
axis.text.y = element_text(face="bold")) +
scale_fill_gradient2(low="blue", high="red", limits=c(-1, 1), midpoint = 0) +
coord_equal() + xlab("") + ylab("")
pdf("testHeatMap.pdf", height=6, width=16)
tmp
dev.off()
tmp <- ggplot(data=cor.vals.sig, aes(x=X1, y=X2,fill=value)) +
geom_tile() +
theme(text=element_text(size=20), axis.text.x = element_text(angle = 45, hjust = 1, face="bold"),
axis.text.y = element_text(face="bold")) +
scale_fill_gradient2(low="blue", high="red", limits=c(-1, 1), midpoint = 0) +
coord_equal() + xlab("") + ylab("")
pdf("testHeatMapSig.pdf", height=6, width=16)
tmp
dev.off()
