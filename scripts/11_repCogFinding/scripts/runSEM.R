## Load library(s)
install_load('blavaan','semPlot','mediation')

## Load data
## Now load the RDS file as our data
data <- readRDS("../../01_dataPrep/scripts/mjPSCogImg.RDS")
data <- data[!is.na(data$marcat),]
data <- data[complete.cases(data$psychosis_ar_4factor),]

## Now create our binary marcat user variable
data <- cbind(data, model.matrix(~ 0 + marcat, data)[,-2])
data$PsyhOU <- data$psychosis_ar_4factor*data$marcatOU
data$PsyhFU <- data$psychosis_ar_4factor*data$marcatFU
## Now create the binary marcat variable

## Now isolate our variables
vars.of.int <- c('psychosis_ar_4factor','marcatOU','marcatFU','PsyhOU','PsyhFU','f1_exec_comp_cog_accuracy','envses','ageatcnb1','ageAtScan1','sex','mprage_jlf_vol_TBV','dti_jlf_tr_MeanWholeBrainTR','marcat')

data.of.int <- data[,vars.of.int]
data.of.int$ageatcnb1 <- scale(data.of.int$ageatcnb1)[,]
data.of.int$ageAtScan1 <- scale(data.of.int$ageAtScan1)[,]
data.of.int$mprage_jlf_vol_TBV <- scale(data.of.int$mprage_jlf_vol_TBV)[,]
data.of.int$dti_jlf_tr_MeanWholeBrainTR <- scale(data.of.int$dti_jlf_tr_MeanWholeBrainTR)[,]

## Now look at lrelationship of MEANTR after TBV regression
data.to.exp <- data.of.int
plot(data.to.exp$dti_jlf_tr_MeanWholeBrainTR, data.to.exp$mprage_jlf_vol_TBV)
data.to.exp[,c(1,6,12)] <- apply(data.to.exp[,c(1,6,12)], 2, function(x) residuals(lm(x~ageatcnb1+sex+mprage_jlf_vol_TBV, data=data.to.exp, na.action=na.exclude)))
plot(data.to.exp$dti_jlf_tr_MeanWholeBrainTR, data.to.exp$f1_exec_comp_cog_accuracy)

## Now regress age and sex out of everything
data.of.int[,c(1,6,11,12)] <- apply(data.of.int[,c(1,6,11,12)], 2, function(x) residuals(lm(x~ageatcnb1+sex+envses, data=data.of.int, na.action=na.exclude)))

## First do the non users volume
set.seed(1234)
Data <- data.frame(X = data.of.int$psychosis_ar_4factor[which(data.of.int$marcat=='NU')], Y = data.of.int$f1_exec_comp_cog_accuracy[which(data.of.int$marcat=='NU')], M = data.of.int$mprage_jlf_vol_TBV[which(data.of.int$marcat=='NU')])
model <- ' # direct effect
Y ~ c*X
# mediator
M ~ a*X
Y ~ b*M
# indirect effect (a*b)
ab := a*b
# total effect
total := c + (a*b)
'
fitNU <- bsem(model, data = Data, jagcontrol=list(method="rjparallel"))
out.vals <- summary(fitNU)
## Now find the indirect path p value
se.val <- as.numeric(out.vals[10,4]) - as.numeric(out.vals[10,4])/ (2*1.96)
z.val <- as.numeric(out.vals[10,2])/se.val
p.val.vol.NU <- pnorm(-abs(z.val))

## Now do the occasional users volume
Data <- data.frame(X = data.of.int$psychosis_ar_4factor[which(data.of.int$marcat=='OU')], Y = data.of.int$f1_exec_comp_cog_accuracy[which(data.of.int$marcat=='OU')], M = data.of.int$mprage_jlf_vol_TBV[which(data.of.int$marcat=='OU')])
fitOU <- bsem(model, data = Data, jagcontrol=list(method="rjparallel"))
summary(fitOU)
out.vals <- summary(fitOU)
## Now find the indirect path p value
se.val <- as.numeric(out.vals[10,4]) - as.numeric(out.vals[10,4])/ (2*1.96)
z.val <- as.numeric(out.vals[10,2])/se.val
p.val.vol.OU <- pnorm(-abs(z.val))

## Now do the frequent users volume
Data <- data.frame(X = data.of.int$psychosis_ar_4factor[which(data.of.int$marcat=='FU')], Y = data.of.int$f1_exec_comp_cog_accuracy[which(data.of.int$marcat=='FU')], M = data.of.int$mprage_jlf_vol_TBV[which(data.of.int$marcat=='FU')])
fitFU <- bsem(model, data = Data, jagcontrol=list(method="rjparallel"))
summary(fitFU)
out.vals <- summary(fitFU)
## Now find the indirect path p value
se.val <- as.numeric(out.vals[10,4]) - as.numeric(out.vals[10,4])/ (2*1.96)
z.val <- as.numeric(out.vals[10,2])/se.val
p.val.vol.FU <- pnorm(-abs(z.val))

## Save the p values
out.p.val.vol <- c(p.val.vol.NU,p.val.vol.OU,p.val.vol.FU)

## Now do the DTI values
## First do the non users volume
Data <- data.frame(X = data.of.int$psychosis_ar_4factor[which(data.of.int$marcat=='NU')], Y = data.of.int$f1_exec_comp_cog_accuracy[which(data.of.int$marcat=='NU')], M = data.of.int$dti_jlf_tr_MeanWholeBrainTR[which(data.of.int$marcat=='NU')])
fitNU <- sem(model, data = Data, se='bootstrap')
summary(fitNU)
out.vals <- summary(fitNU)
## Now find the indirect path p value
se.val <- as.numeric(out.vals[10,4]) - as.numeric(out.vals[10,4])/ (2*1.96)
z.val <- as.numeric(out.vals[10,2])/se.val
p.val.dti.NU <- pnorm(-abs(z.val))

## Now do the occasional users volume
Data <- data.frame(X = data.of.int$psychosis_ar_4factor[which(data.of.int$marcat=='OU')], Y = data.of.int$f1_exec_comp_cog_accuracy[which(data.of.int$marcat=='OU')], M = data.of.int$dti_jlf_tr_MeanWholeBrainTR[which(data.of.int$marcat=='OU')])
fitOU <- sem(model, data = Data, se='bootstrap')
summary(fitOU)
out.vals <- summary(fitOU)
## Now find the indirect path p value
se.val <- as.numeric(out.vals[10,4]) - as.numeric(out.vals[10,4])/ (2*1.96)
z.val <- as.numeric(out.vals[10,2])/se.val
p.val.dti.OU <- pnorm(-abs(z.val))

## Now do the frequent users volume
Data <- data.frame(X = data.of.int$psychosis_ar_4factor[which(data.of.int$marcat=='FU')], Y = data.of.int$f1_exec_comp_cog_accuracy[which(data.of.int$marcat=='FU')], M = data.of.int$dti_jlf_tr_MeanWholeBrainTR[which(data.of.int$marcat=='FU')])
fitFU <- sem(model, data = Data, se='bootstrap')
summary(fitFU)
out.vals <- summary(fitFU)
## Now find the indirect path p value
se.val <- as.numeric(out.vals[10,4]) - as.numeric(out.vals[10,4])/ (2*1.96)
z.val <- as.numeric(out.vals[10,2])/se.val
p.val.dti.FU <- pnorm(-abs(z.val))

out.p.val.dti <- c(p.val.vol.NU,p.val.vol.OU,p.val.vol.FU)

out.p.val <- rbind(out.p.val.vol, out.p.val.dti)
write.csv(out.p.val, "outIndirectPVal.csv", quote=F)
