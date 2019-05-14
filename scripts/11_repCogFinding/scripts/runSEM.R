## Load library(s)
install_load('lavaan','semPlot','mediation')

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
vars.of.int <- c('psychosis_ar_4factor','marcatOU','marcatFU','PsyhOU','PsyhFU','f1_exec_comp_cog_accuracy','envses','ageatcnb1','ageAtScan1','sex','mprage_jlfHiLoLobe_vol_Temporal','dti_jlfHiLoLobe_tr_Occipital','marcat')

data.of.int <- data[,vars.of.int]
data.of.int$ageatcnb1 <- scale(data.of.int$ageatcnb1)[,]
data.of.int$ageAtScan1 <- scale(data.of.int$ageAtScan1)[,]
data.of.int$mprage_jlfHiLoLobe_vol_Temporal <- scale(data.of.int$mprage_jlfHiLoLobe_vol_Temporal)[,]
data.of.int$dti_jlfHiLoLobe_tr_Occipital <- scale(data.of.int$dti_jlfHiLoLobe_tr_Occipital)[,]
## Now create the cross product of the binary marcat variable
f1 <- paste("f1_exec_comp_cog_accuracy~psychosis_ar_4factor+marcatOU+marcatFU+PsyhOU+PsyhFU+envses")
f2 <- paste("dti_jlf_tr_MeanWholeBrainTR~psychosis_ar_4factor+marcatOU+marcatFU+PsyhOU+PsyhFU+envses+sex")
f3 <- paste("mprage_jlf_vol_TBV~psychosis_ar_4factor+marcatOU+marcatFU+PsyhOU+PsyhFU+envses+sex")
f4 <- paste("mprage_jlf_vol_TBV~ageAtScan1+sex")
f5 <- paste("dti_jlf_tr_MeanWholeBrainTR~ageAtScan1+sex")
f6 <- paste("psychosis_ar_4factor~ageatcnb1+sex")
f7 <- paste("f1_exec_comp_cog_accuracy~ageatcnb1+sex")
f8 <- paste("f1_exec_comp_cog_accuracy~mprage_jlf_vol_TBV")
f9 <- paste("f1_exec_comp_cog_accuracy~dti_jlf_tr_MeanWholeBrainTR")
f10 <- paste("ageAtScan1~~ageatcnb1")
f11 <- paste("mprage_jlf_vol_TBV~~dti_jlf_tr_MeanWholeBrainTR")
f12 <- paste("psychosis_ar_4factor~~PsyhOU")
f13 <- paste("psychosis_ar_4factor~~PsyhFU")
f14 <- paste("marcatFU~~PsyhFU")
f15 <- paste("marcatOU~~PsyhOU")
f16 <- paste("marcatOU~~ageAtScan1")
f17 <- paste("marcatFU~~ageAtScan1")
f18 <- paste("psychosis_ar_4factor~~MarcatFU")
f19 <- paste("psychosis_ar_4factor~~MarcatOU")


model <- paste(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,f19,sep="\n")
## Now fit model
fit <- sem(model,data=data.of.int)
semPaths(fit,"std", "hide")


## Now do one following the online tutorial
Mod.Med.Lavaan <- '
#Regressions
#These are the same regression equations from our previous example
#Except in this code we are naming the coefficients that are produced from the regression equations
#E.g., the regression coefficient for the effect of time on pubs is named "a1"
f1_exec_comp_cog_accuracy ~ a1*time.c + a2*alex.c + a3*time.c:alex.c
jobs ~ cdash1*time.c + cdash2*alex.c + cdash3*time.c:alex.c + b1*pubs

#Mean of centered alex (for use in simple slopes)
#This is making a coefficient labeled "alex.c.mean" which equals the intercept because of the "1"
#(Y~1) gives you the intercept, which is the mean for our alex.c variable
alex.c ~ alex.c.mean*1

#Variance of centered alex (for use in simple slopes)
#This is making a coefficient labeled "alex.c.var" which equals the variance because of the "~~"
#Two tildes separating the same variable gives you the variance
alex.c ~~ alex.c.var*alex.c

#Indirect effects conditional on moderator (a1 + a3*ModValue)*b1
indirect.SDbelow := (a1 + a3*(alex.c.mean-sqrt(alex.c.var)))*b1
indirect.SDabove := (a1 + a3*(alex.c.mean+sqrt(alex.c.var)))*b1

#Direct effects conditional on moderator (cdash1 + cdash3*ModValue)
#We have to do it this way because you cannot call the mean and sd functions in lavaan package
direct.SDbelow := cdash1 + cdash3*(alex.c.mean-sqrt(alex.c.var))
direct.SDabove := cdash1 + cdash3*(alex.c.mean+sqrt(alex.c.var))

#Total effects conditional on moderator
total.SDbelow := direct.SDbelow + indirect.SDbelow
total.SDabove := direct.SDabove + indirect.SDabove

#Proportion mediated conditional on moderator
#To match the output of "mediate" package
prop.mediated.SDbelow := indirect.SDbelow / total.SDbelow
prop.mediated.SDabove := indirect.SDabove / total.SDabove

#Index of moderated mediation
#An alternative way of testing if conditional indirect effects are significantly different from each other
index.mod.med := a3*b1
'


## Now try via the mediate package
## First age and sex regress everything
#data.of.int[,c(1,6,11,12)] <- apply(data.of.int[,c(1,6,11,12)], 2, function(x) residuals(lm(x ~ ageatcnb1+sex+envses, data=data.of.int, na.action=na.exclude)))
## Now create the two models
data <- readRDS("../../01_dataPrep/scripts/mjPSCogImg.RDS")
data <- data[!is.na(data$marcat),]
data <- data[complete.cases(data$psychosis_ar_4factor),]

## Now create our binary marcat user variable
data <- cbind(data, model.matrix(~ 0 + marcat, data)[,-2])
data$PsyhOU <- data$psychosis_ar_4factor*data$marcatOU
data$PsyhFU <- data$psychosis_ar_4factor*data$marcatFU
## Now create the binary marcat variable

## Now isolate our variables
vars.of.int <- c('psychosis_ar_4factor','marcatOU','marcatFU','PsyhOU','PsyhFU','f1_exec_comp_cog_accuracy','envses','ageatcnb1','ageAtScan1','sex','mprage_jlfHiLoLobe_vol_Temporal','dti_jlfHiLoLobe_tr_Occipital','marcat')

data.of.int <- data[,vars.of.int]
data.of.int$ageatcnb1 <- scale(data.of.int$ageatcnb1)[,]
data.of.int$ageAtScan1 <- scale(data.of.int$ageAtScan1)[,]
data.of.int$mprage_jlfHiLoLobe_vol_Temporal <- scale(data.of.int$mprage_jlfHiLoLobe_vol_Temporal)[,]
data.of.int$dti_jlfHiLoLobe_tr_Occipital <- scale(data.of.int$dti_jlfHiLoLobe_tr_Occipital)[,]
data.of.int <- data.of.int[complete.cases(data.of.int[,c(9,12)]),]
data.of.int <- data.of.int[-which(data.of.int$marcat=='FU'),]
Mod.Med.Model.1<-lm(f1_exec_comp_cog_accuracy ~ psychosis_ar_4factor*marcat+ageAtScan1+sex+envses, data = data.of.int)
Mod.Med.Model.2 <-lm(f1_exec_comp_cog_accuracy ~ psychosis_ar_4factor*marcat+ageAtScan1+sex+envses+dti_jlfHiLoLobe_tr_Occipital, data = data.of.int)
Mod.Med.LowAlex <- mediate(Mod.Med.Model.1, Mod.Med.Model.2,boot = TRUE,boot.ci.type = "bca", sims = 1000, mediator="dti_jlfHiLoLobe_tr_Occipital",treat='marcat')
summary(Mod.Med.LowAlex)
