#
t1 <- read.csv('alffData.csv')
t2 <- read.csv('cbfData.csv')
t3 <- read.csv('ctData.csv')
t4 <- read.csv('gmdData.csv')
t5 <- read.csv('jhuFALabel.csv')
t6 <- read.csv('jlfTRData.csv')
t7 <- read.csv('rehoData.csv')
t8 <- read.csv('volumeData.csv')

allData <- merge(t1, t2, all=T)
allData <- merge(allData, t3, all=T)
allData <- merge(allData, t4, all=T)
allData <- merge(allData, t5, all=T)
allData <- merge(allData, t6, all=T)
allData <- merge(allData, t7, all=T)
allData <- merge(allData, t8, all=T)

## Load library(s)
install_load('psych', 'ggplot2', 'pROC', 'ggrepel', 'caret', 'randomForest', 'MatchIt', 'glmnet')

# Now start loading the data down here
# This should be run from the dta directory
mjData <- read.csv("n9462_mj_ps_cnb_fortmm.csv")
mjData$dosage <- NA
mjData$dosage[which(mjData$marcat=='MJ Non-User')] <- 0
mjData$dosage[which(mjData$marcat=='MJ User' & mjData$mjpastyr=='')] <- 1
mjData$dosage[which(mjData$mjpastyr=="Less than once a month")] <- 2
mjData$dosage[which(mjData$mjpastyr=="About once a month")] <- 3
mjData$dosage[which(mjData$mjpastyr=="2-3 times a month")] <- 4
mjData$dosage[which(mjData$mjpastyr=="1-2 times a week")] <- 5
mjData$dosage[which(mjData$mjpastyr=="3-4 times a week")] <- 6
mjData$dosage[which(mjData$mjpastyr=="Everyday or nearly every day")] <- 7

# Now load imaging data
img.data <- read.csv('./imagingDataAll.csv')

# Now give us all the values
all.data <- merge(img.data, mjData)
all.data <- all.data[-which(all.data$dosage==1),]

# Now prepare a sex specific values
male.data <- all.data[which(all.data$sex==1),]
#male.data <- male.data[which(male.data$dosage == 0 | male.data$dosage > 5),]
female.data <- all.data[which(all.data$sex==2),]

# Now add a binary matrix for use or no use
male.data$usageBin <- 0
male.data$usageBin[male.data$dosage>1] <- 1

# Now run through every colume and give us a t value for every JLF region
tValsMaleOut <- NULL
colVals <- grep('_jlf_', names(male.data))
colVals <- append(colVals, grep('dtitk_jhulabel', names(male.data)))
for(z in colVals){
  tVal <-t.test(male.data[,z]~usageBin, data=male.data)
  outputRow <- c(colnames(male.data)[z], as.numeric(tVal$statistic), as.numeric(tVal$p.value))
  tValsMaleOut <- rbind(tValsMaleOut, outputRow)
}
tValsMaleOut <- tValsMaleOut[order(tValsMaleOut[,3]),]
tValsMaleSig <- tValsMaleOut[which(tValsMaleOut[,3]<.05),]
rownames(tValsMaleSig) <- NULL
tValsMaleSig <- as.data.frame(tValsMaleSig)

# Now do the same with ROC and find the corellation between roc and
rocValsMaleOut <- NULL
for(z in colVals){
  rocVal <- roc(usageBin~male.data[,z], data=male.data)
  outputRow <- c(colnames(male.data)[z], as.numeric(rocVal$auc))
  rocValsMaleOut <- rbind(rocValsMaleOut, outputRow)
}
rownames(rocValsMaleOut) <- NULL
rocValsMaleOut <- as.data.frame(rocValsMaleOut)

foobar <- merge(tValsMaleSig, rocValsMaleOut, by='V1')
foobar$V2.x <- as.numeric(as.character(foobar$V2.x))
foobar$V3 <- as.numeric(as.character(foobar$V3))
foobar$V2.y <- as.numeric(as.character(foobar$V2.y))
# Now seperate the positive and negative values
foobarPos <- foobar[which(foobar$V2.x>0),]
foobarNeg <- foobar[which(foobar$V2.x<0),]

# Now plot these values in a scatter plot
outPlot <- ggplot(foobarPos, aes(x=V2.y, y=V2.x)) +
  geom_point() +
  geom_smooth(method=lm) +
  geom_label_repel(aes(label=V1,size=2)) +#,box.padding=unit(1,"lines"),point.padding=unit(1,"lines")) +
  xlab("AUC") +
  ylab("t value") +
  #geom_hline(yintercept = 0 , linetype=3) +
  #geom_vline(xintercept = 0 , linetype=3) +
  theme(legend.position="none") +
  geom_text(aes(x=-Inf, y=Inf, hjust=0, vjust=1, label=cor(foobarPos$V2.y,foobarPos$V2.x)))

# Now prepare the t values here
tValsMaleOut <- tValsMaleOut[order(tValsMaleOut[,3]),]

# Now prepare a bar graph with these values
plotData <- tValsMaleOut[which(tValsMaleOut[,3]<.05),]

# Now plot the t values # STart with all positive t values


# Now I want to go through the AUC values and see how adding each additional variable adds
# to our cv performance prediciton - this will be done in a totally n fold cv
# where n is equal to our number of frequent users.

# First build a function which will build a roc model in a forward stepwise manner
# Inputs to the model include the X matrix, and the binary usage outcome
# We are going to have to take a lot of steps to ensure that we use the same rows at each step
buildStepROCModel <- function(x, y, nStep=50, varAdd=1){
    # Prepare the output of this function
    # this will be a matrix with the formula, AUC, significance from previous iteration
    outputData <- NULL
    outPredBase <- rep(NA, length(y))
    #lets first find the highest AUC value
    aucInit <- NULL
    valsToLoop <- seq(1, dim(x)[2])
    for(i in valsToLoop){
      tmp <- roc(y ~ x[,i])
      outVal <- tmp$auc
      aucInit <- append(aucInit, outVal)
    }
    # Now build a base model with these values
    colVal <- which(aucInit==max(aucInit))
    baseModel <- glm(y ~ x[,colVal], family=binomial())
    valsToLoop <- valsToLoop[-colVal]
    outPredBase[as.numeric(names(predict(baseModel)))] <- predict(baseModel)
    # Now prepare the inital row for the output
    initPred <- rep(NA, length(y))
    initPred[as.numeric(names(predict(baseModel)))] <- predict(baseModel)
    rocInit <- roc(y ~ initPred)
    aucInit <- pROC::auc(rocInit)
    modName <- paste(colnames(x)[colVal])
    initRow <- c(0, modName, aucInit, 0)
    outputData <- rbind(outputData, initRow)
    # Now begin the building process
    for(q in 1:length(valsToLoop)){
      # Initialize some variables
      aucVal <- NULL
      for(z in valsToLoop){
          # First build every model and get model performance
          print(z)
          colValNew <- append(colVal, z)
          tmpPredVals <- rep(NA, length(y))
          tmpMod <- glm(y ~ as.matrix(x[,colValNew]), family=binomial())
          tmpPredVals[as.numeric(names(predict(tmpMod)))] <- predict(tmpMod)
          aucVal <- append(aucVal, auc(roc(y ~ tmpPredVals)))
      }
      # Now grab the best performing model and output the model, performance metric, and
      # the difference in model performance compared to the previous iteration
      # then repeat
      oldPred <- rep(NA, length(y))
      newPred <- rep(NA, length(y))
      colVal <- append(colVal, which(floor(rank(aucVal))>=1 & floor(rank(aucVal))<=varAdd))
      #colVal <- append(colVal, which(aucVal==max(aucVal)))
      newModel <- glm(y ~ as.matrix(x[,colVal]), family=binomial())
      oldPred[as.numeric(names(predict(newModel)))] <- predict(baseModel)[as.numeric(names(predict(newModel)))]
      newPred[as.numeric(names(predict(newModel)))] <- predict(newModel)
      rocNew <- roc(y ~ newPred)
      rocOld <- roc(y ~ oldPred)
      modDiff <- roc.test(rocOld, rocNew, alternative='less', method='d')$p.value
      modName <- paste(colnames(x)[colVal], collapse='+')
      modPerf <- pROC::auc(rocNew)
      outputRow <- c(q, modName, modPerf, modDiff)
      # Now prepare for the next loop
      baseModel <- newModel
      valsToLoop <- valsToLoop[-which(floor(rank(aucVal))>=1 & floor(rank(aucVal))<=varAdd)]
      # Now if we don't see an imporvment in AUC break the loop
      if(pROC::auc(rocNew)<=pROC::auc(rocOld)){
        break
      }
      outputData <- rbind(outputData, outputRow)
    }
    return(outputData)
}


# Now lets see how well we can build our model in a cross validated fashion
# This will be done within modality just to explore things
foldsToLoop <- createFolds(male.data$usageBin, table(male.data$usageBin)[2])
foldsToLoop <- createFolds(male.data$usageBin, 5)
cvPredVals <- rep(NA, length(male.data$usageBin))
for(q in seq(1, 10)){
  index <- foldsToLoop[[q]]
  volMod <- buildStepROCModel(y=male.data$usageBin[-index], x=male.data[-index,grep('_jlf_vol_', names(male.data))], varAdd=6)
  outModel <- as.formula(paste('usageBin~', volMod[dim(volMod)[1],2]))
# Now build this model and
  tmpModel <- glm(outModel, data=male.data[-index,], family=binomial())
  cvPredVals[index] <- predict(tmpModel, newdata=male.data[index,], type='response')
  cv.glmnet()
}
plot(roc(male.data$usageBin ~ cvPredVals))
cvPredValsVol <- cvPredVals

cvPredVals <- rep(NA, length(male.data$usageBin))
for(q in seq(1, 10)){
    index <- foldsToLoop[[q]]
    volMod <- buildStepROCModel(y=male.data$usageBin[-index], x=male.data[-index,grep('_jlf_ct_', names(male.data))], varAdd=6)
    outModel <- as.formula(paste('usageBin~', volMod[dim(volMod)[1],2]))
    # Now build this model and
    tmpModel <- glm(outModel, data=male.data[-index,], family=binomial())
    #tmpModel <-  randomForest(formula=outModel,data=male.data[-index,],n.tree=300,mtry=3,na.action="na.omit")
    cvPredVals[index] <- predict(tmpModel, newdata=male.data[index,], type='response')
}
roc(male.data$usageBin ~ cvPredVals)

cvPredVals <- rep(NA, length(male.data$usageBin))
for(q in seq(1, 5)){
    index <- foldsToLoop[[q]]
    volMod <- buildStepROCModel(y=male.data$usageBin[-index], x=male.data[-index,grep('_jlf_', names(male.data))], varAdd=10)
    outModel <- as.formula(paste('usageBin~', volMod[dim(volMod)[1],2]))
    # Now build this model and
    tmpModel <- glm(outModel, data=male.data[-index,], family=binomial())
    tmpModel <-  randomForest(formula=outModel,data=male.data[-index,],n.tree=500,mtry=2,na.action="na.omit")
    cvPredVals[index] <- predict(tmpModel, newdata=male.data[index,], type='response')
}
