# Create a function which will run lasso variable selection mutliple times
runLassoforHiLo <- function(x, y, trainingIterations = 100, nCor=3, nofFolds=10, alphaSequence=seq(0,1,by=0.05)){
    # First thing we have to do is prepare our output
    outputMatrix <- matrix(0, ncol(x), trainingIterations+1)
    outputMatrix[,1] <- colnames(x)
    
    # Now we need to make sure our inputs are the correct formats
    inputData <- as.matrix(x)
    
    # Now set up our parallel environment
    cl <- makeCluster(nCor)
    registerDoParallel(cl)
    
    # Now lets run through each iteration and perform the following steps
    # 1.) Tune to find the optimal alpha for each fold - using tuneAlpha function
    # 2.) create a model using glmnet
    # 3.) return the betas as a numeric vector
    # All of this is going to be done inparallel ... so should be good
    output <- foreach(i=seq(1,trainingIterations), .combine='rbind', .errorhandling=c('remove')) %dopar% {
        # First load required library(s)
        source('/Users/arose/Documents/forCobb/pncMJPS/classPaperScripts/figure6Functions.R')
        install_load('glmnet', 'bootstrap', 'psych', 'caret')
        
        # Now create our BS sample
        sampleVals <- as.numeric(unlist(createDataPartition(y)))
        bootY <- y[sampleVals]
        bootX <- inputData[sampleVals,]
        
        # First find the optimum values to use for glmnet
        optVals <- tuneAlpha(bootX, bootY, alphaSequence, nofFolds)
        
        
        # Now we need to create our final model with our opt vals
        mod <- glmnet(bootX, bootY, standardize=F, alpha=optVals[1], lambda=optVals[2], maxit=100000000, family="binomial")
        
        # Now return our values
        valsToReturn <- as.numeric(mod$beta)
        valsToReturn
    }
    # Kill our cluster
    stopCluster(cl)
    
    # Now return our output
    outputMatrix[,2:ncol(outputMatrix)] <- output
    return(outputMatrix)
}

tuneAlpha <- function(x, y, alphaSequence, nFolds){
    # Prep the outputs
    enet.optparam<-matrix(nrow=length(alphaSequence),ncol=3)
    colnames(enet.optparam)<-c("Alpha","Lambda","CVM")
    
    # Now make sure our data is in the rate format
    inputData <- as.matrix(x)
    
    # Lets iterate through each alpha an print the output to our output data frame
    count=1
    for(a in alphaSequence){
        enet.alphas.cv <- cv.glmnet(inputData, y, alpha=a, standardize=F, nfolds=nFolds, maxit=100000000, family="binomial")
        enet.optparam[count,] <- c(a,
        enet.alphas.cv$lambda.min,
        enet.alphas.cv$cvm[enet.alphas.cv$lambda==enet.alphas.cv$lambda.min])
        count <- count + 1
        
    }
    optval <- enet.optparam[which(enet.optparam[,3]==min(enet.optparam[,3])),]
    optAlpha <- optval[1]
    optLambda <- optval[2]
    output <- cbind(optAlpha, optLambda)
    return(output)
}


## Now create some functions which we can use to trim the fat from the lasso models
returnSelectionCol <- function(outputFromrunLasso){
    sumIndex <- rowSums(abs(sign(apply(outputFromrunLasso[,2:ncol(outputFromrunLasso)],
    2, function(x) as.numeric(as.character(x))))))
    sumIndex <- cbind(rownames(outputFromrunLasso), sumIndex)
    return(sumIndex)
}

rmFat <- function(outputFromrunLasso, imagingData, percentileToApply){
    # Find our cut off from the provided percentile
    cutoffToApply <- quantile(returnSelectionCol(outputFromrunLasso), .5)
    # First create our bool vector
    boo.vec <- rep('FALSE', nrow(outputFromrunLasso))
    sumIndex <- returnSelectionCol(outputFromrunLasso)
    boo.vec[which(sumIndex>=cutoffToApply)] <- 'TRUE'
    index <- which(boo.vec=='TRUE')
    # Now apply our boo vec to the imaging data
    output <- imagingData[,index]
    
    return(output)
}
