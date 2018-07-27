# AFGR June 2018 

## This script will be used to test for equivalence in our marcat groups
## There are three group comparisons we have to make:
##	1. User - Non
##	2. Freq - Non
##	3. Freq - User
## Differences will be calculated by taking the bootstrapped differences in means
## This will give us a confidence interval as well..

## Load library(s)
#source('/home/arosen/adroseHelperScripts/R/afgrHelpFunc.R')
install_load('psych','ggplot2','caret','equivalence', 'mgcv', 'TOSTER')

## Now load the data
all.data <- readRDS('mjAnovaData.RDS')
all.data$marcat <- factor(all.data$marcat, levels=c("MJ Non-User", "MJ Occ User", "MJ Freq User"))
alc.data <- read.csv('alcData.csv')
all.data <- merge(all.data, alc.data)
 
### Begin bootstrapping down here
## First create our resamples
tmp.folds <- createResample(y=all.data$marcat, 1000)

## Now regress out age and sex so we can compare our groups
orig <- all.data
base.model <- paste("s(ageAtScan1) + sex + averageManualRating + marcat + factor(race2) + overall_psychopathology_ar_4factor + alchoholZScore")
vars.of.interest <- c(107:245, 255:352, 353:470,471,1540,1550:1588)
for(v in vars.of.interest){
  name.val <- names(all.data)[v]
  tmp.formula <- as.formula(paste(name.val, "~", base.model))
  tmp.col <- rep(NA, 1504)
  tmp.mod <- gam(tmp.formula, data=all.data)
  index <- which(complete.cases(all.data[,c(name.val, "alchoholZScore")]))
  all.data[index,name.val] <- NA
  all.data[index,name.val] <- scale(residuals(tmp.mod))
}

## Now go through each value and find our bootstrapped mean differences
output.differences <- NULL
for(q in 1:length(tmp.folds)){
  tmpData <- all.data[tmp.folds[[q]],]
  output.vals.tmp <- NULL
  for(v in vars.of.interest){
    name.val <- names(all.data)[v]
    mean.values <- summarySE(data=tmpData, measurevar=name.val, groupvars='marcat', na.rm=T)
    # Now prepare the output data
    output.row <- c(name.val, t(mean.values[,name.val]))
    output.vals.tmp <- rbind(output.vals.tmp, output.row)
  }
  output.differences <- rbind(output.differences, output.vals.tmp)
  print(q)
}
colnames(output.differences) <- c('ROI', 'MJ Frequent', 'Non-User', 'MJ User')
rownames(output.differences) <- NULL
output.differences <- as.data.frame(output.differences)

## Now find our mean differences
output.differences$user.minus.non <- as.numeric(as.character(output.differences[,4])) - as.numeric(as.character(output.differences[,3]))
output.differences$freq.minus.non <- as.numeric(as.character(output.differences[,2])) - as.numeric(as.character(output.differences[,3]))
output.differences$freq.minus.user <- as.numeric(as.character(output.differences[,2])) - as.numeric(as.character(output.differences[,4]))

## Now write this data
write.csv(output.differences, "mjBSValsAlc.csv", quote=F, row.names=F)

## Now get our means and confidence intervals four our differences
mean.vals.1 <- summarySE(data=output.differences, groupvars='ROI', measurevar='user.minus.non')
mean.vals.2 <- summarySE(data=output.differences, groupvars='ROI', measurevar='freq.minus.non')
mean.vals.3 <- summarySE(data=output.differences, groupvars='ROI', measurevar='freq.minus.user')

## And now find those regions that have a min 95% CI > .3
mean.vals.1$flag <- 0
mean.vals.1$flag[which(abs(mean.vals.1[,'user.minus.non']) - mean.vals.1[,'ci']>.3)] <- 1
mean.vals.2$flag <- 0
mean.vals.2$flag[which(abs(mean.vals.2[,'freq.minus.non']) - mean.vals.2[,'ci']>.3)] <- 1
mean.vals.3$flag <- 0
mean.vals.3$flag[which(abs(mean.vals.3[,'freq.minus.user']) - mean.vals.3[,'ci']>.3)] <- 1

## Now grab some p values for our ROI's
## the p value will be obtained by calculating the # of times
## the absolute mean difference is greater than .3
## This proportion will be our p value!
output.differences$user.minus.non.bin <- 0
output.differences$user.minus.non.bin[which(abs(output.differences$user.minus.non)>=.3)] <- 1
output.differences$freq.minus.non.bin <- 0
output.differences$freq.minus.non.bin[which(abs(output.differences$freq.minus.non)>=.3)] <- 1
output.differences$freq.minus.user.bin <- 0
output.differences$freq.minus.user.bin[which(abs(output.differences$freq.minus.user)>=.3)] <- 1

## Now find the bootstrapped p value w/in each ROI
bs.p.vals.out <- NULL
for(r in mean.vals.1$ROI){
  u.m.n <- 1 - dim(output.differences[which(output.differences$ROI==r & output.differences$user.minus.non.bin==1),])[1]/1000
  f.m.n <- 1 - dim(output.differences[which(output.differences$ROI==r & output.differences$freq.minus.non.bin==1),])[1]/1000
  f.m.u <- 1 - dim(output.differences[which(output.differences$ROI==r & output.differences$freq.minus.user.bin==1),])[1]/1000
  out.row <- c(r, u.m.n, f.m.n, f.m.u)
  bs.p.vals.out <- rbind(bs.p.vals.out, out.row)
}
rownames(bs.p.vals.out) <- NULL
colnames(bs.p.vals.out) <- c('ROI', 'umnPVal', 'fmnPVal', 'fmuPVal')
bs.p.vals.out <- as.data.frame(bs.p.vals.out)
bs.p.vals.out[,2:4] <- apply(bs.p.vals.out[,2:4], 2, function(x) as.numeric(as.character(x)))
mean.vals.1 <- merge(mean.vals.1, bs.p.vals.out[,c(1,2)])
mean.vals.2 <- merge(mean.vals.2, bs.p.vals.out[,c(1,3)])
mean.vals.3 <- merge(mean.vals.3, bs.p.vals.out[,c(1,4)])

## Now write the output
output <- merge(mean.vals.1, mean.vals.2, by='ROI')
output <- merge(output, mean.vals.3, by='ROI')

## Now produce our histograms
pdf('user.minus.non.pdf')
for(r in mean.vals.1$ROI){
  ## First grab all of our values
  tmp.dat <- output.differences[which(output.differences$ROI==r),]
  mean.value <- mean.vals.1[which(mean.vals.1$ROI==r),'user.minus.non']
  ci.value <- mean.vals.1[which(mean.vals.1$ROI==r),'ci']
  p.val.string <- paste("p-value = ", mean.vals.1[which(mean.vals.1$ROI==r),'umnPVal'], sep='')
  ## Now plot our histogram for these values
  out.plot <- ggplot(tmp.dat) +
    geom_histogram(aes(user.minus.non), color='red', alpha=.3, bins=100) +
    geom_vline(xintercept=mean.value) +
    geom_vline(xintercept=mean.value+ci.value) +
    geom_vline(xintercept=mean.value-ci.value) +
    geom_vline(xintercept=-.3, linetype="dotted") +
    geom_vline(xintercept=.3, linetype="dotted") +
    ggtitle('') +
    coord_cartesian(xlim=c(-.7, .7), ylim=c(0,50)) +
    theme_bw() +
    annotate("text",  x=Inf, y = Inf, label = p.val.string, vjust=1.5, hjust=1, parse = F, size=10) +
    theme(text = element_text(size=30), axis.title.x = element_blank(), axis.title.y = element_blank())
  print(out.plot)
}
dev.off()
pdf('freq.minus.non.pdf')
for(r in mean.vals.1$ROI){
  ## First grab all of our values
  tmp.dat <- output.differences[which(output.differences$ROI==r),]
  mean.value <- mean.vals.2[which(mean.vals.1$ROI==r),'freq.minus.non']
  ci.value <- mean.vals.2[which(mean.vals.1$ROI==r),'ci']
  p.val.string <- paste("p-value = ", mean.vals.2[which(mean.vals.1$ROI==r),'fmnPVal'], sep='')
  ## Now plot our histogram for these values
  out.plot <- ggplot(tmp.dat) +
    geom_histogram(aes(freq.minus.non), color='red', alpha=.3, bins=100) +
    geom_vline(xintercept=mean.value) +
    geom_vline(xintercept=mean.value+ci.value) +
    geom_vline(xintercept=mean.value-ci.value) +
    geom_vline(xintercept=-.3, linetype="dotted") +
    geom_vline(xintercept=.3, linetype="dotted") +
    ggtitle('') +
    coord_cartesian(xlim=c(-.7, .7), ylim=c(0,50)) +
    theme_bw() +
    annotate("text",  x=Inf, y = Inf, label = p.val.string, vjust=1.5, hjust=1, parse = F, size=10) +
    theme(text = element_text(size=30), axis.title.x = element_blank(), axis.title.y = element_blank())
  print(out.plot)
}
dev.off()
pdf('freq.minus.user.pdf')
for(r in mean.vals.1$ROI){
  ## First grab all of our values
  tmp.dat <- output.differences[which(output.differences$ROI==r),]
  mean.value <- mean.vals.3[which(mean.vals.3$ROI==r),'freq.minus.user']
  ci.value <- mean.vals.3[which(mean.vals.3$ROI==r),'ci']
  p.val.string <- paste("p-value = ", mean.vals.3[which(mean.vals.3$ROI==r),'fmuPVal'], sep='')
  ## Now plot our histogram for these values
  out.plot <- ggplot(tmp.dat) +
    geom_histogram(aes(freq.minus.user), color='red', alpha=.3, bins=100) +
    geom_vline(xintercept=mean.value) +
    geom_vline(xintercept=mean.value+ci.value) +
    geom_vline(xintercept=mean.value-ci.value) +
    geom_vline(xintercept=-.3, linetype="dotted") +
    geom_vline(xintercept=.3, linetype="dotted") +
    ggtitle('') +
    coord_cartesian(xlim=c(-.7, .7), ylim=c(0,50)) +
    theme_bw() +
    annotate("text",  x=Inf, y = Inf, label = p.val.string, vjust=1.5, hjust=1, parse = F, size=10) +
    theme(text = element_text(size=30), axis.title.x = element_blank(), axis.title.y = element_blank())
  print(out.plot)
}
dev.off()

#### Now run equivalence testing down here
equiv.data <- orig
base.model <- paste("s(ageAtScan1) + sex + averageManualRating + marcat + factor(race2) + overall_psychopathology_ar_4factor+alchoholZScore")
vars.of.interest <- c(107:245, 255:352, 353:470,471,1540,1550:1588)
for(v in vars.of.interest){
  name.val <- names(all.data)[v]
  tmp.formula <- as.formula(paste(name.val, "~", base.model))
  tmp.col <- rep(NA, 1504)
  tmp.mod <- gam(tmp.formula, data=equiv.data)
  index <- which(complete.cases(all.data[,c(name.val, "alchoholZScore")]))
  equiv.data[index,name.val] <- NA
  equiv.data[index,name.val] <- scale(residuals(tmp.mod))
}

## Now run through all of our tost tests!
output.tost.vals <- NULL
for(v in vars.of.interest){
  name.val <- names(all.data)[v]
  test.val.one <- tost(x=equiv.data[which(equiv.data$marcat=='MJ Occ User'),v], y=equiv.data[which(equiv.data$marcat=="MJ Non-User"),v], paired=F)
  test.val.two <- tost(x=equiv.data[which(equiv.data$marcat=='MJ Freq User'),v], y=equiv.data[which(equiv.data$marcat=="MJ Non-User"),v], paired=F)
  test.val.three <- tost(x=equiv.data[which(equiv.data$marcat=='MJ Occ User'),v], y=equiv.data[which(equiv.data$marcat=="MJ Freq User"),v], paired=F)
  print(v)
  mean.vals <- summarySE(data=equiv.data, measurevar=v, groupvars='marcat', na.rm=T)
  tost.one <- TOSTtwo(m1=mean.vals[2,3], sd1=mean.vals$sd[2], m2=mean.vals[3,3], sd2=mean.vals$sd[3], n1=mean.vals$N[2], n2=mean.vals$N[3], var.equal=F, low_eqbound=-.6, high_eqbound=.6)
  tost.two <- TOSTtwo(m1=mean.vals[2,3], sd1=mean.vals$sd[2], m2=mean.vals[1,3], sd2=mean.vals$sd[1], n1=mean.vals$N[2], n2=mean.vals$N[1], var.equal=F, low_eqbound=-.6, high_eqbound=.6)
  tost.three <- TOSTtwo(m1=mean.vals[3,3], sd1=mean.vals$sd[3], m2=mean.vals[1,3], sd2=mean.vals$sd[1], n1=mean.vals$N[3], n2=mean.vals$N[1], var.equal=F, low_eqbound=-.6, high_eqbound=.6)
  ## Now prepare our output values
  output.row <- c(name.val, test.val.one$tost.p.value, test.val.two$tost.p.value, test.val.three$tost.p.value, tost.one$TOST_p1, tost.two$TOST_p1, tost.three$TOST_p1)
  output.tost.vals <- rbind(output.tost.vals, output.row)
}
# Now perform fdr correction on our p values
rownames(output.tost.vals) <- NULL
output.tost.vals <- as.data.frame(cbind(output.tost.vals, apply(output.tost.vals[,2:4], 2, p.adjust)))
output.tost.vals[,2:7] <- apply(output.tost.vals[,-1], 2, function(x) as.numeric(as.character(x)))
colnames(output.tost.vals) <- c('ROI','MJUser.MJNonUser','Freq.MJNonUser','MJUser.Freq','MJUser.MJNonUser.Bon','Freq.MJNonUser.Bon','MJUser.Freq.Bon')

## Now combine all of these values
output <- merge(output, output.tost.vals, by='ROI')
write.csv(output, "bsPValResultsAlc.csv", quote=F, row.names=F)
