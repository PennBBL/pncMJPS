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
install_load('psych','ggplot2','caret','equivalence', 'mgcv','TOSTER','foreach','doParallel', 'simpleboot')

## Now load the data
all.data <- readRDS('mjAnovaData.RDS')
 
### Begin bootstrapping down here
## First create our resamples
tmp.folds <- createResample(y=all.data$marcat, 100)

## Now regress out age and sex so we can compare our groups
orig <- all.data
base.model <- paste("s(ageAtScan1) + sex + averageManualRating + marcat + factor(race2) + overall_psychopathology_ar_4factor")
vars.of.interest <- c(107:245, 255:352, 353:470,471,1540,1550:1588)
for(v in vars.of.interest){
  name.val <- names(all.data)[v]
  tmp.formula <- as.formula(paste(name.val, "~", base.model))
  tmp.mod <- gam(tmp.formula, data=all.data)
  index <- which(complete.cases(all.data[,v]))
  all.data[index,name.val] <- NA
  all.data[index,name.val] <- scale(residuals(tmp.mod))
}

## Now I need to loop through each ROI
## and calculate the groupwise differences
cl <- makeCluster(8)
registerDoParallel(cl)
output.ci <- foreach(q=1:length(vars.of.interest), .combine='rbind', .packages='simpleboot') %dopar% {
    # First isolate our values
    v <- vars.of.interest[q]
    name.val <- names(all.data)[v]
    tmp.dat <- all.data[,c(name.val, 'marcat')]
    non.vals <- tmp.dat[which(tmp.dat[,2]=="MJ Non-User"),1]
    occ.vals <- tmp.dat[which(tmp.dat[,2]=="MJ Occ User"),1]
    fre.vals <- tmp.dat[which(tmp.dat[,2]=="MJ Freq User"),1]
    # Now calculate our bootstrapped differences
    n.v.o <- two.boot(non.vals, occ.vals, mean, 10000)
    n.v.f <- two.boot(non.vals, fre.vals, mean, 10000)
    f.v.o <- two.boot(fre.vals, occ.vals, mean, 10000)
    # Now return our confidence intervals
    out.row <- c(name.val,boot.ci(n.v.o)$basic[4:5],boot.ci(n.v.f)$basic[4:5],boot.ci(f.v.o)$basic[4:5])
    out.row
}
colnames(output.ci) <- c("ROI", "nvoMin", "nvoMax", "nvfMin", "nvfMax", "fvoMin", "fvoMax")
## Now find our values that have confidence intervals that do not include 0
## quickly declare a rough function for this
is.between <- function(x, a, b) {
    a <- as.numeric(as.character(a))
    b <- as.numeric(as.character(b))
    (x - a)  *  (b - x) > 0
}
output.ci[which(!is.between(0, output.ci[,2], output.ci[,3])),]
output.ci[which(!is.between(0, output.ci[,4], output.ci[,5])),]
output.ci[which(!is.between(0, output.ci[,6], output.ci[,7])),]

## Now produce our histograms
pdf('user.minus.non.pdf')
for(q in c(358,359,360)){
  ## First grab all of our values
  v <- vars.of.interest[q]
  name.val <- names(all.data)[v]
  tmp.dat <- all.data[,c(name.val, 'marcat')]
  non.vals <- tmp.dat[which(tmp.dat[,2]=="MJ Non-User"),1]
  occ.vals <- tmp.dat[which(tmp.dat[,2]=="MJ Occ User"),1]
  n.v.o <- two.boot(non.vals, occ.vals, mean, 1000, student=T, M=50)
  ci.vals <- round(boot.ci(n.v.o)$student[4:5], digits=3)
  ci.string <- paste("95% CI [", ci.vals[1], ",", ci.vals[2], "] ", sep='')
  ## Now plot our histogram for these values
  out.plot <- ggplot() +
    geom_histogram(aes(n.v.o$t[,1]-n.v.o$t[,2]), color='red', alpha=.3, bins=100) +
    geom_segment(aes(x=ci.vals[1],xend=ci.vals[1], y=0, yend=40)) +
    geom_segment(aes(x=ci.vals[2],xend=ci.vals[2],y=0, yend=40)) +
    ggtitle('') +
    coord_cartesian(xlim=c(-.5, .5), ylim=c(0,50)) +
    theme_bw() +
    annotate("text",  x=Inf, y = Inf, label = ci.string, vjust=3.5, hjust=1, parse = F, size=8) +
    theme(text = element_text(size=30), axis.title.x = element_blank(), axis.title.y = element_blank())
  print(out.plot)
  print(q)
}
dev.off()
pdf('freq.minus.non.pdf')
for(q in c(358,359,360)){
    ## First grab all of our values
    v <- vars.of.interest[q]
    name.val <- names(all.data)[v]
    tmp.dat <- all.data[,c(name.val, 'marcat')]
    non.vals <- tmp.dat[which(tmp.dat[,2]=="MJ Non-User"),1]
    fre.vals <- tmp.dat[which(tmp.dat[,2]=="MJ Freq User"),1]
    n.v.f <- two.boot(non.vals, fre.vals, mean, 1000, student=T, M=50)
    ci.vals <- round(boot.ci(n.v.o)$student[4:5], digits=3)
    ci.string <- paste("95% CI [", ci.vals[1], ",", ci.vals[2], "] ", sep='')
    ## Now plot our histogram for these values
    out.plot <- ggplot() +
    geom_histogram(aes(n.v.f$t[,1]-n.v.o$t[,2]), color='red', alpha=.3, bins=100) +
    geom_segment(aes(x=ci.vals[1],xend=ci.vals[1], y=0, yend=40)) +
    geom_segment(aes(x=ci.vals[2],xend=ci.vals[2],y=0, yend=40)) +
    ggtitle('') +
    #xlab(name.val) +
    coord_cartesian(xlim=c(-.5, .5), ylim=c(0,50)) +
    theme_bw() +
    annotate("text",  x=Inf, y = Inf, label = ci.string, vjust=3.5, hjust=1, parse = F, size=8) +
    theme(text = element_text(size=30), axis.title.x = element_blank(), axis.title.y = element_blank())
    print(out.plot)
}
dev.off()
pdf('freq.minus.user.pdf')
for(q in c(358,359,360)){
    ## First grab all of our values
    v <- vars.of.interest[q]
    name.val <- names(all.data)[v]
    tmp.dat <- all.data[,c(name.val, 'marcat')]
    occ.vals <- tmp.dat[which(tmp.dat[,2]=="MJ Occ User"),1]
    fre.vals <- tmp.dat[which(tmp.dat[,2]=="MJ Freq User"),1]
    f.v.o <- two.boot(fre.vals, occ.vals, mean, 1000, student=T, M=50)
    ci.vals <- round(boot.ci(n.v.o)$student[4:5], digits=3)
    ci.string <- paste("95% CI [", ci.vals[1], ",", ci.vals[2], "] ", sep='')
    ## Now plot our histogram for these values
    out.plot <- ggplot() +
    geom_histogram(aes(f.v.o$t[,1]-n.v.o$t[,2]), color='red', alpha=.3, bins=100) +
    geom_segment(aes(x=ci.vals[1],xend=ci.vals[1], y=0, yend=40)) +
    geom_segment(aes(x=ci.vals[2],xend=ci.vals[2],y=0, yend=40)) +
    ggtitle('') +
    coord_cartesian(xlim=c(-.5, .5), ylim=c(0,50)) +
    theme_bw() +
    annotate("text",  x=Inf, y = Inf, label = ci.string, vjust=3.5, hjust=1, parse = F, size=8) +
    theme(text = element_text(size=30), axis.title.x = element_blank(), axis.title.y = element_blank())
    print(out.plot)
}
dev.off()

#### Now run equivalence testing down here
equiv.data <- orig
base.model <- paste("s(ageAtScan1) + sex + averageManualRating + marcat + factor(race2) + overall_psychopathology_ar_4factor")
vars.of.interest <- c(107:245, 255:352, 353:470,471,1540,1550:1588)
for(v in vars.of.interest){
  name.val <- names(all.data)[v]
  tmp.formula <- as.formula(paste(name.val, "~", base.model))
  tmp.col <- rep(NA, 1504)
  tmp.mod <- gam(tmp.formula, data=equiv.data)
  index <- which(complete.cases(all.data[,v]))
  equiv.data[index,name.val] <- NA
  equiv.data[index,name.val] <- scale(residuals(tmp.mod))
}

## Now run through all of our tost tests!
output.tost.vals <- NULL
for(v in vars.of.interest){
  name.val <- names(all.data)[v]
  test.val.one <- tost(x=equiv.data[which(equiv.data$marcat=='MJ Occ User'),v], y=equiv.data[which(equiv.data$marcat=="MJ Non-User"),v], paired=F, var.equal=F)
  test.val.two <- tost(x=equiv.data[which(equiv.data$marcat=='MJ Freq User'),v], y=equiv.data[which(equiv.data$marcat=="MJ Non-User"),v], paired=F, var.equal=F)
  test.val.three <- tost(x=equiv.data[which(equiv.data$marcat=='MJ Occ User'),v], y=equiv.data[which(equiv.data$marcat=="MJ Freq User"),v], paired=F, var.equal=F)
  print(v)
  test.val.four <- tost(x=equiv.data[which(equiv.data$marcat=='MJ Occ User'),v], y=equiv.data[which(equiv.data$marcat=="MJ Non-User"),v], paired=F, epsilon = .6, var.equal=F)
  test.val.five <- tost(x=equiv.data[which(equiv.data$marcat=='MJ Freq User'),v], y=equiv.data[which(equiv.data$marcat=="MJ Non-User"),v], paired=F, epsilon = .6, var.equal=F)
  test.val.six <- tost(x=equiv.data[which(equiv.data$marcat=='MJ Occ User'),v], y=equiv.data[which(equiv.data$marcat=="MJ Freq User"),v], paired=F, epsilon = .6, var.equal=F)
  ## Now prepare our output values
  output.row <- c(name.val, test.val.one$tost.p.value, test.val.two$tost.p.value, test.val.three$tost.p.value, test.val.four$tost.p.value, test.val.five$tost.p.value, test.val.six$tost.p.value)
  output.tost.vals <- rbind(output.tost.vals, output.row)
}
# Now perform fdr correction on our p values
rownames(output.tost.vals) <- NULL
## Now perform FDR correction
row.grep.patterns <- c("mprage_jlf_vol_","mprage_jlf_ct_","mprage_jlf_gmd_","mprage_jlfLobe_vol_","mprage_jlfLobe_ct_","mprage_jlfLobe_gmd_")
output.tost.corr.vals <- matrix(NA, ncol=dim(output.tost.vals)[2]-1, nrow=dim(output.tost.vals)[1])
for(i in 1:dim(output.tost.corr.vals)[2]){
  orig.vals.col <- i + 1
  for(gpat in row.grep.patterns){
      new.vals <- p.adjust(output.tost.vals[grep(gpat, output.tost.vals[,1]),orig.vals.col], method='fdr')
      output.tost.corr.vals[grep(gpat, output.tost.vals[,1]),i] <- new.vals
  }
}
output.tost.vals <- cbind(output.tost.vals, output.tost.corr.vals)
output.tost.vals <- as.data.frame(output.tost.vals)
output.tost.vals[,2:13] <- apply(output.tost.vals[,-1], 2, function(x) as.numeric(as.character(x)))
colnames(output.tost.vals) <- c('ROI','MJUser.MJNonUser.5','Freq.MJNonUser.5','MJUser.Freq.5','MJUser.MJNonUser.3','Freq.MJNonUser.3','MJUser.Freq.3','MJUser.MJNonUser.5.fdr','Freq.MJNonUser.5.fdr','MJUser.Freq.5.fdr','MJUser.MJNonUser.3.fdr','Freq.MJNonUser.3.fdr','MJUser.Freq.3.fdr')

## Now combine all of these values
output <- merge(output, output.tost.vals, by='ROI')
write.csv(output, "bsPValResults.csv", quote=F, row.names=F)
