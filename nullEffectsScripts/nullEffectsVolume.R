## THis script will be used to prepare the density plots for volume
## with parametric and non parametric variance ratio and group differences
## If running local source this instead
install_load('psych', 'pwr', 'ggplot2', 'caret', 'olsrr', 'mgcv')

## Load the data
all.data <- readRDS('../scripts/07_MJEffects/scripts/mjAnovaData.RDS')
all.data <- all.data[which((all.data$ageAtScan1/12)>=14),]
## Cretae a function which will return a density plot for the input variables
writeDensityPlot <- function(data=all.data, var.of.interest=NULL, covariates="s(ageAtScan1)+sex+averageManualRating+race2+overall_psychopathology_ar_4factor",paraMetricValue=NULL, nonParametricValue=NULL, linMod=FALSE){
    ## First check we have a variable of interest
    if (missing(var.of.interest)) { stop("Gimme a ROI ya idiot; also make it a character string from your DF... ya idiot")}
    ## First thing we need to do is create our value
    ## Lets first regress out our values
    tmp.form <- as.formula(paste(var.of.interest, covariates, sep="~"))
    tmp.mod <- mgcv::gam(tmp.form, data=data)
    if(linMod==TRUE){
        tmp.mod <- lm(tmp.form, data=data)
    }
    
    ## Now create a residual index
    tmp.data <- data
    tmp.data$tmpvals <- NA
    index <- complete.cases(tmp.data[,var.of.interest])
    tmp.data$tmpvals[index] <- scale(as.numeric(residuals(tmp.mod)))
    
    ## Now create our density plot
    out.hist <- ggplot(tmp.data, aes(tmpvals)) +
      geom_density(data=subset(tmp.data,marcat=='MJ Non-User'), fill="#009E73",alpha=.4) +
      geom_density(data=subset(tmp.data,marcat=='MJ User'), fill="#9ad0f3", alpha=.4) +
      geom_density(data=subset(tmp.data,marcat=='MJ Frequent User'), fill="#D55E00",alpha=.4) +
      coord_cartesian(xlim=c(-5, 5), ylim=c(0,.75)) +
      xlab(var.of.interest) +
      theme_bw()
      
    ## Now create our text strings to add to the density plot
    if(!identical(paraMetricValue, NULL)){
      f.string <- paste("F-statistic = ",paraMetricValue)
      out.hist <- out.hist + annotate("text",  x=Inf, y = Inf, label = f.string, vjust=1.5, hjust=1, parse = F)
    }
    if(!identical(nonParametricValue, NULL)){
      chi.string <- paste("Kruskal-Wallis chi^2 = ",nonParametricValue)
      out.hist <- out.hist + annotate("text",  x=Inf, y = Inf, label = chi.string, vjust=3.5, hjust=1, parse=F)
    }
    # Now return the object
    return(out.hist)
}

## Now make a function which will return a legend color scale
returnColorScale <- function(){
    colfunc1 <- colorRampPalette(c("#009E73"))
    colfunc2 <- colorRampPalette(c("#9ad0f3"))
    colfunc3 <- colorRampPalette(c("#D55E00"))
    
    ## Now plot em
    plot(rep(1, 99), col = c(colfunc1(33), colfunc2(33), colfunc3(33)), pch = 15, cex = 10)
}

## Now identify our variables of interest
roi.of.interest <- c(1550:1552,109,110,111,112,114,115,127,128,129,130,131,132, 121, 122)

## Now reduce our data to the
all.data.plot <- all.data[complete.cases(all.data[,roi.of.interest]),]

## Now loop through and create our histogram and values
output.values.nl <- NULL
for(s in roi.of.interest){
    # Get our f value
    mod.one <- gam(all.data.plot[,s] ~ s(ageAtScan1) + sex + averageManualRating + factor(race2) + overall_psychopathology_ar_4factor + marcat, data=all.data.plot)
    aov.mod.one <- anova.gam(mod.one)
    # Now do the non parametric analysis of variance
    # We have to regress out all of the covaraites before we run the test
    all.data.plot$tmpvals <- as.numeric(scale(gam(all.data[,s] ~ s(ageAtScan1) + sex + averageManualRating + factor(race2) + overall_psychopathology_ar_4factor, data=all.data)$residuals))
    mod.two <- kruskal.test(tmpvals ~marcat, data = all.data.plot)
    ## Now do the variance ratio tests
    output.row <- c(names(all.data.plot)[s], round(aov.mod.one$pTerms.table['marcat',c('F')], digits=2), round(mod.two$statistic,digits=2),round(aov.mod.one$pTerms.table['marcat',c('p-value')], digits=2),round(mod.two$p.value, digits=2))
    ## Now attach this to the output
    output.values.nl <- rbind(output.values.nl, output.row)
}

output.values.lin <- NULL
for(s in roi.of.interest){
    # Get our f value
    mod.one <- lm(all.data.plot[,s] ~ ageAtScan1 + sex + averageManualRating + factor(race2) + overall_psychopathology_ar_4factor + marcat, data=all.data.plot)
    aov.mod.one <- aov(mod.one)
    # Now do the non parametric analysis of variance
    # We have to regress out all of the covaraites before we run the test
    all.data.plot$tmpvals <- as.numeric(scale(lm(all.data[,s] ~ ageAtScan1 + sex + averageManualRating + factor(race2) + overall_psychopathology_ar_4factor, data=all.data)$residuals))
    mod.two <- kruskal.test(tmpvals ~marcat, data = all.data.plot)
    ## Now do the variance ratio tests
    output.row <- c(names(all.data.plot)[s], round(summary(aov.mod.one)[[1]]['marcat',c('F value')], digits=2), round(mod.two$statistic,digits=2),round(summary(aov.mod.one)[[1]]['marcat',c('Pr(>F)')], digits=2),round(mod.two$p.value, digits=2))
    ## Now attach this to the output
    output.values.lin <- rbind(output.values.lin, output.row)
}

## Now do a linear model with squared and cubic terms
output.values.nl2 <- NULL
all.data.plot$ageSquared <- scale(all.data.plot$ageAtScan1)^2
all.data.plot$ageCubed <- scale(all.data.plot$ageAtScan1)^3
for(s in roi.of.interest){
    # Get our f value
    mod.one <- lm(all.data.plot[,s] ~ scale(ageAtScan1) + ageSquared + ageCubed + sex + averageManualRating + factor(race2) + overall_psychopathology_ar_4factor + marcat, data=all.data.plot)
    aov.mod.one <- aov(mod.one)
    # Now do the non parametric analysis of variance
    # We have to regress out all of the covaraites before we run the test
    all.data.plot$tmpvals <- as.numeric(scale(lm(all.data[,s] ~ ageAtScan1 + sex + averageManualRating + factor(race2) + overall_psychopathology_ar_4factor, data=all.data)$residuals))
    mod.two <- kruskal.test(tmpvals ~marcat, data = all.data.plot)
    ## Now do the variance ratio tests
    output.row <- c(names(all.data.plot)[s], round(summary(aov.mod.one)[[1]]['marcat',c('F value')], digits=2), round(mod.two$statistic,digits=2),round(summary(aov.mod.one)[[1]]['marcat',c('Pr(>F)')], digits=2),round(mod.two$p.value, digits=2))
    ## Now attach this to the output
    output.values.nl2 <- rbind(output.values.nl2, output.row)
}

## Now produce a density plot for each variable
plot.tbv <- writeDensityPlot(data=all.data.plot, paraMetricValue=output.values.nl[1,2], nonParametricValue=output.values.nl[1,3], var.of.interest='mprage_jlf_vol_TBV')
plot.tbgm <- writeDensityPlot(data=all.data.plot, paraMetricValue=output.values.nl[3,2], nonParametricValue=output.values.nl[3,3], var.of.interest=output.values.nl[2,1],linMod=FALSE)
plot.tbwm <- writeDensityPlot(data=all.data.plot, paraMetricValue=output.values.nl[3,2], nonParametricValue=output.values.nl[3,3], var.of.interest=output.values.nl[3,1],linMod=FALSE)
plot.acc.r <- writeDensityPlot(data=all.data.plot, paraMetricValue=output.values.nl[4,2], nonParametricValue=output.values.nl[4,3], var.of.interest=output.values.nl[4,1],linMod=FALSE)
plot.amy.r <- writeDensityPlot(data=all.data.plot, paraMetricValue=output.values.nl[6,2], nonParametricValue=output.values.nl[6,3], var.of.interest=output.values.nl[6,1],linMod=FALSE)
plot.cau.r <- writeDensityPlot(data=all.data.plot, paraMetricValue=output.values.nl[8,2], nonParametricValue=output.values.nl[8,3], var.of.interest=output.values.nl[8,1],linMod=FALSE)
plot.pal.r <- writeDensityPlot(data=all.data.plot, paraMetricValue=output.values.nl[10,2], nonParametricValue=output.values.nl[10,3], var.of.interest=output.values.nl[10,1],linMod=FALSE)
plot.put.r <- writeDensityPlot(data=all.data.plot, paraMetricValue=output.values.nl[12,2], nonParametricValue=output.values.nl[12,3], var.of.interest=output.values.nl[12,1],linMod=FALSE)
plot.tha.r <- writeDensityPlot(data=all.data.plot, paraMetricValue=output.values.nl[14,2], nonParametricValue=output.values.nl[14,3], var.of.interest=output.values.nl[14,1],linMod=FALSE)
plot.hip.r <- writeDensityPlot(data=all.data.plot, paraMetricValue=output.values.nl[16,2], nonParametricValue=output.values.nl[16,3], var.of.interest=output.values.nl[16,1],linMod=FALSE)
multiplot(plot.tbv, plot.tbgm, plot.tbwm, plot.acc.r, plot.amy.r, plot.cau.r, plot.pal.r, plot.put.r, plot.tha.r, plot.hip.r, cols=4)
