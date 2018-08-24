## THis script will be used to prepare the density plots for volume
## with parametric and non parametric variance ratio and group differences
## If running local source this instead
install_load('psych', 'pwr', 'ggplot2', 'caret', 'olsrr', 'mgcv')

## Load the data
all.data <- readRDS('mjAnovaData.RDS')
## Cretae a function which will return a density plot for the input variables
writeDensityPlot <- function(data=all.data, var.of.interest=NULL, covariates="s(ageAtScan1)+sex+averageManualRating+factor(race2)+overall_psychopathology_ar_4factor",paraMetricValue=NULL, nonParametricValue=NULL, linMod=FALSE){
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
    
    ## Now create a new x axis name
    x.ax.name <- gsub(x=var.of.interest, replacement="", pattern="mprage_jlf_vol_")
    
    ## Now create our density plot
    out.hist <- ggplot(tmp.data, aes(tmpvals)) +
    geom_density(data=subset(tmp.data,marcat=='MJ Non-User'), fill="#009E73",alpha=.4) +
    geom_density(data=subset(tmp.data,marcat=='MJ Occ User'), fill="#9ad0f3", alpha=.4) +
    geom_density(data=subset(tmp.data,marcat=='MJ Freq User'), fill="#D55E00",alpha=.4) +
    coord_cartesian(xlim=c(-5, 5), ylim=c(0,.9)) +
    xlab(x.ax.name) +
    theme_bw() +
    theme(text = element_text(size=30))
    
    ## Now create our text strings to add to the density plot
    if(!identical(paraMetricValue, NULL)){
        f.string <- paste("F-statistic = ",round(as.numeric(paraMetricValue), digits=3), " ")
        out.hist <- out.hist + annotate("text",  x=Inf, y = Inf, label = f.string, vjust=1.5, hjust=1, parse = F, size=10)
    }
    if(!identical(nonParametricValue, NULL)){
        chi.string <- paste("Kruskal-Wallis H = ",round(as.numeric(nonParametricValue), digits=3), " ")
        out.hist <- out.hist + annotate("text",  x=Inf, y = Inf, label = chi.string, vjust=3.5, hjust=1, parse=F, size=10)
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
roi.of.interest <- c(1540,1565:1576)

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

## Now produce a density plot for each variable
pdf('fig2.pdf')#,units='in',res=300,width=10,height=8)
for(i in 1:dim(output.values.nl)[1]){
    tmp.plot <- writeDensityPlot(data=all.data.plot, paraMetricValue=output.values.nl[i,2], nonParametricValue=output.values.nl[i,3], var.of.interest=output.values.nl[i,1])
    ## Now prepare a new x axis name
    tmp.string <- gsub(x=output.values.nl[i,1], pattern="mprage_jlfLobe_ct_", replacement="")
    tmp.string <- gsub(x=tmp.string, pattern="L_", replacement="Left ")
    tmp.string <- gsub(x=tmp.string, pattern="R_", replacement="Right ")
    tmp.string <- gsub(x=tmp.string, pattern="_", replacement=" ")
    if(i ==1){tmp.string <- "Mean Cortical Thickness"}
    # Now remove the "density fdrom the graph"
    tmp.plot <- tmp.plot + theme(axis.title.y = element_text(color="white")) + xlab(tmp.string)
    print(tmp.plot)
}
dev.off()

