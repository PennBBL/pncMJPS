# Load data
alc.data <- read.csv("~/Downloads/go1_alc_tob_14up_may20.csv")

## Now I need to find those that are 1's but have na values
alc.columns <- c(3:14)

# Now go through each  subject
# if they recieved the questions, but have NA values, set all of their values to 0
out.data <- alc.data
for(i in 1:dim(alc.data)[1]){
    tmp.vals <- out.data[i,c(2,alc.columns)]
    if(tmp.vals[1]==1){
        tmp.vals[is.na(tmp.vals)] <- 0
    }
    out.data[i,c(2,alc.columns)] <- tmp.vals
}

## Now z score w/in our alchohol data
z.score.vals <- apply(out.data[,alc.columns], 2, function(x) as.numeric(scale(x)))

## Now create an overall alc score
z.tot <- apply(z.score.vals, 1, function(x) mean(x, na.rm=T))

## Now append our z.tot to our output data
out.data$alchoholZScore <- z.tot

## Now write the csv
write.csv(out.data, "~/alcData.csv", quote=F, row.names=F)
