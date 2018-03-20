# AFGR March 2018 
# This script will be used to run a voxel wise analysis in R 
# for the dosage data - it'll run a binary dosage
# a nonlinear dosage w/ occasiaonl and frequent

## Load library(s)
source("/home/arosen/adroseHelperScripts/R/afgrHelpFunc.R")
install_load('ANTsR','mgcv','utils')

## Load the data
# Start with the imaging data 
input.images <- antsImageRead('/data/joy/BBL/projects/pncMJPS/data/mytbss/n391_pathVals_include_smooth0/fourd.nii.gz',4)
input.mask <- antsImageRead('/data/joy/BBL/projects/pncMJPS/data/mytbss/stats/mean_FA_skeleton_mask.nii.gz',3)
input.matrix.vals <- timeseries2matrix(input.images,input.mask)
# Now do the demographic values
input.demo.vals <- readRDS("/data/joy/BBL/projects/pncMJPS/data/mytbss/subjectData.RDS")

## Now create our data frame to work with
all.values <- cbind(input.matrix.vals, input.demo.vals)

## Find the inital value
check.point.vals <- read.csv("/data/joy/BBL/projects/pncMJPS/data/tmpVoxelVals.csv")
# Now find the last complete.case
max.value <-  max(which(complete.cases(check.point.vals)))

## Now loop through each voxel
output.values <- matrix(NA, dim(input.matrix.vals)[2],9)
output.values <- check.point.vals
rm(check.point.vals)
pb <- txtProgressBar(min=0, dim(input.matrix.vals)[2], initial=max.value, style=3)
for(i in seq(max.value, dim(input.matrix.vals)[2])){
  if(sum(all.values[,i])==0){
    output.values[i,1] <- i
  }
  else{
    tmpMod <- gam(all.values[,i] ~ s(age) + sex + dosage + dti64Tsnr, data=all.values)
    outputVals1 <- as.vector(summary(tmpMod)$s.table[3:4])
    outputVals2 <- as.vector(t(summary(tmpMod)$p.table[2:4,3:4]))
    outputRow <- c(i, outputVals1, outputVals2)
    output.values[i,] <- outputRow
  }
  setTxtProgressBar(pb, i)
  if((i %% 1000)==0){
    write.csv(output.values, "/data/joy/BBL/projects/pncMJPS/data/tmpVoxelVals.csv", quote=F, row.names=F)
  }
}
write.csv(output.values, "/data/joy/BBL/projects/pncMJPS/data/tmpVoxelVals.csv", quote=F, row.names=F)
## Now I need to write out these values

