install_load('gdata', 'WriteXLS', 'xlsx')
is.letter <- function(x) grepl("[:alpha:]", x)
foo <-  read.xls('~/Downloads/mj-meta-database-JJ-for-adon.xls', na.strings=c('999', '99.99'))
bar <-  read.xls('~/Downloads/mj-meta-database-SS-for-adon.xls', na.strings=c('999', '99.99'))
foo1 <-apply(foo, c(1,2), tolower)
bar1 <-apply(bar, c(1,2), tolower)
tmp <-merge(foo1, bar1, by=c('tstname', 'authors', 'pubyear'), suffixes=c(".JJ", ".SS"))
tmp <- tmp[-1,]
# Now create an output matrix - same dimensions as tmp
# which will store boolean values for if things mathc or not
outMat <- matrix(0, ncol=dim(tmp)[2], nrow=dim(tmp)[1])
colnames(outMat) <- names(tmp)

#Now find the dimension of the columns we want to grep over
tmp <- tmp[,-38]
searchVals <- names(tmp)[4:139]
output <- tmp[,c(1,2,3)]
numericValue <- NA
toFix <- NA
for(i in searchVals){
  grepPattern <- gsub(x=i, pattern='.JJ', replacement='')
  print(grepPattern)
  colVals <- grep(grepPattern, names(tmp))
  if(length(colVals)==1){
    break
  }
  if(length(colVals)>3){
    print("Now adjusting grep mismatch")
    colVals <- colVals[c(1,3)]
  }
  # Now go through each row and see if the values are identical after rounding them to 2 sig figs
  outCol <- rep(0, dim(tmp)[1])
  checkVal <- is.letter(tmp[,colVals[1]])
  checkVal <- append(checkVal, is.letter(tmp[,colVals[2]]))
  booTmp <- any(checkVal)
  checkVal <- cbind(colVals[1], booTmp)
  #checkVal <- rbind(checkVal, c(colVals[2], booTmp))
  numericValue <- rbind(numericValue, checkVal)
  for(v in 1:dim(tmp)[1]){
    
    booVal1 <- is.letter(as.character(tmp[v,colVals[1]]))
    booVal2 <- is.letter(as.character(tmp[v,colVals[2]]))
    if(booVal1 == TRUE || booVal2 == TRUE){
          identVal <- identical(as.character(tmp[v,colVals[1]]), as.character(tmp[v,colVals[2]]))
    }  else{
                 identVal <- identical(round(as.numeric(as.character(tmp[v,colVals[1]])), digits=3), round(as.numeric(as.character(tmp[v,colVals[2]])), digits=3))
    }
    if(identVal == "FALSE"){
      print('foobar')
      outVal <- c(as.character(tmp[v,1]), as.character(tmp[v,2]), as.character(tmp[v,3]),names(tmp)[colVals][1], as.character(tmp[v,colVals[1]]), as.character(tmp[v,colVals[2]]))
      toFix <- rbind(toFix, outVal)
    }
    outMat[v,colVals] <- identVal
    print(outMat[v,colVals[1]])
    outCol[v] <- identVal
    
  }
  output <- cbind(output, tmp[,colVals],outCol)
  colnames(output)[dim(output)[2]] <- paste(grepPattern, '.matchVal', sep='')
}
colIndex <- grep('.matchVal', names(output))
output[,colIndex][output[,colIndex] == 1 ] <- 'MATCH'
output[,colIndex][output[,colIndex] == 0 ] <- 'NOMATCH'
write.table(output, file='~/Desktop/test.txt',  sep='\t', col.names=T, row.names=F, quote=F)


# Now merge the fields that agree on numeric values %100
# First find all of the numeric fields
numericFields <- which(numericValue[,2]==0)
mergeNames <- names(bar)[numericFields]
rmField <- grep("_txt", mergeNames)
rmField <- append(rmField, grep("_excl", mergeNames))
rmField <- append(rmField, grep("male_thc", mergeNames))
rmField <- append(rmField, grep("male_nc", mergeNames))
rmField <- append(rmField, grep("thc_sud", mergeNames))
rmField <- append(rmField, grep("thc_dsm_currr", mergeNames))
rmField <- append(rmField, grep("thc_dsm_lt", mergeNames))
rmField <- append(rmField, grep("thc_amph", mergeNames))

mergeNames <- mergeNames[-c(1,2,3,4, rmField)]
for(i in 1:length(mergeNames)){
  allMatch <- merge(foo1, bar1, by=c('tstname', 'authors', mergeNames[1:i]))
  print(mergeNames[1:i])
  print(dim(allMatch))
}
