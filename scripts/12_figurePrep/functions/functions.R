
## Now create a function to mean left and right hemisphere values
## This should only be used for the volume metrics and is a lazy solution
averageLeftAndRightVol <- function(dataFrame){
    # Now get the right data
    dataFrame.right <- dataFrame[, grep('_R_', names(dataFrame))]
    dataFrame.tmp <- dataFrame[, -grep('_R_', names(dataFrame))]
    if(!identical(integer(0),grep('_right', names(dataFrame.tmp)))){
        dataFrame.right <- cbind(dataFrame.right, dataFrame.tmp[, grep('_right', names(dataFrame.tmp))])
        dataFrame.tmp <- dataFrame.tmp[,-grep('_right', names(dataFrame.tmp))]
    }
    if(!identical(integer(0),grep('_rh_', names(dataFrame.tmp)))){
        dataFrame.right <- cbind(dataFrame.right, dataFrame.tmp[, grep('_rh_', names(dataFrame.tmp))])
        dataFrame.tmp <- dataFrame.tmp[,-grep('_rh_', names(dataFrame.tmp))]
    }
    if(dim(dataFrame.tmp)[2] == 0){
        dataFrame.tmp <- dataFrame
        dataFrame.right <- dataFrame.tmp[, grep('_rh_', names(dataFrame.tmp))]
        dataFrame.tmp <- dataFrame.tmp[,-grep('_rh_', names(dataFrame.tmp))]
    }
    # First do the left data
    dataFrame.left <- dataFrame.tmp[, grep('_L_', names(dataFrame.tmp))]
    dataFrame.tmp <- dataFrame.tmp[, -grep('_L_', names(dataFrame.tmp))]
    if(!identical(integer(0),grep('_left', names(dataFrame.tmp)))){
        dataFrame.left <- cbind(dataFrame.left, dataFrame.tmp[, grep('_left', names(dataFrame.tmp))])
        dataFrame.tmp <- dataFrame.tmp[,-grep('_left', names(dataFrame.tmp))]
    }
    if(!identical(integer(0),grep('_lh_', names(dataFrame.tmp)))){
        dataFrame.left <- cbind(dataFrame.left, dataFrame.tmp[, grep('_lh_', names(dataFrame.tmp))])
        dataFrame.tmp <- dataFrame.tmp[,-grep('_lh_', names(dataFrame.tmp))]
    }
    if(dim(dataFrame.tmp)[2] == 0){
        dataFrame.tmp <- dataFrame
        dataFrame.left <- dataFrame.tmp[, grep('_lh', names(dataFrame.tmp))]
        dataFrame.tmp <- dataFrame.tmp[,-grep('_lh_', names(dataFrame.tmp))]
        dataFrame.tmp <- dataFrame.tmp[,-grep('_rh_', names(dataFrame.tmp))]
    }
    # Now combine the data frames
    dataFrame.meaned <- (dataFrame.left + dataFrame.right)
    
    # Now remove the left and right indeices from the names of the meaned data frame
    colnames(dataFrame.meaned) <- gsub(pattern='_L_', replacement = '_', x = colnames(dataFrame.meaned), fixed = TRUE)
    colnames(dataFrame.meaned) <- gsub(pattern='_R_', replacement = '_', x = colnames(dataFrame.meaned), fixed = TRUE)
    colnames(dataFrame.meaned) <- gsub(pattern='_left', replacement = '', x = colnames(dataFrame.meaned), fixed = TRUE)
    colnames(dataFrame.meaned) <- gsub(pattern='_right', replacement = '', x = colnames(dataFrame.meaned), fixed = TRUE)
    colnames(dataFrame.meaned) <- gsub(pattern='_rh_', replacement = '_', x = colnames(dataFrame.meaned), fixed = TRUE)
    colnames(dataFrame.meaned) <- gsub(pattern='_lh_', replacement = '_', x = colnames(dataFrame.meaned), fixed = TRUE)
    # Now rm the left and right values and append the meaned values
    indexToRm <- grep('_L_', names(dataFrame))
    indexToRm <- append(indexToRm, grep('_R_', names(dataFrame)))
    indexToRm <- append(indexToRm, grep('_left', names(dataFrame)))
    indexToRm <- append(indexToRm, grep('_right', names(dataFrame)))
    indexToRm <- append(indexToRm, grep('_rh_', names(dataFrame)))
    indexToRm <- append(indexToRm, grep('_lh_', names(dataFrame)))
    # Now prep the output
    output <- dataFrame[,-indexToRm]
    # Now combine our average values
    output <- cbind(output, dataFrame.meaned)
    # Now return the output
    return(output)
}



findLobe <- function(grepPattern){
    # Declare the rois that we will grep through
    rois<-c("Thal","Putamen","Caudate","Pallidum",  # Basal Ganglia
    "Accumbens", "BasFor", # Basal Ganglia
    "PHG","Hip","PIns","SCA","AIns", # Limbic
    "ACgG","PCgG","Ent","Amygdala","MCgG", # Limbic
    "FO","MFC","MOrG","POrG","OrIFG","TrIFG","AOrG","OpIFG","GRe", # Frontal Orbital
    "FRP", "LOrG", # Frontal Orbital
    "PrG","MSFG","SMC","MFG","SFG", # Frontal Dorsal
    "FuG","PT","PP","ITG","CO","MTG","TMP","STG","TTG", # Temporal
    "PCu","PoG","AnG","PO","SPL","MPrG", # Parietal
    "SMG","MPoG", # Parietal
    "IOG","Cun","LiG","OFuG","MOG","Calc","OCP","SOG", # Occiptal
    "Cerebellum_Exterior", "Cerebellar_Vermal_Lobules_I.V", "Cerebellar_Vermal_Lobules_VI.VII", "Cerebellar_Vermal_Lobules_VIII.X", # Cerebellum
    "Limbic_Lobe_WM", "Insular_Lobe_WM", "Frontal_Lobe_WM", # WM
    "Parietal_Lobe_WM", "Occipital_Lobe_WM", "Temporal_Lobe_WM")  # WM
    
    
    # Declare the index key which corresponds to lobe values
    index.key <- c(1,7,17,28,33,42,50,58,62,68)
    
    # Now find where the pattern matches to the roi list
    for(pattern.match.variable in 1:length(rois)){
        nuclei.to.grep <- rois[pattern.match.variable]
        grep.output <- grep(nuclei.to.grep ,grepPattern)
        if(!identical(integer(0), grep.output)){
            break
        }
        if(pattern.match.variable==67){
            pattern.match.variable <- 68
        }
    }
    OFuGCheck <- grep("OFuG", grepPattern)
    if(!identical(integer(0), OFuGCheck)){
        pattern.match.variable <- 53
    }
    lobe.group <- findInterval(pattern.match.variable, index.key)
    return(lobe.group)
}
