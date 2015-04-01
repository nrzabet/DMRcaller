


# merge.regions - Get rid of small gaps between regions. Optionally, only glue
# together regions which have the same sign. Also merge together regions which
# overlap, regardless of whether they have the same sign.
.mergeRegions <- function(regions, minGap = 1, respectSigns = TRUE) {
  if (!respectSigns) {
    output <- reduce(regions, drop.empty.ranges=TRUE, min.gapwidth=minGap, 
                     ignore.strand=TRUE)
  } else {
    pos  <- which(regions$direction == 1)
    zero <- which(regions$direction == 0)
    neg  <- which(regions$direction == -1)
    
    preg <- reduce(regions[pos], drop.empty.ranges=TRUE, min.gapwidth=minGap, 
                   ignore.strand=TRUE)
    zreg <- reduce(regions[zero], drop.empty.ranges=TRUE, min.gapwidth=minGap, 
                   ignore.strand=TRUE)
    nreg <- reduce(regions[neg], drop.empty.ranges=TRUE, min.gapwidth=minGap, 
                   ignore.strand=TRUE)
    
    if(length(preg)>0){
      preg$direction <- 1 
    }
    if(length(zreg)>0){ 
      zreg$direction <- 0 
    }
    if(length(nreg)>0){
      nreg$direction <- -1
    }
    
    output <- c(preg,zreg,nreg)        
    if(!isDisjoint(output)){
      output <- disjoin(output)
      output$direction <- rep(0, times=length(output))
      output$direction[overlapsAny(output, regions[regions$direction == -1], 
                                   ignore.strand = TRUE)] <- -1
      output$direction[overlapsAny(output, regions[regions$direction == 1], 
                                   ignore.strand = TRUE)] <- 1
    }
    
    output <- output[order(start(output))]
    output$context <- rep(unique(regions$context), length.out=length(output))
  }
  return(output)
}


#' This takes a list of DMRs and attempts to merge DMRs while keeping the new 
#' DMRs statistically significant
.mergeDMRsIteratively <- function(DMRs, 
                                  minGap, 
                                  respectSigns = TRUE, 
                                  methylationData,
                                  minProportionDifference=0.4, 
                                  minReadsPerCytosine = 4, 
                                  pValueThreshold=0.01,
                                  test="fisher",
                                  alternative = "two.sided"){
  
  overlaps <- countOverlaps(DMRs, DMRs, maxgap = minGap, ignore.strand = TRUE)
  notToJoin <- DMRs[overlaps == 1]
  DMRs <- DMRs[overlaps > 1]
  
  joinedAny <- TRUE
  iteration <- 1
  while(joinedAny){
    joinedAny <- FALSE
    bufferDMRs <-GRanges()    
    index <- 1
    localIndex <- index
    joinedCount <- 0
    while(localIndex < length(DMRs)){
      localDMRs <- DMRs[index]
      canJoin <- TRUE
      while(canJoin){
        localIndex <- localIndex + 1
        newDMRs <-.joinDMRs(c(localDMRs, DMRs[localIndex]),
                            minGap = minGap, 
                            respectSigns = respectSigns, 
                            methylationData = methylationData,
                            minProportionDifference=minProportionDifference, 
                            minReadsPerCytosine = minReadsPerCytosine, 
                            pValueThreshold=pValueThreshold,
                            test=test,
                            alternative = alternative)
        if(length(newDMRs) == 1){
          joinedCount <- joinedCount + 1          
          localDMRs <- newDMRs
          joinedAny <- TRUE
          if(localIndex == length(DMRs)){
            canJoin <- FALSE
            bufferDMRs <- c(bufferDMRs, newDMRs)
          }
        } else{
          bufferDMRs <- c(bufferDMRs, localDMRs)
          canJoin <- FALSE
          index <- localIndex
          if(localIndex == length(DMRs)){
            bufferDMRs <- c(bufferDMRs, DMRs[localIndex])
          }
        }
      }
    }
    DMRs <- bufferDMRs
    iteration <- iteration + 1
  }
  
  DMRs <- c(DMRs, notToJoin)
  
  DMRs <- DMRs[order(DMRs)]
  
  return(DMRs)
}

# join if possible a set of DMRs
.joinDMRs <- function(DMRs,
                      minGap = minGap, 
                      respectSigns = TRUE, 
                      methylationData = methylationData,
                      minProportionDifference=0.4, 
                      minReadsPerCytosine = 4, 
                      pValueThreshold=0.01,
                      test="fisher",
                      alternative = "two.sided"){

  #are within the requeired distance
  if(length(reduce(DMRs, drop.empty.ranges=TRUE, min.gapwidth=minGap, 
                   ignore.strand=TRUE)) == 1){
    #all have the same direction
    if(!respectSigns | length(unique(DMRs$direction)) == 1){
      direction <- unique(DMRs$direction)
      if(length(unique(strand(DMRs))) == 1){
        localDMR <- DMRs[1]
        end(localDMR) <- max(end(DMRs))
        start(localDMR) <- min(start(DMRs))
        localDMR <- .analyseReadsInsideRegions(methylationData, localDMR)
        localDMR$pValue <-.computeaAjustedPValuesInDMRs(test, localDMR, 
                                                        alternative = alternative)
        if(abs(localDMR$proportion1 - localDMR$proportion2) >= minProportionDifference & 
             localDMR$pValue <= pValueThreshold &
             (localDMR$sumReadsN1 / localDMR$cytosinesCount) >= minReadsPerCytosine &
             (localDMR$sumReadsN2 / localDMR$cytosinesCount) >= minReadsPerCytosine){
          DMRs <- localDMR
        }
      }
    }
  }
  return(DMRs)
}

.getLongestDMRs <- function(DMRs,
                            minGap = minGap, 
                            respectSigns = TRUE, 
                            methylationData = methylationData,
                            minProportionDifference=0.4, 
                            minReadsPerCytosine = 4, 
                            pValueThreshold=0.01,
                            test="fisher",
                            alternative = "two.sided"){
  newDMRs <-.joinDMRs(DMRs,
                      minGap = minGap, 
                      respectSigns = respectSigns, 
                      methylationData = methylationData,
                      minProportionDifference=minProportionDifference, 
                      minReadsPerCytosine = minReadsPerCytosine, 
                      pValueThreshold=pValueThreshold,
                      test=test,
                      alternative = alternative)
  if(length(newDMRs) == 1){
    result <- newDMRs
  } else{
    result <- .mergeDMRsIteratively(DMRs, 
                                    minGap = minGap, 
                                    respectSigns = respectSigns, 
                                    methylationData = methylationData,
                                    minProportionDifference=minProportionDifference, 
                                    minReadsPerCytosine = minReadsPerCytosine, 
                                    pValueThreshold=pValueThreshold,
                                    test=test,
                                    alternative = alternative)
      
  }
  return(result)
}

#' This takes a list of DMRs and attempts to merge DMRs while keeping the new 
#' DMRs statistically significant
.smartMergeDMRs <- function(DMRs, 
                            minGap, 
                            respectSigns = TRUE, 
                            methylationData,
                            minProportionDifference=0.4, 
                            minReadsPerCytosine = 4, 
                            pValueThreshold=0.01,
                            test="fisher",
                            alternative = "two.sided",
                            cores = 1){  
  overlaps <- countOverlaps(DMRs, DMRs, maxgap = minGap, ignore.strand = TRUE)
  notToJoin <- DMRs[overlaps == 1]
  
  DMRs <- DMRs[overlaps > 1]
  

  if(length(DMRs) > 0){
    overlaps <- findOverlaps(DMRs, 
                             reduce(DMRs, min.gapwidth = minGap, 
                                    ignore.strand=TRUE), 
                             maxgap = minGap, ignore.strand = TRUE)
    DMRsList <- IRanges::splitAsList(DMRs[queryHits(overlaps)],  
                                     subjectHits(overlaps))
    
    if(cores > 1){
      bufferDMRs <- parallel::mclapply(1:length(DMRsList), function(i){ .getLongestDMRs(DMRsList[[i]],
                                                                                        minGap = minGap, 
                                                                                        respectSigns = respectSigns, 
                                                                                        methylationData = methylationData,
                                                                                        minProportionDifference=minProportionDifference, 
                                                                                        minReadsPerCytosine = minReadsPerCytosine, 
                                                                                        pValueThreshold=pValueThreshold,
                                                                                        test=test,
                                                                                        alternative = alternative)}, 
                                       mc.cores = cores)
      bufferDMRs <- unlist(GRangesList(bufferDMRs))
     } else{
       bufferDMRs <- GRanges()
       for(i in 1:length(DMRsList)){
         bufferDMRs <- c(bufferDMRs, .getLongestDMRs(DMRsList[[i]],
                                                     minGap = minGap, 
                                                     respectSigns = respectSigns, 
                                                     methylationData = methylationData,
                                                     minProportionDifference=minProportionDifference, 
                                                     minReadsPerCytosine = minReadsPerCytosine, 
                                                     pValueThreshold=pValueThreshold,
                                                     test=test,
                                                     alternative = alternative))
       }
    } 

    DMRs <- c(bufferDMRs, notToJoin)
  } else{
    DMRs <- notToJoin
  }
  
  DMRs <- DMRs[order(DMRs)]
  
  return(DMRs)
}








#returns the sum of values in a range
.getSumInRange <- function(vector, regions, start){
  sums <- rep(0, length(regions))
  for(i in 1:length(regions)){
    sums[i] <- sum(vector[(start(regions[i]) - start + 1):(end(regions[i]) - start + 1)])
  }
  return(sums)
}



#' Performs the analysis in all regions in a \code{\link{GRanges}} object
#'
#' @title Analyse reads inside regions
#' @param methylationData a \code{\link{GRanges}} object with five metadata 
#' columns see \code{\link{methylationDataList}}
#' @param regions a \code{\link{GRanges}} object with the identified regions
#' @return a \code{\link{GRanges}} object with eaual sized tiles of the regions. 
#' The object consists of the following metadata 
#' \describe{
#'  \item{sumReadsM1}{the number of methylated reads in condition 1}
#'  \item{sumReadsN1}{the total number of reads in condition 1}
#'  \item{proportion1}{the proportion of methylated reads in condition 1}
#'  \item{sumReadsM2}{the number of methylated reads in condition 2} 
#'  \item{sumReadsN2}{the total number of reads in condition 2} 
#'  \item{proportion2}{the proportion of methylated reads in condition 2} 
#'  \item{cytosinesCount}{the number of cytosines in the correct context} 
#' }  
#'       
#' @author Radu Zabet
.analyseReadsInsideRegions <- function(methylationData, regions){
  
  overlaps <- findOverlaps(methylationData, regions, ignore.strand = TRUE)
  methylationDataContextList <- IRanges::splitAsList(methylationData[queryHits(overlaps)],  subjectHits(overlaps))
  regionsIndexes <- as.integer(names(methylationDataContextList))
  
  regions$sumReadsM1 <- rep(0, times=length(regions))
  regions$sumReadsN1 <- rep(0, times=length(regions))    
  regions$proportion1 <- rep(0, times=length(regions))        
  regions$sumReadsM2 <- rep(0, times=length(regions))            
  regions$sumReadsN2 <- rep(0, times=length(regions))            
  regions$proportion2 <- rep(0, times=length(regions))        
  regions$cytosinesCount <- rep(0, times=length(regions))        
  
 
  if(length(regionsIndexes) > 0){  
    regions$sumReadsM1[regionsIndexes] <- sapply(methylationDataContextList,.sumReadsM1)
    regions$sumReadsN1[regionsIndexes] <- sapply(methylationDataContextList,.sumReadsN1)               
    regions$sumReadsM2[regionsIndexes] <- sapply(methylationDataContextList,.sumReadsM2)
    regions$sumReadsN2[regionsIndexes] <- sapply(methylationDataContextList,.sumReadsN2)
    regions$cytosinesCount[regionsIndexes] <- sapply(methylationDataContextList,length)
    
    valid <- regions$cytosinesCount[regionsIndexes] > 0
    regions$proportion1[regionsIndexes[valid]] <- regions$sumReadsM1[regionsIndexes[valid]]/regions$sumReadsN1[regionsIndexes[valid]]      
    regions$proportion2[regionsIndexes[valid]] <- regions$sumReadsM2[regionsIndexes[valid]]/regions$sumReadsN2[regionsIndexes[valid]]
  }
  return(regions)
}


.sumReadsM1 <- function(methylationData){
  return(sum(methylationData$readsM1))
}

.sumReadsN1 <- function(methylationData){
  return(sum(methylationData$readsN1))
}

.sumReadsM2 <- function(methylationData){
  return(sum(methylationData$readsM2))
}

.sumReadsN2 <- function(methylationData){
  return(sum(methylationData$readsN2))
}



#' Performs the analysis in all regions in a \code{\link{GRanges}} object
#'
#' @title Analyse reads inside regions
#' @param methylationData a \code{\link{GRanges}} object with five metadata 
#' columns see \code{\link{methylationDataList}}
#' @param regions a \code{\link{GRanges}} object with the identified regions
#' @return a \code{\link{GRanges}} object with eaual sized tiles of the regions. 
#' The object consists of the following metadata 
#' \describe{
#'  \item{sumReadsM}{the number of methylated reads in condition 1}
#'  \item{sumReadsN}{the total number of reads in condition 1}
#'  \item{proportion}{the proportion of methylated reads in condition 1}
#'  \item{cytosinesCount}{the number of cytosines in the correct context} 
#' }  
#'       
#' @author Radu Zabet
.analyseReadsInsideRegionsForCondition <- function(regions, methylationData, label="", context=""){
  
  overlaps <- findOverlaps(methylationData, regions, ignore.strand = TRUE)
  methylationDataContextList <- IRanges::splitAsList(methylationData[queryHits(overlaps)],  subjectHits(overlaps))
  regionsIndexes <- as.integer(names(methylationDataContextList))
  
  sumReadsM <- rep(0, times=length(regions))
  sumReadsN <- rep(0, times=length(regions))    
  proportion <- rep(0, times=length(regions))        
  cytosinesCount <- rep(0, times=length(regions))        
    
  
  sumReadsM[regionsIndexes] <- sapply(methylationDataContextList,.sumReadsM)
  sumReadsN[regionsIndexes] <- sapply(methylationDataContextList,.sumReadsN)               
  cytosinesCount[regionsIndexes] <- sapply(methylationDataContextList,length)
  
  proportion[regionsIndexes] <- sumReadsM[regionsIndexes]/sumReadsN[regionsIndexes]      
  
  
  done <- FALSE
  metadata <- DataFrame("sumReadsM" = sumReadsM, 
                         "sumReadsN" = sumReadsN, 
                         "proportion" = proportion, 
                         "cytosinesCount" = cytosinesCount)
  if(!is.null(label)){
    if(length(label) == 1){
      colnames(metadata) <- c(paste0("sumReadsM", label, context),
                              paste0("sumReadsN", label, context),
                              paste0("proportion", label, context),
                              paste0("cytosinesCount", context))
      done <- TRUE
    }
  }
  
  if(!done){
    colnames(metadata) <- c(paste0("sumReadsM", context),
                            paste0("sumReadsN", context),
                            paste0("proportion", context),
                            paste0("cytosinesCount", context))
  }
  values(regions) <- cbind(values(regions), metadata)
  return(regions)
}


.sumReadsM <- function(methylationData){
  return(sum(methylationData$readsM))
}

.sumReadsN <- function(methylationData){
  return(sum(methylationData$readsN))
}



#' This function counts the number of cytosines in each DMR
#'
#' @title Count cytosines inside
#' @param methylationData a \code{\link{GRanges}} object with five metadata 
#' columns; see \code{\link{methylationDataList}}
#' @param regions a \code{\link{GRanges}} with the computed regions
#' @return a \code{vector} with the number of cytosines in each DMR
#'  
#' @author Jonathan Michael Foonlan Tsang
.countCytosinesInside <- function(methylationData, regions) {
  return(countSubjectHits(findOverlaps(methylationData, regions, ignore.strand = TRUE)))
}


#' Performs the analysis in equal width regions of an \code{\link{GRanges}} 
#' object
#'
#' @title Analyse reads inside regions
#' @param methylationData a \code{\link{GRanges}} object with five metadata 
#' columns; see \code{\link{methylationDataList}}
#' @param regions a \code{\link{GRanges}} object with the identified regions
#' @return a \code{\link{GRanges}} object with eaual sized tiles of the regions. 
#' The object consists of the following metadata 
#' \describe{
#'  \item{sumReadsM1}{the number of methylated reads in condition 1}
#'  \item{sumReadsN1}{the total number of reads in condition 1}
#'  \item{proportion1}{the proportion of methylated reads in condition 1}
#'  \item{sumReadsM2}{the number of methylated reads in condition 2} 
#'  \item{sumReadsN2}{the total number of reads in condition 2} 
#'  \item{proportion2}{the proportion of methylated reads in condition 2} 
#'  \item{cytosinesCount}{the number of cytosines in the correct context} 
#' }  
#'       
#' @author Radu Zabet
.analyseReadsInsideBins <- function(methylationData, bins, currentRegion){
  
  binSize <- min(unique(width(bins)))
  #Rcpp
  readsM1 <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), methylationData$readsM1, windowSize = binSize)
  sumReadsM1 <- readsM1[seq(1,length(readsM1)-binSize, by=binSize)]
  
  readsN1 <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), methylationData$readsN1, windowSize = binSize)
  sumReadsN1 <- readsN1[seq(1,length(readsN1)-binSize, by=binSize)]
  
  proportion1 <- sumReadsM1/sumReadsN1
  
  readsM2 <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), methylationData$readsM2, windowSize = binSize)
  sumReadsM2 <- readsM2[seq(1,length(readsM2)-binSize, by=binSize)]
  
  readsN2 <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), methylationData$readsN2, windowSize = binSize)
  sumReadsN2 <- readsN2[seq(1,length(readsN2)-binSize, by=binSize)]
  
  proportion2 <- sumReadsM2/sumReadsN2
  
  
  cytosines <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), rep(1, length(start(methylationData))), windowSize = binSize)
  cytosinesCount <- cytosines[seq(1,length(cytosines)-binSize, by=binSize)]
  
  bins$sumReadsM1 <- sumReadsM1
  bins$sumReadsN1 <- sumReadsN1    
  bins$proportion1 <- proportion1        
  bins$sumReadsM2 <- sumReadsM2           
  bins$sumReadsN2 <- sumReadsN2            
  bins$proportion2 <- proportion2       
  bins$cytosinesCount <- cytosinesCount 

  return(bins)
}


#' This function splits a set of GRanges into a list of equally width GRanges 
#' objects
#'
#' @title .splitGRangesEqualy
#' @param regions a \code{\link{GRanges}} object 
#' @param breaks number of elements to break the GRanges
#' @return a \code{list} consisting of GRanges objects of same total width
#'  
#' @author Radu Zabet
.splitGRangesEqualy <- function(regions, breaks=1){
  result <- lapply(split(unlist(tile(regions, n=breaks)), rep(1:breaks, each=length(regions))), reduce, ignore.strand=TRUE)
  return(result)
}




#' Performs the analysis in equal width regions of an \code{\link{GRanges}} 
#' object
#'
#' @title Analyse reads inside regions for one sample
#' @param methylationData a \code{\link{GRanges}} object with five metadata 
#' columns; see \code{\link{methylationDataList}}.
#' @param bins a \code{\link{GRanges}} object with the bins
#' @return a \code{\link{GRanges}} object with eaual sized tiles of the regions. 
#' The object consists of the following metadata 
#' \describe{
#'  \item{sumReadsM}{the number of methylated reads}
#'  \item{sumReadsN}{the total number of reads}
#'  \item{Proportion}{the proportion of methylated reads}
#' }  
#'       
#' @author Radu Zabet
.analyseReadsInsideBinsOneSample <- function(methylationData, bins, currentRegion){
  
  binSize <- min(unique(width(bins)))
  #Rcpp
  readsM <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), methylationData$readsM, windowSize = binSize)
  sumReadsM <- readsM[seq(1,length(readsM)-binSize, by=binSize)]
  
  readsN <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), methylationData$readsN, windowSize = binSize)
  sumReadsN <- readsN[seq(1,length(readsN)-binSize, by=binSize)]
  
  Proportion <- sumReadsM/sumReadsN

  
  bins$sumReadsM <- sumReadsM
  bins$sumReadsN <- sumReadsN    
  bins$Proportion <- Proportion        
  
  return(bins)
}




#' Counts the number of cytosines in all the regions
#'
#' @title Count total number of Cytosines
#' @param methylationData a \code{\link{GRanges}} object with five metadata 
#' columns: see \code{\link{methylationDataList}}
#' @param regions a \code{\link{GRanges}} object with the identified regions
#' @return a the total number of Cytosines
#'       
#' @author Radu Zabet
.countTotalNumberOfCytosines <- function(methylationData, regions){
  
  total_cytocines <- 0
  for (index in 1:length(regions)) {
    currentRegion <- regions[index]
    total_cytocines <- total_cytocines + length(methylationData[queryHits(findOverlaps(methylationData, currentRegion, ignore.strand = TRUE))])
  }
  return(total_cytocines)
}


#' Performs the analysis in all regions in a \code{\link{GRanges}} object
#'
#' @title Analyse reads inside regions
#' @param methylationData a \code{\link{GRanges}} object with four metadata 
#' columns; see \code{\link{loadMethylationDataList}}.
#' @param regions a \code{\link{GRanges}} object with the identified regions
#' @return a \code{\link{GRanges}} object with eaual sized tiles of the regions. 
#' The object consists of the following metadata.
#' \describe{
#'  \item{sumReadsM}{the number of methylated reads in condition 1}
#'  \item{sumReadsN}{the total number of reads in condition 1}
#'  \item{Proportion}{the proportion of methylated reads in condition 1}
#' }  
#'       
#' @author Radu Zabet
.analyseReadsInsideRegionsOneSample <- function(methylationData, regions){
  overlaps <- findOverlaps(methylationData, regions, ignore.strand = TRUE)
  methylationDataContextList <- IRanges::splitAsList(methylationData[queryHits(overlaps)],  subjectHits(overlaps))
  regionsIndexes <- as.integer(names(methylationDataContextList))
  
  regions$sumReadsM <- rep(0, times=length(regions))
  regions$sumReadsN <- rep(0, times=length(regions))    
  regions$Proportion <- rep(0, times=length(regions))        

  
  regions$sumReadsM[regionsIndexes] <- sapply(methylationDataContextList,.sumReadsM)
  regions$sumReadsN[regionsIndexes] <- sapply(methylationDataContextList,.sumReadsN)  
  
  regions$Proportion[regionsIndexes] <- regions$sumReadsM[regionsIndexes]/regions$sumReadsN[regionsIndexes]      
  return(regions)
}

.sumReadsM <- function(methylationData){
  return(sum(methylationData$readsM))
}

.sumReadsN <- function(methylationData){
  return(sum(methylationData$readsN))
}


