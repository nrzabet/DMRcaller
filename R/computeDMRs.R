#' This function computes the differentially methylated regions between two 
#' conditions.  
#'
#' @title Compute DMRs
#' @param methylationData1 the methylation data in condition 1
#' (see \code{\link{methylationDataList}}).
#' @param methylationData2 the methylation data in condition 2
#' (see \code{\link{methylationDataList}}).
#' @param regions a \code{\link{GRanges}} object with the regions where to 
#' compute the DMRs. If \code{NULL}, the DMRs are computed genome-wide.
#' @param context the context in which the DMRs are computed (\code{"CG"}, 
#' \code{"CHG"} or \code{"CHH"}).
#' @param method the method used to compute the DMRs (\code{"noise_filter"}, 
#' \code{"neighbourhood"} or \code{"bins"}). The \code{"noise_filter"} method 
#' uses a triangular kernel to smooth the number of reads and then performs a 
#' statistical test to determine which regions dispay different levels of 
#' methylation in the two conditions. The \code{"neighbourhood"} method 
#' computates differentially methylated cytosines. Finally, the \code{"bins"} 
#' method partiones the genome into equal sized tilling bins and performs the 
#' statistical test between the two conditions in each bin. For all three 
#' methods, the cytosines or bins are then merged into DMRs without affecting 
#' the inital parameters used when calling the differentiall methylated 
#' cytosines/bins (p-value, difference in methylation levels, minimum number of 
#' reads per cytosine).
#' @param windowSize the size of the triangle base measured in nucleotides. 
#' This parameter is required only if the selected method is 
#' \code{"noise_filter"}. 
#' @param kernelFunction a \code{character} indicating which kernel function to 
#' be used. Can be one of \code{"uniform"}, \code{"triangular"}, 
#' \code{"gaussian"} or \code{"epanechnicov"}. This is required only if the 
#' selected method is \code{"noise_filter"}. 
#' @param lambda numeric value required for the Gaussian filter 
#' (\code{K(x) = exp(-lambda*x^2)}). This is required only if the selected 
#' method is \code{"noise_filter"} and the selected kernel function is 
#' \code{"gaussian"}. 
#' @param binSize the size of the tiling bins in nucleotides. This parameter is 
#' required only if the selected method is \code{"bins"}.
#' @param test the statistical test used to call DMRs (\code{"fisher"} for 
#' Fisher's exact test or \code{"score"} for Score test). 
#' @param pValueThreshold DMRs with p-values (when performing the statistical 
#' test; see \code{test}) higher or equal than \code{pValueThreshold} are 
#' discarded. Note that we adjust the p-values using the Benjamini and 
#' Hochberg's method to control the false discovery rate.
#' @param minCytosinesCount DMRs with less cytosines in the specified context 
#' than \code{minCytosinesCount} will be discarded.
#' @param minProportionDifference DMRs where the difference in methylation 
#' proportion between the two conditions is lower than 
#' \code{minProportionDifference} are discarded.
#' @param minGap DMRs separated by a gap of at least \code{minGap} are not 
#' merged. Note that only DMRs where the change in methylation is in the same 
#' direction are joined.
#' @param minSize DMRs with a size smaller than \code{minSize} are discarded.
#' @param minReadsPerCytosine  DMRs with the average number of reads lower than 
#' \code{minReadsPerCytosine} are discarded. 
#' @param cores the number of cores used to compute the DMRs. 
#' @return the DMRs stored as a \code{\link{GRanges}} object with the following 
#' metadata columns:
#' \describe{
#'  \item{direction}{a number indicating whether the region lost (-1)  or gain 
#'  (+1) methylation in condition 2 compared to condition 1.}
#'  \item{context}{the context in which the DMRs was computed (\code{"CG"}, 
#'  \code{"CHG"} or \code{"CHH"}).}
#'  \item{sumReadsM1}{the number of methylated reads in condition 1.}
#'  \item{sumReadsN1}{the total number of reads in condition 1.} 
#'  \item{proportion1}{the proportion methylated reads in condition 1.} 
#'  \item{sumReadsM2}{the number of methylated reads in condition 2.} 
#'  \item{sumReadsN2}{the total number reads in condition 2.}
#'  \item{proportion2}{the proportion methylated reads in condition 2.} 
#'  \item{cytosinesCount}{the number of cytosines in the DMR.} 
#'  \item{regionType}{a string indicating whether the region lost (\code{"loss"})  
#'  or gained (\code{"gain"}) methylation in condition 2 compared to condition 1.} 
#'  \item{pValue}{the p-value (adjusted to control the false discovery rate with 
#'  the Benjamini and Hochberg's method) of the statistical test when the DMR was 
#'  called.} 
#' }
#' @seealso \code{\link{filterDMRs}}, \code{\link{mergeDMRsIteratively}}, 
#' \code{\link{analyseReadsInsideRegionsForCondition}} and 
#' \code{\link{DMRsNoiseFilterCG}}
#' @examples
#' 
#' # load the methylation data
#' data(methylationDataList)
#' 
#' # the regions where to compute the DMRs
#' regions <- GRanges(seqnames = Rle("Chr3"), ranges = IRanges(1,1E5))
#' 
#' # compute the DMRs in CG context with noise_filter method
#' DMRsNoiseFilterCG <- computeDMRs(methylationDataList[["WT"]], 
#'                      methylationDataList[["met1-3"]], regions = regions, 
#'                      context = "CG", method = "noise_filter", 
#'                      windowSize = 100, kernelFunction = "triangular",  
#'                      test = "score", pValueThreshold = 0.01, 
#'                      minCytosinesCount = 4, minProportionDifference = 0.4, 
#'                      minGap = 200, minSize = 50, minReadsPerCytosine = 4, 
#'                      cores = 1)
#' 
#' \dontrun{
#' # compute the DMRs in CG context with neighbourhood method
#' DMRsNeighbourhoodCG <- computeDMRs(methylationDataList[["WT"]], 
#'                        methylationDataList[["met1-3"]], regions = regions, 
#'                        context = "CG", method = "neighbourhood", 
#'                        test = "score", pValueThreshold = 0.01, 
#'                        minCytosinesCount = 4, minProportionDifference = 0.4, 
#'                        minGap = 200, minSize = 50, minReadsPerCytosine = 4, 
#'                        cores = 1)
#' 
#' # compute the DMRs in CG context with bins method
#' DMRsBinsCG <- computeDMRs(methylationDataList[["WT"]], 
#'                methylationDataList[["met1-3"]], regions = regions, 
#'                context = "CG", method = "bins", binSize = 100, 
#'                test = "score", pValueThreshold = 0.01, minCytosinesCount = 4, 
#'                minProportionDifference = 0.4, minGap = 200, minSize = 50, 
#'                minReadsPerCytosine = 4, cores = 1)
#' 
#' }
#' @author Nicolae Radu Zabet and Jonathan Michael Foonlan Tsang
#' @export
computeDMRs <- function(methylationData1, 
                        methylationData2, 
                        regions = NULL, 
                        context = "CG", 
                        method="noise_filter",
                        windowSize = 100,
                        kernelFunction = "triangular", 
                        lambda = 0.5,
                        binSize = 100,
                        test = "fisher", 
                        pValueThreshold = 0.01, 
                        minCytosinesCount = 4, 
                        minProportionDifference = 0.4,
                        minGap = 200, 
                        minSize = 50, 
                        minReadsPerCytosine = 4, 
                        cores = 1) {
  ##Parameters checking
  cat("Parameters checking ...\n")
  
  .validateMethylationData(methylationData1, variableName="methylationData1")
  .validateMethylationData(methylationData2, variableName="methylationData2")
  
  
  regions <- union(.validateGRanges(regions, methylationData1), .validateGRanges(regions, methylationData2))
  
  .validateContext(context)
  
  .stopIfNotAll(c(!is.null(method), 
                  all(is.character(method)),
                  length(method) == 1,
                  all(method %in% c("noise_filter","neighbourhood","bins"))),
                " method can be only noise_filter, neighbourhood or bins")  
  
  if(method == "noise_filter"){
    .stopIfNotAll(c(.isInteger(windowSize, positive=TRUE)), 
                  " the window size used by the interpolation method is an integer higher than 0")
    
    .stopIfNotAll(c(!is.null(kernelFunction), 
                    kernelFunction%in%c("uniform", "triangular", "gaussian", "epanechnicov")), 
                  paste("Unknown kernel function: ", kernelFunction, ". 
                        kernelFunction should be one of \"uniform\", \"triangular\", 
                        \"gaussian\", \"epanechnicov\"",sep=""))
    
    if(kernelFunction == "gaussian"){
      .stopIfNotAll(c(!is.null(lambda),
                      is.numeric(lambda)), 
                    " lambda needs to be a numeric value.")
    }
    
  }
  
  if(method == "bins"){
    .stopIfNotAll(c(.isInteger(binSize, positive=TRUE)), 
                  " the bin size used by the method is an integer higher than 0")
    
  }
  
  .validateStatisticalTest(test)
  
  .stopIfNotAll(c(!is.null(pValueThreshold), 
                  is.numeric(pValueThreshold), 
                  pValueThreshold > 0, 
                  pValueThreshold < 1),
                " the p-value threshold needs to be in the interval (0,1)")
  
  .stopIfNotAll(c(.isInteger(minCytosinesCount, positive=TRUE)), 
                " the minCytosinesCount is an integer higher or equal to 0")
  
  .stopIfNotAll(c(!is.null(minProportionDifference), 
                  is.numeric(minProportionDifference), 
                  minProportionDifference > 0, 
                  minProportionDifference < 1), 
                " the minimum difference in methylation needs to be in the interval (0,1)")
  
  .stopIfNotAll(c(.isInteger(minGap, positive=TRUE)),
                " the minimum gap between DMRs is an integer higher or equal to 0")
  
  .stopIfNotAll(c(.isInteger(minSize, positive=TRUE)),
                " the minimum size of a DMR is an integer higher or equal to 0")
  
  .stopIfNotAll(c(.isInteger(minReadsPerCytosine, positive=TRUE)), 
                " the minimum average number of reads in a DMR is an integer higher or equal to 0")
  
  .stopIfNotAll(c(.isInteger(cores, positive=TRUE)), 
                " the number of cores to used when computing the DMRs needs to be an integer  hirger or equal to 1.")
  
  computedDMRs <- GRanges()
  
  if(method == "noise_filter"){
    computedDMRs <- .computeDMRsNoiseFilter(methylationData1 = methylationData1, 
                                            methylationData2 = methylationData2, 
                                            regions = regions, 
                                            context = context,
                                            windowSize = windowSize,
                                            kernelFunction = kernelFunction, 
                                            lambda = lambda,
                                            test = test, 
                                            pValueThreshold = pValueThreshold, 
                                            minCytosinesCount = minCytosinesCount, 
                                            minProportionDifference =minProportionDifference, 
                                            minGap = minGap, 
                                            minSize = minSize, 
                                            minReadsPerCytosine = minReadsPerCytosine, 
                                            cores = cores)
  } else if(method == "neighbourhood"){
    computedDMRs <- .computeDMRsNeighbourhood(methylationData1 = methylationData1, 
                                              methylationData2 = methylationData2, 
                                              regions = regions, 
                                              context = context, 
                                              test = test, 
                                              pValueThreshold = pValueThreshold, 
                                              minCytosinesCount = minCytosinesCount, 
                                              minProportionDifference =minProportionDifference, 
                                              minGap = minGap, 
                                              minSize = minSize, 
                                              minReadsPerCytosine = minReadsPerCytosine, 
                                              cores = cores)
  } else if(method == "bins"){
    computedDMRs <- .computeDMRsBins(methylationData1 = methylationData1, 
                                     methylationData2 = methylationData2, 
                                     regions = regions, 
                                     context = context, 
                                     binSize = binSize,
                                     test = test, 
                                     pValueThreshold = pValueThreshold, 
                                     minCytosinesCount = minCytosinesCount, 
                                     minProportionDifference =minProportionDifference, 
                                     minGap = minGap, 
                                     minSize = minSize, 
                                     minReadsPerCytosine = minReadsPerCytosine, 
                                     cores = cores)
    
  } else{
    cat("Unknown method: ",method," \n")
  } 
  
  return(computedDMRs)  
}

#' This function computes the differentially methylated regions between two conditions
#' using the noise filter method.  
.computeDMRsNoiseFilter <- function(methylationData1, 
                                    methylationData2, 
                                    regions = NULL, 
                                    context = "CG",
                                    windowSize = 100,
                                    kernelFunction="triangular", 
                                    lambda=0.5,                                    
                                    test = "fisher", 
                                    pValueThreshold = 0.01, 
                                    minCytosinesCount = 4, 
                                    minProportionDifference = 0.4, 
                                    minGap = 200, 
                                    minSize = 50, 
                                    minReadsPerCytosine = 4, 
                                    cores = 1) {
  
  regions <- reduce(regions)
  
  # extract the methylation data in the correct context
  cat("Extract methylation in the corresponding context \n")
  
  
  contextMethylationData1 <- methylationData1[methylationData1$context%in%context]
  contextMethylationData2 <- methylationData2[methylationData2$context%in%context]
  
  
  localContextMethylationData1 <- contextMethylationData1[queryHits(findOverlaps(contextMethylationData1, regions))]
  localContextMethylationData2 <- contextMethylationData2[queryHits(findOverlaps(contextMethylationData2, regions))]
  
  
  localContextMethylationData <- .joinMethylationData(localContextMethylationData1, localContextMethylationData2)
  
  
  regionsList <- .splitGRangesEqualy(regions, cores)
  
  # inner loop function for parallel::mclapply
  .computeDMRsInterpolationLoop = function(i){
    computedDMRs <- GRanges()    
    for(index in 1:length(regionsList[[i]])){
      currentRegion <- regionsList[[i]][index]
      
      cat("Computing DMRs at ",.printGenomicRanges(currentRegion),"\n")
      
      # Select the points in methylationData that we're interested in. These are the 
      # points that lie within 'regions', as well as any that lie within 
      # window.size of them. 
      windowSizeHalf <- floor((windowSize - 1)/2)
      extendedRegion <- currentRegion
      start(extendedRegion) <- start(currentRegion) - windowSizeHalf
      end(extendedRegion) <- end(currentRegion) + windowSizeHalf
      
      overlaps <- findOverlaps(localContextMethylationData, extendedRegion)
      if(length(overlaps) > 0){
        localMethylationData <- localContextMethylationData[queryHits(overlaps)]
        
        cat("Calculating interpolations...\n")
        
        #Rcpp
        movingAverageMethylReads1 <- round(.movingAverage(start(currentRegion), 
                                                          end(currentRegion), 
                                                          start(localMethylationData), 
                                                          localMethylationData$readsM1, 
                                                          windowSizeHalf = windowSizeHalf))
        movingAverageTotalReads1 <- round(.movingAverage(start(currentRegion), 
                                                         end(currentRegion), 
                                                         start(localMethylationData), 
                                                         localMethylationData$readsN1, 
                                                         windowSizeHalf = windowSizeHalf))
        movingAverageProportion1 <- movingAverageMethylReads1 / movingAverageTotalReads1
        
        #Rcpp
        movingAverageMethylReads2 <- round(.movingAverage(start(currentRegion), 
                                                          end(currentRegion), 
                                                          start(localMethylationData), 
                                                          localMethylationData$readsM2, 
                                                          windowSizeHalf = windowSizeHalf))
        movingAverageTotalReads2 <- round(.movingAverage(start(currentRegion), 
                                                         end(currentRegion), 
                                                         start(localMethylationData), 
                                                         localMethylationData$readsN2, 
                                                         windowSizeHalf = windowSizeHalf))
        movingAverageProportion2 <- movingAverageMethylReads2 / movingAverageTotalReads2
        
        cat("Identifying DMRs...\n")    
        pValue <- .computeAdjuestedPValues(test, 
                                           movingAverageMethylReads1, 
                                           movingAverageTotalReads1, 
                                           movingAverageMethylReads2, 
                                           movingAverageTotalReads2, 
                                           alternative = "two.sided")
        
        
        # compute the differentially methylated cytosines
        DMPs <- rep(0, times=width(currentRegion))
        DMPs[is.na(pValue)] <- -2
        bufferIndex <- !is.na(pValue) & 
          pValue < pValueThreshold & 
          abs(movingAverageProportion1 - movingAverageProportion2) >= minProportionDifference & 
          movingAverageTotalReads1 >=minReadsPerCytosine &
          movingAverageTotalReads2 >=minReadsPerCytosine
        DMPs[bufferIndex] <- sign(movingAverageProportion2[bufferIndex] - movingAverageProportion1[bufferIndex])
        
        #join the differentially methylated cytosines into regions
        rle <- rle(DMPs)
        rle$cumulative <- cumsum(rle$lengths)
        endOfRuns <- rle$cumulative + start(currentRegion) - 1
        
        DMRs <- GRanges(
          seqnames    = seqnames(currentRegion),
          ranges      = IRanges(endOfRuns - rle$lengths + 1, endOfRuns),
          strand      = strand(currentRegion),
          direction   = rle$values,
          context     = paste(context, collapse = "_")
        )
        
        
        
        DMRs$direction[DMRs$direction == -2] <- NA
        
        # Select the crude list of DMRs
        DMRs <- DMRs[!is.na(DMRs$direction) & (DMRs$direction == 1 | DMRs$direction == -1)]
        
        
        # append current DMRs to the global list of DMRs 
        if(length(computedDMRs) == 0){
          computedDMRs <- DMRs
        } else{
          computedDMRs <- c(computedDMRs,DMRs)
        }
        
      }
    }
    
    return(computedDMRs)
    
  }
  
  # compute the DMRs
  if(cores > 1){
    cat("Compute the DMRs using ", cores, "cores\n")
    computedDMRs <- parallel::mclapply(1:length(regionsList), .computeDMRsInterpolationLoop, mc.cores = cores)
  } else {
    computedDMRs <- lapply(1:length(regionsList), .computeDMRsInterpolationLoop)
  }
  
  
  computedDMRs <- unlist(GRangesList(computedDMRs))
  
  
  if(length(computedDMRs) > 0){
    computedDMRs <- computedDMRs[order(computedDMRs)]
    cat("Analysed reads inside DMRs\n")
    overlaps <- countOverlaps(localContextMethylationData, computedDMRs)
    localContextMethylationDataDMRs <- localContextMethylationData[overlaps > 0]
    if(cores > 1){
      computedDMRsList <- IRanges::splitAsList(computedDMRs,  rep(1:cores, length.out=length(computedDMRs)))
      bufferComputedDMRsList <- parallel::mclapply(1:length(computedDMRsList), function(i){ 
        .analyseReadsInsideRegions(localContextMethylationDataDMRs, computedDMRsList[[i]])}, 
        mc.cores = cores)
      computedDMRs <- unlist(GRangesList(bufferComputedDMRsList))
      computedDMRs <- computedDMRs[order(computedDMRs)]
      
    } else{
      computedDMRs <- .analyseReadsInsideRegions(localContextMethylationDataDMRs, computedDMRs)
    }
    
    computedDMRs$pValue <- .computeaAjustedPValuesInDMRs(test, computedDMRs ,alternative = "two.sided")
    
    
    computedDMRs <- computedDMRs[computedDMRs$pValue < pValueThreshold & 
                                   abs(computedDMRs$proportion1 - computedDMRs$proportion2) >= minProportionDifference & 
                                   computedDMRs$sumReadsN1/computedDMRs$cytosinesCount >=minReadsPerCytosine &
                                   computedDMRs$sumReadsN2/computedDMRs$cytosinesCount >=minReadsPerCytosine]
    
    cat("Merge DMRs iteratively\n")    
    # Get rid of small gaps between DMRs.
    #computedDMRs <- .mergeRegions(computedDMRs, minGap = minGap, respectSigns = TRUE)
    if(minGap > 0){
      computedDMRs <- .smartMergeDMRs(computedDMRs, 
                                      minGap = minGap, 
                                      respectSigns = TRUE, 
                                      methylationData = localContextMethylationData,
                                      minProportionDifference=minProportionDifference, 
                                      minReadsPerCytosine = minReadsPerCytosine, 
                                      pValueThreshold=pValueThreshold,
                                      test=test, 
                                      alternative = "two.sided",
                                      cores = cores)
    }  
    
    cat("Filter DMRs \n")    
    if(length(computedDMRs) > 0){
      #remove small DMRs 
      computedDMRs <- computedDMRs[width(computedDMRs) >= minSize]
      if(length(computedDMRs) > 0){
        #remove DMRswith few cytosines
        computedDMRs <- computedDMRs[computedDMRs$cytosinesCount >= minCytosinesCount]
        #recompute the adjusted p-values
        if(length(computedDMRs) > 0){
          computedDMRs$pValue <- .computeaAjustedPValuesInDMRs(test, computedDMRs ,alternative = "two.sided")
          computedDMRs$regionType <- rep("loss", length(computedDMRs))
          computedDMRs$regionType[which(computedDMRs$proportion1 < computedDMRs$proportion2)] <- "gain"
          
        }
      }
      
    }
  }  
  
  return(computedDMRs)  
}

#' This function computes the differentially methylated regions between two conditions
#' using the neighbourhood method. This assumes the computation of differentially methylated 
#' cytosines followed by smart merging of these cytosines while keeping the new DMRs 
#' statistically significant.   
.computeDMRsNeighbourhood <- function(methylationData1, 
                                      methylationData2, 
                                      regions = NULL, 
                                      context = "CG",
                                      test = "fisher", 
                                      pValueThreshold = 0.01, 
                                      minCytosinesCount = 4, 
                                      minProportionDifference = 0.4,
                                      minGap = 200, 
                                      minSize = 50, 
                                      minReadsPerCytosine = 4, 
                                      cores = 1) {  
  
  regions <- reduce(regions)
  
  # extract the methylation data in the correct context
  cat("Extract methylation in the corresponding context \n")
  
  
  contextMethylationData1 <- methylationData1[methylationData1$context%in%context]
  contextMethylationData2 <- methylationData2[methylationData2$context%in%context]
  
  
  localContextMethylationData1 <- contextMethylationData1[queryHits(findOverlaps(contextMethylationData1, regions))]
  localContextMethylationData2 <- contextMethylationData2[queryHits(findOverlaps(contextMethylationData2, regions))]
  
  
  localContextMethylationData <- .joinMethylationData(localContextMethylationData1, localContextMethylationData2)
  
  
  cat("Computing DMRs \n")
  DMPs <- GRanges()
  if(length(localContextMethylationData) > 0){
    DMPs <- localContextMethylationData
    DMPs$pValue <- .computeAdjuestedPValues(test, 
                                            DMPs$readsM1, 
                                            DMPs$readsN1, 
                                            DMPs$readsM2, 
                                            DMPs$readsN2, 
                                            alternative = "two.sided")
    DMPs <- DMPs[!is.na(DMPs$pValue)]
    DMPs$sumReadsM1 <- DMPs$readsM1
    DMPs$sumReadsN1 <- DMPs$readsN1
    DMPs$proportion1 <- DMPs$readsM1 / DMPs$readsN1
    DMPs$sumReadsM2 <- DMPs$readsM2
    DMPs$sumReadsN2 <- DMPs$readsN2
    DMPs$proportion2 <- DMPs$readsM2 / DMPs$readsN2
    DMPs$cytosinesCount <- 1
    DMPs$direction <- sign(DMPs$proportion2 - DMPs$proportion1)
    
    bufferIndex <- DMPs$pValue < pValueThreshold & 
      abs(DMPs$proportion2 - DMPs$proportion1) >= minProportionDifference & 
      DMPs$sumReadsN1 >=minReadsPerCytosine &
      DMPs$sumReadsN2 >=minReadsPerCytosine
    DMPs <- DMPs[bufferIndex]
    strand(DMPs) <- "*" 
  }
  
  computedDMRs <- GRanges()
  if(length(DMPs) > 0){    
    cat("Merge DMRs iteratively\n")    
    # Get rid of small gaps between DMRs.
    if(minGap > 0){
      computedDMRs <- .smartMergeDMRs(DMPs, 
                                      minGap = minGap, 
                                      respectSigns = TRUE, 
                                      methylationData = localContextMethylationData,
                                      minProportionDifference=minProportionDifference, 
                                      minReadsPerCytosine = minReadsPerCytosine, 
                                      pValueThreshold=pValueThreshold,
                                      test=test, 
                                      alternative = "two.sided",
                                      cores = cores)
    } else{
      computedDMRs <- DMPs
    } 
    
    cat("Filter DMRs \n")    
    if(length(computedDMRs) > 0){
      #remove small DMRs 
      computedDMRs <- computedDMRs[width(computedDMRs) >= minSize]
      if(length(computedDMRs) > 0){
        #remove DMRs with few cytosines
        computedDMRs <- computedDMRs[computedDMRs$cytosinesCount >= minCytosinesCount]
        #recompute the adjusted p-values
        if(length(computedDMRs) > 0){
          computedDMRs$pValue <- .computeaAjustedPValuesInDMRs(test, computedDMRs ,alternative = "two.sided")
          computedDMRs$regionType <- rep("loss", length(computedDMRs))
          computedDMRs$regionType[which(computedDMRs$proportion1 < computedDMRs$proportion2)] <- "gain"
          
        }
      }
      
    }
  }  
  
  return(computedDMRs)  
}


#' This function computes the differentially methylated regions between two conditions
#' using the bins method.  
.computeDMRsBins <- function(methylationData1, 
                             methylationData2, 
                             regions = NULL, 
                             context = "CG", 
                             binSize = 100,
                             test = "fisher", 
                             pValueThreshold = 0.01, 
                             minCytosinesCount = 4,
                             minProportionDifference = 0.4, 
                             minGap = 200, 
                             minSize = 50, 
                             minReadsPerCytosine = 4, 
                             cores = 1) {
  
  regions <- reduce(regions)
  
  # extract the methylation data in the correct context
  cat("Extract methylation in the corresponding context \n")
  
  contextMethylationData1 <- methylationData1[methylationData1$context%in%context]
  contextMethylationData2 <- methylationData2[methylationData2$context%in%context]
  
  
  localContextMethylationData1 <- contextMethylationData1[queryHits(findOverlaps(contextMethylationData1, regions))]
  localContextMethylationData2 <- contextMethylationData2[queryHits(findOverlaps(contextMethylationData2, regions))]
  
  
  localContextMethylationData <- .joinMethylationData(localContextMethylationData1, localContextMethylationData2)
  
  
  
  regionsList <- .splitGRangesEqualy(regions, cores)
  
  # inner loop function for parallel::mclapply
  .computeDMRsBinsLoop = function(i){
    computedDMRs <- GRanges()
    for(index in 1:length(regionsList[[i]])){
      currentRegion <- regionsList[[i]][index]
      
      
      cat("Computing DMRs at ",.printGenomicRanges(currentRegion),"\n")
      
      seqs <- seq(start(currentRegion), (end(currentRegion)-binSize), by = binSize);
      
      bins <- GRanges(seqnames(currentRegion), IRanges(seqs, (seqs+binSize-1)))
      
      overlapsBins <- findOverlaps(localContextMethylationData, currentRegion)
      
      if(length(overlapsBins) > 0){
        localMethylationData <- localContextMethylationData[queryHits(overlapsBins)]
        
        cat("Count inside each bin...\n")
        #bins <- .analyseReadsInsideRegions(localMethylationData, bins, context, cores)
        bins <- .analyseReadsInsideBins(localMethylationData, bins, currentRegion)
        
        cat("Filter the bins...\n")
        # Get rid of the bins with fewer than minCytosinesCount cytosines inside.  
        bins  <- bins[bins$cytosinesCount >= minCytosinesCount]
        
        # Get rid of the bins with fewer than minReadsPerCytosine reads per cytosine.  
        bins  <- bins[(bins$sumReadsN1/bins$cytosinesCount >= minReadsPerCytosine) & 
                        (bins$sumReadsN2/bins$cytosinesCount >= minReadsPerCytosine)]
        
        # Get rid of the bins with small difference in proportion of methylation 
        bins  <- bins[(abs(bins$proportion1 - bins$proportion2) >= minProportionDifference)]
        
        
        cat("Identifying DMRs...\n")    
        pValue <- .computeAdjuestedPValues(test, bins$sumReadsM1, bins$sumReadsN1, bins$sumReadsM2, bins$sumReadsN2, alternative = "two.sided")
        
        bins <- bins[!is.na(pValue) & pValue < pValueThreshold ]  
        bins$context <- rep(paste(context, collapse = "_"), length(bins))
        bins$direction <- rep(NA, length(bins))
        bins$direction <- sign(bins$proportion2 - bins$proportion1)
        
        # Select the crude list of DMRs
        DMRs <- bins[!is.na(bins$direction) & (bins$direction == 1 | bins$direction == -1)]
        
        # append current DMRs to the global list of DMRs 
        if(length(computedDMRs) == 0){
          computedDMRs <- DMRs
        } else{
          computedDMRs <- c(computedDMRs,DMRs)
        }
      }
    }
    return(computedDMRs)
  }
  
  # compute the DMRs
  if(cores > 1){
    cat("Compute the DMRs using ", cores, "cores\n")
    computedDMRs <- parallel::mclapply(1:length(regionsList), .computeDMRsBinsLoop, mc.cores = cores)
  } else {
    computedDMRs <- lapply(1:length(regionsList), .computeDMRsBinsLoop)
  }
  
  
  computedDMRs <- unlist(GRangesList(computedDMRs))
  
  if(length(computedDMRs) > 0){
    
    cat("Merge adjacent DMRs\n")    
    computedDMRs <- computedDMRs[order(computedDMRs)]
    computedDMRs$pValue <- .computeaAjustedPValuesInDMRs(test, computedDMRs ,alternative = "two.sided")
    
    cat("Merge DMRs iteratively\n")    
    # Get rid of small gaps between DMRs.
    if(minGap > 0){
      computedDMRs <- .smartMergeDMRs(computedDMRs, 
                                      minGap = minGap, 
                                      respectSigns = TRUE, 
                                      methylationData = localContextMethylationData,
                                      minProportionDifference=minProportionDifference, 
                                      minReadsPerCytosine = minReadsPerCytosine, 
                                      pValueThreshold=pValueThreshold,
                                      test=test, 
                                      alternative = "two.sided",
                                      cores = cores)
    }  
    
    cat("Filter DMRs \n")    
    if(length(computedDMRs) > 0){
      #remove small DMRs 
      computedDMRs <- computedDMRs[width(computedDMRs) >= minSize]
      
      if(length(computedDMRs) > 0){
        #remove DMRs with few cytosines
        computedDMRs <- computedDMRs[computedDMRs$cytosinesCount >= minCytosinesCount]
        #recompute the adjusted p-values
        if(length(computedDMRs) > 0){
          computedDMRs$pValue <- .computeaAjustedPValuesInDMRs(test, computedDMRs ,alternative = "two.sided")
          computedDMRs$regionType <- rep("loss", length(computedDMRs))
          computedDMRs$regionType[which(computedDMRs$proportion1 < computedDMRs$proportion2)] <- "gain"
        }
      }
    }
  }  
  
  return(computedDMRs)
}




#' This function verifies whether a set of pottential DMRs (e.g. genes, 
#' transposons, CpG islands) are differentially methylated or not.
#'
#' @title Filter DMRs 
#' @param methylationData1 the methylation data in condition 1
#' (see \code{\link{methylationDataList}}).
#' @param methylationData2 the methylation data in condition 2
#' (see \code{\link{methylationDataList}}).
#' @param potentialDMRs a \code{\link{GRanges}} object with potential DMRs 
#' where to compute the DMRs. This can be a a list of gene and/or transposable 
#' elements coordinates.
#' @param context the context in which the DMRs are computed (\code{"CG"}, 
#' \code{"CHG"} or \code{"CHH"}).
#' @param test the statistical test used to call DMRs (\code{"fisher"} for 
#' Fisher's exact test or \code{"score"} for Score test). 
#' @param pValueThreshold DMRs with p-values (when performing the statistical 
#' test; see \code{test}) higher or equal than \code{pValueThreshold} are 
#' discarded. Note that we adjust the p-values using the Benjamini and 
#' Hochberg's method to control the false discovery rate.
#' @param minCytosinesCount DMRs with less cytosines in the specified context 
#' than \code{minCytosinesCount} will be discarded.
#' @param minProportionDifference DMRs where the difference in methylation 
#' proportion between the two conditions is lower than 
#' \code{minProportionDifference} are discarded.
#' @param minReadsPerCytosine  DMRs with the average number of reads lower than 
#' \code{minReadsPerCytosine} are discarded. 
#' @param cores the number of cores used to compute the DMRs. 
#' @return a \code{\link{GRanges}} object with 11 metadata columns that contain 
#' the DMRs; see \code{\link{computeDMRs}}.
#' @seealso \code{\link{DMRsNoiseFilterCG}}, \code{\link{computeDMRs}}, 
#' \code{\link{analyseReadsInsideRegionsForCondition}}  
#' and \code{\link{mergeDMRsIteratively}}
#' @examples
#' # load the methylation data
#' data(methylationDataList)
#' # load the gene annotation data
#' data(GEs)
#' 
#' #select the genes
#' genes <- GEs[which(GEs$type == "gene")]
#' 
#' # the regions where to compute the DMRs
#' regions <- GRanges(seqnames = Rle("Chr3"), ranges = IRanges(1,1E5))
#' genes <- genes[overlapsAny(genes, regions)]
#' 
#' # filter genes that are differntially methylated in the two conditions
#' DMRsGenesCG <- filterDMRs(methylationDataList[["WT"]], 
#'                methylationDataList[["met1-3"]], potentialDMRs = genes, 
#'                context = "CG", test = "score", pValueThreshold = 0.01, 
#'                minCytosinesCount = 4, minProportionDifference = 0.4, 
#'                minReadsPerCytosine = 3, cores = 1)
#'
#' @author Nicolae Radu Zabet
#' @export
filterDMRs <- function(methylationData1, 
                       methylationData2, 
                       potentialDMRs, 
                       context = "CG",   
                       test = "fisher", 
                       pValueThreshold = 0.01, 
                       minCytosinesCount = 4, 
                       minProportionDifference = 0.4, 
                       minReadsPerCytosine = 3, 
                       cores = 1) {
  ##Parameters checking
  cat("Parameters checking ...\n")
  
  .validateMethylationData(methylationData1, variableName="methylationData1")
  .validateMethylationData(methylationData2, variableName="methylationData2")
  
  regions <- union(getWholeChromosomes(methylationData1), 
                   getWholeChromosomes(methylationData2))
  
  .validateContext(context)
  
  
  .validateGRanges(potentialDMRs, generateGenomeWide=FALSE, variableName="potentialDMRs", minLength=NULL)
  
  .validateStatisticalTest(test)
  
  .stopIfNotAll(c(!is.null(pValueThreshold), is.numeric(pValueThreshold), pValueThreshold > 0, pValueThreshold < 1),
                " the p-value threshold needs to be in the interval (0,1)")
  
  .stopIfNotAll(c(.isInteger(minCytosinesCount, positive=TRUE)), 
                " the minCytosinesCount is an integer higher or equal to 0")
  
  .stopIfNotAll(c(!is.null(minProportionDifference), is.numeric(minProportionDifference), minProportionDifference > 0, minProportionDifference < 1), 
                " the minimum difference in methylation needs to be in the interval (0,1)")
  
  .stopIfNotAll(c(.isInteger(minReadsPerCytosine, positive=TRUE)), 
                " the minimum number of reads in a bin is an integer higher or equal to 0")
  
  .stopIfNotAll(c(.isInteger(cores, positive=TRUE)), 
                " the number of cores to use when computing the DMRs.")
  
  regions <- reduce(regions)
  
  if(length(potentialDMRs) > 0){
    
    
    
    # extract the methylation data in the correct context
    cat("Extract methylation in the corresponding context \n")
    
    contextMethylationData1 <- methylationData1[methylationData1$context%in%context]
    contextMethylationData2 <- methylationData2[methylationData2$context%in%context]
    
    
    localContextMethylationData1 <- contextMethylationData1[queryHits(findOverlaps(contextMethylationData1, regions))]
    localContextMethylationData2 <- contextMethylationData2[queryHits(findOverlaps(contextMethylationData2, regions))]
    
    
    localContextMethylationData <- .joinMethylationData(localContextMethylationData1, localContextMethylationData2)
    
    
    regionsList <- .splitGRangesEqualy(regions, cores)
    
    # inner loop function for parallel::mclapply
    .filterDMRsLoop = function(i){
      computedDMRs <- GRanges()  
      for(index in 1:length(regionsList[[i]])){
        currentRegion <- regionsList[[i]][index]
        
        
        cat("Computing DMRs at ",.printGenomicRanges(currentRegion),"\n")
        
        cat("Selecting data...\n")
        
        # Select the points in methylationData that we're interested in. These are the 
        # points that lie within 'regions', as well as any that lie within 
        # window.size of them. 
        
        overlapsPotentialDMRs <- findOverlaps(potentialDMRs, currentRegion)
        if(length(overlapsPotentialDMRs) > 0){
          potentialDMRsLocal <- potentialDMRs[queryHits(overlapsPotentialDMRs)]
          
          localMethylationData <- localContextMethylationData[queryHits(findOverlaps(localContextMethylationData, currentRegion))]
          potentialDMRsLocal <- .analyseReadsInsideRegions(localMethylationData, potentialDMRsLocal)
          
          if(length(computedDMRs) == 0){
            computedDMRs <- potentialDMRsLocal
          } else{
            computedDMRs <- c(computedDMRs,potentialDMRsLocal)
          }
        } 
      }
      return(computedDMRs)
    }
    
    # compute the DMRs
    if(cores > 1){
      cat("Compute the DMRs using ", cores, "cores\n")
      computedDMRs <- parallel::mclapply(1:length(regionsList), .filterDMRsLoop, mc.cores = cores)
    } else {
      computedDMRs <- lapply(1:length(regionsList), .filterDMRsLoop)
    }
    
    
    computedDMRs <-  unlist(GRangesList(computedDMRs))
    
    if(length(computedDMRs) > 0){
      cat("Identifying DMRs...\n")    
      pValue <- .computeAdjuestedPValues(test, computedDMRs$sumReadsM1, computedDMRs$sumReadsN1, computedDMRs$sumReadsM2, computedDMRs$sumReadsN2, alternative = "two.sided")
      
      
      bufferIndex <- !is.na(pValue) & 
        pValue < pValueThreshold & 
        abs(computedDMRs$proportion1 - computedDMRs$proportion2) >= minProportionDifference & 
        computedDMRs$sumReadsN1/computedDMRs$cytosinesCount >= minReadsPerCytosine &
        computedDMRs$sumReadsN2/computedDMRs$cytosinesCount >= minReadsPerCytosine &
        computedDMRs$cytosinesCount >= minCytosinesCount
      
      computedDMRs <- computedDMRs[bufferIndex]  
    }
  } else{
    computedDMRs <- GRanges() 
  }
  
  if(length(computedDMRs) > 0){
    computedDMRs$pValue <- .computeaAjustedPValuesInDMRs(test, computedDMRs ,alternative = "two.sided")
    computedDMRs$regionType <- rep("loss", length(computedDMRs))
    computedDMRs$regionType[which(computedDMRs$proportion1 < computedDMRs$proportion2)] <- "gain"
    computedDMRs$direction <- rep(-1, length(computedDMRs))
    computedDMRs$direction[which(computedDMRs$proportion1 < computedDMRs$proportion2)] <- 1
    computedDMRs <- computedDMRs[order(computedDMRs)]
  }
  
  return(computedDMRs)
  
}



#' This function takes a list of DMRs and attempts to merge DMRs while keeping 
#' the new DMRs statistically significant.
#'
#' @title Merge DMRs iteratively
#' @param DMRs the list of DMRs as a \code{\link{GRanges}} object; 
#' e.g. see \code{\link{computeDMRs}}
#' @param minGap DMRs separated by a gap of at least \code{minGap} are not 
#' merged.
#' @param respectSigns logical value indicating whether to respect the sign when 
#' joining DMRs.
#' @param methylationData1 the methylation data in condition 1
#' (see \code{\link{methylationDataList}}).
#' @param methylationData2 the methylation data in condition 2
#' (see \code{\link{methylationDataList}}).
#' @param context the context in which the DMRs are computed (\code{"CG"}, 
#' \code{"CHG"} 
#' or \code{"CHH"}).
#' @param minProportionDifference two adjacent DMRs are merged only if the 
#' difference in methylation proportion of the new DMR is higher than 
#' \code{minProportionDifference}.
#' @param minReadsPerCytosine two adjacent DMRs are merged only if the number of 
#' reads per cytosine of the new DMR is higher than \code{minReadsPerCytosine}.
#' @param pValueThreshold two adjacent DMRs are merged only if the p-value of 
#' the new DMR (see \code{test} below) is lower than \code{pValueThreshold}. 
#' Note that we adjust the p-values using the Benjamini and Hochberg's method to 
#' control the false discovery rate.
#' @param test the statistical test used to call DMRs (\code{"fisher"} for 
#' Fisher's exact test or \code{"score"} for Score test). 
#' @param alternative indicates the alternative hypothesis and must be one of 
#' \code{"two.sided"}, \code{"greater"} or \code{"less"}.
#' @param cores the number of cores used to compute the DMRs. 
#' @return the reduced list of DMRs as a \code{\link{GRanges}} object; 
#' e.g. see \code{\link{computeDMRs}}
#' @seealso \code{\link{filterDMRs}}, \code{\link{computeDMRs}}, 
#' \code{\link{analyseReadsInsideRegionsForCondition}} and 
#' \code{\link{DMRsNoiseFilterCG}}
#' @examples
#' # load the methylation data
#' data(methylationDataList)
#' 
#' #load the DMRs in CG context they were computed with minGap = 200
#' data(DMRsNoiseFilterCG)
#' 
#'
#' #merge the DMRs 
#' DMRsNoiseFilterCGLarger <- mergeDMRsIteratively(DMRsNoiseFilterCG[1:100], 
#'                            minGap = 500, respectSigns = TRUE, 
#'                            methylationDataList[["WT"]], 
#'                            methylationDataList[["met1-3"]],
#'                            context = "CG", minProportionDifference=0.4, 
#'                            minReadsPerCytosine = 1, pValueThreshold=0.01, 
#'                            test="score",alternative = "two.sided")
#' 
#' 
#' \dontrun{
#' #set genomic coordinates where to compute DMRs
#' regions <- GRanges(seqnames = Rle("Chr3"), ranges = IRanges(1,1E5))
#' 
#' # compute DMRs and remove gaps smaller than 200 bp
#' DMRsNoiseFilterCG200 <- computeDMRs(methylationDataList[["WT"]], 
#'                        methylationDataList[["met1-3"]], regions = regions, 
#'                        context = "CG", method = "noise_filter", 
#'                        windowSize = 100, kernelFunction = "triangular",  
#'                        test = "score", pValueThreshold = 0.01, 
#'                        minCytosinesCount = 1, minProportionDifference = 0.4, 
#'                        minGap = 200, minSize = 0, minReadsPerCytosine = 1, 
#'                        cores = 1)
#'                        
#' DMRsNoiseFilterCG0 <- computeDMRs(methylationDataList[["WT"]], 
#'                        methylationDataList[["met1-3"]], regions = regions, 
#'                        context = "CG", method = "noise_filter", 
#'                        windowSize = 100, kernelFunction = "triangular",  
#'                        test = "score", pValueThreshold = 0.01, 
#'                        minCytosinesCount = 1, minProportionDifference = 0.4, 
#'                        minGap = 0, minSize = 0, minReadsPerCytosine = 1, 
#'                        cores = 1)
#' DMRsNoiseFilterCG0Merged200 <- mergeDMRsIteratively(DMRsNoiseFilterCG0, 
#'                              minGap = 200, respectSigns = TRUE, 
#'                              methylationDataList[["WT"]], 
#'                              methylationDataList[["met1-3"]],
#'                              context = "CG", minProportionDifference=0.4, 
#'                              minReadsPerCytosine = 1, pValueThreshold=0.01, 
#'                              test="score",alternative = "two.sided")                      
#' 
#' #check that all newley computed DMRs are identical
#' print(all(DMRsNoiseFilterCG200 == DMRsNoiseFilterCG0Merged200))
#' 
#' }
#' 
#' @author Nicolae Radu Zabet 
#' 
#' @export
mergeDMRsIteratively <- function(DMRs, 
                                 minGap, 
                                 respectSigns = TRUE, 
                                 methylationData1,
                                 methylationData2,
                                 context = "CG",
                                 minProportionDifference=0.4, 
                                 minReadsPerCytosine = 4, 
                                 pValueThreshold=0.01,
                                 test="fisher",
                                 alternative = "two.sided",
                                 cores = 1){
  
  ##Parameters checking
  cat("Parameters checking ...\n")
  
  .validateMethylationData(methylationData1, variableName="methylationData1")
  .validateMethylationData(methylationData2, variableName="methylationData2")
  .validateContext(context)
  .validateGRanges(DMRs, generateGenomeWide=FALSE, variableName="DMRs", minLength=NULL)
  .validateStatisticalTest(test)
  .stopIfNotAll(c(!is.null(pValueThreshold), 
                  is.numeric(pValueThreshold), 
                  pValueThreshold > 0, 
                  pValueThreshold < 1),
                " the p-value threshold needs to be in the interval (0,1)")
  .stopIfNotAll(c(!is.null(minProportionDifference), 
                  is.numeric(minProportionDifference), 
                  minProportionDifference > 0, 
                  minProportionDifference < 1), 
                " the minimum difference in methylation needs to be in the interval (0,1)")
  
  .stopIfNotAll(c(.isInteger(minGap, positive=TRUE)),
                " the minimum gap between DMRs is an integer higher or equal to 0")
  
  .stopIfNotAll(c(.isInteger(minReadsPerCytosine, positive=TRUE)), 
                " the minimum average number of reads in a DMR is an integer higher or equal to 0")
  
  .stopIfNotAll(c(.isInteger(cores, positive=TRUE)), 
                " the number of cores to use when computing the DMRs.")
  
  contextMethylationData1 <- methylationData1[methylationData1$context%in%context]
  contextMethylationData2 <- methylationData2[methylationData2$context%in%context]
  
  totalRegion <- reduce(DMRs, drop.empty.ranges=FALSE, min.gapwidth=minGap, ignore.strand=TRUE)
  localContextMethylationData1 <- contextMethylationData1[queryHits(findOverlaps(contextMethylationData1, totalRegion))]
  localContextMethylationData2 <- contextMethylationData2[queryHits(findOverlaps(contextMethylationData2, totalRegion))]
  localContextMethylationData <- .joinMethylationData(localContextMethylationData1, localContextMethylationData2)
  
  cat("Merge DMRs iteratively ...\n")
  
  
  return(.smartMergeDMRs(DMRs, 
                         minGap = minGap, 
                         respectSigns = respectSigns, 
                         methylationData = localContextMethylationData,
                         minProportionDifference=minProportionDifference, 
                         minReadsPerCytosine = minReadsPerCytosine, 
                         pValueThreshold=pValueThreshold,
                         test=test, 
                         alternative = alternative, 
                         cores = cores))                                               
}


#' This function extracts from the methylation data the total number of reads, 
#' the number of methylated reads and the number of cytosines in the specific 
#' context from a region (e.g. DMRs)
#'
#' @title Analyse reads inside regions for condition
#' @param regions a \code{\link{GRanges}} object with a list of regions on the 
#' genome; e.g. could be a list of DMRs
#' @param methylationData the methylation data in one condition
#' (see \code{\link{methylationDataList}}).
#' @param context the context in which to extract the reads (\code{"CG"}, 
#' \code{"CHG"} or \code{"CHH"}).
#' @param label a string to be added to the columns to identify the condition
#' @param cores the number of cores used to compute the DMRs. 
#' @return a \code{\link{GRanges}} object with additional four metadata columns
#' \describe{
#'  \item{sumReadsM}{the number of methylated reads}
#'  \item{sumReadsN}{the total number of reads} 
#'  \item{proportion}{the proportion methylated reads} 
#'  \item{cytosinesCount}{the number of cytosines in the regions} 
#' }
#' @seealso \code{\link{filterDMRs}}, \code{\link{computeDMRs}}, 
#' \code{\link{DMRsNoiseFilterCG}}, and \code{\link{mergeDMRsIteratively}}
#' @examples
#' 
#' # load the methylation data
#' data(methylationDataList)
#'  
#' #load the DMRs in CG context. These DMRs were computed with minGap = 200.
#' data(DMRsNoiseFilterCG)
#' 
#' #retrive the number of reads in CHH context in WT
#' DMRsNoiseFilterCGreadsCHH <- analyseReadsInsideRegionsForCondition(
#'                              DMRsNoiseFilterCG[1:10], 
#'                              methylationDataList[["WT"]], context = "CHH", 
#'                              label = "WT")
#' 
#' 
#' @author Nicolae Radu Zabet 
#' 
#' @export
analyseReadsInsideRegionsForCondition <- function(regions,
                                                  methylationData,
                                                  context,
                                                  label = "",
                                                  cores = 1){
  
  ##Parameters checking
  cat("Parameters checking ...\n")
  .validateGRanges(regions, generateGenomeWide=FALSE, variableName="regions", minLength=NULL)
  .validateMethylationData(methylationData, variableName="methylationData")
  .validateContext(context)
  .stopIfNotAll(c(.isInteger(cores, positive=TRUE)), 
                " the number of cores to use when computing the DMRs.")
  
  cat("Extract methylation levels in corresponding context ...\n")
  contextMethylationData <- methylationData[methylationData$context%in%context]
  
  cat("Compute reads inside each region ...\n")
  if(length(regions) > 0){
    if(cores > 1){
      cat("Analyse reads in regions using ", cores, "cores\n")
      regionsList <- split(regions, rep(1:cores, length.out=length(regions)))
      .analyseReadsInsideRegionsForConditionLoop = function(i){
        regionsLocal <- .analyseReadsInsideRegionsForCondition(regionsList[[i]], 
                                                               contextMethylationData, 
                                                               label = label,
                                                               context = context)
        return(regionsLocal)
      }
      regions <- parallel::mclapply(1:length(regionsList), 
                                    .analyseReadsInsideRegionsForConditionLoop, mc.cores = cores)
      regions <- unlist(GRangesList(regions))
    } else{
      regions <- .analyseReadsInsideRegionsForCondition(regions, 
                                                        contextMethylationData, 
                                                        label = label,
                                                        context = context)
    }
    regions <- regions[order(regions)]
  }
  return(regions)
}  
