#' This function computes the number of cytosines with more reads than a 
#' threshold in the methylation data
#' @title Count methylation data coverage
#' @param methylationData the methylation data stored as a \code{\link{GRanges}} 
#' object with four metadata columns (see \code{\link{methylationDataList}}).
#' @param regions a \code{\link{GRanges}} object with the regions where to 
#' compute the coverage. If NULL, the coverage is computed genome-wide.
#' @param threshold the threshold used when counting cytosines that have at 
#' least a certain number of reads
#' 
#' @return a the number of cytosines with more reads than the threshold in the 
#' specified region 
#' 
#' @author Radu Zabet
.countMethylationDataCoverage <- function(methylationData, regions, threshold){
  count <- 0
  for (index in 1:length(regions)) {
    currentRegion <- regions[index]
    localMethylationData <- methylationData[queryHits(findOverlaps(methylationData, currentRegion))]
    count <- count + length(which(localMethylationData$readsN > threshold))
  }
  return(count)
}


#' This function computes the coverage for bisulfite sequencing data. It 
#' returns a \code{vector} with the proportion (or raw count) of cytosines that 
#' have the number of reads higher or equal than a \code{vector} of specified 
#' thresholds.
#'
#' @title Compute methylation data coverage
#' @param methylationData the methylation data stored as a \code{\link{GRanges}} 
#' object with four metadata columns (see \code{\link{methylationDataList}}).
#' @param regions a \code{\link{GRanges}} object with the regions where to 
#' compute the coverage. If \code{NULL}, the coverage is computed genome-wide.
#' @param context the context in which the DMRs are computed (\code{"CG"}, 
#' \code{"CHG"} or \code{"CHH"}).
#' @param breaks a \code{numeric} vector specifing the different values for the 
#' thresholds when computing the coverage.
#' @param proportion a \code{logical} value indicating whether to compute the 
#' proportion (\code{TRUE}) or raw counts (\code{FALSE}).
#'  
#' @return a \code{vector} with the proportion (or raw count) of cytosines that 
#' have the number of reads higher or equal than the threshold values specified 
#' in the \code{breaks} \code{vector}.
#' @seealso \code{\link{plotMethylationDataCoverage}}, 
#' \code{\link{methylationDataList}}
#' @examples
#' 
#' # load the methylation data
#' data(methylationDataList)
#' 
#' # compute coverage in CG context
#' breaks <- c(1,5,10,15)
#' coverage_CG_wt <- computeMethylationDataCoverage(methylationDataList[["WT"]], 
#'                  context="CG", breaks=breaks)
#' 
#' 
#' @author Nicolae Radu Zabet and Jonathan Michael Foonlan Tsang
#' @export
computeMethylationDataCoverage <- function(methylationData, 
                                           regions = NULL, 
                                           context = "CG", 
                                           breaks=NULL, 
                                           proportion = TRUE) {
  
  .validateMethylationData(methylationData)
  
  .stopIfNotAll(c(all(sapply(breaks,.isInteger, positive=TRUE)),
                  (length(breaks) > 1)), 
                " the breaks is a vector of integers with the number of reads used 
                as thresholds when computing the coverage.")
  
  regions <- .validateGRanges(regions, methylationData)
  
  .validateContext(context)
  
  contextMethylationData <- methylationData[methylationData$context%in%context];
  
  
  counts <- c()
  for(i in breaks){
    counts <- c(counts,
                .countMethylationDataCoverage(contextMethylationData,regions,i))
  }
  
  
  if(proportion){
    counts <- counts/.countTotalNumberOfCytosines(contextMethylationData, regions)
  }
  
  return(counts)
}



#' This function plots the coverage for the bisulfite sequencing data.
#'
#' @title Plot methylation data coverage
#' @param methylationData1 the methylation data in condition 1
#' (see \code{\link{methylationDataList}}).
#' @param methylationData2 the methylation data in condition 2
#' (see \code{\link{methylationDataList}}). This is optional.
#' @param breaks a \code{numeric} vector specifing the different values for the 
#' thresholds when computing the coverage.
#' @param regions a \code{\link{GRanges}} object with the regions where to 
#' compute the coverage. If \code{NULL}, the coverage is computed genome-wide.
#' @param conditionsNames a vector of character with the names of the conditions 
#' for \code{methylationData1} and \code{methylationData2}. 
#' @param context the context in which the DMRs are computed (\code{"CG"}, 
#' \code{"CHG"} or \code{"CHH"}).
#' @param proportion a \code{logical} value indicating whether proportion or 
#' counts will be plotted.
#' @param labels a \code{vector} of \code{character} used to add a subfigure 
#' character to the plot. If \code{NULL} nothing is added. 
#' @param col a \code{character} vector with the colors. It needs to contain a 
#' minimum of 2 colors per condition. If not or if \code{NULL}, the defalut 
#' colors will be used.
#' @param pch the R symbols used to plot the data. It needs to contain a minimum 
#' of 2 symbols per condition. If not or if \code{NULL}, the defalut symbols 
#' will be used.
#' @param lty the line types used to plot the data. It needs to contain a 
#' minimum of 2 line types per condition. If not or if \code{NULL}, the defalut 
#' line types will be used.
#' @param contextPerRow a \code{logical} value indicating if the each row 
#' represents an individual context. If \code{FALSE}, each column will represent 
#' an individual context.
#' @details This function plots the proportion of cytosines in a specific 
#' context that have at least a certain number of reads (x-axis)
#' @seealso \code{\link{computeMethylationDataCoverage}}, 
#' \code{\link{methylationDataList}}
#' @return Invisibly returns \code{NULL}
#' @examples
#'  
#' # load the methylation data
#' data(methylationDataList)
#' 
#' # plot the coverage in CG context
#' par(mar=c(4, 4, 3, 1)+0.1)
#' plotMethylationDataCoverage(methylationDataList[["WT"]], 
#'                            methylationDataList[["met1-3"]], 
#'                            breaks = c(1,5,10,15), regions = NULL, 
#'                            conditionsNames = c("WT","met1-3"), 
#'                            context = c("CG"), proportion = TRUE, 
#'                            labels = LETTERS, col = NULL,  
#'                            pch = c(1,0,16,2,15,17), lty = c(4,1,3,2,6,5), 
#'                            contextPerRow = FALSE)
#' 
#' \dontrun{
#' # plot the coverage in all three contexts
#' plotMethylationDataCoverage(methylationDataList[["WT"]], 
#'                            methylationDataList[["met1-3"]], 
#'                            breaks = 1:15, regions = NULL, 
#'                            conditionsNames = c("WT","met1-3"), 
#'                            context = c("CG", "CHG", "CHH"), 
#'                            proportion = TRUE, labels = LETTERS, col = NULL,  
#'                            pch = c(1,0,16,2,15,17), lty = c(4,1,3,2,6,5), 
#'                            contextPerRow = FALSE)
#' }
#' 
#' @author Nicolae Radu Zabet
#' @export
plotMethylationDataCoverage <- function(methylationData1, 
                                        methylationData2 = NULL, 
                                        breaks, 
                                        regions = NULL, 
                                        conditionsNames = NULL, 
                                        context = "CG", 
                                        proportion = TRUE, 
                                        labels=NULL,
                                        col=NULL,  
                                        pch = c(1,0,16,2,15,17), 
                                        lty = c(4,1,3,2,6,5), 
                                        contextPerRow = FALSE) {
  .validateContext(context)
  .validateMethylationData(methylationData1, variableName="methylationData1")
  numberOfConditions <- 1
  if(!is.null(methylationData2)){
    .validateMethylationData(methylationData2, variableName="methylationData2")
    numberOfConditions <- 2
  }
  
  if(is.null(conditionsNames) | length(conditionsNames) < numberOfConditions){
    conditionsNames <- paste("condition ",(1:numberOfConditions),sep="")
  }
  
  number_of_symbols <- numberOfConditions*length(context)
  
  if(!.isColor(col, minLength = number_of_symbols )){
    if(number_of_symbols <= 6){
      col <- c("#D55E00","#E69F00", "#0072B2", "#56B4E9", "#F0E442", "#009E73")
    } else if(number_of_symbols <= 8){
      col <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    } else{
      col <- rainbow(number_of_symbols)
    }
  }
  
  .stopIfNotAll(c(!is.null(pch), 
                  all(as.integer(pch) == pch),
                  all(pch >= 0),
                  length(pch) >= number_of_symbols), 
                paste(" pch is a vector of positive integers of size ",number_of_symbols, sep=""))
  
  .stopIfNotAll(c(!is.null(lty), 
                  all(as.integer(lty) == lty),
                  all(lty >= 0),
                  length(lty) >= number_of_symbols), 
                paste(" lty is a vector of positive integers of size ",number_of_symbols, sep=""))
  
  
  if(!is.null(labels) & (length(labels) < length(coverage) | !is.character(labels))){
    labels <- LETTERS[1:length(coverage)]  
  }
  
  
  
  par(mar=c(4, 4, 3, 1)+0.1)
  if(contextPerRow){
    par(mfrow = c(length(context),1))  
  } else{
    par(mfrow = c(1,length(context)))
  }
  
  for(i in 1:length(context)){
    buffer_coverage <- matrix(0, nrow = length(breaks), ncol=numberOfConditions)
    buffer_coverage[,1] <- computeMethylationDataCoverage(methylationData1, 
                                                          regions = regions, 
                                                          context = context[i], 
                                                          breaks=breaks, 
                                                          proportion = proportion)
    
    if(numberOfConditions > 1){
      buffer_coverage[,2] <- computeMethylationDataCoverage(methylationData2, 
                                                            regions = regions, 
                                                            context = context[i], 
                                                            breaks=breaks, 
                                                            proportion = proportion)
    }  
    colnames(buffer_coverage) <- conditionsNames
    
    symbols_ind <- (numberOfConditions*(i-1) + 1):(numberOfConditions*i)
    
    plot(breaks, buffer_coverage[, 1], type = "o", ylim = c(0,1), main = paste("Coverage in ",context[i]," context",sep=""), 
         xlab="minimum number of reads", ylab="coverage", col=col[symbols_ind[1]], xaxt="n", yaxt="n", 
         lty = lty[symbols_ind[1]], pch = pch[symbols_ind[1]])
    
    if(numberOfConditions > 1){
      lines(breaks, buffer_coverage[, 2], type = "o", col=col[symbols_ind[2]], lty = lty[symbols_ind[2]], pch = pch[symbols_ind[2]])
    }
    
    axis(1, at=breaks,labels=breaks, las=1)
    axis(2, at=seq(0,1,0.1),labels=seq(0,1,0.1), las=1)
    
    
    legend("topright", legend = colnames(buffer_coverage), lty = lty[symbols_ind], col = col[symbols_ind], pch = pch[symbols_ind], bty="n")
    
    
    mtext(labels[i], line = 0.7, adj = 0, cex=1.0);  
    
  }
  invisible(NULL)
}
