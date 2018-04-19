#' This function computes the correlation of the methylation levels for
#' Cytosines located at a certain distance apart. The function returns a the
#' correlation of methylation levels at distance equal to the specified
#' thresholds.
#'
#' @title Compute methylation data spatial correlation
#' @param methylationData the methylation data stored as a \code{\link{GRanges}}
#' object with four metadata columns (see \code{\link{methylationDataList}}).
#' @param regions a \code{\link{GRanges}} object with the regions where to
#' compute the correlation If NULL, the correlation is computed genome-wide.
#' @param distance the distance used when computing the correlation of
#' methylation levels
#'
#' @return correlation of the methylation levels for Cytosines located at a
#' certain distance apart.
#'
#' @author Radu Zabet
.computeMethylationDataSpatialCorrelation <- function(methylationData, regions, distance){
  values1 <- c()
  values2 <- c()
  for(i in 1:length(regions)){
    localMethylationData <- methylationData[queryHits(findOverlaps(methylationData, regions[i]))]
    localMethylationData <- localMethylationData[localMethylationData$readsN > 0]
    localMethylationDataValues <- localMethylationData$readsM/localMethylationData$readsN
    localMethylationDataPositions <- start(localMethylationData)

    localMethylationDatamatch <- match(localMethylationDataPositions,
                                       (localMethylationDataPositions + distance))


    ids <- !is.na(localMethylationDatamatch)

    values1 <- c(values1, localMethylationDataValues[ids])
    values2 <- c(values2, localMethylationDataValues[localMethylationDatamatch[ids]])

  }

  return(cor(values1, values2))
}

#' This function computes the correlation of the methylation levels as a
#' function of the distances between the Cytosines. The function returns a
#' \code{vector} with the correlation of methylation levels at distance equal to
#' a \code{vector} of specified thresholds.
#'
#' @title Compute methylation data spatial correlation
#' @param methylationData the methylation data stored as a \code{\link{GRanges}}
#' object with four metadata columns (see \code{\link{methylationDataList}}).
#' @param regions a \code{\link{GRanges}} object with the regions where to
#' compute the correlation. If \code{NULL}, the correlation is computed
#' genome-wide.
#' @param context the context in which the correlation is computed (\code{"CG"},
#' \code{"CHG"} or \code{"CHH"}).
#' @param distances a \code{numeric} vector specifing the different values for
#' the distances when computing the correlation.
#'
#' @return a \code{vector} with the correlation of the methylation levels for
#' Cytosines located at distances specified in the \code{distances}
#' \code{vector}.
#'
#' @seealso \code{\link{plotMethylationDataSpatialCorrelation}},
#' \code{\link{methylationDataList}}
#' @examples
#'
#' \dontrun{
#' # load the methylation data
#' data(methylationDataList)
#'
#' # compute spatial correlation in CG context
#' distances <- c(1,5,10,15)
#' correlation_CG_wt <- computeMethylationDataSpatialCorrelation(methylationDataList[["WT"]],
#'                  context="CG", distances=distances)
#'
#'
#' }
#' @author Nicolae Radu Zabet
#' @export
computeMethylationDataSpatialCorrelation <- function(methylationData,
                                                      regions = NULL,
                                                      context = "CG",
                                                      distances=NULL) {

  .validateMethylationData(methylationData)

  .stopIfNotAll(c(all(sapply(distances,.isInteger, positive=TRUE)),
                  (length(distances) > 1)),
                " the distances is a vector of integers with the distance between the positions .")

  regions <- .validateGRanges(regions, methylationData)

  .validateContext(context)

  contextMethylationData <- methylationData[methylationData$context%in%context];
  rm(methylationData)
  localMethylationData <- contextMethylationData[queryHits(findOverlaps(contextMethylationData, regions))]
  rm(contextMethylationData)

  correlation <- c()

  for(k in distances){
    cat("Computing methylation levels correlation for distances of ",k," bp\n")
    correlation <- c(correlation,
                     .computeMethylationDataSpatialCorrelation(localMethylationData,
                                                                regions,
                                                                k))
  }
  names(correlation) <- distances

  return(correlation)
}



#' This function plots the correlation of methylation levels for Cytosines
#' located at a certain distance apart.
#'
#' @title Plot methylation data spatial correlation
#' @param methylationData1 the methylation data in condition 1
#' (see \code{\link{methylationDataList}}).
#' @param methylationData2 the methylation data in condition 2
#' (see \code{\link{methylationDataList}}). This is optional.
#' @param distances a \code{numeric} vector specifing the different values for
#' the distances when computing the correlation.
#' @param regions a \code{\link{GRanges}} object with the regions where to
#' compute the correlation. If \code{NULL}, the coverage is computed genome-wide.
#' @param conditionsNames a vector of character with the names of the conditions
#' for \code{methylationData1} and \code{methylationData2}.
#' @param context the context in which the DMRs are computed (\code{"CG"},
#' \code{"CHG"} or \code{"CHH"}).
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
#' @param log a \code{character} indicating if any of the axes will be displayed
#' on log scale. This argument will be passed to \code{\link{plot}} function.
#' @details This function plots the proportion of cytosines in a specific
#' context that have at least a certain number of reads (x-axis)
#' @seealso \code{\link{computeMethylationDataSpatialCorrelation}},
#' \code{\link{methylationDataList}}
#' @return Invisibly returns \code{NULL}
#' @examples
#'
#' \dontrun{
#' # load the methylation data
#' data(methylationDataList)
#'
#' # plot the spatial correlation in CG context
#' par(mar=c(4, 4, 3, 1)+0.1)
#' plotMethylationDataSpatialCorrelation(methylationDataList[["WT"]],
#'                            distances = c(1,5,10,15), regions = NULL,
#'                            conditionsNames = c("WT","met1-3"),
#'                            context = c("CG"),
#'                            labels = LETTERS, col = NULL,
#'                            pch = c(1,0,16,2,15,17), lty = c(4,1,3,2,6,5),
#'                            contextPerRow = FALSE)
#'
#' # plot the spatial correlation in all three contexts
#' plotMethylationDataSpatialCorrelation(methylationDataList[["WT"]],
#'                            methylationDataList[["met1-3"]],
#'                            distances = c(1,5,10,15,20,50,100,150,200,500,1000),
#'                            regions = NULL, conditionsNames = c("WT","met1-3"),
#'                            context = c("CG", "CHG", "CHH"),
#'                            labels = LETTERS, col = NULL,
#'                            pch = c(1,0,16,2,15,17), lty = c(4,1,3,2,6,5),
#'                            contextPerRow = FALSE, log="x")
#' }
#'
#' @author Nicolae Radu Zabet
#' @export
plotMethylationDataSpatialCorrelation <- function(methylationData1,
                                                   methylationData2 = NULL,
                                                   distances,
                                                   regions = NULL,
                                                   conditionsNames = NULL,
                                                   context = "CG",
                                                   labels=NULL,
                                                   col=NULL,
                                                   pch = c(1,0,16,2,15,17),
                                                   lty = c(4,1,3,2,6,5),
                                                   contextPerRow = FALSE,
                                                   log="") {
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
    buffer_correlation <- matrix(0, nrow = length(distances), ncol=numberOfConditions)
    buffer_correlation[,1] <- computeMethylationDataSpatialCorrelation(methylationData1,
                                                                     regions = regions,
                                                                     context = context[i],
                                                                     distances=distances)

    if(numberOfConditions > 1){
      buffer_correlation[,2] <- computeMethylationDataSpatialCorrelation(methylationData2,
                                                                       regions = regions,
                                                                       context = context[i],
                                                                       distances=distances)
    }
    colnames(buffer_correlation) <- conditionsNames

    symbols_ind <- (numberOfConditions*(i-1) + 1):(numberOfConditions*i)

    plot(distances, buffer_correlation[, 1], type = "o", ylim = c(0,1), main = paste("Correlation of methylation levels in ",context[i]," context",sep=""),
         xlab="distance between cytosines (bp)", ylab="correlation", col=col[symbols_ind[1]], xaxt="n", yaxt="n",
         lty = lty[symbols_ind[1]], pch = pch[symbols_ind[1]], log=log)

    if(numberOfConditions > 1){
      lines(distances, buffer_correlation[, 2], type = "o", col=col[symbols_ind[2]], lty = lty[symbols_ind[2]], pch = pch[symbols_ind[2]])
    }

    axis(1, at=distances,labels=distances, las=1)
    axis(2, at=seq(0,1,0.1),labels=seq(0,1,0.1), las=1)


    legend("topright", legend = colnames(buffer_correlation), lty = lty[symbols_ind], col = col[symbols_ind], pch = pch[symbols_ind], bty="n")


    mtext(labels[i], line = 0.7, adj = 0, cex=1.0);

  }
  invisible(NULL)
}
