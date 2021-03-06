
#' This function computes the low resolution profiles for the bisulfite
#' sequencing data.
#'
#' @title Compute methylation profile
#' @param methylationData the methylation data stored as a \code{\link{GRanges}}
#' object with four metadata columns (see \code{\link{methylationDataList}}).
#' @param region a \code{\link{GRanges}} object with the regions where to
#' compute the DMRs.
#' @param windowSize a \code{numeric} value indicating the size of the window in
#' which methylation is averaged.
#' @param context the context in which the DMRs are computed (\code{"CG"},
#' \code{"CHG"} or \code{"CHH"}).
#' @return a \code{\link{GRanges}} object with equal sized tiles of the
#' \code{region}. The object consists of the following metadata
#' \describe{
#'  \item{sumReadsM}{the number of methylated reads.}
#'  \item{sumReadsN}{the total number of reads.}
#'  \item{Proportion}{the proportion of methylated reads.}
#'  \item{context}{the context (\code{"CG"}, \code{"CHG"} or \code{"CHH"}).}
#' }
#' @seealso \code{\link{plotMethylationProfileFromData}},
#' \code{\link{plotMethylationProfile}}, \code{\link{methylationDataList}}
#' @examples
#'
#' # load the methylation data
#' data(methylationDataList)
#'
#' # the region where to compute the profile
#' region <- GRanges(seqnames = Rle("Chr3"), ranges = IRanges(1,1E6))
#'
#' # compute low resolution profile in 20 Kb windows
#' lowResProfileWTCHH <- computeMethylationProfile(methylationDataList[["WT"]],
#'                      region, windowSize = 20000, context = "CHH")
#'
#' \dontrun{
#' # compute low resolution profile in 10 Kb windows
#' lowResProfileWTCG <- computeMethylationProfile(methylationDataList[["WT"]],
#'                      region, windowSize = 10000, context = "CG")
#'
#' lowResProfileMet13CG <- computeMethylationProfile(
#'                      methylationDataList[["met1-3"]],  region,
#'                      windowSize = 10000, context = "CG")
#' }
#'
#' @author Nicolae Radu Zabet and Jonathan Michael Foonlan Tsang
#' @export
computeMethylationProfile <- function(methylationData,
                                      region,
                                      windowSize = floor(width(region)/500),
                                      context = "CG") {
  .validateMethylationData(methylationData)

  .validateGRanges(region, methylationData, length=1, generateGenomeWide=FALSE)

  .validateContext(context)

  .stopIfNotAll(c(.isInteger(windowSize, positive=TRUE)),
                " the window size used to compute the methylation profile is an integer higher or equal to 0")


  seqname <- seqnames(region)
  minPos <- start(region)
  maxPos <- end(region)

  hits <- findOverlaps(methylationData, region)



  localMethylationData <- methylationData[queryHits(hits)];
  rm(methylationData)
  contextMethylationData <- localMethylationData[localMethylationData$context%in%context];
  rm(localMethylationData)

  cat("Calculating methylation profile for ", .printGenomicRanges(region), " using a window of ", windowSize, " bp \n")
  seqs = seq(minPos, maxPos - windowSize, windowSize);

  ranges <- GRanges(seqname, IRanges(seqs, seqs+windowSize-1))

  #for windowSize <= 2000 Bins method is faster, but for anything else, regions method is faster
  if(windowSize <= 2000){
    ranges <- .analyseReadsInsideBinsOneSample(contextMethylationData, ranges, region)
  } else{
    ranges <- .analyseReadsInsideRegionsOneSample(contextMethylationData, ranges)
  }

  ranges$context <- paste(context, collapse = "_")

  return(ranges)
}


#' This function plots the low resolution profiles for the bisulfite sequencing
#' data.
#'
#' @title Plot Methylation Profile
#' @param methylationProfiles a \code{GRangesList} object. Each
#' \code{\link{GRanges}} object in the list is generated by calling  the
#' function \code{\link{computeMethylationProfile}}.
#' @param autoscale a \code{logical} value indicating whether the values are
#' autoscalled for each context or not.
#' @param labels a \code{vector} of \code{character} used to add a subfigure
#' characters to the plot. If \code{NULL} nothing is added.
#' @param title the plot title.
#' @param col a \code{character} vector with the colours. It needs to contain a
#' minimum of 2 colours per context. If not or if NULL, the defalut colours will
#' be used.
#' @param pch the R symbols used to plot the data.
#' @param lty the line types used to plot the data.
#' @return Invisibly returns \code{NULL}
#' @seealso \code{\link{plotMethylationProfileFromData}},
#' \code{\link{computeMethylationProfile}} and \code{\link{methylationDataList}}
#' @examples
#'
#'
#' # load the methylation data
#' data(methylationDataList)
#'
#' # the region where to compute the profile
#' region <- GRanges(seqnames = Rle("Chr3"), ranges = IRanges(1,1E6))
#'
#' # compute low resolution profile in 20 Kb windows
#' lowResProfileWTCG <- computeMethylationProfile(methylationDataList[["WT"]],
#'                      region, windowSize = 20000, context = "CG")
#'
#' lowResProfilsCG <- GRangesList("WT" = lowResProfileWTCG)
#'
#' #plot the low resolution profile
#' par(mar=c(4, 4, 3, 1)+0.1)
#' par(mfrow=c(1,1))
#' plotMethylationProfile(lowResProfilsCG, autoscale = FALSE,
#'                        title="CG methylation on Chromosome 3",
#'                        col=c("#D55E00","#E69F00"),  pch = c(1,0),
#'                        lty = c(4,1))
#'
#' \dontrun{
#' # compute low resolution profile in 10 Kb windows in CG context
#' lowResProfileWTCG <- computeMethylationProfile(methylationDataList[["WT"]],
#'                      region, windowSize = 10000, context = "CG")
#'
#' lowResProfileMet13CG <- computeMethylationProfile(
#'                      methylationDataList[["met1-3"]], region,
#'                      windowSize = 10000, context = "CG")
#'
#' lowResProfileCG <- GRangesList("WT" = lowResProfileWTCG,
#'                    "met1-3" = lowResProfileMet13CG)
#'
#' # compute low resolution profile in 10 Kb windows in CHG context
#' lowResProfileWTCHG <- computeMethylationProfile(methylationDataList[["WT"]],
#'                      region, windowSize = 10000, context = "CHG")
#'
#' lowResProfileMet13CHG <- computeMethylationProfile(
#'                      methylationDataList[["met1-3"]], region,
#'                      windowSize = 10000, context = "CHG")
#'
#' lowResProfileCHG <- GRangesList("WT" = lowResProfileWTCHG,
#'                    "met1-3" = lowResProfileMet13CHG)
#'
#' # plot the low resolution profile
#' par(mar=c(4, 4, 3, 1)+0.1)
#' par(mfrow=c(2,1))
#' plotMethylationProfile(lowResProfileCG, autoscale = FALSE,
#'                        labels = LETTERS[1],
#'                        title="CG methylation on Chromosome 3",
#'                        col=c("#D55E00","#E69F00"),  pch = c(1,0),
#'                        lty = c(4,1))
#' plotMethylationProfile(lowResProfileCHG, autoscale = FALSE,
#'                        labels = LETTERS[2],
#'                        title="CHG methylation on Chromosome 3",
#'                        col=c("#0072B2", "#56B4E9"),  pch = c(16,2),
#'                        lty = c(3,2))
#' }
#'
#' @author Nicolae Radu Zabet
#' @export
plotMethylationProfile <- function(methylationProfiles,
                                   autoscale = FALSE,
                                   labels=NULL,
                                   title = "",
                                   col = NULL,
                                   pch = c(1,0,16,2,15,17),
                                   lty = c(4,1,3,2,6,5)) {
  .validateMethylationProfile(methylationProfiles)



  defaultConditionsNames <- paste("condition ", (1:length(methylationProfiles)), sep="")
  if(is.null(names(methylationProfiles))){
    names(methylationProfiles) <- defaultConditionsNames
  }

  if(!.isColor(col, minLength = length(methylationProfiles) )){
    if(length(methylationProfiles) <= 6){
      col <- c("#D55E00","#E69F00", "#0072B2", "#56B4E9", "#F0E442", "#009E73")
    } else if(length(methylationProfiles) <= 8){
      col <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    } else{
      col <- rainbow(length(methylationProfiles))
    }
  }

  .stopIfNotAll(c(!is.null(autoscale),
                  is.logical(autoscale),
                  length(autoscale) == 1),
                " autoscale is a logical value")
  if(is.null(title)){
    title=""
  } else if(is.character(title[1]) > 1){
    title[1] = as.character(title[1])
  }

  .stopIfNotAll(c(!is.null(pch),
                  all(as.integer(pch) == pch),
                  all(pch >= 0),
                  length(pch) >= length(methylationProfiles)),
                paste(" pch is a vector of positive integers of size ",length(methylationProfiles), sep=""))

  .stopIfNotAll(c(!is.null(lty),
                  all(as.integer(lty) == lty),
                  all(lty >= 0),
                  length(lty) >= length(methylationProfiles)),
                paste(" lty is a vector of positive integers of size ",length(methylationProfiles), sep=""))


  ymax <- 1
  if(autoscale){
    ymax <- 0
    for(i in 1:length(methylationProfiles)){
      if(ymax < max(methylationProfiles[[i]]$Proportion)){
        ymax = max(methylationProfiles[[i]]$Proportion)
      }
    }
  }

  options(scipen=10)
  pos <- (start(methylationProfiles[[1]])+end(methylationProfiles[[1]]))/2
  plot(pos, methylationProfiles[[1]]$Proportion, type = "o", ylim = c(0,ymax),
       xlab="genomic coordinate", ylab="methylation", col=col[1], yaxt="n",
       pch = pch[1], lty = lty[1], main = title[1])
  if(length(methylationProfiles) > 1){
    for(i in 2:length(methylationProfiles)){
      lines(pos, methylationProfiles[[i]]$Proportion, type = "o", col=col[i], lty = lty[i], pch = pch[i])
    }
  }

  legend("topright", legend = names(methylationProfiles),
         lty = lty[1:length(methylationProfiles)],
         col = col[1:length(methylationProfiles)],
         pch = pch[1:length(methylationProfiles)], bty="n")
  if(!is.null(labels)){
    mtext(labels[1], line = 0.7, adj = 0, cex=1.4);
  }
  if(autoscale){
    axis(2,c(0,signif(ymax/2, 1),signif(ymax, 1)));
  } else{
    axis(2,c(0,.5,1));
  }
  invisible(NULL)
}

#' This function plots the low resolution profiles for all bisulfite sequencing
#' data.
#'
#' @title Plot methylation profile from data
#' @param methylationData1 the methylation data in condition 1
#' (see \code{\link{methylationDataList}}).
#' @param methylationData2 the methylation data in condition 2
#' (see \code{\link{methylationDataList}}). This is optional.
#' @param regions a \code{\link{GRanges}} object with the regions where to plot
#' the profiles.
#' @param conditionsNames the names of the two conditions. This will be used to
#' plot the legend.
#' @param context a \code{vector} with all contexts in which the DMRs are
#' computed (\code{"CG"}, \code{"CHG"} or \code{"CHH"}).
#' @param windowSize a \code{numeric} value indicating the size of the window in
#' which methylation is averaged.
#' @param autoscale a \code{logical} value indicating whether the values are
#' autoscalled for each context or not.
#' @param labels a \code{vector} of \code{character} used to add a subfigure
#' character to the plot. If \code{NULL} nothing is added.
#' @param col a \code{character} vector with the colours. It needs to contain a
#' minimum of 2 colours per condition. If not or if \code{NULL}, the defalut
#' colours will be used.
#' @param pch the R symbols used to plot the data It needs to contain a minimum
#' of 2 symbols per condition. If not or if \code{NULL}, the defalut symbols
#' will be used.
#' @param lty the line types used to plot the data. It needs to contain a
#' minimum of 2 line types per condition. If not or if \code{NULL}, the defalut
#' line types will be used.
#' @param contextPerRow a \code{logical} value indicating if the each row
#' represents an individual context. If \code{FALSE}, each column will represent
#' an individual context.
#' @return Invisibly returns \code{NULL}
#' @seealso \code{\link{plotMethylationProfile}},
#' \code{\link{computeMethylationProfile}} and \code{\link{methylationDataList}}
#' @examples
#'
#'
#' # load the methylation data
#' data(methylationDataList)
#'
#' #plot the low resolution profile at 10 Kb resolution
#' par(mar=c(4, 4, 3, 1)+0.1)
#' plotMethylationProfileFromData(methylationDataList[["WT"]],
#'                                methylationDataList[["met1-3"]],
#'                                conditionsNames=c("WT", "met1-3"),
#'                                windowSize = 20000, autoscale = TRUE,
#'                                context = c("CHG"))
#'
#' \dontrun{
#' #plot the low resolution profile at 5 Kb resolution
#' par(mar=c(4, 4, 3, 1)+0.1)
#' plotMethylationProfileFromData(methylationDataList[["WT"]],
#'                                methylationDataList[["met1-3"]],
#'                                conditionsNames=c("WT", "met1-3"),
#'                                windowSize = 5000, autoscale = TRUE,
#'                                context = c("CG", "CHG", "CHH"),
#'                                labels = LETTERS)
#' }
#'
#' @author Nicolae Radu Zabet
#' @export
plotMethylationProfileFromData <- function(methylationData1,
                                           methylationData2 = NULL,
                                           regions = NULL,
                                           conditionsNames = NULL,
                                           context = "CG",
                                           windowSize = NULL,
                                           autoscale = FALSE,
                                           labels=NULL,
                                           col=NULL,
                                           pch = c(1,0,16,2,15,17),
                                           lty = c(4,1,3,2,6,5),
                                           contextPerRow = TRUE) {
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


  number_of_symbols <-numberOfConditions*length(context)

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



  #compute all genomic regions for which data is available
  recomputeRegions <- FALSE
  if(is.null(regions)){
    recomputeRegions <- TRUE
  }
  if(!recomputeRegions & (!is(regions, "GRanges") | length(regions) < 1)){
    recomputeRegions <- TRUE
  }


  if(recomputeRegions){
    cat("Recompute regions... \n")
    regions <- GRanges()
    regions <- union(regions, getWholeChromosomes(methylationData1))
    if(numberOfConditions > 1){
      regions <- union(regions, getWholeChromosomes(methylationData2))
    }
  }
  if(is.null(windowSize)){
    windowSize = floor(min(width(regions))/500)
  }


  index <- 1;
  cat("Computing low resolution profiles... \n")


  par(mar=c(4, 4, 3, 1)+0.1)
  if(contextPerRow){
    par(mfrow=c(length(context),length(regions)))
    for(i in 1:length(context)){
      for(j in 1:length(regions)){
        buffer_profile <- GRangesList(computeMethylationProfile(methylationData1, regions[j], windowSize = windowSize, context = context[i]))
        if(numberOfConditions > 1){
          buffer_profile <- c(buffer_profile, GRangesList(computeMethylationProfile(methylationData2, regions[j], windowSize = windowSize, context = context[i])))
        }
        names(buffer_profile) <- conditionsNames

        symbols_ind <- (numberOfConditions*(i-1) + 1):(numberOfConditions*i)

        plotMethylationProfile(buffer_profile, autoscale = autoscale, labels = labels[index],
                               title=paste(context[i], " methylation on chromosome ", as.character(seqnames(regions)[j]),sep=""), col=col[symbols_ind],
                               pch = pch[symbols_ind], lty = lty[symbols_ind])
        index <- index + 1
      }
    }
  } else{
    par(mfrow=c(length(regions),length(context)))
    for(j in 1:length(regions)){
      for(i in 1:length(context)){
        buffer_profile <- GRangesList(computeMethylationProfile(methylationData1, regions[j], windowSize = windowSize, context = context[i]))
        if(numberOfConditions > 1){
          buffer_profile <- c(buffer_profile, GRangesList(computeMethylationProfile(methylationData2, regions[j], windowSize = windowSize, context = context[i])))
        }
        names(buffer_profile) <- conditionsNames

        symbols_ind <- (numberOfConditions*(i-1) + 1):(numberOfConditions*i)
        plotMethylationProfile(buffer_profile, autoscale = autoscale, labels = labels[index],
                               title=paste(context[i], " methylation on chromosome ", as.character(seqnames(regions)[j]),sep=""), col=col[symbols_ind],
                               pch = pch[symbols_ind], lty = lty[symbols_ind])
        index <- index + 1
      }
    }
  }
  invisible(NULL)
}




.countOverlaps <- function(region, subRegions){
  return(sum(countOverlaps(subRegions, region)))
}


.measureOverlaps <- function(region, subRegions){
  return(sum(width(intersect(subRegions, region))))
}



#' This function computes the distribution of a subset of regions
#' (\code{\link{GRanges}} object)  over a large region (\code{\link{GRanges}}
#' object)
#'
#' @title Compute Overlaps Profile
#' @param subRegions a \code{\link{GRanges}} object with the sub regions; e.g.
#' can be the DMRs.
#' @param largeRegion a \code{\link{GRanges}} object with the region where to
#' compute the overlaps; e.g. a chromosome
#' @param windowSize The \code{largeRegion} is partitioned into equal sized
#' tiles of width \code{windowSize}.
#' @param binary a value indicating whether to count 1 for each overlap or to
#' compute the width of the overlap
#' @param cores the number of cores used to compute the DMRs.
#' @return a \code{\link{GRanges}} object with equal sized tiles of the regions.
#' The object has one metadata file \code{score} which represents: the number of
#' subRegions overlapping with the tile, in the case of \code{binary = TRUE},
#' and the width of the subRegions overlapping with the tile , in the case of
#' \code{binary = FALSE}.
#' @seealso \code{\link{plotOverlapProfile}}, \code{\link{filterDMRs}},
#' \code{\link{computeDMRs}} and \code{\link{mergeDMRsIteratively}}
#' @examples
#' # load the methylation data
#' data(methylationDataList)
#'
#' # load the DMRs in CG context
#' data(DMRsNoiseFilterCG)
#'
#' # the coordinates of the area to be plotted
#' largeRegion <- GRanges(seqnames = Rle("Chr3"), ranges = IRanges(1,1E5))
#'
#' # compute overlaps distribution
#' hotspots <- computeOverlapProfile(DMRsNoiseFilterCG, largeRegion,
#'            windowSize = 10000, binary = FALSE)
#' @author Nicolae Radu Zabet
#' @export
computeOverlapProfile <- function(subRegions,
                                  largeRegion,
                                  windowSize = floor(width(largeRegion)/500),
                                  binary=TRUE,
                                  cores = 1){
  .stopIfNotAll(c(!is.null(subRegions),
                  is(subRegions, "GRanges")),
                paste(" subRegions neads to be a GRanges object.",sep=""))

  .stopIfNotAll(c(!is.null(largeRegion),
                  is(largeRegion, "GRanges"),
                  length(largeRegion) == 1),
                paste(" largeRegion neads to be a GRanges object with one location (e.g. entire chromosome).",sep=""))

  .stopIfNotAll(c(!is.null(binary),
                  is.logical(binary),
                  length(binary) == 1),
                paste(" binary neads to be a logical variable.",sep=""))

  .stopIfNotAll(c(.isInteger(windowSize, positive=TRUE)),
                " the window size used to compute the methylation profile is an integer higher or equal to 0")

  .stopIfNotAll(c(.isInteger(cores, positive=TRUE)),
                " the number of cores to used when computing the DMRs needs to be an integer  hirger or equal to 1.")

  seqname <- seqnames(largeRegion)
  minPos <- start(largeRegion)
  maxPos <- end(largeRegion)


  cat("Calculating overlaps for ", .printGenomicRanges(largeRegion), " using a window of ", windowSize, " bp \n")
  seqs = seq(minPos, maxPos - windowSize, windowSize);

  ranges <- GRanges(seqname, IRanges(seqs, seqs+windowSize-1))
  ranges$score <- 0


  rangesList <- tile(ranges, cores)
  if(binary){
    if(cores > 1){
      scores <- parallel::mclapply(rangesList, .countOverlaps, subRegions = subRegions, mc.cores = cores)
    } else {
      scores <- lapply(rangesList, .countOverlaps, subRegions = subRegions)
    }
  } else{
    if(cores > 1){
      scores <- parallel::mclapply(rangesList, .measureOverlaps, subRegions = subRegions, mc.cores = cores)
    } else {
      scores <- lapply(rangesList, .measureOverlaps, subRegions = subRegions)
    }
  }
  ranges$score <- unlist(scores)

  return(ranges)
}



#' This function plots the distribution of a set of subregions on a large
#' region.
#'
#' @title Plot overlap profile
#' @param overlapsProfiles1 a  \code{\link{GRanges}} object with the overlaps
#' profile; see \code{\link{computeOverlapProfile}}.
#' @param overlapsProfiles2 a  \code{\link{GRanges}} object with the overlaps
#' profile; see \code{\link{computeOverlapProfile}}. This is optional. For
#' example, one can be use \code{overlapsProfiles1} to display hypomethylated
#' regions and \code{overlapsProfiles2} the hypermethylated regions.
#' @param names a \code{vector} of \code{character} to add labels for the two
#' overlapsProfiles. This is an optinal parameter.
#' @param labels a \code{vector} of \code{character} used to add a subfigure
#' character to the plot. If \code{NULL} nothing is added.
#' @param col a \code{character} vector with the colours. It needs to contain 2
#' colours. If not or if \code{NULL}, the defalut colours will be used.
#' @param title the title of the plot.
#' @param logscale a \code{logical} value indicating if the colours are on
#' logscale or not.
#' @param maxValue a maximum value in a region. Used for the colour scheme.
#' @return Invisibly returns \code{NULL}.
#' @seealso \code{\link{computeOverlapProfile}}, \code{\link{filterDMRs}},
#' \code{\link{computeDMRs}} and \code{\link{mergeDMRsIteratively}}
#' @examples
#'
#' # load the methylation data
#' data(methylationDataList)
#'
#' # load the DMRs in CG context
#' data(DMRsNoiseFilterCG)
#'
#' # the coordinates of the area to be plotted
#' largeRegion <- GRanges(seqnames = Rle("Chr3"), ranges = IRanges(1,1E5))
#'
#' # compute overlaps distribution
#' hotspotsHypo <- computeOverlapProfile(DMRsNoiseFilterCG, largeRegion,
#'                  windowSize = 10000, binary = FALSE)
#'
#' plotOverlapProfile(GRangesList("Chr3"=hotspotsHypo),
#'                    names = c("hypomethylated"), title = "CG methylation")
#'
#' \dontrun{
#'
#' largeRegion <- GRanges(seqnames = Rle("Chr3"), ranges = IRanges(1,1E6))
#'
#' hotspotsHypo <- computeOverlapProfile(
#'                DMRsNoiseFilterCG[(DMRsNoiseFilterCG$regionType == "loss")],
#'                largeRegion, windowSize=2000, binary=TRUE, cores=1)
#'
#' hotspotsHyper <- computeOverlapProfile(
#'                DMRsNoiseFilterCG[(DMRsNoiseFilterCG$regionType == "gain")],
#'                largeRegion, windowSize=2000, binary=TRUE, cores=1)
#'
#' plotOverlapProfile(GRangesList("Chr3"=hotspotsHypo),
#'                    GRangesList("Chr3"=hotspotsHyper),
#'                    names=c("loss", "gain"), title="CG methylation")
#' }
#'
#' @author Nicolae Radu Zabet
#' @export
plotOverlapProfile <- function(overlapsProfiles1, overlapsProfiles2=NULL, names=NULL,
                               labels=NULL, col=NULL, title="", logscale=FALSE, maxValue=NULL) {

  .stopIfNotAll(c(!is.null(overlapsProfiles1),
                  is(overlapsProfiles1, "GRangesList"),
                  length(overlapsProfiles1) > 0),
                " overlapsProfiles1 needs to be a GRangesList object")
  unlisted_overlapsProfiles1 <- unlist(overlapsProfiles1, use.names=FALSE)
  score <- mcols(unlisted_overlapsProfiles1)$score
  if(is.null(score)) stop("elements of the overlapsProfiles1 are incorrect")
  maxScore <- max(score)
  minScore <- min(score)

  numberOfConditions <- 1

  if(!is.null(overlapsProfiles2)){
    .stopIfNotAll(c(is(overlapsProfiles2, "GRangesList"),
                    length(overlapsProfiles2) == length(overlapsProfiles1)),
                  " overlapsProfiles2 needs to be a GRangesList object with the same number of elements as overlapsProfiles1")
    unlisted_overlapsProfiles2 <- unlist(overlapsProfiles2, use.names=FALSE)
    score <- mcols(unlisted_overlapsProfiles2)$score
    if(is.null(score)) stop("elements of the overlapsProfiles2 are incorrect")
    maxScore <- max(maxScore, score)
    minScore <- min(minScore, score)

    numberOfConditions <- 2
  }




  .stopIfNotAll(c(!is.null(logscale),
                  is.logical(logscale),
                  length(logscale) == 1),
                paste(" logscale neads to be a logical variable.",sep=""))


  if(!.isColor(col) | length(col) < numberOfConditions){
    col <- c("#0072B2", "#D55E00")
  }

  if(!is.null(maxValue)){
    if(is.numeric(maxValue)){
      if(round(maxValue) > 0){
        maxScore <- round(maxValue)
      }
    }
  }

  if(logscale){
    minScore <- log10(1)
    maxScore <- log10(maxScore)
  }

  minScore <- round(minScore)
  maxScore <- round(maxScore)


  cols1 <- c(colorRampPalette(c("white",col[1]))(1 + maxScore - minScore))
  cols2 <- c(colorRampPalette(c("white",col[2]))(1 + maxScore - minScore))

  numberOfProfiles<-length(overlapsProfiles1)

  xmin <- min(unlist(lapply(overlapsProfiles1, start)))
  xmax <- max(unlist(lapply(overlapsProfiles1, end)))

  par(mar=c(4, 4, 3, 1)+0.1)
  plot(c(xmin,xmax),c(0, 1*numberOfProfiles), type="n",xlab="position (bp)", ylab="", yaxt="n", bty="n", main=title)

  for(i in 1:length(overlapsProfiles1)){
    if(!is.null(overlapsProfiles2)){
      #plot two conditions
      for(j in 1:length(overlapsProfiles1[[i]])){
        colId <- overlapsProfiles1[[i]][j]$score
        if(colId > 0){
          if(logscale){
            colId <- log10(colId)
          }
          rect(start(overlapsProfiles1[[i]][j]), (length(overlapsProfiles1) - i + 0.4), end(overlapsProfiles1[[i]][j]), (length(overlapsProfiles1) - i + 0.7), col=cols1[round(colId)+1], border = NA)
        }
      }

      for(j in 1:length(overlapsProfiles2[[i]])){
        colId <- overlapsProfiles2[[i]][j]$score
        if(colId > 0){
          if(logscale){
            colId <- log10(colId)
          }
          rect(start(overlapsProfiles2[[i]][j]), (length(overlapsProfiles2) - i + 0.1), end(overlapsProfiles2[[i]][j]), (length(overlapsProfiles2) - i + 0.4), col=cols2[round(colId)+1], border = NA)
        }
      }

      if(!is.null(names)){
        if(all(is.character(names)) & length(names) >= numberOfConditions){
          at_pos <- c(0.55,0.25)
          for(j in 1:numberOfConditions){
            mtext(names[j], side = 2, line=-0.5, at=(length(overlapsProfiles1) - i + at_pos[j]), col="black", las=2)
          }
        }
      }
    } else{
      #plot one conditions
      for(j in 1:length(overlapsProfiles1[[i]])){
        colId <- overlapsProfiles1[[i]][j]$score
        if(colId > 0){
          if(logscale){
            colId <- log10(colId)
          }
          rect(start(overlapsProfiles1[[i]][j]), (length(overlapsProfiles1) - i + 0.1), end(overlapsProfiles1[[i]][j]), (length(overlapsProfiles1) - i + 0.7), col=cols1[round(colId)+1], border = NA)
        }
      }
      if(!is.null(names)){
        if(all(is.character(names)) & length(names) >= numberOfConditions){
          at_pos <- c(0.4)
          mtext(names[1], side = 2, line=-0.5, at=(length(overlapsProfiles1) - i + at_pos[1]), col="black")
        }
      }
    }


    rect(xmin, (length(overlapsProfiles1) - i + 0.1), max(end(overlapsProfiles1[[i]])), (length(overlapsProfiles1) - i + 0.7), col=NA, border="#999999")
    text(0, (length(overlapsProfiles1) - i + 0.80), names(overlapsProfiles1)[i], pos=4, col="black")
  }

  if(!is.null(labels)){
    mtext(labels[1], line = 0.7, adj = 0, cex=1.4);
  }
  invisible(NULL)

}
