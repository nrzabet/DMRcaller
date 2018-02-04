
#' This function draws the genetic elements over a plot
#'
#' @title Plot genetic elements
#' @param gff a \code{\link{GRanges}} object with all elements usually imported
#' from a GFF3 file
#' @param region the location on the genome where to draw the genetic elements
#' @param col the colors used to draw the exons and the transposable elemnets
#'
#' @author Radu Zabet
.plotGeneticElements <- function(gff, region, col){

  seqname <- seqnames(region)
  minPos <- start(region)
  maxPos <- end(region)


  # Select the genes that lie in the region of interest.
  gff <- gff[queryHits(findOverlaps(gff, region))]
  # Chop off the ends of anything sticking out...
  start(gff) <- pmax(start(gff), minPos)
  end(gff)   <- pmin(end(gff), maxPos)

  genes <- gff[gff$type == 'gene']
  genesPos <- genes[strand(genes) == '+' | strand(genes) == '*']
  genesNeg <- genes[strand(genes) == '-' | strand(genes) == '*']
  exons <- gff[gff$type == 'exon']
  exons <- exons[overlapsAny(exons, genes)];
  exonsPos <- exons[strand(exons) == '+' | strand(exons) == '*']
  exonsNeg <- exons[strand(exons) == '-' | strand(exons) == '*']

  transposons <- gff[gff$type == 'transposable_element'];
  transposonsPos <- transposons[strand(transposons) == '+' | strand(transposons) == '*'];
  transposonsNeg <- transposons[strand(transposons) == '-' | strand(transposons) == '*'];

  negativeStrandPosition <- -0.175
  positiveStrandPosition <- -0.075


  text(maxPos + (maxPos-minPos)/100, positiveStrandPosition, '+');
  text(maxPos + (maxPos-minPos)/100, negativeStrandPosition, '-');
  lines(c(minPos, maxPos),c(-0.14,-0.14), lty=1, lwd=0.75, col="black")

  if(length(genesPos)>0) {
    segments(start(genesPos), positiveStrandPosition, end(genesPos), positiveStrandPosition);
    text(start(genesPos), -0.115, genesPos$ID, pos=4, cex=0.5);
  }
  if(length(genesNeg)>0) {
    segments(start(genesNeg), -0.175, end(genesNeg), negativeStrandPosition);
    text(start(genesNeg), -0.23, genesNeg$ID, pos=4, cex=0.5);
  }
  if(length(exonsPos)>0) {
    rect(start(exonsPos), -0.05, end(exonsPos), -0.09, col=col[1], border = NA);
  }
  if(length(exonsNeg)>0) {
    rect(start(exonsNeg), -0.16, end(exonsNeg), -0.2, col=col[1], border = NA);
  }
  if(length(transposonsPos)>0) {
    rect(start(transposonsPos), -0.05, end(transposonsPos), -0.09, col=col[2], border = col[2], density=30, angle=30);
    text(start(transposonsPos), -0.115, transposonsPos$ID, pos=4, cex=0.5);
  }
  if(length(transposonsNeg)>0) {
    rect(start(transposonsNeg), -0.16, end(transposonsNeg), -0.2, col=col[2], border = col[2], density=30, angle=30);
    text(start(transposonsNeg), -0.23, transposonsNeg$ID, pos=4, cex=0.5);
  }

}


#' This function plots the methylation profile at one locus for the bisulfite
#' sequencing data.The points on the graph represent methylation proportion of
#' individual cytosines, their colour which sample they belong to and the
#' intesity of the the colour how many reads that particular cytosine had. This
#' means that darker colors indicate stronger evidence that the corresponding
#' cytosine has the corresponding methylation proportion, while lighter colors
#' indicate a weaker evidence. The solid lines represent the smoothed profiles
#' and the intensity of the line the coverage at the corresponding position
#' (darker colors indicate more reads while lighter ones less reads). The boxes
#' on top represent the DMRs, where a filled box will represent a DMR which
#' gained methylation while a box with a pattern represent a DMR that lost
#' methylation. The DMRs need to have a metadafield \code{"regionType"} which
#' can be either \code{"gain"} (where there is more methylation in condition 2
#' compared to condition 1) or \code{"loss"} (where there is less methylation in
#' condition 2 compared to condition 1). In case this metadafield is missing all
#' DMRs are drawn using a filled box. Finally, we also allow annotation of the
#' DNA sequence. We represent by a black boxes all the exons, which are joined
#' by a horizontal black line, thus, marking the full body of the gene. With
#' grey boxes we mark the transposable elements. Both for genes and transposable
#' elements we plot them over a mid line if they are on the positive strand and
#' under the mid line if they are on the negative strand.
#'
#' @title Plot local methylation profile
#' @param methylationData1 the methylation data in condition 1
#' (see \code{\link{methylationDataList}}).
#' @param methylationData2 the methylation data in condition 2
#' (see \code{\link{methylationDataList}}).
#' @param region a \code{\link{GRanges}} object with the region where to plot
#' the high resolution profile.
#' @param DMRs a \code{\link{GRangesList}} object or a list with the list of
#' DMRs (see \code{\link{computeDMRs}} or \code{\link{filterDMRs}}.
#' @param conditionsNames the names of the two conditions. This will be used to
#' plot the legend.
#' @param gff a \code{\link{GRanges}} object with all elements usually imported
#' from a GFF3 file. The gff file needs to have an metafield \code{"type"}. Only
#' the elements of type  \code{"gene"}, \code{"exon"} and
#' \code{"transposable_element"} are plotted. Genes are represented as
#' horizontal black lines, exons as a black rectangle and transposable elements
#' as a grey rectangle. The elements are plotted on the corresponding strand
#' (\code{+} or \code{-}).
#' @param windowSize the size of the triangle base used to smooth the average
#' methylation profile.
#' @param context the context in which the DMRs are computed (\code{"CG"},
#' \code{"CHG"} or \code{"CHH"}).
#' @param labels a \code{vector} of \code{character} used to add a subfigure
#' characters to the plot. If \code{NULL} nothing is added.
#' @param col a \code{character} vector with the colors. It needs to contain a
#' minimum of \code{4 length(DMRs)} colors. If not or if \code{NULL}, the
#' defalut colors will be used.
#' @param main a \code{character} with the title of the plot
#' @param plotMeanLines a \code{logical} value indicating whether to plot the
#' mean lines or not.
#' @param plotPoints a \code{logical} value indicating whether to plot the
#' points or not.
#' @return Invisibly returns \code{NULL}
#' @examples
#'
#' # load the methylation data
#' data(methylationDataList)
#' # load the gene annotation data
#' data(GEs)
#'
#' #select the genes
#' genes <- GEs[which(GEs$type == "gene")]
#'
#' # the coordinates of the area to be plotted
#' chr3Reg <- GRanges(seqnames = Rle("Chr3"), ranges = IRanges(510000,530000))
#'
#' # load the DMRs in CG context
#' data(DMRsNoiseFilterCG)
#'
#' DMRsCGlist <- list("noise filter"=DMRsNoiseFilterCG)
#'
#'
#' # plot the CG methylation
#' par(mar=c(4, 4, 3, 1)+0.1)
#' par(mfrow=c(1,1))
#' plotLocalMethylationProfile(methylationDataList[["WT"]],
#'                            methylationDataList[["met1-3"]], chr3Reg,
#'                            DMRsCGlist, c("WT", "met1-3"), GEs,
#'                            windowSize=100, main="CG methylation")
#'
#' @author Nicolae Radu Zabet
#' @export
plotLocalMethylationProfile <- function(methylationData1,
                                        methylationData2,
                                        region,
                                        DMRs = NULL,
                                        conditionsNames = NULL,
                                        gff = NULL,
                                        windowSize = 150,
                                        context = "CG",
                                        labels=NULL,
                                        col=NULL,
                                        main="",
                                        plotMeanLines=TRUE,
                                        plotPoints=TRUE) {

  .validateMethylationData(methylationData1, variableName="methylationData1")
  .validateMethylationData(methylationData2, variableName="methylationData1")

  .stopIfNotAll(c(!is.null(plotMeanLines),
                  is.logical(plotMeanLines)),
                " plotMeanLines is a logical parameter and can be only TRUE or FALSE.")
  .stopIfNotAll(c(!is.null(plotPoints),
                  is.logical(plotPoints)),
                " plotPoints is a logical parameter and can be only TRUE or FALSE.")

  numberOfConditions <- 2

  if(is.null(conditionsNames) | length(conditionsNames) < numberOfConditions | !all(is.character(conditionsNames))){
    conditionsNames <- paste("condition ",(1:numberOfConditions),sep="")
  }

  .validateGRanges(region, length=1, generateGenomeWide=FALSE)

  .stopIfNotAll(c(.isInteger(windowSize, positive=TRUE)),
                " the window size used by the interpolation method is an integer higher or equal to 0")
  .validateContext(context)



  numberOfDMRs <- 0
  if(!is.null(DMRs)){
    numberOfDMRs <- length(DMRs)
  }

  if(!is.null(labels) & (length(labels) < 1 | !is.character(labels))){
    labels <- LETTERS[1:length(labels)]
  }

  if(!.isColor(col, minLength = (4+numberOfDMRs))){
    col <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

    cond1Color <- col[7]
    cond2Color <- col[6]
    geneColor <- col[1]
    TEColor <- col[9]
    DMRsColor <- col[c(4, 5, 8, 2, 3)]

  } else{
    cond1Color <- col[1]
    cond2Color <- col[2]
    geneColor <- col[3]
    TEColor <- col[4]
    DMRsColor <- col[5:length(col)]
  }

  seqname <- seqnames(region)
  minPos <- start(region)
  maxPos <- end(region)

  methylationData <- .joinMethylationData(methylationData1, methylationData2)

  hits <- findOverlaps(methylationData, region)

  localMethylationData <- methylationData[queryHits(hits)];
  contextMethylationData <- localMethylationData[localMethylationData$context%in%context];


  ramp1 <- colorRampPalette(c("white", cond1Color))
  colramp1 <- ramp1(100)

  ramp2 <- colorRampPalette(c("white", cond2Color))
  colramp2 <- ramp2(100)



  # plot the points
  proportion1 <- rep(0, length(contextMethylationData))
  index <- which(contextMethylationData$readsM1 >= 0 & contextMethylationData$readsN1 > 0)
  proportion1[index] <- contextMethylationData$readsM1[index] / contextMethylationData$readsN1[index]

  proportion2 <- rep(0, length(contextMethylationData))
  index <- which(contextMethylationData$readsM2 >= 0 & contextMethylationData$readsN2 > 0)
  proportion2[index] <- contextMethylationData$readsM2[index] / contextMethylationData$readsN2[index]

  maxColor <- max(c(contextMethylationData$readsN1[!is.na(contextMethylationData$readsN1)],
                    contextMethylationData$readsN2[!is.na(contextMethylationData$readsN2)]))



  plot(start(contextMethylationData), proportion1,
       col=colramp1[round(99 * log(contextMethylationData$readsN1) / log(maxColor))+1],
       pch = 16, cex=0.7,
       xlim = c(minPos, maxPos),
       ylim = c(-.2, 1.2 + numberOfDMRs*.1),
       xlab = paste("genomic coordinate on chromosome ",seqname, sep="") ,
       ylab = "methylation proportion",
       yaxt = "n",
       main = main,
       mar  = c(4,4,2,2)+.1,
       type = "n")
  axis(2,c(0,.5,1));


  if(plotPoints){
    points(start(contextMethylationData), proportion1,
           col=colramp1[round(99 * log(contextMethylationData$readsN1) / log(maxColor))+1],
           pch = 16, cex=0.7)
    points(start(contextMethylationData), proportion2,
           col=colramp2[round(99 * log(contextMethylationData$readsN2) / log(maxColor))+1],
           pch = 15, cex=0.7)
    if(!plotMeanLines){
      legend("topright", bty="n", col=c(cond1Color, cond2Color), legend=conditionsNames, pch=c(15,16), horiz=TRUE)
    }
  }

  #axis

  windowSizeHalf <- floor((windowSize - 1)/2)



  if(plotMeanLines){
    positions <- start(region) : end(region)


    movingAverageMethylReads1 <- .movingAverage(start(region),
                                                end(region),
                                                start(contextMethylationData),
                                                contextMethylationData$readsM1,
                                                windowSizeHalf = windowSizeHalf)
    movingAverageTotalReads1 <- .movingAverage(start(region),
                                               end(region),
                                               start(contextMethylationData),
                                               contextMethylationData$readsN1,
                                               windowSizeHalf = windowSizeHalf)
    movingAverageProportion1 <- movingAverageMethylReads1 / movingAverageTotalReads1

    movingAverageMethylReads2 <- .movingAverage(start(region),
                                                end(region),
                                                start(contextMethylationData),
                                                contextMethylationData$readsM2,
                                                windowSizeHalf = windowSizeHalf)
    movingAverageTotalReads2 <- .movingAverage(start(region),
                                               end(region),
                                               start(contextMethylationData),
                                               contextMethylationData$readsN2,
                                               windowSizeHalf = windowSizeHalf)
    movingAverageProportion2 <- movingAverageMethylReads2 / movingAverageTotalReads2

    maxColor <- max(c(movingAverageTotalReads1[!is.na(movingAverageTotalReads1)],
                      movingAverageTotalReads2[!is.na(movingAverageTotalReads2)]))


    #lines(positions, movingAverageProportion1, lty=1, lwd=2, col=cond1Color)
    #lines(positions, movingAverageProportion2, lty=1, lwd=2, col=cond2Color)
    #compute coordinates for mean lines condition 1
    bufferIndex <- !is.na(movingAverageProportion1) & !is.na(movingAverageTotalReads1)
    if(any(bufferIndex)){
      minID <- min(which(bufferIndex))
      maxID <- max(which(bufferIndex))
      colorId <- round(99 * log(round(movingAverageTotalReads1)) / log(maxColor))+1
      colorId <- round((colorId[minID:(maxID-1)] +  colorId[(minID+1):(maxID)])/2)
      colorId[is.na(colorId)] <- 1
      segments(positions[minID:(maxID-1)],
               movingAverageProportion1[minID:(maxID-1)],
               positions[(minID+1):(maxID)],
               movingAverageProportion1[(minID+1):(maxID)],
               col=colramp1[colorId],
               lty=1, lwd=2)


    }
    #compute coordinates for mean lines condition 2
    bufferIndex <- !is.na(movingAverageProportion2) & !is.na(movingAverageTotalReads2)
    if(any(bufferIndex)){
      minID <- min(which(bufferIndex))
      maxID <- max(which(bufferIndex))
      colorId <- round(99 * log(round(movingAverageTotalReads2)) / log(maxColor))+1
      colorId <- round((colorId[minID:(maxID-1)] +  colorId[(minID+1):(maxID)])/2)
      colorId[is.na(colorId)] <- 1
      segments(positions[minID:(maxID-1)],
               movingAverageProportion2[minID:(maxID-1)],
               positions[(minID+1):(maxID)],
               movingAverageProportion2[(minID+1):(maxID)],
               col=colramp2[colorId],
               lty=1, lwd=2)
    }

    legend("topright", bty="n", col=c(cond1Color, cond2Color), legend=conditionsNames, lty=1, lwd = 2, horiz=TRUE)
  }




  if(numberOfDMRs > 0) {
    for(i in 1:numberOfDMRs) {
      range <- DMRs[[i]][queryHits(findOverlaps(DMRs[[i]],region))]

      if(length(range) > 0 ){
        start(range) <- pmax(start(range), start(region))
        end(range)   <- pmin(end(range), end(region))

        if(length(which(range$regionType == "gain"))>0){
          rect(start(range)[range$regionType == "gain"], 1.025 + (length(DMRs) - i)*.1,
               end(range)[range$regionType == "gain" ], 1.075 + (length(DMRs)-i)*.1,
               col=DMRsColor[i], border=DMRsColor[i])
        }
        if(length(which(range$regionType == "loss"))>0){
          rect(start(range)[range$regionType == "loss"], 1.025 + (length(DMRs) - i)*.1,
               end(range)[range$regionType == "loss"], 1.075 + (length(DMRs)-i)*.1,
               col=DMRsColor[i], border=DMRsColor[i], density =30)
        }

        if(length(range) > 0 & (is.null(range$regionType))){
          rect(start(range), (1.025 + (length(DMRs) - i)*.1), end(range), (1.075 + (length(DMRs)-i)*.1),
               col=DMRsColor[i], border=DMRsColor[i])
        }

      }
      text(par('usr')[1], 1.050 + (length(DMRs)-i)*.1, pos=4, names(DMRs)[i])
    }
  }

  if(!is.null(gff)){
    if(is(gff, "GRanges")){
      .plotGeneticElements(gff,region, c(geneColor, TEColor))
    }
  }

  if(!is.null(labels)){
    if(length(labels) >= 1){
      mtext(labels[1], line = 0.7, adj = 0, cex=1.0);
    }
  }
  invisible(NULL)
}






