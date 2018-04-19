#' This function asses an expression and if this is FALSE stops the execution
#' and prints a message
#'
#' @title .stopIfNotAll
#' @param expr an array of logical expressions
#' @param errorMsg the error message to be printed
#'
#' @author Radu Zabet
.stopIfNotAll <- function(exprs, errorMsg) {
  for(expr in exprs){
    if (! expr)
      stop(errorMsg, call. = FALSE)
  }
}

#' This function checks is a variable is an integer number
#'
#' @title .is.integer
#' @param x the variable
#' @param positive a logical variable indicating whether to check if the number
#' is a positive integer
#' @return a logical number whether the variable x is an integer or not
#'
#' @author Radu Zabet
.isInteger <- function(x, positive = FALSE){
  isInteger <- TRUE
  if (is.null(x)){
    isInteger <- FALSE
  } else{
    if (!is.numeric(x)){
      isInteger <- FALSE
    } else{
      if(x%%1!=0){
        isInteger <- FALSE
      } else{
        if(positive & x < 0){
          isInteger <- FALSE
        }
      }
    }
  }
  return(isInteger)
}


#' This function that prints a GenomicRanges object
#'
#' @title Print GenomicRanges
#' @param gr the GenomicRanges object
#' @return a vector of type \code{character} with template chr:start..end
#'
#' @author Radu Zabet
.printGenomicRanges <- function(gr){

  .stopIfNotAll(c(!is.null(gr),
                  typeof(gr) == "S4",
                  class(gr)[1] == "GRanges"),
                " gr is a GenomicRanges object");
  result=c();
  for(index in 1:length(gr)){
    result <- c(result,
                paste(seqnames(gr)[index],":",start(gr)[index],"..",end(gr)[index],sep=""))
  }
  return(result)
}

#' This function checks whether the argument is a vector containing colours
#'
#' @title Is color
#' @param x the vector to be validated
#' @param minLength the  minimum length of the vector. If NULL the minimum
#' length is 1
#' @return a \code{logical} value indicating whether \code{x} is a vector
#' containing only colors
#'
#' @author Radu Zabet
.isColor <- function(x, minLength=NULL){
  isColor <- TRUE
  if(is.null(x)){
    isColor <- FALSE
  }

  if(is.null(minLength)){
    minLength = 1
  }

  if(isColor & length(x) < minLength){
    isColor <- FALSE
  }

  if(isColor){
    for(i in 1:length(x)){
      if(!(x[i]%in%colors()) & length(grep("^#[0-9A-Fa-f]{6}$", x[i])) < 1){
        isColor <- FALSE
      }
    }
  }
  return(isColor)
}


#' Returns a \code{\link{GRanges}} object spanning from the first cytocine until
#' the last one on each chromosome
#'
#' @title Get whole chromosomes from methylation data
#' @param methylationData the methylation data stored as a \code{\link{GRanges}}
#' object with four metadata columns (see \code{\link{methylationDataList}}).
#' @return a \code{\link{GRanges}} object will all chromosomes.
#'
#' @examples
#' # load the methylation data
#' data(methylationDataList)
#'
#' # get all chromosomes
#' chromosomes <- getWholeChromosomes(methylationDataList[["WT"]])
#'
#' @author Nicolae Radu Zabet
#' @export
getWholeChromosomes <- function(methylationData){
  max <- c()
  min <- c()
  seqnames <- c()
  for(chr in levels(seqnames(methylationData))){
    indexes=which(as.character(seqnames(methylationData)) == chr)
    if(length(indexes) > 0){
      seqnames <- c(seqnames, chr)
      max <- c(max, max(start(ranges(methylationData))[indexes]))
      min <- c(min, min(start(ranges(methylationData))[indexes]))
    }
  }

  regions <- GRanges(seqnames = Rle(seqnames), ranges   = IRanges(min,max))

  return(regions)
}

#' Checks whether the passed parameter has the correct format for methylation data
#'
#' @title Validate methylation data
#' @param methylationData the methylation data stored as a \code{\link{GRanges}}
#' object containing all the replicates.
#'
#' @author Alessandro Pio Greco and Nicolae Radu Zabet
.validateMethylationData <- function(methylationData, variableName="methylationData",
                                     manageDuplicates = "mean"){
  .stopIfNotAll(c(!is.null(methylationData),
                  typeof(methylationData) == "S4",
                  class(methylationData)[1] == "GRanges"),
                " methylationData needs to be a GRanges object")
  .stopIfNotAll(c(ncol(mcols(methylationData)[grepl("reads",
                  names(mcols(methylationData)))])%%2 == 0,
                  length(methylationData) > 0,
                  any(grepl("context", names(mcols(methylationData)))) == TRUE,
                  any(grepl("trinucleotide_context", names(mcols(methylationData)))) == TRUE),
                paste(" ",variableName," the object does not contain the correct metadata columns", sep=" "))
  if(any(duplicated(methylationData)) == TRUE){
    indexesDuplicated <- which(ranges(methylationData) ==
                               ranges(methylationData)[duplicated(methylationData)])
    checkMetadataEqual <- all(mcols(methylationData[indexesDuplicated[1:(length(indexesDuplicated)/2)]]) ==
                                mcols(methylationData[indexesDuplicated[((length(indexesDuplicated)/2)+1):
                                                              length(indexesDuplicated)]]))
    if(all(checkMetadataEqual) == TRUE){
      if(manageDuplicates == "mean"){
      cat("Cytosines that were duplicated and had the different metadata columns were merged by meaning readings \n")


    } else if(manageDuplicates == "sum"){


    } else if(manageDuplicates == "discard"){
      cat("Cytosines that were duplicated (",indexesDuplicated, ") and had the same metadata columns were discarded \n", sep = " ")
      methylationData <- unique(methylationData)
    }


    } else{
      stop(" context or trinucleotide context on duplicated cytosines (", indexesDuplicated, ") are not equal")
    }

  }
}


#' Checks whether the passed parameter has the correct format for methylation data
#'
#' @title Validate methylation data
#' @param methylationData the methylation data stored as a \code{\link{GRanges}}
#' object with six metadata columns (see \code{\link{methylationData}}).
#'
#' @author Radu Zabet
.validateMethylationDataList <- function(methylationDataList){
  .stopIfNotAll(c(!is.null(methylationDataList),
                  typeof(methylationDataList) == "S4",
                  class(methylationDataList)[1] == "GRangesList",
                  length(methylationDataList) > 0),
                " methylationDataList needs to be a GRangesList object")

  for(i in 1:length(methylationDataList)){
    .stopIfNotAll(c(.validateMethylationData(methylationDataList[[i]])),
                 paste(" element ", i," of the methylationDataList is incorrect", sep=""))
  }

}

#' Checks whether the passed parameter is the context
#'
#' @title Validate context
#' @param context the context in which the DMRs are computed (\code{"CG"},
#' \code{"CHG"} or \code{"CHH"})
#' @param length the expected length of the vector. If NULL any length is
#' allowed
#'
#' @author Radu Zabet
.validateContext <- function(context, length=NULL){
  .stopIfNotAll(c(!is.null(context), all(is.character(context)),
                  all(context %in% c("CG","CHG","CHH"))),
                " context can be only CG,CHG or CHH")
  if(!is.null(length)){
    .stopIfNotAll(c(is.numeric(length), length(context) == length),
                  paste(" context needs to contain exactly ", length," elements", sep=""))
  }
}



#' Checks whether the passed parameter is the statistical test
#'
#' @title Validate statistial test
#' @param test the statistical test used to call DMRs (\code{"fisher"} for
#' Fisher's exact test or \code{"score"} for Score test).
#'
#' @author Radu Zabet
.validateStatisticalTest <- function(test){
  .stopIfNotAll(c(!is.null(test), is.character(test), length(test) == 1, test %in% c("fisher","score")),
                " test needs to be one of the following \"fisher\" for Fisher's exact test or \"score\" for Score test")
}


#' Checks whether the passed parameter is a \code{\link{GRanges}} object
#'
#' @title Validate GRanges
#' @param regions a \code{\link{GRanges}} object. If \code{NULL} and
#' \code{generateGenomeWide} is true it uses the \code{methylationData} to
#' compute the regions and returns this \code{\link{GRanges}} object
#' @param methylationData the methylation data stored as a \code{\link{GRanges}}
#' object with six metadata columns (see \code{\link{methylationData}}).
#' @param length the expected length of the vector. If \code{NULL} any length is
#' allowed.
#' @param generateGenomeWide logical value to indicate whether to compute the
#' regions that span over all the \code{methylationData}
#' @return a \code{\link{GRanges}} object representing the regions
#'
#' @author Radu Zabet
.validateGRanges <- function(regions, methylationData, length=NULL, generateGenomeWide=TRUE, variableName="regions", minLength=0){

  if(is.null(regions) & generateGenomeWide){
    if(typeof(methylationData) == "S4" & class(methylationData)[1] == "GRangesList" & length(methylationData) > 0){
      regions <- NULL
      for(i in 1:length(methylationData)){
        if(is.null(regions)){
          regions <- getWholeChromosomes(methylationData[[i]])
        } else{
          regions <- union(regions, getWholeChromosomes(methylationData[[i]]))
        }
      }
    } else if(typeof(methylationData) == "S4" & class(methylationData)[1] == "GRanges"){
      regions <- getWholeChromosomes(methylationData)
    }
  }

  .stopIfNotAll(c(!is.null(regions),
                  typeof(regions) == "S4",
                  class(regions)[1] == "GRanges"),
                paste(" ",variableName," neads to be a GRanges object. If NULL, the DMRs are computed genome-wide.",sep=""))

  if(!is.null(length)){
    .stopIfNotAll(c(is.numeric(length), length(regions) == length),
                  paste(" ",variableName," needs to contain exactly ", length," elements", sep=""))
  }
  if(!is.null(minLength)){
    .stopIfNotAll(c(is.numeric(minLength), length(regions) >= minLength),
                  paste(" ",variableName," needs to contain more than ", minLength," elements", sep=""))
  }
  return(regions)
}

#' Checks whether the passed parameter is a \code{GRangesList} object that
#' represents the methylation profile
#'
#' @title Validate methylation profile
#' @param methylationProfile a \code{GRangesList} object
#' (see \code{\link{computeMethylationProfile}}).
#'
#' @author Radu Zabet
.validateMethylationProfile <- function(methylationProfile){
  .stopIfNotAll(c(!is.null(methylationProfile),
                  typeof(methylationProfile) == "S4",
                  class(methylationProfile)[1] == "GRangesList"),
                " methylationProfile needs to be a GRangesList")


  for(i in 1:length(methylationProfile)){
    .stopIfNotAll(c(!is.null(methylationProfile[[i]]),
                    typeof(methylationProfile[[i]]) == "S4",
                    class(methylationProfile[[i]])[1] == "GRanges"),
                  paste(" element ",i," of the methylationProfile is not a GRanges object", sep=""))
    .stopIfNotAll(c(ncol(mcols(methylationProfile[[i]])) == 4,
                    length(methylationProfile[[i]]) > 0),
                  paste(" element ",i," of the methylationProfile is not a GRanges object with four metadata columns (see computeMethylationProfile function).", sep=""))
  }

}
