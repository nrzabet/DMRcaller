
#' This function extracts GC sites in the genome 
#'
#' @title Extract GC
#' @param methylationData the methylation data stored as a \code{\link{GRanges}}
#' object with four metadata columns (see \code{\link{methylationDataList}}).
#' 
#' @param genome a BSgenome with the DNA sequence of the organism
#' 
#' @param contexts the context in which the DMRs are computed (\code{"ALL"}, 
#' \code{"CG"}, \code{"CHG"} or \code{"CHH"}).
#' @return the a subset of \code{methylationData} consisting of all GC sites. 
#'
#' @examples
#' 
#' \dontrun{
#' # load the genome sequence
#' if(!require("BSgenome.Athaliana.TAIR.TAIR9", character.only = TRUE)){
#'   source("https://bioconductor.org/biocLite.R")
#'   biocLite("BSgenome.Athaliana.TAIR.TAIR9")
#' }
#' library(BSgenome.Athaliana.TAIR.TAIR9)
#' 
#' # load the methylation data
#' data(methylationDataList)
#' 
#' methylationDataWTGpCpG <- extractGC(methylationDataList[["WT"]], 
#'                                     BSgenome.Athaliana.TAIR.TAIR9,
#'                                     "CG")
#' 
#' }
#' 
#' @author Ryan Merritt
extractGC <- function(methylationData, 
                      genome, 
                      contexts=c("ALL","CG","CHG","CHH")){
  # Splitting data according to the context specified
  message("Splitting data according to context \n")
  contexts <- match.arg(contexts)
  switch(contexts,
         "ALL" = metData <- methylationData,
         "CG" = metData <- methylationData[which(methylationData$context=="CG")],
         "CHH" = metData <- methylationData[which(methylationData$context=="CHH")],
         "CHG" = metData <- methylationData[which(methylationData$context=="CHG")]
  )
  # Shifting strands to get the base before the cytosine
  message("Shifting Strands \n")
  start(metData[strand(metData)=="+"]) <- start(metData[strand(metData)=="+"])-1
  end(metData[strand(metData)=="-"]) <- end(metData[strand(metData)=="-"])+1
  # Removing cytosines that are at the ends of each chromosome
  message("Removing data outside of chromosome length ranges \n")
  chromlength <- width(getSeq(genome))
  names(chromlength) <- seqlevels(metData)
  chromlength <- chromlength[which(seqlevels(metData)%in%seqlevels(genome))]
  for(i in 1:length(chromlength)){
    metData <- metData[!(seqnames(metData)==names(chromlength)[i] &
                           end(metData)>chromlength[i])]
  }
  for(i in 1:length(chromlength)){
    metData <- metData[!(seqnames(metData)==names(chromlength)[i] &
                           start(metData)<=0)]
  }
  message("Extracting GC \n")
  metSeq <- getSeq(genome,metData)
  metDataGC <- which(metSeq =="GC")
  metDataGC <- metData[metDataGC]
  message("Unshifting strands \n")
  start(metDataGC[strand(metDataGC)=="+"]) <- start(metDataGC[strand(metDataGC)
                                                              =="+"])+1
  end(metDataGC[strand(metDataGC)=="-"]) <- end(metDataGC[strand(metDataGC)
                                                          =="-"])-1
  if(contexts=="ALL"){
    metDataGC <- split(metDataGC,metDataGC$context)
  }
  return(metDataGC)
}
