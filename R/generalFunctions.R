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
                  is(gr, "GRanges")),
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
                  is(methylationData, "GRanges")),
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
                is(methylationDataList, "GRangesList"),
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
    if(is(methylationData, "GRangesList") & length(methylationData) > 0){
      regions <- NULL
      for(i in 1:length(methylationData)){
        if(is.null(regions)){
          regions <- getWholeChromosomes(methylationData[[i]])
        } else{
          regions <- union(regions, getWholeChromosomes(methylationData[[i]]))
        }
      }
    } else if(is(methylationData, "GRanges")){
      regions <- getWholeChromosomes(methylationData)
    }
  }

  .stopIfNotAll(c(!is.null(regions),
                  is(regions, "GRanges")),
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
                is(methylationProfile, "GRangesList")),
                " methylationProfile needs to be a GRangesList")


  for(i in 1:length(methylationProfile)){
    .stopIfNotAll(c(!is.null(methylationProfile[[i]]),
                  is(methylationProfile[[i]], "GRanges")),
                  paste(" element ",i," of the methylationProfile is not a GRanges object", sep=""))
    .stopIfNotAll(c(ncol(mcols(methylationProfile[[i]])) == 5,
                    length(methylationProfile[[i]]) > 0),
                  paste(" element ",i," of the methylationProfile is not a GRanges object with five metadata columns (see computeMethylationProfile function).", sep=""))
  }

}

#' Checks whether the passed parameter is the statistical test for detecting Co-methylation
#'
#' @title Validate statistial test for Co-methylation
#' @param test the statistical test used to call Co-methylation (\code{"fisher"} for
#' Fisher's exact test, \code{"score"} for Score test or \code{"permutation"} for Permutation test).
#' @keywords internal
#' @author Radu Zabet and Young Jun Kim
.validateCoMethylationStatTest <- function(test){
  .stopIfNotAll(c(!is.null(test), is.character(test), length(test) == 1, test %in% c("fisher","score","permutation","binom")),
                " Test needs to be one of the following \"fisher\" for Fisher's exact test, \"score\" for Score test or \"permutation\" for Permutation test" )
}

#' Checks whether the passed parameter is the alternative for detecting Co-methylation
#'
#' @title Validate alternative value for Co-methylation
#' @param test the alternative valuet used to Fisher exact's test must be one of
#' \code{"two.sided"}, \code{"greater"} or \code{"less"}. You can specify just
#' the initial letter.
#' @keywords internal
#' @author Radu Zabet and Young Jun Kim
.validateCoMethylationAlternative <- function(alternative){
  .stopIfNotAll(c(!is.null(alternative), is.character(alternative), length(alternative) == 1, alternative %in% c("two.sided","greater","less")),
                " Alternative needs to be one of the following \"two.sided\", \"greater\" or \"less\" for Fisher's exact test" )
}

#' Checks whether the passed parameter is the modified context 
#'
#' @title Validate modified context 
#' @param test the  modified context used to call methylation information from bamfile must be one of
#' sequence context for \code{selectCytosine()} (e.g. \code{"CG"}, \code{"CHG"}, \code{"CHH"})
#' @keywords internal
#' @author Radu Zabet and Young Jun Kim
.validateContext <- function(context){
  .stopIfNotAll(c(!is.null(context), is.character(context), context %in% c("CG","CHG","CHH")),
                " Modified context needs to be one of the following \"CG\", \"CHG\" or \"CHH\"" )
}


#' @title Validate BAM Filename
#' @description
#' Checks that `bamfile` exists and has a “.bam” extension.
#'
#' @param bamfile Character scalar. Path to a BAM file.
#' @return Invisibly `TRUE` if file exists and extension is “.bam”; otherwise errors.
#' @keywords internal
#' @author Radu Zabet and Young Jun Kim
.validateBamfile <- function(bamfile){
  # stop if bamfile directory is not include the .bam file format  
  .stopIfNotAll(c(!is.null(bamfile),is.character(bamfile), file.exists(bamfile), endsWith(tolower(bamfile), ".bam")),
                "bamfile directory name should need .bam extension")
}

#' @title Validate BSgenome Object
#' @description
#' Ensures `genome` inherits from BSgenome and that its package is installed.
#'
#' @param genome A BSgenome object (e.g. BSgenome.Hsapiens.UCSC.hg38).
#' @return Invisibly `TRUE` if valid; otherwise throws an error.
#' @keywords internal
#' @author Radu Zabet and Young Jun Kim
.validateGenome <- function(genome) {
  # 1) must be a BSgenome object
  if (!inherits(genome, "BSgenome")) {
    stop(
      "`genome` must be supplied as a BSgenome object, not as a string literal;\n",
      "please call e.g. `.validateGenome(BSgenome.Hsapiens.UCSC.hg38)` without quotes.",
      call. = FALSE
    )
  }
  
  ## 2) extract its package name
  pkg <- genome@pkgname
  
  ## 3) make sure BSgenome machinery is there
  if (!requireNamespace("BSgenome", quietly=TRUE)) {
    stop("please install the BSgenome package before specifying a genome", call. = FALSE)
  }
  
  ## 4) get the official list of BSgenome packages
  avail <- tryCatch(
    {suppressMessages(suppressPackageStartupMessages(BSgenome::available.genomes()))},
    error = function(e) NULL
  )
  
  ## 5) fallback: any installed BSgenome.* packages
  if (is.null(avail) || !is.character(avail)) {
    avail <- grep("^BSgenome\\.", rownames(utils::installed.packages()), value=TRUE)
  }
  
  ## 6) error if this genome isn’t known
  if (!pkg %in% avail) {
    stop(
      sprintf(
        "Loaded BSgenome object has package name '%s',\n  but that package is not among the installed/available BSgenome.* packages.\n  Known packages: %s%s",
        pkg,
        paste(head(avail, 3), collapse=", "),
        if (length(avail)>3) ", …" else ""
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @title Validate Chromosome Names
#' @description
#' Checks that `chr` entries are present in the supplied BSgenome.
#'
#' @param chr Character vector of chromosome names.
#' @param genome A validated BSgenome object.
#' @return Invisibly `TRUE` if all `chr` are in `seqnames(genome)`; otherwise errors.
#' @keywords internal
#' @author Radu Zabet and Young Jun Kim
.validateChromosome <- function(chr, genome){
  # check the Chromosome is included in the genome
  .stopIfNotAll(c(is.character(chr), chr %in% seqnames(genome)),
                " Chromosome should select same code as the BSgenome provided" )
}

#' @title Validate Modification Tag in BAM
#' @description
#' Scans the BAM for MM/ML tags in the specified region (or whole file)
#' and checks that `modif` is one of the observed modification codes.
#'
#' @param modif Character scalar. The MM tag code to validate (e.g. “C+m?”).
#' @param bamfile Character path to BAM (must exist plus .bai index).
#' @param chr Optional character vector of chromosome(s) to restrict scan.
#' @param genome A BSgenome, used to determine seqlengths if `chr` is set.
#' @return Invisibly `TRUE` if `modif` is found among the MM tags; otherwise errors.
#' @author Radu Zabet and Young Jun Kim
.validateModif <- function(modif, bamfile, chr, genome){
  # stop if modif is not included in the tags in the bam
  if (!file.exists(paste0(bamfile,".bai")))
    Rsamtools::indexBam(bamfile)
  if (!is.null(chr) && length(chr) > 0) {
    which_region <- GenomicRanges::GRanges(seqnames = chr,
      ranges = IRanges::IRanges(1, seqlengths(genome)[chr]))
    param <- Rsamtools::ScanBamParam(tag = c("MM","ML"), which = which_region)
  } else {
    param <- Rsamtools::ScanBamParam(tag = c("MM","ML"))
  }
  
  bam <- Rsamtools::scanBam(bamfile, param = param)[[1]]
  bam_tag <- unique(bam$tag$MM)
  bam_tag <- bam_tag[!is.na(bam_tag)]
  
  if (length(bam_tag) == 0) {
    stop("No MM tag found in the BAM file for the selected region/chromosome.")
  }
  
  # flatten all tags
  bam_tag_split <- unlist(strsplit(as.character(bam_tag), ";", fixed = TRUE))
  bam_tag_split <- bam_tag_split[nzchar(bam_tag_split)]  # remove empty strings
  

  code_list <- character()
  for (i in seq_along(bam_tag_split)) {
    parts <- strsplit(bam_tag_split[[i]], ",", fixed = TRUE)[[1]]
    code_list <- c(code_list, parts[1])
  }
  
  codes <- unique(code_list)
  .stopIfNotAll(c(!is.null(modif), is.character(modif), length(modif) == 1, modif %in% codes),
                c("Modified context must be one of: ", paste(shQuote(codes), collapse=", " )))
}

#' @title Validate or Auto‐Select a BiocParallelParam
#' @description
#' If `BPPARAM` is `NULL`, picks `SnowParam` or `MulticoreParam` based on env vars
#' or OS; otherwise checks class and fork compatibility.
#'
#' @param BPPARAM A \code{BiocParallelParam} or `NULL`.
#' @param progressbar Logical; whether to show a progress bar.
#' @return A valid `BiocParallelParam` object.
#' @keywords internal
#' @author Radu Zabet and Young Jun Kim
.validateBPPARAM <- function(BPPARAM = NULL, progressbar = FALSE) {
  if (is.null(BPPARAM)) {
    BPPARAM <- .chooseBPPARAM(progressbar = progressbar)
    return(BPPARAM)
  }
  # must inherit from BiocParallelParam
  if (!inherits(BPPARAM, "BiocParallelParam")) {
    stop("`BPPARAM` must be a BiocParallelParam object, got: ",
         paste(class(BPPARAM), collapse = "/"))
  }
  # no fork on Windows
  if (.Platform$OS.type == "windows" &&
      methods::is(BPPARAM, "MulticoreParam")) {
    stop("MulticoreParam() (fork) not supported on Windows; please use SnowParam().")
  }
  invisible(BPPARAM)
}

#' @title Auto‐choose BiocParallelParam Backend
#' @description
#' 1) If `NSLOTS` or `SLURM_CPUS_ON_NODE` are set, uses `SnowParam()` with that many workers.  
#' 2) On Windows, `SnowParam(detectCores())`.  
#' 3) Otherwise on Unix, `MulticoreParam(detectCores())`.
#'
#' @param workers Optional integer to override detected cores.
#' @param progressbar Logical; show progress bar.
#' @param cluster_type Type passed to `SnowParam()`, defaults to "SOCK".
#' @param ... Further args passed to `SnowParam` or `MulticoreParam`.
#' @return A `SnowParam` or `MulticoreParam` configured with the desired workers.
#' @keywords internal
#' @author Radu Zabet and Young Jun Kim
.chooseBPPARAM <- function(workers = NULL, progressbar = FALSE, cluster_type = "SOCK", ...) {
  # 1) HPC scheduler slots
  nslots <- as.integer(Sys.getenv("NSLOTS", unset = NA))
  if (!is.na(nslots) && nslots > 1) {
    cat("Detected SGE/PBS (NSLOTS=", nslots, ") -> using SnowParam()")
    return(SnowParam(workers = nslots, type = cluster_type,
                     progressbar = progressbar, ...))
  }
  # SLURM
  cpus <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE", unset = NA))
  if (!is.na(cpus) && cpus > 1) {
    cat("Detected SLURM (SLURM_CPUS_ON_NODE=", cpus, ") -> using SnowParam()")
    return(SnowParam(workers = cpus, type = cluster_type,
                     progressbar = progressbar, ...))
  }
  
  # 2) local fallback
  os <- .Platform$OS.type
  ncore <- if (is.null(workers)) parallel::detectCores(logical = FALSE) else workers
  if (os == "windows") {
    cat("Detected Windows -> using SnowParam() with ", ncore, " workers")
    return(SnowParam(workers = ncore, progressbar = progressbar, ...))
  }
  # 3) Unix
  cat("Local Unix -> using MulticoreParam() with forking and ", ncore, " workers")
  MulticoreParam(workers = ncore, progressbar = progressbar, ...)
}