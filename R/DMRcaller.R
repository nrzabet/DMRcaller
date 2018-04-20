#' Call Differentially Methylated Regions (DMRs) between two samples
#'
#' Uses bisulfite sequencing data in two conditions and identifies
#' differentially methylated regions between the conditions in CG and
#' non-CG context. The input is the CX report files produced by Bismark
#' and the output is a list of DMRs stored as GRanges objects.
#'
#' The most important functions in the \pkg{DMRcaller} are:
#' \describe{
#'  \item{\code{\link{readBismark}}}{reads the Bismark CX report files in a
#'        \code{\link{GRanges}} object.}
#'  \item{\code{\link{readBismarkPool}}}{Reads multiple CX report files and
#'        pools them together.}
#'  \item{\code{\link{saveBismark}}}{saves the methylation data stored in a
#'        \code{\link{GRanges}} object into a Bismark CX report file.}
#'  \item{\code{\link{poolMethylationDatasets}}}{pools together multiple
#'        methylation datasets.}
#'  \item{\code{\link{poolTwoMethylationDatasets}}}{pools together two
#'        methylation datasets.}
#'  \item{\code{\link{computeMethylationDataCoverage}}}{Computes the coverage
#'        for the bisulfite sequencing data.}
#'  \item{\code{\link{plotMethylationDataCoverage}}}{Plots the coverage for the
#'        bisulfite sequencing data.}
#'  \item{\code{\link{computeMethylationDataSpatialCorrelation}}}{Computes the
#'        correlation between methylation levels  as a function of the distances
#'        between the Cytosines.}
#'  \item{\code{\link{plotMethylationDataSpatialCorrelation}}}{Plots the
#'        correlation of methylation levels for Cytosines located at a certain
#'        distance apart.}
#'  \item{\code{\link{computeMethylationProfile}}}{Computes the low resolution
#'        profiles for the bisulfite sequencing data at certain locations.}
#'  \item{\code{\link{plotMethylationProfile}}}{Plots the low resolution
#'        profiles for the bisulfite sequencing data at certain locations.}
#'  \item{\code{\link{plotMethylationProfileFromData}}}{Plots the low resolution
#'        profiles for the loaded bisulfite sequencing data.}
#'  \item{\code{\link{computeDMRs}}}{Computes the differentially methylated
#'        regions between two conditions.}
#'  \item{\code{\link{filterDMRs}}}{Filters a list of (potential) differentially
#'        methylated regions.}
#'  \item{\code{\link{mergeDMRsIteratively}}}{Merge DMRs iteratively.}
#'  \item{\code{\link{analyseReadsInsideRegionsForCondition}}}{Analyse reads
#'        inside regions for condition.}
#'  \item{\code{\link{plotLocalMethylationProfile}}}{Plots the methylation
#'        profile at one locus for the bisulfite sequencing data.}
#'  \item{\code{\link{computeOverlapProfile}}}{Computes the distribution of a
#'        set of subregions on a large region.}
#'  \item{\code{\link{plotOverlapProfile}}}{Plots the distribution of a set of
#'        subregions on a large region.}
#'  \item{\code{\link{getWholeChromosomes}}}{Computes the GRanges objects with
#'        each chromosome as an element from the methylationData.}
#'  \item{\code{\link{joinReplicates}}}{Merges two GRanges objects with single
#'        reads columns. It is necessary to start the analysis of DMRs with
#'        biological replicates.}
#'  \item{\code{\link{computeDMRsReplicates}}}{Computes the differentially
#'                    methylated regions between two conditions with multiple
#'                    biological replicates.}
#' }
#'
#' @author
#' Nicolae Radu Zabet \email{n.r.zabet@@gen.cam.ac.uk},
#' Jonathan Michael Foonlan Tsang \email{jmft2@@cam.ac.uk}
#' Alessandro Pio Greco \email{apgrec@@essex.ac.uk}
#'
#' Maintainer: Nicolae Radu Zabet \email{n.r.zabet@@gen.cam.ac.uk}
#' @name DMRcaller
#' @docType package
#' @import GenomicRanges
#' @import IRanges
#' @import S4Vectors
#' @import Rcpp
#' @importFrom RcppRoll roll_sum
#' @importFrom parallel mclapply
#' @examples
#' \dontrun{
#' # load the methylation data
#' data(methylationDataList)
#'
#' #plot the low resolution profile at 5 Kb resolution
#' par(mar=c(4, 4, 3, 1)+0.1)
#' plotMethylationProfileFromData(methylationDataList[["WT"]],
#'                                methylationDataList[["met1-3"]],
#'                                conditionsNames=c("WT", "met1-3"),
#'                                windowSize = 5000, autoscale = TRUE,
#'                                context = c("CG", "CHG", "CHH"),
#'                                labels = LETTERS)
#'
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
#'
#' # plot the coverage in all three contexts
#' plotMethylationDataCoverage(methylationDataList[["WT"]],
#'                            methylationDataList[["met1-3"]],
#'                            breaks = 1:15, regions = NULL,
#'                            conditionsNames = c("WT","met1-3"),
#'                            context = c("CG", "CHG", "CHH"),
#'                            proportion = TRUE, labels = LETTERS, col = NULL,
#'                            pch = c(1,0,16,2,15,17), lty = c(4,1,3,2,6,5),
#'                            contextPerRow = FALSE)
#'
#' # plot the correlation of methylation levels as a function of distance
#' plotMethylationDataSpatialCorrelation(methylationDataList[["WT"]],
#'                            distances = c(1,15), regions = NULL,
#'                            conditionsNames = c("WT","met1-3"),
#'                            context = c("CG"),
#'                            labels = LETTERS, col = NULL,
#'                            pch = c(1,0,16,2,15,17), lty = c(4,1,3,2,6,5),
#'                            contextPerRow = FALSE)
#'
#' # the regions where to compute the DMRs
#' regions <- GRanges(seqnames = Rle("Chr3"), ranges = IRanges(1,1E6))
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
#' # load the gene annotation data
#' data(GEs)
#'
#' #select the genes
#' genes <- GEs[which(GEs$type == "gene")]
#'
#' # the regions where to compute the DMRs
#' genes <- genes[overlapsAny(genes, regions)]
#'
#' # filter genes that are differntially methylated in the two conditions
#' DMRsGenesCG <- filterDMRs(methylationDataList[["WT"]],
#'                methylationDataList[["met1-3"]], potentialDMRs = genes,
#'                context = "CG", test = "score", pValueThreshold = 0.01,
#'                minCytosinesCount = 4, minProportionDifference = 0.4,
#'                minReadsPerCytosine = 3, cores = 1)
#'
#' #merge the DMRs
#' DMRsNoiseFilterCGLarger <- mergeDMRsIteratively(DMRsNoiseFilterCG,
#'                            minGap = 500, respectSigns = TRUE,
#'                            methylationDataList[["WT"]],
#'                            methylationDataList[["met1-3"]],
#'                            context = "CG", minProportionDifference=0.4,
#'                            minReadsPerCytosine = 1, pValueThreshold=0.01,
#'                            test="score",alternative = "two.sided")
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
#' DMRsCGlist <- list("noise filter"=DMRsNoiseFilterCG,
#'                    "neighbourhood"=DMRsNeighbourhoodCG,
#'                    "bins"=DMRsBinsCG,
#'                    "genes"=DMRsGenesCG)
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
#'
#' hotspotsHypo <- computeOverlapProfile(
#'                DMRsNoiseFilterCG[(DMRsNoiseFilterCG$regionType == "loss")],
#'                region, windowSize=2000, binary=TRUE, cores=1)
#'
#' hotspotsHyper <- computeOverlapProfile(
#'                DMRsNoiseFilterCG[(DMRsNoiseFilterCG$regionType == "gain")],
#'                region, windowSize=2000, binary=TRUE, cores=1)
#'
#' plotOverlapProfile(GRangesList("Chr3"=hotspotsHypo),
#'                    GRangesList("Chr3"=hotspotsHyper),
#'                    names=c("loss", "gain"), title="CG methylation")
#' }
#'
#' # loading synthetic data
#' data("syntheticDataReplicates")
#'
#' # creating condition vector
#' condition <- c("a", "a", "b", "b")
#'
#' # computing DMRs using the neighbourhood method
#' DMRsReplicatesBins <- computeDMRsReplicates(methylationData = methylationData,
#'                                                      condition = condition,
#'                                                      regions = NULL,
#'                                                      context = "CHH",
#'                                                      method = "bins",
#'                                                      binSize = 100,
#'                                                      test = "betareg",
#'                                                      pseudocountM = 1,
#'                                                      pseudocountN = 2,
#'                                                      pValueThreshold = 0.01,
#'                                                      minCytosinesCount = 4,
#'                                                      minProportionDifference = 0.4,
#'                                                      minGap = 200,
#'                                                      minSize = 50,
#'                                                      minReadsPerCytosine = 4,
#'                                                      cores = 1)
#'
#' @seealso See \code{vignette("rd", package = "DMRcaller")} for an overview
#' of the package.
NULL


#' @name syntheticDataReplicates
#' @title Simulated data for biological replicates
#' @docType data
#' @description
#' A \code{GRanges} object containing simulated date for methylation in four
#' samples. The conditions assciated witch each sample are a, a, b and b.
#'
#' @format A \code{\link{GRanges}} object containing multiple metadata
#' columns with the reads from each object passed as parameter
#'
#' @source The object was created by calling \code{\link{joinReplicates}}
#'  function.
NULL

#' @name GEs
#' @title The genetic elements data
#' @docType data
#' @description
#' A \code{GRanges} object containing the annotation of the Arabidopsis thaliana
#'
#' @format A \code{GRanges} object
#'
#' @source The object was created by calling \code{import.gff3} function
#' from \code{rtracklayer} package for
#' \url{ftp://ftp.arabidopsis.org/Maps/gbrowse_data/TAIR10/TAIR10_GFF3_genes_transposons.gff}
NULL


#' @name methylationDataList
#' @title The methylation data list
#' @description
#' A \code{GRangesList} object containing the methylation data at each cytosine
#' location in the genome in  Wild Type (WT) and  met1-3 mutant (met1-3) in
#' Arabidopsis thaliana. The data only contains the first 1 Mbp from Chromosome 3.
#'
#' @format The \code{GRanges} elements contain four metadata columns
#' \describe{
#'  \item{context}{the context in which the DMRs are computed (\code{"CG"},
#'                \code{"CHG"} or \code{"CHH"}).}
#'  \item{readsM}{the number of methylated reads.}
#'  \item{readsN}{the total number of reads.}
#'  \item{trinucleotide_context}{the specific context of the cytosine (H is
#'                              replaced by the actual nucleotide).}
#' }
#'
#' @source Each element was created by by calling \code{\link{readBismark}}
#' function on the CX report files generated by Bismark
#' \url{http://www.bioinformatics.babraham.ac.uk/projects/bismark/}
#' for \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM980986} dataset
#' in the case of Wild Type (WT) and
#' \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM981032}
#' dataset in the case of met1-3 mutant (met1-3).
NULL


#' @name DMRsNoiseFilterCG
#' @title The DMRs between WT and met1-3 in CG context
#' @description
#' A \code{GRangesList} object containing the DMRs between  Wild Type (WT) and
#' met1-3 mutant (met1-3) in Arabidopsis thaliana
#' (see \code{\link{methylationDataList}}). The DMRs were computed on the first
#' 1 Mbp from Chromosome 3 with noise filter method using a triangular kernel
#' and a windowSize of 100 bp
#'
#' @format The \code{\link{GRanges}} element contain 11 metadata columns;
#' see \code{\link{computeDMRs}}
#' @seealso \code{\link{filterDMRs}}, \code{\link{computeDMRs}},
#' \code{\link{analyseReadsInsideRegionsForCondition}}
#' and \code{\link{mergeDMRsIteratively}}
#'
NULL
