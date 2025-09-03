#' Call Differentially Methylated Regions (DMRs) between two samples
#'
#' Supports both Bisulfite Sequencing (Bismark CX reports) and Oxford Nanopore
#' Sequencing (MM/ML tags) for per-site methylation calling. Identifies
#' differentially methylated regions between two samples in CG and non-CG
#' contexts.
#'
#' For bisulfite data, the input is Bismark CX report files and the output
#' is a list of DMRs stored as a \code{GRanges} object.
#' \itemize{
#'   \item \code{readsM} — count of modified reads per site
#'   \item \code{readsN} — total same-strand coverage per site
#' }
#' 
#' For Nanopore data, the input is an indexed ONT BAM with calling 
#' \code{\link{readONTbam}} function in this package and the output is a 
#' \code{GRanges} augmented with metadata columns:
#' \itemize{
#'   \item \code{ONT_Cm} — comma-delimited read-indices called “modified”
#'   \item \code{ONT_C}  — comma-delimited read-indices covering but not modified
#'   \item \code{readsM} — count of modified reads per site
#'   \item \code{readsN} — total same-strand coverage per site
#' }
#'
#' The most important functions in the \pkg{DMRcaller} are:
#' \describe{
#'  \item{\code{\link{readBismark}}}{reads the Bismark CX report files in a
#'        \code{\link{GRanges}} object.}
#'  \item{\code{\link{readBismarkPool}}}{Reads multiple CX report files and
#'        pools them together.}
#'  \item{\code{\link{saveBismark}}}{saves the methylation data stored in a
#'        \code{\link{GRanges}} object into a Bismark CX report file.}
#'  \item{\code{\link{selectCytosine}}}{Enumerates cytosine positions in a 
#'        BSgenome reference, optionally filtering by methylation context 
#'        (CG/CHG/CHH), chromosome and genomic region.}
#'  \item{\code{\link{readONTbam}}}{Loads an Oxford Nanopore BAM 
#'        (with MM/ML tags), decodes per‐C modification probabilities and 
#'        counts modified vs. unmodified reads per site.}
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
#'        methylated regions between two conditions with multiple
#'        biological replicates.}
#'  \item{\code{\link{selectCytosine}}}{Enumerates cytosine positions in a 
#'        BSgenome reference.}  
#'  \item{\code{\link{readONTbam}}}{Loads an ONT BAM (MM/ML tags), 
#'        decodes per‐C modification probabilities,  
#'        and counts modified vs. unmodified reads per site.}  
#'  \item{\code{\link{computePMDs}}}{Partitions the genome into PMDs via 
#'        three methods ("noise_filter", "neighbourhood", "bins").}
#'  \item{\code{\link{filterPMDs}}}{Filters a set of PMDs by methylation level 
#'        and read depth.}
#'  \item{\code{\link{mergePMDsIteratively}}}{Merge PMDs while preserving 
#'        statistical significance.}
#'  \item{\code{\link{analyseReadsInsideRegionsForConditionPMD}}}{Counts reads 
#'        in each PMD for one condition.}
#'  \item{\code{\link{computeCoMethylation}}}{Computes pairwise co‐methylation 
#'        between CpGs within regions.}
#' }
#'
#' @author
#' Nicolae Radu Zabet \email{n.r.zabet@@gen.cam.ac.uk},
#' Jonathan Michael Foonlan Tsang \email{jmft2@@cam.ac.uk},
#' Alessandro Pio Greco \email{apgrec@@essex.ac.uk},
#' Young Jun Kim \email{qc25039@@qmul.ac.uk}
#'
#' Maintainer: Nicolae Radu Zabet \email{n.r.zabet@@gen.cam.ac.uk}
#' @name DMRcaller
#' @docType package
#' @import GenomicRanges
#' @import IRanges
#' @import S4Vectors
#' @import Rcpp
#' @import Rsamtools
#' @importFrom RcppRoll roll_sum
#' @importFrom parallel mclapply
#' @useDynLib DMRcaller, .registration = TRUE
#' @examples
#' \dontrun{
#' # load the methylation data
#' data(methylationDataList)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' 
#' # All cytosines in hg38:
#' gr_all <- selectCytosine()
#'
#' # Only CpG sites on chr1 and chr2:
#' gr_chr1_2 <- selectCytosine(context="CG", chr=c("chr1","chr2"))
#'
#' # CHH sites in a specific region on chr3:
#' my_region <- GRanges("chr3", IRanges(1e6, 1e6 + 1e5))
#' gr_region <- selectCytosine(context="CHH", chr="chr3", region=my_region)
#' 
#' # set the bam file directory
#' bam_path <- system.file("extdata", "scanBamChr1Random5.bam", package="DMRcaller")
#' 
#' # read ONTbam file (chromosome 1 only) in CG context with BSgenome.Hsapiens.UCSC.hg38
#' ONTSampleGRanges <- readONTbam(bamfile = bam_path, ref_gr = NULL, modif = "C+m?",
#'                          prob_thresh = 0.50,genome = BSgenome.Hsapiens.UCSC.hg38,
#'                          context = "CG", chr = "chr1", region = NULL,
#'                          synonymous = FALSE, parallel = FALSE, BPPARAM = NULL)
#'                    
#' # plot the low resolution profile at 5 Kb resolution
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
#'                            distances = c(1,5,10,15), regions = NULL,
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
#' # select the genes
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
#' # merge the DMRs
#' DMRsNoiseFilterCGLarger <- mergeDMRsIteratively(DMRsNoiseFilterCG,
#'                            minGap = 500, respectSigns = TRUE,
#'                            methylationDataList[["WT"]],
#'                            methylationDataList[["met1-3"]],
#'                            context = "CG", minProportionDifference=0.4,
#'                            minReadsPerCytosine = 1, pValueThreshold=0.01,
#'                            test="score",alternative = "two.sided")
#'
#' # select the genes
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
#' # loading synthetic data
#' data("syntheticDataReplicates")
#'
#' # creating condition vector
#' condition <- c("a", "a", "b", "b")
#'
#' # computing DMRs using the neighbourhood method
#' DMRsReplicatesNeighbourhood <- computeDMRsReplicates(methylationData = methylationData,
#'                                                      condition = condition,
#'                                                      regions = NULL,
#'                                                      context = "CHH",
#'                                                      method = "neighbourhood",
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
#' # load the ONT methylation data 
#' data(ontSampleGRangesList)
#' 
#' # the regions where to compute the PMDs
#' chr1_ranges <- GRanges(seqnames = Rle("chr1"), ranges = IRanges(1E6+5E5,2E6))
#' 
#' # compute the PMDs in CG context with noise_filter method
#' PMDsNoiseFilterCG <- computePMDs(ontSampleGRangesList[["GM18501"]],
#'                                  regions = chr1_ranges,
#'                                  context = "CG",
#'                                  windowSize = 100,
#'                                  method = "noise_filter",
#'                                  kernelFunction = "triangular", 
#'                                  lambda = 0.5,
#'                                  minCytosinesCount = 4, 
#'                                  minMethylation = 0.4, 
#'                                  maxMethylation = 0.6, 
#'                                  minGap = 200, 
#'                                  minSize = 50, 
#'                                  minReadsPerCytosine = 4, 
#'                                  cores = 1,
#'                                  parallel = FALSE)
#'                                  
#' # compute the PMDs in CG context with neighbourhood method
#' PMDsNeighbourhoodCG <- computePMDs(ontSampleGRangesList[["GM18501"]],
#'                                    regions = chr1_ranges,
#'                                    context = "CG",
#'                                    method = "neighbourhood"
#'                                    minCytosinesCount = 4, 
#'                                    minMethylation = 0.4, 
#'                                    maxMethylation = 0.6, 
#'                                    minGap = 200, 
#'                                    minSize = 50, 
#'                                    minReadsPerCytosine = 4, 
#'                                    cores = 1,
#'                                    parallel = FALSE)
#'                                    
#' # compute the PMDs in CG context with bins method
#' PMDsBinsCG <- computePMDs(ontSampleGRangesList[["GM18501"]],
#'                           regions = chr1_ranges,
#'                           context = "CG",
#'                           method = "bins",
#'                           binSize = 100,
#'                           minCytosinesCount = 4, 
#'                           minMethylation = 0.4, 
#'                           maxMethylation = 0.6, 
#'                           minGap = 200, 
#'                           minSize = 50, 
#'                           minReadsPerCytosine = 4, 
#'                           cores = 1,
#'                           parallel = FALSE)  
#'  
#' # load the gene annotation data
#' data(GEs_hg38)
#' 
#' # select the transcript
#' transcript <- GEs_hg38[which(GEs_hg38$type == "transcript")]
#' 
#' # the regions where to compute the PMDs
#' regions <- GRanges(seqnames = Rle("chr1"), ranges = IRanges(1E6+5E5,2E6))
#' transcript <- transcript[overlapsAny(transcript, regions)]
#' 
#' # filter genes that are partially methylated in the two conditions
#' PMDsGenesCG <- filterPMDs(ontSampleGRangesList[["GM18501"]], 
#'                potentialPMDs = transcript, 
#'                context = "CG", minMethylation = 0.4, maxMethylation = 0.6,
#'                minCytosinesCount = 4, minReadsPerCytosine = 3, cores = 1)
#'                
#' # load the PMDs in CG context they were computed with minGap = 200
#' data(PMDsNoiseFilterCG)
#' 
#'
#' # merge the PMDs 
#' PMDsNoiseFilterCGLarger <- mergePMDsIteratively(PMDsNoiseFilterCG[1:100], 
#'                            minGap = 500, respectSigns = TRUE, 
#'                            ontSampleGRangesList[["GM18501"]], context = "CG", 
#'                            minReadsPerCytosine = 4, minMethylation = 0.4, 
#'                            maxMethylation = 0.6, cores = 1)
#' 
#' # set genomic coordinates where to compute PMDs
#' chr1_ranges <- GRanges(seqnames = Rle("chr1"), ranges = IRanges(1E6+5E5,2E6))
#' 
#' # compute PMDs and remove gaps smaller than 200 bp
#' PMDsNoiseFilterCG200 <- computePMDs(ontSampleGRangesList[["GM18501"]], 
#'                        regions = chr1_ranges, context = "CG", method = "noise_filter", 
#'                        windowSize = 100, kernelFunction = "triangular",  
#'                        minCytosinesCount = 1, minMethylation = 0.4, 
#'                        maxMethylation = 0.6, minGap = 0, minSize = 200, 
#'                        minReadsPerCytosine = 1, cores = 1)
#' PMDsNoiseFilterCG0 <- computePMDs(ontSampleGRangesList[["GM18501"]], 
#'                        regions = chr1_ranges, context = "CG", method = "noise_filter", 
#'                        windowSize = 100, kernelFunction = "triangular", 
#'                        minCytosinesCount = 1, minMethylation = 0.4, 
#'                        maxMethylation = 0.6, minGap = 0, minSize = 0, 
#'                        minReadsPerCytosine = 1, cores = 1)
#' PMDsNoiseFilterCG0Merged200 <- mergePMDsIteratively(PMDsNoiseFilterCG0, 
#'                              minGap = 200, respectSigns = TRUE, 
#'                              ontSampleGRangesList[["GM18501"]], context = "CG",  
#'                              minReadsPerCytosine = 4, minMethylation = 0.4, 
#'                              maxMethylation = 0.6, cores = 1)                      
#' 
#' #check that all newley computed PMDs are identical
#' print(all(PMDsNoiseFilterCG200 == PMDsNoiseFilterCG0Merged200))   
#' 
#' #retrive the number of reads in CG context in GM18501
#' PMDsNoiseFilterCGreadsCG <- analyseReadsInsideRegionsForConditionPMD(
#'                              PMDsNoiseFilterCG[1:10], 
#'                              ontSampleGRangesList[["GM18501"]], context = "CG", 
#'                              label = "GM18501")     
#'                              
#' # load the PMD data
#' data(PMDsBinsCG)
#' 
#' # compute the co-methylations with Fisher's exact test
#' coMetylationFisher <- computeCoMethylation(
#'   ontSampleGRangesList[[1]],
#'   regions = PMDsBinsCG,
#'   minDistance = 150,
#'   maxDistance = 1000,
#'   minCoverage = 4,
#'   pValueThreshold = 0.01,
#'   test = "fisher",
#'   parallel = FALSE)
#'   
#' # compute the co-methylations with Permuation test
#' coMetylationPermutation <- computeCoMethylation(
#'   ontSampleGRangesList[[1]],
#'   regions = PMDsBinsCG,
#'   minDistance = 150,
#'   maxDistance = 1000,
#'   minCoverage = 4,
#'   pValueThreshold = 0.01,
#'   test = "permutation",
#'   parallel = FALSE) # highly recommended to set as TRUE   
#'   
#' # select the transcript
#' transcript <- GEs_hg38[which(GEs_hg38$type == "transcript")]
#' 
#' # the regions where to compute the PMDs
#' regions <- GRanges(seqnames = Rle("chr1"), ranges = IRanges(1E6+5E5,2E6))
#' transcript <- transcript[overlapsAny(transcript, regions)]
#' 
#' # filter genes that are differntially methylated in the two conditions
#' VMRsGenesCG <- filterVMRsONT(ontSampleGRangesList[["GM18501"]], 
#'                ontSampleGRangesList[["GM18876"]], potentialVMRs = transcript, 
#'                context = "CG", pValueThreshold = 0.01, 
#'                minCytosinesCount = 4, minProportionDifference = 0.01, 
#'                minReadsPerCytosine = 3, ciExcludesOne = TRUE,
#'                varRatioFc = NULL, parallel = TRUE) # parallel recommended             
#'
#' }
#'
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


#' @name scanBamChr1Random5
#' @title The bam file from ONT nanopore .pod5 files
#' @docType data
#' @description
#' A \code{.bam} file containing the basecalling result of 5 sequences 
#' containing with the MM, ML tag from \code{dorado}.
#'
#' @format A \code{.bam} object
#'
#' @source The object was created by base-calling and alignment by \code{dorado, ver.0.9.6}
#' \url{https://github.com/nanoporetech/dorado?tab=readme-ov-file#dna-models}
#' to generate the \code{bam} files with \code{dna_r10.4.1_e8.2_400bps_hac@@v5.2.0} as basecalling model
#' from GM18501 cell line \url{https://s3.amazonaws.com/1000g-ont/index.html?prefix=pod5_data/GM18501_R9/}.
#' To subset, call the \code{bam} file by \code{\link{scanBam}} function from 
#' \code{\link{Rsamtools}} package and randomly select the 5 sequences. 
#' and repackaging the subsetted \code{list} to the \code{bam} file for test running.
#' 
NULL


#' @name ontSampleGRangesList
#' @title The ONT methylation data list
#' @description
#' A \code{GRangesList} object containing the methylation data at each cytosine
#' location in the genome in  GM18501 and  GM18876 B-Lymphocyte cell lines in
#' Homo sapiens. The data only contains the 1.5 ~ 2 Mbp from Chromosome 1.
#'
#' @format The \code{GRanges} elements contain six metadata columns
#' \describe{
#'  \item{context}{the context in which the DMRs are computed \code{"CG"}.}
#'  \item{trinucleotide_context}{the specific context of the cytosine (H is
#'                              replaced by the actual nucleotide).}
#'  \item{ONT_Cm}{comma-delimited read‐indices called modified}
#'  \item{ONT_C}{comma-delimited read‐indices covering but unmodified}
#'  \item{readsM}{the number of methylated reads.}
#'  \item{readsN}{the total number of reads.}
#' }
#'
#' @source Each element was created by calling \code{bam} files with \code{\link{readONTbam}}
#' function which in \code{\link{DMRcaller}} package.
#' The sample pod5 files were from the nanopore dataset from 1000 genome project 
#' \url{https://pmc.ncbi.nlm.nih.gov/articles/PMC10942501/}.
#' \url{https://s3.amazonaws.com/1000g-ont/index.html?prefix=pod5_data/GM18501_R9/} 
#' \code{.pod5} files in the case of GM18501 and
#' \url{https://s3.amazonaws.com/1000g-ont/index.html?prefix=pod5_data/GM18876_R9/}
#' \code{.pod5} files in the case of GM18876 cell line.
#' For base-calling and alignment, run \code{dorado, ver.0.9.6} 
#' \url{https://github.com/nanoporetech/dorado?tab=readme-ov-file#dna-models}
#' to generate the \code{bam} files with \code{dna_r10.4.1_e8.2_400bps_hac@@v5.2.0} as basecalling model.
NULL


#' @name GEs_hg38
#' @title The genetic elements data of GRCh38 Genome Reference
#' @docType data
#' @description
#' A \code{GRanges} object containing the annotation of the Homo sapiens (hg38)
#'
#' @format A \code{GRanges} object
#'
#' @source The object was created by loading the \code{gtf} file from 
#' \url{https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz}.
#' and concatenated by range in 1.5 ~ 2 Mbp from Chromosome 1
NULL


#' @name PMDsBinsCG
#' @title The PMDs between GM18501 and GM18876 using Bins method
#' @description
#' A \code{GRangesList} object containing the PMDs between GM18501 and
#' GM18876 B-Lymphocyte cell lines in Homo sapiens
#' (see \code{\link{ontSampleGRangesList}}). The PMDs were computed on the
#' 1.5 ~ 2 Mbp from Chromosome 1 with bins method using binSize of 1 kbp
#'
#' @format The \code{\link{GRanges}} element contain 5 metadata columns;
#' see \code{\link{computePMDs}}
#' @seealso \code{\link{filterPMDs}}, \code{\link{computePMDs}},
#' \code{\link{analyseReadsInsideRegionsForConditionPMD}}
#' and \code{\link{mergePMDsIteratively}}
#'
NULL


#' @name PMDsNoiseFilterCG
#' @title The PMDs between GM18501 and GM18876 using Noise_filter method
#' @description
#' A \code{GRangesList} object containing the PMDs between GM18501 and
#' GM18876 B-Lymphocyte cell lines in Homo sapiens
#' (see \code{\link{ontSampleGRangesList}}). The PMDs were computed on the
#' 1.5 ~ 2 Mbp from Chromosome 1 with noise filter method using a triangular kernel
#' and a windowSize of 100 bp
#'
#' @format The \code{\link{GRanges}} element contain 5 metadata columns;
#' see \code{\link{computePMDs}}
#' @seealso \code{\link{filterPMDs}}, \code{\link{computePMDs}},
#' \code{\link{analyseReadsInsideRegionsForConditionPMD}}
#' and \code{\link{mergePMDsIteratively}}
#'
NULL





