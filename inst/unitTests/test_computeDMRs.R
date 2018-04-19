library(DMRcaller)
library(RUnit)

# create synthetic data
step <- 25
numberOfCytosines <- 200
startPos <- seq(from = 1, by=step, length.out=numberOfCytosines)

syntheticMethylationData1 <- GRanges(
  seqnames              = factor(rep("Chr1", length(startPos))),
  ranges                = IRanges(startPos, startPos),
  strand                = rep("+", length(startPos)),
  context               = rep("CG", length(startPos)),
  readsM                = rep(0, length(startPos)),
  readsN                = rep(0, length(startPos)),
  trinucleotide_context = rep("CGH", length(startPos)))

syntheticMethylationData1$readsN <- sample(50:150,
                                           length(syntheticMethylationData1$readsN),
                                           replace = TRUE)
syntheticMethylationData1$readsN <- c(rep(50,50),rep(100, 100), rep(50,50))


proportion1 <- c(rep(0.1, 10), rep(0.5, 30), rep(0.2, 3), rep(0.45, 7),
                 rep(0.2, 50), rep(0.8, 20), rep(0.3, 10), rep(0.045, 20),
                 rep(0.3, 20), rep(0.7, 20), rep(0.3, 10))
syntheticMethylationData1$readsM <- syntheticMethylationData1$readsN*proportion1

syntheticMethylationData2 <- syntheticMethylationData1
proportion2 <- c(rep(0.045, 130), rep(0.6, 20), rep(0.045, 50))
syntheticMethylationData2$readsM <- round(proportion2*syntheticMethylationData2$readsN)


#DMRs
vals <- rep(0, length(proportion1))
vals[(syntheticMethylationData1$readsM/syntheticMethylationData1$readsN -
        syntheticMethylationData2$readsM/syntheticMethylationData2$readsN)>0.4] <- 1
vals[(syntheticMethylationData2$readsM/syntheticMethylationData2$readsN -
        syntheticMethylationData1$readsM/syntheticMethylationData1$readsN)>0.4] <- -1



#join the differentially methylated cytosines into regions
rle <- rle(vals)
rle$cumulative <- cumsum(rle$lengths)
endOfRuns <- rle$cumulative

DMRs <- GRanges(
  seqnames    = "Chr1",
  ranges      = IRanges(((endOfRuns - rle$lengths + 1)*step), (endOfRuns*step)),
  strand      = "+",
  direction   = rle$values
)
DMRs <- DMRs[DMRs$direction!=0]


genes <- GRanges(
  seqnames    = "Chr1",
  ranges      = IRanges(c(300, 1800),c(1200,2100)),
  strand      = "+"
)



# check if compute coverage works
test_computeMethylationDataCoverage <- function() {
  breaks <- c(1,10,40,90,100)
  answer <- c(1,1,1,0.5,0)
  checkEqualsNumeric(computeMethylationDataCoverage(syntheticMethylationData1,
                                                    context="CG",
                                                    breaks=breaks),
                     answer, tolerance=1.0e-8)
  checkEqualsNumeric(computeMethylationDataCoverage(syntheticMethylationData2,
                                                    context="CG",
                                                    breaks=breaks),
                     answer, tolerance=1.0e-8)
}


# check if compute coverage works
test_computeDMRs <- function() {
  DMRs_CG_noise_filter <-  computeDMRs(syntheticMethylationData1,
                                       syntheticMethylationData2,
                                       NULL,
                                       method = "noise_filter",
                                       context = "CG",
                                       pValueThreshold = 0.01,
                                       minReadsPerCytosine = 3,
                                       minProportionDifference = 0.4,
                                       minGap = 50,
                                       minSize = 50,
                                       test = "score",
                                       windowSize=50)

  checkTrue(all(overlapsAny(DMRs_CG_noise_filter, DMRs)))
  checkTrue(all(overlapsAny(DMRs, DMRs_CG_noise_filter)))

  DMRs_CG_bins <-  computeDMRs(syntheticMethylationData1,
                               syntheticMethylationData2,
                               NULL,
                               method = "bins",
                               context = "CG",
                               binSize = 100,
                               test = "score",
                               pValueThreshold = 0.01,
                               minCytosinesCount = 3,
                               minProportionDifference = 0.4,
                               minGap = 100,
                               minReadsPerCytosine = 3,
                               cores = 1)
  checkTrue(all(overlapsAny(DMRs_CG_bins, DMRs)))
  checkTrue(all(overlapsAny(DMRs, DMRs_CG_bins)))
}

test_computeDMRs <- function() {
}


test_filterDMRs <- function() {
  DMRs_CG_genes <-  filterDMRs(syntheticMethylationData1,
                               syntheticMethylationData2,
                               genes,
                               context = "CG",
                               test = "score",
                               pValueThreshold=0.01,
                               minCytosinesCount = 1,
                               minProportionDifference = 0.4,
                               minReadsPerCytosine = 3)
  checkTrue(DMRs_CG_genes == genes[overlapsAny(genes, DMRs)])
}

test_filterDMRs <- function() {
  DMRs_CG_noise_filter <-  computeDMRs(syntheticMethylationData1,
                                       syntheticMethylationData2,
                                       NULL,
                                       method = "noise_filter",
                                       context="CG",
                                       pValueThreshold = 0.01,
                                       minReadsPerCytosine = 3,
                                       minProportionDifference = 0.4,
                                       minGap = 50,
                                       minSize = 50,
                                       test = "score",
                                       windowSize=50)

  DMRs_CG_noise_filter_filtered <- filterDMRs(syntheticMethylationData1,
                                              syntheticMethylationData2,
                                              DMRs_CG_noise_filter,
                                              context = "CG",
                                              minCytosinesCount = 10,
                                              minProportionDifference = 0.55,
                                              minReadsPerCytosine = 20,
                                              pValueThreshold = 0.01)

  checkTrue(all(overlapsAny(DMRs_CG_noise_filter_filtered, DMRs[3:5])))
  checkTrue(all(overlapsAny(DMRs[3:5], DMRs_CG_noise_filter_filtered)))
}

test_mergeDMRsIteratively <- function() {
  DMRs_CG_bins <-  computeDMRs(syntheticMethylationData1,
                               syntheticMethylationData2,
                               NULL,
                               context = "CG",
                               binSize = 100,
                               method = "bins",
                               test = "score",
                               pValueThreshold = 0.01,
                               minCytosinesCount = 3,
                               minProportionDifference = 0.4,
                               minGap = 0,
                               minReadsPerCytosine = 3,
                               cores = 1)

  DMRs_CG_bins_merged <- mergeDMRsIteratively(DMRs_CG_bins,
                                              minGap = 200,
                                              respectSigns = TRUE,
                                              syntheticMethylationData1,
                                              syntheticMethylationData2,
                                              context = "CG",
                                              minProportionDifference = 0.4,
                                              minReadsPerCytosine = 3,
                                              pValueThreshold = 0.01,
                                              test = "score",
                                              alternative = "two.sided")
  checkTrue(length(DMRs_CG_bins_merged) < length(DMRs_CG_bins))
  checkTrue(all(overlapsAny(DMRs_CG_bins_merged, DMRs_CG_bins)))
  checkTrue(all(overlapsAny(DMRs_CG_bins, DMRs_CG_bins_merged)))
}

# testing replicates
data("methylationDataList")

# testing joinReplicates
test_joinReplicates <- function(){
  methylationDataReplicates <- joinReplicates(methylationData1 = methylationDataList$WT,
                 methylationData2 = methylationDataList$`met1-3`,
                 usecomplete = FALSE)
  checkTrue(length(methylationDataList$WT) == length(methylationDataReplicates))
  checkTrue(length(methylationDataList$`met1-3`) == length(methylationDataReplicates))
  # checks that the number of reads columns is even
  checkTrue(length(grep("reads", colnames(mcols(methylationDataReplicates)))) %% 2 == 0)
}

# testing computeDMRsReplicates
data("syntheticDataReplicates")

condition <- c("a", "a", "b", "b")

test_computeDMRsReplicates <- function(){
  DMRsReplicatesNeighbourhood <- computeDMRsReplicates(methylationData = methylationData,
                                                       condition = condition,
                                                       regions = NULL,
                                                       context = "CHH",
                                                       method = "neighbourhood",
                                                       test = "betareg",
                                                       pseudocountM = 1,
                                                       pseudocountN = 2,
                                                       pValueThreshold = 0.01,
                                                       minCytosinesCount = 4,
                                                       minProportionDifference = 0.4,
                                                       minGap = 200,
                                                       minSize = 50,
                                                       minReadsPerCytosine = 4,
                                                       cores = 1)
  checkTrue(DMRsReplicatesNeighbourhood$pValue < 0.01)
  checkTrue(DMRsReplicatesNeighbourhood$proportion2 -
            DMRsReplicatesNeighbourhood$proportion1 < 0.4)
  checkTrue(DMRsReplicatesNeighbourhood$context == "CHH")
  checkTrue(overlapsAny(ranges(DMRsReplicatesNeighbourhood),
                        ranges(methylationData)))

}
