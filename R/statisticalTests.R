

#' This function computes the p-values of the Score test.  
#'
#' @title Score test
#' @param m1 the number of methylated reads in condition 1
#' @param n1 the total number of methylated reads in condition 1
#' @param m2 the number of methylated reads in condition 2
#' @param n2 the total number of methylated reads in condition 2
#' @return The p-values of the Score test
#'
#' @author Radu Zabet
.scoreTest <-  function(m1, n1, m2, n2){
  p <- (m1 + m2)/(n1 + n2)

  zScore <- rep(NA, times=length(m1))
  bufferIndexes <- which((p * (1-p))!=0)
  p1 <- m1/n1
  p2 <- m2/n2
  zScore[bufferIndexes]  <- (p1[bufferIndexes] - p2[bufferIndexes]) / (sqrt(p[bufferIndexes] * (1-p[bufferIndexes])*(1/n1[bufferIndexes]+1/n2[bufferIndexes])))
  bufferIndexes <- which((p * (1-p))==0)
  zScore[bufferIndexes] <- 0

  pValue <- rep(NA, times=length(zScore))
  bufferIndexes <- which(!is.na(zScore))
  pValue[bufferIndexes] <- 2*pnorm(-abs(zScore[bufferIndexes]))

  return(pValue)
}


#' This function computes the p-values of the Fisher's exact test.
#'
#' @title Fisher test
#' @param m1 the number of methylated reads in condition 1
#' @param n1 the total number of methylated reads in condition 1
#' @param m2 the number of methylated reads in condition 2
#' @param n2 the total number of methylated reads in condition 2
#' @param alternative indicates the alternative hypothesis and must be one of
#' \code{"two.sided"}, \code{"greater"} or \code{"less"}. You can specify just
#' the initial letter. Only used in the 2 by 2 case.
#' @return The p-values of Fisher's exact test.
#'
#' @author Radu Zabet
.fisherTest <-  function(m1, n1, m2, n2, alternative = c("two.sided", "less", "greater")){

  u1 <- n1 - m1
  u2 <- n2 - m2
  dataM <- cbind(m1, u1, m2, u2)

  pValue <- rep(NA, times=length(m1))

  bufferIndexes <- which(!(is.na(u1) | is.na(u2) | n1 <= 0 | n2 <= 0 | m1 < 0 | m2 < 0 | m1 > n1 | m2 > n2))

  if(length(bufferIndexes) > 0){
    pValue[bufferIndexes] <- unlist(lapply( apply(matrix(dataM[bufferIndexes, ], ncol=4), 1, list),
                                            .fisherTestPValue,
                                            alternative = alternative))
  }
  #for(i in bufferIndexes){
  #    pValue[i] <- fisher.test(matrix(dataM[i, ],nrow=2, byrow=TRUE), alternative = alternative)$p.value
  #}

  #pValue[bufferIndexes] <- apply(dataM[bufferIndexes, ],1, .fisherTestPValue, alternative = alternative)

  pValue[pValue > 1] <- 1
  pValue[pValue < 0] <- 0

  return(pValue)
}

#' This function computes the p-values of the Fisher's exact test.
#'
#' @title Fisher test
#' @param x a vector with 4 values
#' @return The p-values of Fisher's exact test.
#'
#' @author Radu Zabet
.fisherTestPValue <-  function(x, alternative = c("two.sided", "less", "greater")){
  return(fisher.test(matrix(unlist(x),nrow=2, byrow=TRUE), alternative = alternative)$p.value)
}

#' This function computes the adjusted p-values (using Benjamini & Hochberg
#' method)
#'
#'
#' @title Compute adjusted p-values
#' @param test a \code{character} indicating the statistical test used.
#' (\code{"fisher"} for Fisher's exact test or \code{"score"} for Score test).
#' @param m1 the number of methylated reads in condition 1.
#' @param n1 the total number of methylated reads in condition 1.
#' @param m2 the number of methylated reads in condition 2.
#' @param n2 the total number of methylated reads in condition 2.
#' @param alternative indicates the alternative hypothesis and must be one of
#' \code{"two.sided"}, \code{"greater"} or \code{"less"}. You can specify just
#' the initial letter. Only used in the 2 by 2 case. This is used only for
#' Fisher's test
#' @return The adjusted p-values of the statistical test.
#' @author Radu Zabet
.computeAdjuestedPValues <- function(test, m1, n1, m2, n2, alternative = "two.sided"){
  if(test == "score"){
    pValue <- .scoreTest(m1,n1,m2,n2)
  } else{
    pValue <- .fisherTest(m1,n1,m2,n2, alternative="two.sided")
  }

  # convert p-values to FDR
  adjPValue <- rep(NA, times=length(pValue))

  adjPValue[which(!is.na(pValue))] <- p.adjust(pValue[which(!is.na(pValue))], method="fdr")

  return(adjPValue)
}



#' This function computes the adjusted p-values (using Benjamini & Hochberg
#' method) using one the statistical tests
#'
#' @title Compute adjusted p-values in DMRs
#' @param test a \code{character} indicating the statistical test used .
#' (\code{"fisher"} for Fisher's exact test or \code{"score"} for Score test).
#' @param DMRs a \code{\link{GRanges}} object with the DMRs.
#' @param alternative indicates the alternative hypothesis and must be one of
#' \code{"two.sided"}, \code{"greater"} or \code{"less"}. You can specify just
#' the initial letter. Only used in the 2 by 2 case. This is used only for
#' Fisher's test.
#' @return The adjuested p-values of the statistical test.
#'
#' @author Radu Zabet
.computeaAjustedPValuesInDMRs <- function(test, DMRs ,alternative = "two.sided"){


  m1 <- DMRs$sumReadsM1
  n1 <- DMRs$sumReadsN1
  m2 <- DMRs$sumReadsM2
  n2 <- DMRs$sumReadsN2

  if(test == "score"){
    pValue <- .scoreTest(m1,n1,m2,n2)
  } else{
    pValue <- .fisherTest(m1,n1,m2,n2, alternative="two.sided")
  }



  # convert p-values to FDR
  adjPValue <- rep(NA, times=length(pValue))

  adjPValue[which(!is.na(pValue))] <- p.adjust(pValue[which(!is.na(pValue))], method="fdr")

  return(adjPValue)
}
