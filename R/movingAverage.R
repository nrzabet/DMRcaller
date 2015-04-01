#' This function computes the moving average for a vector that has missing 
#' values. 
#'
#' @title Moving average 
#' @param minPos the first position of the input fragment
#' @param maxPos the last position of the input fragment
#' @param pos the positions where the vector has values
#' @param val the values of the input vector
#' @param weights the weights associated with the moving average 
#' @param windowSize the size of the moving average
#' @param normalize \code{logical} value that indicates whether the moving 
#' average is normalised
#' @param kernelFunction a \code{character} indicating which kernel function to 
#' be used. Can be one of \code{"uniform"}, \code{"triangular"}, 
#' \code{"gaussian"} or \code{"epanechnicov"}.
#' @param lambda numeric value required for the Gaussian filter 
#' (\code{K(x) = exp(-lambda*x^2)})
#' @return A numeric vector with the values computed by the moving average
#'  
#' @author Radu Zabet
.movingAverage <- function(minPos, maxPos, pos, val, weights=1, windowSizeHalf=150, normalize = FALSE, kernelFunction="triangular", lambda=0.5) {
  
  .stopIfNotAll(c(length(pos)==length(val)), "pos and val vectors need to have the same length")
  
  if(length(weights) < length(val)){
    weights <- rep(weights, length.out=length(val))
  } else if(length(weights) > length(val)){
    weights <- weights[1:length(val)]
  }
  
  # Filter out NAs
  keepIndexes <- which(!is.na(pos) & !is.na(val) & !is.na(weights))
  pos <- pos[keepIndexes]
  weights <- weights[keepIndexes]  
  val <- val[keepIndexes]
    
  #set the values
  rawVector <- rep(0, maxPos - minPos + 2*windowSizeHalf + 1)
  rawVector[(pos - minPos + windowSizeHalf + 1)] <-weights*val
  
  normVector <- rep(0, maxPos - minPos + 2*windowSizeHalf + 1)
  normVector[(pos - minPos + windowSizeHalf + 1)] <-weights
  
  
  
  # Define the (triangular) kernel.
  if(kernelFunction =="uniform"){
    kernel <- rep(1, times=(2*windowSizeHalf + 1))  
  } else if(kernelFunction =="triangular"){
    kernel <- 1 - abs(-windowSizeHalf:windowSizeHalf)/windowSizeHalf 
  } else if(kernelFunction =="gaussian"){
    kernel <- .gaussianKernel(lambda, - windowSizeHalf:windowSizeHalf)
  } else if(kernelFunction =="epanechnicov"){
    kernel <- .epanechnicovKernel(-windowSizeHalf:windowSizeHalf)
  } else{
    stop(paste("Unknown kernel function: ", kernelFunction, ". 
               It should be one of \"uniform\", \"triangular\", \"gaussian\", \"epanechnicov\"",sep=""))
  }
  
  kernel <- kernel / sum(kernel)
  
  if(windowSizeHalf >= 1){
    smoothedVector <- RcppRoll::roll_sum(rawVector, length(kernel), weights = kernel, normalize = normalize) / RcppRoll::roll_sum(normVector, length(kernel), weights = kernel, normalize = normalize)
  } else{
    smoothedVector <- rawVector 
  }
  
  
  return(smoothedVector);
}


.gaussianKernel <- function(lambda,x){
  return(exp(-lambda*x^2))
}

.epanechnicovKernel <- function(x){
  return((3*(1-x^2))/4)
}

#' This function computes the moving sum for a vector that has missing values. 
#'
#' @title Moving sum 
#' @param minPos the first position of the input fragment
#' @param maxPos the last position of the input fragment
#' @param pos the positions where the vector has values
#' @param val the values of the input vector
#' @param weights the weights associated with the moving sum 
#' @param windowSize the size of the moving sum
#' @param normalize \code{logical} value that indicates whether the moving sum is normalised
#' @return A numeric vector with the values computed by the moving sum
#'  
#' @author Radu Zabet
.movingSum <- function(minPos, maxPos, pos, val, weights=1, windowSize=150, normalize = FALSE) {
  
  .stopIfNotAll(c(length(pos)==length(val)), "pos and val vectors need to have the same length")
  
  if(length(weights) < length(val)){
    weights <- rep(weights, length.out=length(val))
  } else if(length(weights) > length(val)){
    weights <- weights[1:length(val)]
  }
  
  # Filter out NAs
  keepIndexes <- which(!is.na(pos) & !is.na(val) & !is.na(weights))
  pos <- pos[keepIndexes]
  weights <- weights[keepIndexes]  
  val <- val[keepIndexes]
  
  
  #set the values
  rawVector <- rep(0, (maxPos - minPos + 1))
  rawVector[pos - minPos + 1] <-weights*val
  rawVector <- c(rawVector,rep(0,(windowSize-1)))
  
  # Define the (triangular) kernel.
  kernel <- c(rep(1, (windowSize)))
  
  
  smoothedVector <- RcppRoll::roll_sum(rawVector, length(kernel), weights = kernel, normalize = normalize) 
  
  #smoothedVector <- smoothedVector[seq(1,length(smoothedVector), by=windowSize)]
  
  return(smoothedVector);
}

