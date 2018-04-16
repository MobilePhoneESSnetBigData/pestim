#' @title Confluent hypergeometric or Kummer function
#'
#' @description Partial implementation of the confluent hypergeometric function
#' \eqn{{}_{1}F_{1}(x ; a; b)}
#'
#' @param x, a, b numeric vectors of the same length
#'
#' @param relTol relative tolerance (default value \code{1e-6}) understood as the ratio of each term
#' in the series relative to the sum
#'
#' @param nThreads number (default the number of all cores, including logical cores) to use for computation
#'
#' @return Return a numeric vector with the values of the function
#'
#' @details This function is implemented in C++. It is based on Pearson et al (2016). It only
#' implements the Taylor series method together with an asymtoptic expansion based on Watson's lemma
#'
#' @author Luis Sanguiao Bogdan Oancea
#'
#' @useDynLib pestim
#' @importFrom Rcpp sourceCpp
#' @export
kummer <- function(x, a, b, relTol = 1e-8, nThreads = RcppParallel::defaultNumThreads()){
  
  n <- length(x)
  if (n > 1 & length(a) == 1) a <- rep(a, n)
  if (n > 1 & length(b) == 1) b <- rep(b, n)
  if (length(a) != n) stop('[pestim::kummer] a must have the same length as x.')
  if (length(b) != n) stop('[pestim::kummer] b must have the same length as x.')
  
  err <- as.integer(rep(0, n))
  if(nThreads > 1) {
    
    RcppParallel::setThreadOptions(numThreads = nThreads)
    output <- pKummer(x, a, b, relTol, err)
  
  } else {
    
    output <- Kummer(x, a, b, relTol, err)
  
  }

  output[is.infinite(output)] <- .Machine$double.xmax
  return(output)
}