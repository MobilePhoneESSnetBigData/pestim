#' @title Confluent hypergeometric or Kummer function
#'
#' @description A parallel implementation of the confluent hypergeometric function
#' \eqn{{}_{1}F_{1}(x ; a; b)} using RcppParallel package
#'
#' @param x, a, b numeric vectors of the same length
#'
#' @param relTol relative tolerance (default value \code{1e-6}) understood as the ratio of each term
#' in the series relative to the sum
#'
#' @return Return a numeric vector with the values of the function
#'
#' @details This function is implemented in C++. It is based on Pearson et al (2016). It
#' implements the Taylor series method together with an asymtoptic expansion based on Watson's lemma. The
#' actual implementation uses parallelFor from RcppParallel package together with a functor that computes
#'confluent hypergeometric function for chunks of vectors x, a, b
#'
#' @author Bogdan Oancea
#'
#' @useDynLib pestim
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel defaultNumThreads setThreadOptions
#' @export
pkummer <- function(x, a, b, relTol = 1e-6, nThreads = RcppParallel::defaultNumThreads()){
  RcppParallel::setThreadOptions(numThreads = nThreads)
  output <- pKummer(x, a, b, relTol)
  output[is.infinite(output)] <- .Machine$double.xmax
  return(output)
}
