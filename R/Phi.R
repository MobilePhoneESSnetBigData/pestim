#' @title The product of ratioBeta and Kummer functions
#'
#' @description Compute the product of \code{\link{ratioBeta}} and
#' \code{\link{kummer}} functions with a specific set of arguments
#'
#' @param alpha, beta non-negative numeric vectors
#'
#' @param lambda numeric vector
#'
#' @param n non-negative integer vector
#'
#' @param relTol relative tolerance in the computation of the \code{\link{kummer}} function. Default
#' value is \code{1e-6}
#'
#' @param nThreads number (default the number of all cores, including logical cores) to use for computation
#'
#' @return \code{Phi} returns
#' \eqn{\frac{B(alpha + m, beta + n)}{B(alpha, beta)}\cdot {}_{1}F_{1}(lambda; alpha; beta)}, where
#' \eqn{{}_{1}F_{1}} stands for the confluent hypergeometric function
#'
#' The lengths of the input vectors must be all equal except when their length is 1, which are
#' recycled. Otherwise \code{NA}s are produced.
#'
#' @seealso \code{\link{ratioBeta}}, \code{\link{kummer}} for related functions.
#'
#' @examples
#' Phi(1, 1, 0.5, 10)
#' Phi(1:10, 10:1, seq(0, 1, length.out = 10), 3)
#' Phi(1:4, 4:1, c(2, 3), c(4, 3, 1))
#'
#' @include kummer.R ratioBeta.R
#'
#' @export
Phi <- function(alpha, beta, lambda, n, relTol = 1e-6, nThreads = RcppParallel::defaultNumThreads()){
  output <- ratioBeta(alpha, beta, n, 0) * kummer(lambda, beta, alpha + beta + n, relTol, nThreads)
  return(output)

}

