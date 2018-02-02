#' @title The ratio of two beta functions.
#'
#' @description Ratio of two beta functions whose arguments differ by integer numbers
#'
#' @param alpha, beta non-negative numeric vectors
#'
#' @param m, n non-negative integer vectors
#'
#' @return \code{ratioBeta} gives \eqn{\frac{B(alpha + m, beta + n)}{B(alpha, beta)}}
#'
#' The lengths of the input vectors must be all equal except when their length is 1, which are
#' recycled. Otherwise \code{NA}s are produced.
#'
#' @seealso \code{\link{beta}}, \code{\link{lbeta}} for related functions.
#'
#' @examples
#' ratioBeta(10, 13, 2, 3)
#' ratioBeta(1:10, 10:1, 2, 3)
#' ratioBeta(1:3, 3:1, c(2, 3), 4)
#'
#' @export
#'
ratioBeta <- function(alpha, beta, m, n){

  argExp <- numeric(length(alpha))
  if (length(m) == 1) m <- rep(m, length(alpha))
  if (length(n) == 1) n <- rep(n, length(alpha))

  nonZero <- (alpha > .Machine$double.eps & beta > .Machine$double.eps)
  argExp[nonZero] <- lbeta( alpha[nonZero] + m[nonZero], beta[nonZero] + n[nonZero] ) -
    lbeta( alpha[nonZero], beta[nonZero] )
  output <- exp(argExp)
  return(output)
}
