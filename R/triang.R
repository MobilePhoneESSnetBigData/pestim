#' @name triang
#' @aliases dtriang
#' @aliases ptriang
#' @aliases qtriang
#' @aliases rtriang
#'
#' @title The Triangular Distribution.
#'
#' @description Density, distribution funtion, quantile function and random generation for the
#' triangular distribution
#'
#' @param x, q vector of quantiles
#'
#' @param p vector pf probabilities
#'
#' @param n number of observations
#'
#' @param xMin  the minimum values of the range of the random variable
#'
#' @param xMax the maximum values of the range of the random variable
#'
#' @param xMode the mode of the random variable
#'
#' @return \code{dtriang} gives the density, \code{ptriang} gives the distribution function,
#' \code{qtriang} gives the quantile function, and \code{rtriang} generates random deviates.
#'
#'
#' @seealso \code{\link{Distributions}} for other distributions
#'
#' @examples
#' curve(dtriang(x, 0, 3, 1), xlim = c(0, 3))
#' curve(ptriang(x, 0, 3, 1), xlim = c(0, 3))
#' curve(qtriang(x, 0, 3, 1), xlim = c(0, 1))
#' hist(rtriang(1e6, 0, 3, 1), breaks = seq(0, 3, by = 0.01))
#'
#' @rdname triang
#' @export
#'
#'
#' @rdname triang
#' @export
dtriang <- function(x, xMin, xMax, xMode){

  if (xMin > xMax) stop('xMax must be greater than xMin.')
  if (!((xMode >= xMin) & (xMode <= xMax))) stop('xMode must be between xMin and xMax.')


  output <- x
  output[x <= xMin | x >= xMax] <- 0
  range1 <- (x > xMin & x <= xMode)
  output[range1] <- (2 * (x[range1] - xMin)) / ((xMax - xMin) * (xMode - xMin))
  range2 <- (x >= xMode & x < xMax)
  output[range2] <- (2 * (xMax - x[range2])) / ((xMax - xMin) * (xMax - xMode))
  return(output)
}


#' @rdname triang
#' @export
ptriang <- function(q, xMin, xMax, xMode){

  if (xMin > xMax) stop('xMax must be greater than xMin.')
  if (!((xMode >= xMin) & (xMode <= xMax))) stop('xMode must be between xMin and xMax.')

  output <- q
  output[q <= xMin] <- 0
  output[q >= xMax] <- 1
  range1 <- (q > xMin & q <= xMode)
  output[range1] <- ((output[range1] - xMin)^2) / ((xMax - xMin) * (xMode - xMin))
  range2 <- (q > xMode & q < xMax)
  output[range2] <- 1 - ((output[range2] - xMax)^2) / ((xMax - xMin) * (xMax - xMode))
  return(output)
}


#' @rdname triang
#' @export
qtriang <- function(q, xMin, xMax, xMode){

  if (xMin > xMax) stop('xMax must be greater than xMin.')
  if (!((xMode >= xMin) & (xMode <= xMax))) stop('xMode must be between xMin and xMax.')


  output <- q
  range1 <- (q < (xMode - xMin) / (xMax - xMin))
  output[range1] <- xMin + sqrt(q[range1] * (xMax - xMin) * (xMode - xMin))
  range2 <- ( q > (xMode - xMin) / (xMax - xMin))
  output[range2] <- xMax - sqrt((1 - q[range2]) * (xMax - xMin) * (xMax - xMode))
  return(output)
}

#' @rdname triang
#' @export
rtriang <- function(n, xMin, xMax, xMode){
  if (xMin > xMax) stop('xMax must be greater than xMin.')
  if (!((xMode >= xMin) & (xMode <= xMax))) stop('xMode must be between xMin and xMax.')
  u <- runif(n)
  mc <- match.call()
  mc[[1L]] <- qtriang
  mc[['n']] <- NULL
  mc[['q']] <- u
  output <- eval(mc)
  return(output)
}

