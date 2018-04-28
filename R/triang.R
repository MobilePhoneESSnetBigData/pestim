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
dtriang <- function(x, xMin, xMax, xMode){

  if (any(xMin > xMax)) stop('xMax must be greater than xMin.')
  if (!(all(xMode >= xMin) & all(xMode <= xMax))) stop('xMode must be between xMin and xMax.')


  n <- length(x)
  if (length(xMin) == 1) xMin <- rep(xMin, n)
  if (length(xMax) == 1) xMax <- rep(xMax, n)
  if (length(xMode) == 1) xMode <- rep(xMode, n)

  output <- x
  output[x <= xMin | x >= xMax] <- 0
  range1 <- (x > xMin & x <= xMode)
  output[range1] <- (2 * (x[range1] - xMin[range1])) / ((xMax[range1] - xMin[range1]) * (xMode[range1] - xMin[range1]))
  range2 <- (x >= xMode & x < xMax)
  output[range2] <- (2 * (xMax[range2] - x[range2])) / ((xMax[range2] - xMin[range2]) * (xMax[range2] - xMode[range2]))
  return(output)
}

#' @rdname triang
#' @export
ptriang <- function(q, xMin, xMax, xMode){

  if (any(xMin > xMax)) stop('xMax must be greater than xMin.')
  if (!(all(xMode >= xMin) & all(xMode <= xMax))) stop('xMode must be between xMin and xMax.')

  n <- length(q)
  if (length(xMin) == 1) xMin <- rep(xMin, n)
  if (length(xMax) == 1) xMax <- rep(xMax, n)
  if (length(xMode) == 1) xMode <- rep(xMode, n)

  output <- q
  output[q <= xMin] <- 0
  output[q >= xMax] <- 1
  range1 <- (q > xMin & q <= xMode)
  output[range1] <- ((output[range1] - xMin[range1])^2) / ((xMax[range1] - xMin[range1]) * (xMode[range1] - xMin[range1]))
  range2 <- (q > xMode & q < xMax)
  output[range2] <- 1 - ((output[range2] - xMax[range2])^2) / ((xMax[range2] - xMin[range2]) * (xMax[range2] - xMode[range2]))
  return(output)
}

#' @rdname triang
#' @export
qtriang <- function(q, xMin, xMax, xMode){

  if (any(xMin > xMax)) stop('xMax must be greater than xMin.')
  if (!(all(xMode >= xMin) & all(xMode <= xMax))) stop('xMode must be between xMin and xMax.')

  n <- length(q)
  if (length(xMin) == 1) xMin <- rep(xMin, n)
  if (length(xMax) == 1) xMax <- rep(xMax, n)
  if (length(xMode) == 1) xMode <- rep(xMode, n)

  output <- q
  range1 <- (q < (xMode - xMin) / (xMax - xMin))
  output[range1] <- xMin[range1] + sqrt(q[range1] * (xMax[range1] - xMin[range1]) * (xMode[range1] - xMin[range1]))
  range2 <- ( q > (xMode - xMin) / (xMax - xMin))
  output[range2] <- xMax[range2] - sqrt((1 - q[range2]) * (xMax[range2] - xMin[range2]) * (xMax[range2] - xMode[range2]))
  return(output)
}

#' @rdname triang
#' @export
rtriang <- function(n, xMin, xMax, xMode){

  if (any(xMin > xMax)) stop('xMax must be greater than xMin.')
  if (!(all(xMode >= xMin) & all(xMode <= xMax))) stop('xMode must be between xMin and xMax.')

  u <- runif(n)
  mc <- match.call()
  mc[[1L]] <- qtriang
  mc[['n']] <- NULL
  mc[['q']] <- u
  
  output <- eval(mc, sys.frame(sys.parent()))
  return(output)
}
