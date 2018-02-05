#' @title Time evolution of the posterior mean, median, and mode for the number of individuals unconditional on their estimated initial values
#' @description Compute the time evolution of the posterior mean, median, and mode for the number of individuals
#' unconditional on their estimated initial values, but starting from the input data themselves.
#' @param nMNOmat ??
#'
#' @param nReg ??
#'
#' @param fu, fv named lists with the prior marginal distributions of the two-dimensional points
#' for the Monte Carlo integration
#'
#' @param flambda named list with the prior distribution of the lambda parameter
#'
#' @param n number of points to generate in the posterior distribution for the computation. Default
#' value is 1e3
#'
#' @param scale numeric vector with the scale to count the number of individuals. Default value is 1
#'
#' @param relTol relative tolerance in the computation of the \code{\link{kummer}} function. Default
#' value is \code{1e-5}
#'
#' @param nSim number of two-dimensional points to generate to compute the integral. Default value
#' is \code{1e3}
#'
#' @param nStrata integer vector of length 2 with the number of strata in each dimension. Default
#' value is \code{c(1, 1e2)}
#'
#' @param verbose logical (default \code{FALSE}) to report progress of the computation
#'
#' @return \code{postNt} computes the posterior mean, median, and mode of the posterior distribution
#' for each cell. The function returns a matrix with the estimates in columns and the cells in rows.
#'
#' @details The prior distributions are specified as named lists where the first component of each list must
#' be the name of distribution ('unif', 'triang', 'degen', 'gamma') and the rest components must be
#' named according to the name of the parameters of the random generator of the corresponding
#' distribution according to:
#'
#'   \itemize{
#'
#'     \item unif: \code{xMin}, \code{xMax} for the minimum, maximum of the sampled interval.
#'     \item degen: \code{x0} for the degenerate value of the random variable.
#'     \item triang: \code{xMin}, \code{xMax}, \code{xMode} for minimum, maximum and mode (see
#'     \code{\link{qtriang}}).
#'     \item gamma: \code{scale} and \code{shape} with the same meaning as in \code{\link{rgamma}}.
#'   }

#'
#' @seealso \code{\link{rNt}}
#'
#' @examples
#'
#' @include rNt.R
#'
#' @author David Salgado
#' @export
#'
postNt <- function(nMNOmat, nReg, fu, fv, flambda, distNames, variation, scale = 1, n = 1e3, relTol = 1e-5, nSim = 1e3, nStrata = c(1, 1e2), verbose = FALSE){

  Ntmat <- rNt(n, nMNOmat, nReg, fu, fv, flambda, distNames, variation, scale, relTol, nSim, nStrata, verbose)
  postMean <- Ntmat[, round(mean(N)), by = c('cellID')]
  setnames(postMean, 'V1', 'value')
  postMean[, variable := 'postMean']
  postMedian <- Ntmat[, round(median(N)), by = c('cellID')]
  setnames(postMedian, 'V1', 'value')
  postMedian[, variable := 'postMedian']
  fmode <- function(N){N[which.max(names(table(N)))]}
  postMode <- Ntmat[, fmode(N), by = c('cellID')]
  setnames(postMode, 'V1', 'value')
  postMode[, variable := 'postMode']
  DT <- rbindlist(list(postMean, postMedian, postMode))
  output <- dcast(DT, cellID ~ variable, value.var = 'value')
  output[, cellID := NULL]
  output <- as.matrix(output)
  return(output)

}
