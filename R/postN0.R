#' @title Posterior mean, median, and mode for the number of individuals at the initial time.
#'
#' @description Compute the posterior mean, median, and mode for the number of individuals
#' generating posterior distribution according to the hierarchical model at the initial time instant
#'
#' @param nMNO, nReg non-negative integer vectors with the number of individuals detected in each
#' cell according to the network operator and the register
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
#' value is \code{1e-6}
#'
#' @param nSim number of two-dimensional points to generate to compute the integral. Default value
#' is \code{1e4}
#'
#' @param nStrata integer vector of length 2 with the number of strata in each dimension. Default
#' value is \code{c(1, 1e2)}
#'
#' @param verbose logical (default \code{FALSE}) to report progress of the computation
#'
#' @param nThreads number (default the number of all cores, including logical cores) to use for computation
#'
#' @param alpha the significance level for accuracy measures. Default value is 0.05
#'
#' @return \code{postN0} computes the posterior mean, median, and mode of the posterior distribution
#' for each cell. The function returns a matrix with the estimates in columns and the cells in rows.
#'
#' @details The prior distributions are specified as named lists where the first component of each
#' list must be the name of distribution ('unif', 'triang', 'degen', 'gamma') and the rest of
#' components must be named according to the name of the parameters of the random generator of the
#' corresponding distribution according to:
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
#' @seealso \code{\link{rN0}}
#'
#' @examples
#' # It takes a couple of minutes
#' postN0(nMNO = 20, nReg = 115, fu = list('unif', xMin = 0.3, xMax = 0.5),
#'         fv = list('unif', xMin = 100, xMax = 120),
#'         flambda = list('gamma', shape = 11, scale = 12))
#'
#' @include rN0.R
#' @include utils.R
#' @import HDInterval
#'
#' @export
postN0 <- function(nMNO, nReg, fu, fv, flambda, n = 1e3, scale = 1, relTol = 1e-8, nSim = 1e3,
                   nStrata = c(1, 1e2), verbose = FALSE, nThreads = RcppParallel::defaultNumThreads(), alpha = 0.05){

  nCells <- length(nMNO)
  if (length(nReg) != nCells) stop('nReg and nMNO must have the same length.')

    Nvalues <- rN0(n, nMNO, nReg, fu, fv, flambda, scale, relTol, nSim, nStrata, verbose, nThreads)
    setDT(Nvalues, key  = 'cellID')
    #postMean <- round(mean(Nvalues))
    postMean<-Nvalues[, .SD[, round(mean(N0))], by = cellID][[2]]
    postSD<-Nvalues[, .SD[, round(sd(N0),2)], by = cellID][[2]]
    postCV<-Nvalues[, .SD[, round(sd(N0) / mean(N0) * 100, 2)], by = cellID][[2]]

    #postMedian <- round(median(Nvalues))
    postMedian<-Nvalues[, .SD[, round(median(N0))], by = cellID][[2]]
    postMedian_CILB<-Nvalues[, .SD[, round(equalTailedInt(N0, alpha))]['lower'], by = cellID][[2]]
    postMedian_CIUB<-Nvalues[, .SD[, round(equalTailedInt(N0, alpha))]['upper'], by = cellID][[2]]
    postMedianQuantileCV<-Nvalues[, .SD[, round( IQR(N0) / median(N0) * 100, 2)], by = cellID][[2]]

    #postMode <- Nvalues[which.max(names(table(Nvalues)))]
    postMode<-Nvalues[, .SD[, Mode(N0)], by = cellID][[2]]
    postMode_CILB<-Nvalues[, .SD[, hdi(N0, 1-alpha)]['lower'], by = cellID][[2]]
    postMode_CIUB<-Nvalues[, .SD[, hdi(N0, 1-alpha)]['upper'], by = cellID][[2]]
    postModeQuantileCV<-Nvalues[, .SD[, round( IQR(N0) / Mode(N0) * 100, 2)], by = cellID][[2]]


    output <- cbind(postMean, postSD, postCV, postMedian, postMedian_CILB, postMedian_CIUB, postMedianQuantileCV, postMode, postMode_CILB, postMode_CIUB, postModeQuantileCV)
    colnames(output)<-c("postMean", "StDev", "CV","postMedian", "Median_CI_LB", "Median_CI_UB","MedianQuantileCV", "postMode", "Mode_CI_LB", "Mode_CI_UB", "ModeQuantileCV")
    return(output)
}

