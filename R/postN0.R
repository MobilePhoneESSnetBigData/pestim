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
#'
#' @export
postN0 <- function(nMNO, nReg, fu, fv, flambda, n = 1e3, scale = 1, relTol = 1e-8, nSim = 1e3,
                   nStrata = c(1, 1e2), verbose = FALSE, nThreads = RcppParallel::defaultNumThreads()){

  nCells <- length(nMNO)
  if (length(nReg) != nCells) stop('nReg and nMNO must have the same length.')

  if (nCells == 1) {

    Nvalues <- rN0(n, nMNO, nReg, fu, fv, flambda, scale, relTol, nSim, nStrata, verbose, nThreads)[['N0']]
    postMean <- round(mean(Nvalues))
    postMedian <- round(median(Nvalues))
    postMode <- Nvalues[which.max(names(table(Nvalues)))]
    output <- c(postMean = postMean, postMedian = postMedian, postMode = postMode)
    return(output)

  } else {

    output <- lapply(seq(along = nMNO), function(i){

      postN0(nMNO[i], nReg[i], fu[[i]], fv[[i]], flambda[[i]],
             n, scale, relTol, nSim, nStrata, verbose, nThreads)

    })
    output <- Reduce(rbind, lapply(output, rbind))
    return(output)
  }
}
