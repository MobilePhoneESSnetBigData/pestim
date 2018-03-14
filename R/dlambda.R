#' @title Posterior density function of the lambda parameter.
#'
#' @description Compute the unnormalized posterior density function of the parameter \eqn{\lambda}
#' in the hierarchical model to estimate population counts given by
#' \deqn{f(\lambda\big | N^{\textrm{MNO}}; N^{\textrm{Nreg}})\propto f(\lambda)\cdot
#' \textrm{dpois}(N^{\textrm{MNO}}; \lambda)\cdot S(\lambda; N^{\textrm{MNO}}, N^{\textrm{Nreg}}), }
#'  where \code{\link{dpois}} is the probability density function of a Poisson distribution and
#'  \eqn{S} is defined in the bibliographic reference.
#'
#' @param lambda numeric vector
#'
#' @param nMNO, nReg non-negative integer vectors with the number of individuals
#' detected in each cell according to the network operator and the register
#'
#' @param fu, fv named lists with the prior marginal distributions of the
#' two-dimensional points for the Monte Carlo integration
#'
#' @param flambda named list with the prior distribution of the lambda parameter
#'
#' @param relTol relative tolerance in the computation of the \code{\link{kummer}} function. Default
#' value is \code{1e-6}
#'
#' @param nSim number of two-dimensional points to generate to compute the integral. Default value
#' is \code{1e3}
#'
#' @param nStrata integer vector of length 2 with the number of strata in each dimension. Default
#' value is \code{c(1, 1e2)}
#'
#' @param verbose logical (default \code{FALSE}) to report progress of the computation
#'
#' @param nThreads number (default the number of all cores, including logical cores) to use for computation
#'
#' @return \code{dlambda} returns a \linkS4class{data.table} with the values of the density function
#' (column \code{probLambda}) for each value of lambda together with additional variables:
#'
#'  \itemize{
#'
#'    \item The common length of \code{nMNO} and \code{nReg} identifies the number of territorial
#'    cells in which the number of individuals detected by the telecommunication network and
#'    official data. The column \code{cellID} identifies these territorial cells.
#'
#'    \item The length of \code{lambda} identifies the number of parameters upon which the integral
#'    will be computed in each cell. The column \code{parID} identifies each of these input
#'    parameters.
#'
#'    \item The inputs \code{nMNO} and \code{nReg} are also included in the output
#'    \linkS4class{data.table} in columns under the same name.
#'
#'    \item The value on the integral times the Poisson density function ifalso included under the
#'    column \code{integral}
#'
#'  }
#'
#' @details The lengths of the input vectors \code{nMNO} and \code{nReg} must be both equal to 1 and
#' independent of the length of the input vector \code{lambda}. The integral is computed using with
#' Monte Carlo techniques using \code{nSim} points for each of the values \code{lambda} specified so
#' that the final \linkS4class{data.table} has \code{length(lambda)} rows.
#'
#' The prior distributions are specified as named lists where the first component of each list must
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
#' It is important to know that currently this function accepts only parameters for a single cell at
#' a time. In case of interest for the density function values for a set of cells, the user should
#' program his/her own routine to apply this function to every cell.
#'
#' @seealso \code{\link{genUV}}, \code{\link{Phi}} for related functions.
#'
#' @references \url{https://github.com/MobilePhoneESSnetBigData}
#'
#' @examples
#' # This data.table must have 5x3= 15 rows
#' dlambda(seq(0, 1, length.out = 5),
#'         nMNO = c(20, 17, 25), nReg = c(115, 123, 119),
#'         fu = list('unif', xMin = 0.3, xMax = 0.5), fv = list('gamma', shape = 11, scale = 12),
#'         flambda = list('gamma', shape = 11, scale = 12))
#'
#' # Easily, a function to draw conditioned on the parameters:
#' f <- function(x){
#'   dlambda(x, nMNO = 20, nReg = 115,
#'           fu = list('unif', xMin = 0.3, xMax = 0.5), fv = list('unif', xMin = 100, xMax = 120),
#'           flambda = list('gamma', shape = 11, scale = 12))$probLambda
#' }
#' curve(f, xlim = c(0, 150))
#'
#' @include triang.R
#'
#' @import data.table
#'
#' @export
dlambda <- function(lambda, nMNO, nReg, fu, fv, flambda, relTol = 1e-6, nSim = 1e3, nStrata = c(1, 1e2), verbose = FALSE, nThreads = RcppParallel::defaultNumThreads()){

  if (any(lambda < 0)) stop('lambda must be nonnegative.')
  if (any(nMNO < 0)) stop('nMNO must be nonnegative.')
  if (any(nReg < 0)) stop('nReg must be nonnegative.')

  if (verbose) cat(paste0('Generating ', nSim, ' variables u and v in ', paste0(nStrata, collapse = ', '), ' strata...'))
  uv.DT <- genUV(nSim, nStrata, fu, fv, lambda, nMNO, nReg)
  if (verbose) cat(' ok.\n')
  uv.DT[, alpha := u * v]
  uv.DT[, beta := (1 - u) * v]
  if (verbose) cat(paste0('Computing Poisson factors...'))
  uv.DT[, factorPoisson := exp(-lambda + nMNO * log(lambda) - lfactorial(nMNO))]
  if (verbose) cat(paste0(' ok.\n'))
  if (verbose) cat(paste0('Computing phi factors...'))
  uv.DT[factorPoisson < .Machine$double.xmin, phiValues := 0]
  uv.DT[factorPoisson >= .Machine$double.xmin, phiValues := Phi(alpha, beta, lambda, nMNO, relTol, nThreads)]
  uv.DT[, prob := factorPoisson * phiValues]
  if (verbose) cat(paste0(' ok.\n'))
  if (verbose) cat('Computing the integral(s) by Monte Carlo approximations...')
  integral <- uv.DT[, sum(prob) / nSim, by = c('parID', 'cellID')]
  setnames(integral, 'V1', 'integral')
  uv.DT <- uv.DT[, c('cellID', 'parID', 'lambda', 'nMNO', 'nReg'), with = FALSE]
  setkeyv(uv.DT, c('cellID', 'parID'))
  uv.DT <- uv.DT[!duplicated(uv.DT, by = key(uv.DT))]
  integral <- merge(uv.DT, integral, by = c('cellID', 'parID'))

  if (verbose) cat(' ok.\n')

  if (flambda[[1L]] == 'gamma'){

    integral[, scale := flambda[['scale']]]
    integral[, shape := flambda[['shape']]]
    integral[, flambda := dgamma(lambda, scale = scale, shape = shape)]
    integral[, probLambda := integral * flambda]
    integral[, scale := NULL]
    integral[, shape := NULL]
    integral[, flambda := NULL]
    return(integral[])
  }
}
