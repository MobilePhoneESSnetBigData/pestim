#' @title Estimates unconditioned on the initial estimate of the population in each cell
#' @description estimates unconditioned on the initial estimate of the population in each cell
#' @author David Salgado
#' #'
#' @param n number of values to generate
#'
#' @param nMNOmat
#'
#' @param fu, fv named lists with the prior marginal distributions of the two-dimensional points
#' for the Monte Carlo integration
#'
#' @param flambda named list with the prior distribution of the lambda parameter
#'
#' @param distNames a list of prior distribution names. It could contain the following names: 'unif', 'degen', 'triang' or 'gamma'
#'
#' @param variation ??
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
#' @return \code{rNt} generates \code{n} points according to the posterior distribution. The
#' function returns a \linkS4class{data.table} with these points (under the column \code{N0})
#' together with the additional variables:
#'
#'  \itemize{
#'
#'    \item The common length of \code{nMNO} and \code{nReg} identifies the number of territorial
#'    cells in which the number of individuals detected by the telecommunication network and
#'    official data. The column \code{cellID} identifies these territorial cells.
#'
#'    \item The different values of the generated values of lambda are returned under the column
#'     \code{lambda}.
#'
#'    \item The inputs \code{nMNO} and \code{nReg} are also included in the output
#'    \linkS4class{data.table} in columns under the same name.
#'
#'  }
#'
#' @details The posterior distribution is a Poisson distribution with parameter \code{lambda *
#' scale}, where the values of \code{lambda} are generated with the function \code{\link{rlambda}}.
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
#' @seealso \code{\link{dlambda}}, \code{\link{rg}} for related functions.
#'
#' @examples
#'
#' @include rN0.R 
#' @include rNtcondN0.R
#' 
#' @export
#'

rNt <- function(n, nMNOmat, nReg, fu, fv, flambda, distNames, variation, scale = 1, relTol = 1e-5, nSim = 1e3, nStrata = c(1, 1e2), verbose = FALSE){

  if (!is.matrix(nMNOmat)) stop('nMNOmat must be a square matrix.')
  nCells <- length(nReg)
  if (dim(nMNOmat)[1] != nCells) stop('The number of rows of nMNOmat must coincide with the number of cells.')
  if (dim(nMNOmat)[2] != nCells) stop('The number of columns of nMNOmat must coincide with the number of cells.')
  nMNO <- rowSums(nMNOmat)
  if (length(fu) != 1 && length(fu) != nCells) stop('The length of fu must 1 or coincide with the numbers of cells.')
  if (length(fv) != 1 && length(fv) != nCells) stop('The length of fv must 1 or coincide with the numbers of cells.')
  if (length(flambda) != 1 && length(flambda) != nCells) stop('The length of flambda must 1 or coincide with the numbers of cells.')

  DT <- rN0(n, nMNO, nReg, fu, fv, flambda, scale,  relTol, nSim, nStrata, verbose)
  DT[, lambda := NULL]
  DT[, n := rep(1:n, nCells)]
  DT[, N := rNtcondN0(1, N0, nMNOmat, distNames, variation), by = 'n']
  DT[, c('nMNO', 'nReg', 'N0') := NULL]
  setcolorder(DT, c('n', 'cellID', 'N'))
  return(DT[])
}
