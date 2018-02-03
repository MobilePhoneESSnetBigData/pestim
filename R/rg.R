#' @title Generation of random deviates of the candidate distribution.
#'
#' @description Generate random points according to the candidate probability distribution in the
#' accept-reject method of generation of random variables applied to the distribution of the lambda
#' parameter
#'
#' @param n number of values to generate
#'
#' @param nMNO, nReg non-negative integer vectors with the number of individuals detected in each
#' cell according to the network operator and the register
#'
#' @param fu, fv named lists with the prior marginal distributions of the two-dimensional points
#' for the Monte Carlo integration
#'
#' @param flambda named list with the prior distribution of the lambda parameter
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
#' @return \code{rg} generates \code{n} points according to the candidate distribution.
#'
#' @details The candidate distribution is a gamma distribution with parameters shape = \code{nMNO} +
#' 1 and scale = \deq{\lambda^{*}} / \code{nMNO}, where \deq{\lambda^{*}} stands for the mode of the
#' posterior distribution of the lambda parameter.
#'
#' It is important to know that currently this function accepts only parameters for a single cell at
#' a time. In case of interest for the candidate density function values for a set of cells, the
#' user should program his/her own routine to apply this function to every cell.
#'
#' @seealso \code{\link{modeLambda}}, \code{\link{dlambda}} for related functions.
#'
#' @examples
#' hist(rg(1e6, nMNO = 20, nReg = 115, fu = list('unif', xMin = 0.3, xMax = 0.5),
#'         fv = list('unif', xMin = 100, xMax = 120),
#'         flambda = list('gamma', shape = 11, scale = 12)), breaks = seq(1, 200, by = 2), main ='')
#'
#' @export
#'
rg <- function(n, nMNO, nReg, fu, fv, flambda, relTol = 1e-6, nSim = 1e4, nStrata = c(1, 1e2),
               verbose = FALSE){

  nCells <- length(nMNO)
  if (nCells != 1) stop('Only one cell at a time.')
  if (length(nReg) != nCells) stop('nReg and nMNO must have the same length.')

  ######  Computing lambdaOpt
  lambdaOpt <- modeLambda(nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose)

  if (verbose) cat('Generating values... ')
  x <- rgamma(n, shape = nMNO + 1, scale = lambdaOpt / nMNO)
  if (verbose) cat(' ok.\n')
  return(x)

}
