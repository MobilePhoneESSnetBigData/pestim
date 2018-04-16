#' @title Density function of a candidate distribution in the accept-reject method.
#'
#' @description Generate values of a candidate distribution density function in the accept-reject
#' method of generation of random variables applied to the distribution of the lambda parameter
#'
#' @param lambda numeric vector with the lambda parameter values
#'
#' @param nMNO non-negative integer vectors with the number of individuals detected in each cell 
#' according to the network operator
#'
#' @param nReg non-negative integer vectors with the number of individuals detected in each cell 
#' according to the population register
#'
#' @param fu named list with the prior marginal distribution for the parameter \eqn{u}
#'
#' @param fv named list with the prior marginal distribution for the parameter \eqn{v}
#'
#' @param flambda named list with the prior distribution of the parameter \eqn{\lambda}
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
#' @return \code{dg} generates \code{length(lambda)} values of the density probability function of
#' the candidate distribution in the accept-reject method.
#'
#' @details The candidate distribution is a gamma distribution with parameters shape = \code{nMNO} +
#' 1 and scale = \eqn{\lambda^{*}} / \code{nMNO}, where \eqn{\lambda^{*}} stands for the mode of the
#' posterior distribution of the lambda parameter.
#'
#' It is important to know that currently this function accepts only parameters for a single cell at
#' a time. In case of interest for the candidate density function values for a set of cells, the
#' user should program his/her own routine to apply this function to every cell.
#'
#' @seealso \code{\link{modeLambda}}, \code{\link{dlambda}} for related functions.
#'
#' @examples
#' curve(dg(x, nMNO = 20, nReg = 115, fu = list('unif', xMin = 0.3, xMax = 0.5),
#'         fv = list('unif', xMin = 100, xMax = 120),
#'         flambda = list('gamma', shape = 11, scale = 12)), xlim = c(0, 150),
#'         main = '', ylab = 'density', xlab = 'lambda')
#' 
#' @include modeLambda.R
#' 
#' @export
dg <- function(lambda, nMNO, nReg, fu, fv, flambda, relTol = 1e-6, nSim = 1e4, nStrata = c(1, 1e2),
               verbose = FALSE){

  nCells <- length(nMNO)
  if (nCells != 1) stop('Only one cell at a time.')
  if (length(nReg) != nCells) stop('nReg and nMNO must have the same length.')

  lambdaOpt <- modeLambda(nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose)
  output <- dgamma(lambda, shape = nMNO + 1, scale = lambdaOpt / nMNO)
  return(output)

}
