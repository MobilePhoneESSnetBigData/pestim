#' @title rg
#' @description
#' @author David Salgado
#' @export
#'
rg <- function(n, nMNO, nReg, fu, fv, flambda, relTol = 1e-2, nSim = 1e3, nStrata = c(1, 1e2), verbose = FALSE){

  nCells <- length(nMNO)
  if (nCells != 1) stop('Only one cell at a time.')
  if (length(nReg) != nCells) stop('nReg and nMNO must have the same length.')

  ######  Computing betaOpt
  lambdaOpt <- modeLambda(nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose)

  if (verbose) cat('Generating values... ')
  x <- rgamma(n, shape = nMNO + 1, scale = lambdaOpt / nMNO)
  if (verbose) cat(' ok.\n')
  return(x)

}
