#' @title dg
#' @description
#' @author David Salgado
#' @export
#'
dg <- function(lambda, nMNO, nReg, fu, fv, flambda, relTol = 1e-2, nSim = 1e3, nStrata = c(1, 1e2), verbose = FALSE){

  nCells <- length(nMNO)
  if (nCells != 1) stop('Only one cell at a time.')
  if (length(nReg) != nCells) stop('nReg and nMNO must have the same length.')

  lambdaOpt <- modeLambda(nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose)
  output <- dgamma(lambda, shape = nMNO + 1, scale = lambdaOpt / nMNO)
  return(output)

}
