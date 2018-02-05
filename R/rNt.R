#' @title kummmer
#' @description postestimates of
#' @author David Salgado
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
