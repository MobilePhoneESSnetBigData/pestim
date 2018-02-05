#' @title estimates of the population conditional on their initial values
#' @description estimates of the population conditional on their initial values
#' @author David Salgado
#'
#'@param n number of values to generate
#'
#'@param N0 a numeric vector with the initial population estimates for each cell
#'
#'@param nMNOmat
#'
#'@param distNames a list of prior distribution names. It could contain the following names: 'unif', 'degen', 'triang' or 'gamma'
#' 
#'@param variation
#'
#'@return
#'
#'@import rmatProb.R
#'
#' @export
#'

rNtcondN0 <- function(n, N0, nMNOmat,  distNames, variation){

  nCells <- length(N0)
  if (!is.matrix(nMNOmat)) stop('nMNOmat must be a square matrix.')
  if (dim(nMNOmat)[1] != nCells) stop('The number of rows of nMNOmat must coincide with the number of cells.')
  if (dim(nMNOmat)[2] != nCells) stop('The number of columns of nMNOmat must coincide with the number of cells.')
  if (length(distNames) != nCells) stop('The length of distNames must coincide with the number of cells.')
  if (length(variation) != nCells) stop('The length of variation must coincide with the number of cells.')

  matProb <- rmatProb(n, nMNOmat, distNames, variation)

  output <- lapply(matProb, function(matP){

    outLocal <- N0 %*% matP
    outLocal <- round(outLocal)
    return(outLocal)

  })

  output <- Reduce(rbind, output)
  dimnames(output) <- NULL
  return(output)

}
