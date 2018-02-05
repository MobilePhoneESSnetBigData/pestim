#' @title rmatProb
#' @description postestimates of
#' @author David Salgado
#' 
#' @param n
#' 
#' @param nMNOmat
#' 
#' @param distNames a list of prior distribution names. It could contain the following names: 'unif', 'degen', 'triang' or 'gamma'
#' 
#' @param variation
#' 
#' @return 
#' 
#' @include alphaPrior.R
#' 
#' @export
#'

rmatProb <- function(n, nMNOmat, distNames, variation){

  if (!is.matrix(nMNOmat)) stop('nMNOmat must be a square matrix.')
  nCells <- dim(nMNOmat)[1]
  if (dim(nMNOmat)[2] != nCells) stop('nMNOmat must be a square matrix.')
  if (!all(nMNOmat >= 0)) stop('nMNOmat must have nonnegative values.')
  if (length(distNames) != nCells) stop('The length of distNames must coincide with the number of cells.')
  if (length(variation) != nCells) stop('The length of variation must coincide with the number of cells.')



  simList <- lapply(1:nCells, function(i){

    flist <- alphaPrior(nMNOmat[i,], distNames, variation)
    outLocal <- data.table(rp(n, flist))
    outLocal[, sim := 1:n]
    return(outLocal)

  })

  output <- rbindlist(simList)
  output <- split(output, by = 'sim', keep.by = FALSE)
  output <- lapply(output, as.matrix)
  dimnames(output) <- NULL
  return(output)

}
