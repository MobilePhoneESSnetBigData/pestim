#' @title postNestimates
#' @description
#' @author David Salgado
#' @export
#'
postNestimates <- function(nMNO, nReg, fu, fv, flambda, n = 1e3, relTol = 1e-4, nSim = 1e3, nStrata = c(1, 1e2), verbose = FALSE){

  nCells <- length(nMNO)
  if (length(nReg) != nCells) stop('nReg and nMNO must have the same length.')

  if (nCells == 1) {

    Nvalues <- rN(n, nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose)[['N']]
    postMean <- round(mean(Nvalues))
    postMedian <- round(median(Nvalues))
    postMode <- Nvalues[which.max(names(table(Nvalues)))]
    output <- list(postMean = postMean, postMedian = postMedian, postMode = postMode)
    return(output)

  } else {

    output <- sapply(seq(along = nMNO), function(i){
      postNestimates(nMNO[i], nReg[i], fu, fv, flambda, n, relTol, nSim, nStrata, verbose)
    })
    return(output)
  }
}
