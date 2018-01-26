#' @title rN
#' @description
#' @author David Salgado
#' @export
#'
rN <- function(n, nMNO, nReg, fu, fv, flambda, relTol = 1e-4, nSim = 1e3, nStrata = c(1, 1e2), verbose = FALSE){

  nCell <- length(nMNO)
  if (nCell != length(nReg)) stop('nMNO and nReg must have the same length.')
  lambda <- rlambda(n, nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose)
  DT <- data.table(lambda = lambda)
  DT[ , lambdaID := 1:.N]
  nMNO <- DT[, list(nMNO = rep(nMNO, each = .N), nReg = rep(nReg, each = .N), cellID = rep(seq(along = nMNO))), by = 'lambdaID']
  DT <- merge(DT, nMNO, by = 'lambdaID', allow.cartesian = TRUE)
  DT[, N := rpois(1, lambda), by = c('cellID', 'lambdaID')]
  DT[, lambdaID := NULL]
  DT[, cellID := NULL]
  setcolorder(DT, c('nMNO', 'nReg', 'lambda', 'N'))
  return(DT)
}
