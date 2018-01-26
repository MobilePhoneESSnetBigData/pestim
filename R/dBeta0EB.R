#' @title dBeta0EB
#' @description
#' @author David Salgado
#' @export
#'
dBeta0EB <- function(beta0, nMNO, nReg, uMode, vMode, relTol = 1e-2){

  DT <- data.table(beta0 = beta0)
  DT[ , betaID := 1:.N]
  nMNO <- DT[, list(nMNO = rep(nMNO, each = .N), nReg = rep(nReg, each = .N), cellID = rep(seq(along = nMNO))), by = 'betaID']
  DT <- merge(DT, nMNO, by = 'betaID', allow.cartesian = TRUE)
  DT[, factorPoisson := exp(-beta0 * nReg + nMNO * log(beta0 * nReg) - lfactorial(nMNO))]
  DT[, factorKummer :=  kummer(beta0 * nReg, (1 - uMode) * vMode, vMode + nMNO, relTol)]
  DT[, factorGamma := exp(lgamma(uMode * vMode + nMNO) - lgamma(uMode * vMode))]
  DT[, probBeta0 := factorPoisson * factorKummer * factorGamma]
  DT <- DT[, c('beta0', 'nMNO', 'nReg', 'probBeta0'), with = FALSE]
  return(DT)

}
