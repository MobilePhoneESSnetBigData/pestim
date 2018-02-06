#' @title factorBayes
#' @description used for testing purposes only
#' @author David Salgado
#' @export
#'
factorBayes <- function(beta0, nMNO, nReg, uMode, vMode, relTol = 1e-2){

  kummer(beta0 * nReg, (1 - uMode) * vMode, vMode + nMNO) * exp(lgamma(uMode * vMode + nMNO) - lgamma(uMode * vMode))

}
