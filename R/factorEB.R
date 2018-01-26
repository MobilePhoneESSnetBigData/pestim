#' @title factorEB
#' @description
#' @author David Salgado
#' @export
#'
factorEB <- function(beta0, nMNO, nReg, uMode, vMode){

  ((uMode * vMode - 1) / (vMode - 2))^nMNO * exp(beta0 * nReg *( (1 - uMode) * vMode ) / (vMode - 2))

}
