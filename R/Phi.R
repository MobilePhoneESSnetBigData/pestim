#' @title Phi
#' @description
#' @author David Salgado
#' @export
#'
Phi <- function(alpha, beta, lambda, nMNO, nReg, relTol = 1e-3){

  output <- ratioBeta(alpha, beta, nMNO, 0) * kummer(lambda, beta, alpha + beta + nMNO, relTol)
  return(output)

}

