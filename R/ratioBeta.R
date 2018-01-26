#' @title ratioBeta
#' @description
#' @author David Salgado
#' @export
#'
ratioBeta <- function(alpha, beta, nMNO, n){

  argExp <- numeric(length(alpha))
  if (length(n) == 1) n <- rep(n, length(alpha))
  if (length(nMNO) == 1) nMNO <- rep(nMNO, length(alpha))
  nonZero <- (alpha > .Machine$double.eps & beta > .Machine$double.eps)
  argExp[nonZero] <- lbeta( alpha[nonZero] + nMNO[nonZero], beta[nonZero] + n[nonZero] ) -
    lbeta( alpha[nonZero], beta[nonZero] )
  output <- exp(argExp)
  return(output)
}
