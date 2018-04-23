#' @title Generate random deviates of the number of individuals detected by the MNO.
#'
#' @description Generate random deviates of the number of individuals detected by the MNO according
#' to the hierarchical model.
#'
#' @param n number of points to generate
#'
#' @param lambda values of the parameter \eqn{\lambda}
#' 
#' @param u values of the parameter \eqn{u}
#' 
#' @param v values of the parameter \eqn{v}
#'
#' @return \code{rNMNO} returns the random deviates of the number of indviduals detected by the 
#' mobile network operator according to the hierarchical model.
#'
#'
#' @seealso \code{\link{rlambda}}, \code{\link{ruvlambda}}
#'
#' @examples
#' rNMNO(n = 1e3, lambda = c(10.3, 10.9), u = c(0.35, 0.41), v = c(100.3, 98.9))
#'
#' @import data.table 
#' 
#' @export
rNMNO <- function(n, lambda, u, v){
  
  nlambda <- length(lambda)
  nu <- length(u)
  nv <- length(v)
  if (nlambda != nu) stop('lambda and u must have the same length.')
  if (nlambda != nv) stop('lambda and v must have the same length.')
  
  nPars <- max(nlambda, nu, nv)
  if (nPars > 1 & nlambda == 1) lambda <- rep(lambda, nPars)
  if (nPars > 1 & nu == 1) u <- rep(u, nPars)
  if (nPars > 1 & nv == 1) v <- rep(v, nPars)

  output <- lapply(1:nPars, function(nPar){
    
    ui <- u[nPar]
    vi <- v[nPar]
    lambdai <- lambda[nPar]
    alpha <- ui * vi
    beta <- (1 - ui) * vi
    p <- rbeta(n, alpha, beta)
    N <- rpois(n, lambdai)
    NMNOout <- rbinom(n, N, p)
    locOutput <- data.table(lambda = lambdai, u = ui, v = vi, nMNO = NMNOout)
    return(locOutput)
  })
  output <- rbindlist(output)
  setcolorder(output, c('u', 'v', 'lambda', 'nMNO'))
  return(output)  
  
}

