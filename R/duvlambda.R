#' @title Posterior density function of the parameters u, v, and lambda.
#'
#' @description Compute the unnormalized posterior density function of the parameter \eqn{\lambda}
#' in the hierarchical model.
#'
#' @param lambda values of the parameter \eqn{\lambda}
#' 
#' @param u values of the parameter \eqn{u}
#' 
#' @param v values of the parameter \eqn{v}
#'
#' @param nMNO non-negative integer vectors with the number of individuals detected according to the
#'  network operator
#' 
#' @param fu named list with the prior marginal distributions of the hyperparameter \eqn{u}
#' 
#' @param fv named list with the prior marginal distributions of the hyperparameter \eqn{v}
#'
#' @param flambda named list with the prior marginal distribution of the parameter \eqn{\lambda}
#' 
#' @param relTol relative tolerance in the computation of the \code{\link{kummer}} function. Default
#' value is \code{1e-6}
#'
#' @param nThreads number (default the number of all cores, including logical cores) to use for 
#' computation
#'
#' @return \code{duvlambda} returns the probability mass function of the number of indviduals 
#' detected by the mobile network operator according to the hierarchical model. It depends on priors
#' for the parameters \deqn{u}, \deqn{v}, \deqn{\lambda}.
#'
#' @details The prior distributions are specified as named lists where the first component of each 
#' list must be the name of distribution ('unif', 'triang', 'gamma') and the rest components must be
#' named according to the name of the parameters of the density/probability function of the 
#' corresponding distribution according to:
#'
#'   \itemize{
#'
#'     \item unif: \code{xMin}, \code{xMax} for the minimum, maximum of the sampled interval.
#'     \item triang: \code{xMin}, \code{xMax}, \code{xMode} for minimum, maximum and mode (see
#'     \code{\link{qtriang}}).
#'     \item gamma: \code{scale} and \code{shape} with the same meaning as in \code{\link{rgamma}}.
#'   }

#'
#' @seealso \code{\link{dlambda}}
#'
#' @examples
#' f <- function(x){duvlambda(x, 0.346, 97, 19, 
#'            fu = list('unif', xMin = 0.3, xMax = 0.4),
#'            fv = list('unif', xMin = 90, xMax = 105),
#'            flambda = list('gamma', shape = 11, scale = 97 / 10))$prob}
#' curve(f, 0, 150)
#' 
#' f <- function(x){duvlambda(70, x, 97, 19, 
#'            fu = list('unif', xMin = 0.3, xMax = 0.4),
#'            fv = list('unif', xMin = 90, xMax = 105),
#'            flambda = list('gamma', shape = 11, scale = 97 / 10))$prob}
#' curve(f, 0, 1)
#' 
#' f <- function(x){duvlambda(70, 0.35, x, 19, 
#'            fu = list('unif', xMin = 0.3, xMax = 0.4),
#'            fv = list('unif', xMin = 90, xMax = 105),
#'            flambda = list('gamma', shape = 11, scale = 97 / 10))$prob}
#' curve(f, 80, 115)

#' 
#' @include Phi.R triang.R
#' @import doParallel
#' @import foreach
#' @export
duvlambda <- function(lambda, u, v, nMNO, fu, fv, flambda, relTol = 1e-6, 
                       nThreads = RcppParallel::defaultNumThreads()){
  
  alpha <- u * v
  beta <- (1 - u) * v
  if (fu[[1]] == 'unif') {
    
    fuValues <- dunif(u, min = fu[['xMin']], max = fu[['xMax']])
      
  }
  
  if (fu[[1]] == 'triang') {
    
    fuValues <- dtriang(u, xMin = fu[['xMin']], xMax = fu[['xMax']], xMode = fu[['xMode']])
    
  }
  
  if (fv[[1]] == 'unif') {
    
    fvValues <- dunif(v, min = fv[['xMin']], max = fv[['xMax']])
    
  }
  
  if (fv[[1]] == 'triang') {
    
    fvValues <- dtriang(v, xMin = fv[['xMin']], xMax = fv[['xMax']], xMode = fv[['xMode']])
    
  }
  
  if (fv[[1]] == 'gamma') {
    
    fvValues <- dgamma(v, shape = fv[['shape']], scale = fv[['scale']])
    
  }
  
  if (flambda[[1]] == 'gamma') {
    
    flambdaValues <- dgamma(lambda, shape = flambda[['shape']], scale = flambda[['scale']])
    
  }
  
  DT <- data.table(lambda = lambda, nMNO = nMNO, u = u, v= v, alpha = alpha, beta = beta, 
                   fuValues = fuValues, fvValues = fvValues, flambdaValues = flambdaValues)
  DT[, PoissonValues := exp(-lambda + nMNO * log(lambda) - lfactorial(nMNO))]
  DT[PoissonValues < .Machine$double.xmin, phiValues := 0]
  DT[PoissonValues >= .Machine$double.xmin, 
     phiValues := Phi(alpha, beta, lambda, nMNO, relTol, nThreads)]
  DT[ , phiValues := phiValues * PoissonValues]
  DT[, prob := phiValues * fuValues * fvValues * flambdaValues]
  DT <- DT[, c('lambda', 'u', 'v', 'nMNO', 'prob'), with = FALSE]
  return(DT[])  
  
}

