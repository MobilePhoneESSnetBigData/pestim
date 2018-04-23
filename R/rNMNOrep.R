#' @title Generation random deviates of the replicated number of individuals according to the MNO.
#'
#' @description Generate random deviates of the number \eqn{N^{\textrm{MNO, rep}}} of individuals 
#' according to the MNO according to the original number of individuals \eqn{N^{\textrm{MNO}}}
#' detected by the MNO
#'
#' @param n number of points to generate
#'
#' @param nMNO non-negative integer vectors with the number of individuals detected according to the
#'  network operator
#' 
#' @param nReg non-negative integer vectors with the number of individuals detected according to the
#'  population register
#' 
#' @param fu named list with the prior marginal distribution of the parameter \code{u}
#'
#' @param fv named list with the prior marginal distributions of the parameter \code{v}
#'
#' @param flambda named list with the prior distribution of the parameter \eqn{\lambda}
#'
#' @param relTol relative tolerance in the computation of the \code{\link{kummer}} function. Default
#' value is \code{1e-6}
#'
#' @param nSim number of two-dimensional points to generate to compute the integral. Default value
#' is \code{1e4}
#'
#' @param nStrata integer vector of length 3 with the number of strata in each dimension. Default
#' value is \code{c(1, 1e2, 1e2)}
#'
#' @param verbose logical (default \code{FALSE}) to report progress of the computation
#'
#' @param nThreads number (default the number of all cores, including logical cores) to use for computation
#'
#' @return \code{rMNOrep} returns a \linkS4class{data.table} with the values of generated points
#' (column \code{nMNO}) for each value of the parameters \eqn{u, v, \lambda} together with priors 
#' for the hierarchical model.
#'
#'
#' @details The prior distributions are specified as named lists where the first component of each 
#' list must be the name of distribution ('unif', 'triang', 'gamma') and the rest components must be
#' named according to the name of the parameters of the random generator of the corresponding
#' distribution according to:
#'
#'   \itemize{
#'
#'     \item unif: \code{xMin}, \code{xMax} for the minimum, maximum of the sampled interval.
#'     \item degen: \code{x0} for the degenerate value of the random variable.
#'     \item triang: \code{xMin}, \code{xMax}, \code{xMode} for minimum, maximum and mode (see
#'     \code{\link{qtriang}}).
#'     \item gamma: \code{scale} and \code{shape} with the same meaning as in \code{\link{rgamma}}.
#'   }
#'
#'
#' @seealso \code{\link{rNMNO}}
#'
#' @references \url{https://github.com/MobilePhoneESSnetBigData}
#'
#' @examples
#'
#' rNMNOrep(10, nMNO = 15, nReg = 110, 
#'          fu = list('unif', xMin = 0.1, xMax = 0.2),
#'          fv = list('unif', xMin = 100, xMax = 120),
#'          flambda = list('gamma', shape = 11, scale = 110 / 10))
#'
#' @include ruvlambda.R rNMNO.R
#'
#' @import data.table
#'
#' @export
rNMNOrep <- function(n, nMNO, nReg, fu, fv, flambda, relTol = 1e-6, nSim = 1e6, 
                    nStrata = c(1, 1e2, 1e2), verbose = FALSE, 
                    nThreads = RcppParallel::defaultNumThreads()){
  
  nCells <- length(nMNO)
  if (nCells == 1) {
  
    uvlambda <- ruvlambda(n, nMNO = nMNO, nReg = nReg, 
                          fu = fu, fv = fv, flambda = flambda,
                          relTol = relTol, nSim = nSim, 
                          nStrata = nStrata, verbose = verbose, 
                          nThreads = nThreads)
    
    output <- lapply(1:n, function(i){
      
      locOut <- rNMNO(1, lambda = uvlambda[i][['lambda']], 
                      u = uvlambda[i][['u']], v = uvlambda[i][['v']])
      return(locOut)
      
    })

    output <- rbindlist(output)    
    return(output)
    
  } else {
    
    
  }  

}
