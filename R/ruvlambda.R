#' @title Generation of random deviates of the posterior distribution of parameters u,v, lambda.
#'
#' @description Generate random points according to the posterior probability distribution of the
#' parameters u,v, and lambda in the hierarchical model
#'
#' @param n number of values to generate
#'
#' @param nMNO non-negative integer vector with the number of individuals detected in each cell 
#' according to the network operator
#' 
#' @param nReg non-negative integer vector with the number of individuals detected in each cell 
#' according to the register
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
#' @return \code{ruvlambda} generates \code{n} points according to the posterior distribution of
#' the parameters \eqn{u, v, \lambda}. The function returns a \linkS4class{data.table} with these 
#' points.
#'
#' @details The points are generated according to the accept-reject method using as candidate
#' distribution the unnormalised distribution given by 
#' \deqn{\mathbb{P}(u, v, \lambda|N^{\textrm{MNO}})\propto\mathbb{P}(u,v|\lambda, N^{\textrm{MNO}})
#' \cdot\mathbb{P}(\lambda|N^{\textrm{MNO}})}.
#'
#' The prior distributions are specified as named lists where the first component of each list must
#' be the name of distribution ('unif', 'triang', 'gamma') and the rest components must be
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
#' @seealso \code{\link{rlambda}}, \code{\link{duvlambda}}
#'
#' @examples
#' ruvlambda(10, nMNO = 20, nReg = 115, 
#'           fu = list('unif', xMin = 0.10, xMax = 0.20),
#'           fv = list('unif', xMin = 100, xMax = 120),
#'           flambda = list('gamma', shape = 11, scale = 12))
#'           
#'           
#' ruvlambda(10, nMNO = c(19, 20), nReg = c(115, 117), 
#'           fu = list(list('unif', xMin = 0.10, xMax = 0.20),
#'                     list('unif', xMin = 0.11, xMax = 0.22)),
#'           fv = list(list('unif', xMin = 100, xMax = 120),
#'                     list('unif', xMin = 97, xMax = 123)),
#'           flambda = list(list('gamma', shape = 11, scale = 105 / 10),
#'                          list('gamma', shape = 12, scale = 1107 / 11)))
#'           
#' @include rlambda.R duvlambda.R triang.R
#'
#' @import data.table
#'
#' @export
ruvlambda <- function(n, nMNO, nReg, fu, fv, flambda, relTol = 1e-6, nSim = 1e6,
                      nStrata = c(1, 1e2, 1e2), verbose = FALSE, 
                      nThreads = RcppParallel::defaultNumThreads()){
  
  nCells <- length(nMNO)
  if (length(nReg) != nCells) stop('nReg and nMNO must have the same length.')
  
  if (nCells == 1) {
    
    if (verbose) cat('Generating random values of lambda... ')    
    lambda0s <- rlambda(n, nMNO, nReg, fu, fv, flambda, relTol, nSim = 1e4, nStrata[1:2], verbose, nThreads)
    if (verbose) cat(' ok.\n')
    
    flist <- lapply(lambda0s, function(lambda0){
      
      function(u, v){ duvlambda(lambda0, u, v, nMNO, fu, fv, flambda, relTol, nThreads)$prob }
      
    })

    if (fu[[1]] == 'unif') {
      
      uMin <- fu[['xMin']]
      uMax <- fu[['xMax']]
      uMode <- (uMin + uMax) / 2
      gu <- function(u){ dunif(u, min = uMin, max = uMax) }
      
    
    }
    
    if (fu[[1]] == 'triang') {
      
      uMin <- fu[['xMin']]
      uMax <- fu[['xMax']]
      uMode <- fu[['xMode']]
      gu <- function(u){ dtriang(u, xMin = uMin, xMax = uMax, xMode = uMode) }
      
    }
    
    if (fv[[1]] == 'unif') {
      
      vMin <- fv[['xMin']]
      vMax <- fv[['xMax']]
      vMode <- (vMin + vMax) / 2
      gv <- function(v){ dunif(v, min = vMin, max = vMax) }
      
    }
    
    if (fv[[1]] == 'triang') {
   
      vMin <- fv[['xMin']]
      vMax <- fv[['xMax']]
      vMode <- fv[['xMode']]
      gv <- function(v){ dtriang(v, xMin = vMin, xMax = vMax, xMode = vMode) }
      
    }
    
    if (fv[[1]] == 'gamma') {
      
      vshape <- fv[['shape']]
      vscale <- fv[['scale']]
      vMode <- (vshape - 1) * vscale
      vMin <- vMode * 0.85
      vMax <- vMode * 1.15
      gv <- function(v){ dgamma(v, shape = vshape, scale = vscale) }
      
    }
    
    ###### Computing rejection rate  #####################
    g <- function(u, v){ gu(u) * gv(v) }
        

    optimC <- sapply(seq(along = flist), function(i){
      
      fun <- function(x){
      
        output <- numeric(1)
        guv <- g(x[1], x[2])
        if (guv >= .Machine$double.xmin) output <- g(x[1], x[2]) / flist[[i]](x[1], x[2])  
        return(-output)
      }
      output <- -1 / optim(c(uMode, vMode), fun)$value
      return(output)
      
    })
    if (verbose) cat(' ok.\n')

    if (verbose) cat('Generating and accepting/rejecting values...\n')

    if (fu[[1]] == 'unif') {
      
      ru <- runif(length(lambda0s), min = uMin, max = uMax)
      
    }
    
    if (fu[[1]] == 'triang') {
      
      ru <- rtriang(length(lambda0s), xMin = uMin, xMax = uMax, xMode = uMode)
      
    }
    
    if (fv[[1]] == 'unif') {
      
      rv <- runif(length(lambda0s), min = vMin, max = vMax)
      
    }
    
    if (fv[[1]] == 'triang') {

      rv <- rtriang(length(lambda0s), xMin = vMin, xMax = vMax, xMode = vMode)
      
    }
    
    if (fv[[1]] == 'gamma') {
      
      rv <- rgamma(length(lambda0s), shape = vshape, scale = vscale)
      
    }
    
    ruv0 <- cbind(u = ru, v = rv, lambda = lambda0s)
    v <- runif(n)
    if (verbose) cat('   of target distribution...\n')
    fx <- sapply(seq(along = flist), function(i) do.call(flist[[i]], list(u = ruv0[i, 1], v = ruv0[i, 2])))
    if (verbose) cat('   ok.\n')
    if (verbose) cat('   of candidate distribution...\n')
    #gx <- dcauchy(x, location = location, scale = scale)
    gx <- sapply(seq(along = flist), function(i) do.call(g, list(u = ruv0[i, 1], v = ruv0[i, 2]))) 
    if (verbose) cat('   ok.\n')
    indexOut <- (v <= fx / (optimC * gx))
    
    if (nrow(ruv0) == 1) {
      
      output <- t(as.matrix(ruv0[indexOut, ]))
    
    } else {
      
      output <- ruv0[indexOut, ]
      
    }

    while (nrow(output) < n) {
      
      lambda0s <- rlambda(n - nrow(output), nMNO, nReg, fu, fv, flambda, relTol, nSim = 1e4, nStrata[1:2], verbose, nThreads)
      
      flist <- lapply(lambda0s, function(lambda0){
        
        function(u, v){ duvlambda(lambda0, u, v, nMNO, fu, fv, flambda, relTol, nThreads)$prob }
        
      })
      
      optimC <- sapply(seq(along = flist), function(i){
        
        fun <- function(x){
          
          output <- numeric(1)
          guv <- g(x[1], x[2])
          if (guv >= .Machine$double.xmin) output <- g(x[1], x[2]) / flist[[i]](x[1], x[2])  
          return(-output)
        }
        output <- -1 / optim(c(uMode, vMode), fun)$value
        return(output)
        
      })
      
      if (verbose) cat(' ok.\n')
      
      if (verbose) cat('Generating and accepting/rejecting values...\n')
      
      if (fu[[1]] == 'unif') {
        
        ru <- runif(length(lambda0s), min = uMin, max = uMax)
        
      }
      
      if (fu[[1]] == 'triang') {
        
        ru <- rtriang(length(lambda0s), xMin = uMin, xMax = uMax, xMode = uMode)
        
      }
      
      if (fv[[1]] == 'unif') {
        
        rv <- runif(length(lambda0s), min = vMin, max = vMax)
        
      }
      
      if (fv[[1]] == 'triang') {
        
        rv <- rtriang(length(lambda0s), xMin = vMin, xMax = vMax, xMode = vMode)
        
      }
      
      if (fv[[1]] == 'gamma') {
        
        rv <- rgamma(length(lambda0s), shape = vshape, scale = vscale)
        
      }
      
      ruv0 <- cbind(u = ru, v = rv, lambda = lambda0s)
      v <- runif(n - nrow(output))
      if (verbose) cat('   of target distribution...\n')
      fx <- sapply(seq(along = flist), function(i) do.call(flist[[i]], list(u = ruv0[i, 1], v = ruv0[i, 2])))
      if (verbose) cat('   ok.\n')
      if (verbose) cat('   of candidate distribution...\n')
      gx <- sapply(seq(along = flist), function(i) do.call(g, list(u = ruv0[i, 1], v = ruv0[i, 2]))) 
      if (verbose) cat('   ok.\n')
      indexOut <- (v <= fx / (optimC * gx))
      
      if (nrow(ruv0) == 1) {
        
        aux <- t(as.matrix(ruv0[indexOut, ]))
        
      } else {
        
        aux <- ruv0[indexOut, ]
        
      }
      
      output <- rbind(output, aux)
      
    }
    
    if (verbose) cat(paste0(length(output), ' points selected.\n'))
    if (verbose) cat(' ok.\n')
    output <- as.data.table(output)
    output[, nMNO := rep(nMNO, each = n)]
    output[, nReg := rep(nReg, each = n)]
    setcolorder(output, c('nMNO', 'nReg', 'u', 'v', 'lambda'))
    return(output[])
    
  } else {
    
    output <- lapply(seq(along = nMNO), function(i){
      
      ruvlambda(n, nMNO[i], nReg[i], fu[[i]], fv[[i]], flambda[[i]], relTol, nSim, nStrata, verbose,  nThreads)
      
    })
    
    output <- rbindlist(output)
    return(output[])
  }
}
