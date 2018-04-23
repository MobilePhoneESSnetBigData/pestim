#' @title Generation of random deviates of the posterior distribution of initial population counts.
#'
#' @description Generate random points according to the posterior probability distribution of the
#' number of individuals in the hierarchical model.
#'
#' @param n number of values to generate
#'
#' @param nMNO, nReg non-negative integer vectors with the number of individuals detected in each
#' cell according to the network operator and the register
#'
#' @param fu, fv named lists with the prior marginal distributions of the two-dimensional points
#' for the Monte Carlo integration
#'
#' @param flambda named list with the prior distribution of the lambda parameter
#'
#' @param scale numeric vector with the scale to count the number of individuals. Default value is 1
#'
#' @param relTol relative tolerance in the computation of the \code{\link{kummer}} function. Default
#' value is \code{1e-6}
#'
#' @param nSim number of two-dimensional points to generate to compute the integral. Default value
#' is \code{1e4}
#'
#' @param nStrata integer vector of length 2 with the number of strata in each dimension. Default
#' value is \code{c(1, 1e2)}
#'
#' @param verbose logical (default \code{FALSE}) to report progress of the computation
#'
#' @param nThreads number (default the number of all cores, including logical cores) to use for computation
#'
#' @return \code{rN0} generates \code{n} points according to the posterior distribution. The
#' function returns a \linkS4class{data.table} with these points (under the column \code{N0})
#' together with the additional variables:
#'
#'  \itemize{
#'
#'    \item The common length of \code{nMNO} and \code{nReg} identifies the number of territorial
#'    cells in which the number of individuals detected by the telecommunication network and
#'    official data. The column \code{cellID} identifies these territorial cells.
#'
#'    \item The different values of the generated values of lambda are returned under the column
#'     \code{lambda}.
#'
#'    \item The inputs \code{nMNO} and \code{nReg} are also included in the output
#'    \linkS4class{data.table} in columns under the same name.
#'
#'  }
#'
#' @details The posterior distribution is a Poisson distribution with parameter \code{lambda *
#' scale}, where the values of \code{lambda} are generated with the function \code{\link{rlambda}}.
#'
#' The prior distributions are specified as named lists where the first component of each list must
#' be the name of distribution ('unif', 'triang', 'degen', 'gamma') and the rest components must be
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
#' @seealso \code{\link{rlambda}}, \code{\link{rg}}, \code{\link{rNt}} for related functions.
#'
#' @examples
#' # It takes a couple of minutes
#' hist(rN0(500, nMNO = 20, nReg = 115, fu = list('unif', xMin = 0.3, xMax = 0.5),
#'         fv = list('unif', xMin = 100, xMax = 120),
#'         flambda = list('gamma', shape = 11, scale = 12))$N0,
#'         breaks = seq(1, 200, by = 1), main ='', xlab = 'number of individuals')
#'
#' @include rlambda.R
#' @import parallel
#' @import doParallel
#' @import foreach
#' @export
rN0 <- function(n, nMNO, nReg, fu, fv, flambda, scale = 1, relTol = 1e-6, nSim = 1e4,
                nStrata = c(1, 1e2), verbose = FALSE, nThreads = RcppParallel::defaultNumThreads()){

  nCell <- length(nMNO)
  if (nCell != length(nReg)) stop('nMNO and nReg must have the same length.')

  #  lambda <- rlambda(n, nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose)
  #  DT <- data.table(lambda = as.vector(lambda))
  #  DT[, lambdaID := rep(1:nCell, each = n)]
  #  nMNO <- DT[, list(nMNO = rep(nMNO, each = n), nReg = rep(nReg, each = n), cellID = rep(seq(along = nMNO), each = n)), by = 'lambdaID']
  #  DT <- merge(DT, nMNO, by = 'lambdaID', allow.cartesian = TRUE)
  #  DT[, N := rpois(1, lambda), by = c('cellID', 'lambdaID')]
  #  return(DT[])


  if (nCell == 1){

    lambda <- rlambda(n, nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose, nThreads)
    DT <- data.table(lambda = lambda)
    DT[ , lambdaID := 1:.N]
    nMNO <- DT[, list(nMNO = rep(nMNO, each = .N), nReg = rep(nReg, each = .N), cellID = rep(seq(along = nMNO))), by = 'lambdaID']
    DT <- merge(DT, nMNO, by = 'lambdaID', allow.cartesian = TRUE)
    DT[, N0 := rpois(1, lambda * scale), by = c('cellID', 'lambdaID')]
    DT[, nMNO := nMNO * scale]
    DT[, nReg := nReg * scale]
    DT[, lambdaID := NULL]
    #DT[, cellID := NULL]
    setcolorder(DT, c('cellID', 'nMNO', 'nReg', 'lambda', 'N0'))

    return(DT[])

  }
  else {
      if(nCell < nThreads) {
        cl <- makeCluster(nCell)
      } else {
        cl <- makeCluster(nThreads)
      }
      registerDoParallel(cl)
      output <- foreach(i=1:nCell, .combine = rbind, .options.snow = list(preschedule = TRUE)) %dopar% {
        outLocal <- rN0(n, nMNO[i], nReg[i], fu[[i]], fv[[i]], flambda[[i]], scale, relTol, nSim, nStrata, verbose, nThreads)
        outLocal[, cellID := i]
        #setDT(outLocal, key  = 'cellID')
      }
      stopCluster(cl)
      setDT(output, key  = 'cellID')
      #output <- rbindlist(output)
      return(output)
    }
}

