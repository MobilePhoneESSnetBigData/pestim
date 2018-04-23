#' @title Indicators for model checking.
#'
#' @description Compute the indicators for model checking given by ...
#'
#' @param nSimPar number of simulations to compute the underlying integrals
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
#' is \code{1e6}
#'
#' @param nStrata integer vector of length 3 with the number of strata in each dimension. Default
#' value is \code{c(1, 1e2, 1e2)}
#'
#' @param verbose logical (default \code{FALSE}) to report progress of the computation
#'
#' @param nThreads number (default the number of all cores, including logical cores) to use for computation
#'
#' @return \code{modelCheckInd} returns a \linkS4class{data.table} with the values of the indicators
#' for model checking:
#'
#'  \itemize{
#'
#'    \item nMNO:
#'
#'    \item nReg:
#'
#'    \item B:
#'
#'    \item relB:
#'
#'    \item V:
#'
#'    \item relV:
#'
#'    \item MSE:
#'
#'    \item relMSE:
#'
#'  }
#'
#' @details The underlying integrals are computed using with Monte Carlo techniques using
#' \code{nSimPar} points for each of the generated random deviate values.
#'
#' The prior distributions are specified as named lists where the first component of each list must
#' be the name of distribution ('unif', 'triang', 'degen', 'gamma') and the rest components must be
#' named according to the name of the parameters of the random generator of the corresponding
#' distribution according to:
#'
#'   \itemize{
#'
#'     \item unif: \code{xMin}, \code{xMax} for the minimum, maximum of the sampled interval.
#'     \item triang: \code{xMin}, \code{xMax}, \code{xMode} for minimum, maximum and mode (see
#'     \code{\link{qtriang}}).
#'     \item gamma: \code{scale} and \code{shape} with the same meaning as in \code{\link{rgamma}}.
#'   }
#'
#'
#' @seealso \code{\link{rNMNOrep}}
#'
#' @references \url{https://github.com/MobilePhoneESSnetBigData}
#'
#' @examples
#' # Easily, a function to draw conditioned on the parameters:
#' modelCheckInd(nSimPar = 10, nMNO = 29, nReg = 123,
#'          fu = list('unif', xMin = 0.2, xMax = 0.25),
#'          fv = list('unif', xMin = 115, xMax = 130),
#'          flambda = list('gamma', shape = 21, scale = 123 / 20))
#'
#' modelCheckInd(nSimPar = 10, nMNO = c(29, 31), nReg = c(123, 119),
#'          fu = list(list('unif', xMin = 0.2, xMax = 0.25),
#'                    list('unif', xMin = 0.21, xMax = 0.26)),
#'          fv = list(list('unif', xMin = 115, xMax = 130),
#'                    list('unif', xMin = 114, xMax = 131)),
#'          flambda = list(list('gamma', shape = 21, scale = 123 / 20),
#'                         list('gamma', shape = 11, scale = 124 / 10)))
#'
#'
#' @include rNMNO.R
#' @include ruvlambda.R
#' @import data.table
#' @import parallel
#' @import doParallel
#' @import foreach
#' @export
modelCheckInd <- function(nSimPar, nMNO, nReg, fu, fv, flambda, relTol = 1e-6, nSim = 1e6,
                           nStrata = c(1, 1e2, 1e2), verbose = FALSE,
                           nThreads = RcppParallel::defaultNumThreads()){
  nCell <- length(nMNO)
  if (nCell == 1){
    uvlambda <- ruvlambda(nSimPar, nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose, nThreads)
    nMNOrep <- lapply(1:nSimPar, function(i){
      output <- rNMNO(nSimPar,
                      lambda = uvlambda[i][['lambda']],
                      u = uvlambda[i][['u']],
                      v = uvlambda[i][['v']])
      return(output)
    })

    nSim2 = nSimPar**2
    DT <- rbindlist(nMNOrep)
    setnames(DT, 'nMNO', 'nMNOrep')
    DT[, nMNOrep2 := nMNOrep ** 2]
    DT[, difnMNOrep1 := (nMNOrep - nMNO)]
    DT[, relDifnMNOrep1 := difnMNOrep1 / nMNO]
    DT[, difnMNOrep2 := difnMNOrep1 ** 2]
    DT[, relDifnMNOrep2 := (relDifnMNOrep1) ** 2]
    indicators <- DT[, list(B = sum(difnMNOrep1) / nSim2,
                            relB = sum(relDifnMNOrep1) / nSim2,
                            m2 = sum(nMNOrep2) / nSim2,
                            e2 = ( sum(nMNOrep) / nSim2 ) ** 2,
                            rele2 = ( sum(relDifnMNOrep1) / nSim2 ) ** 2,
                            MSE = sum(difnMNOrep2) / nSim2,
                            relMSE = sum(relDifnMNOrep2) / nSim2)]

    indicators[, V := m2 - e2]
    indicators[, relV := relMSE - rele2]
    indicators[, nMNO := nMNO]
    indicators[, nReg := nReg]
    indicators <- indicators[, c('nMNO', 'nReg', 'B', 'relB', 'V', 'relV', 'MSE', 'relMSE'), with = FALSE]
    return(indicators[])

  } else {

    if(nCell < nThreads) {
      cl <- makeCluster(nCell)
    } else {
      cl <- makeCluster(nThreads)
    }
    registerDoParallel(cl)
    output <- foreach(i=1:nCell, .combine = rbind, .options.snow = list(preschedule = TRUE)) %dopar% {
      outLocal <- modelCheckInd(nSimPar, nMNO[i], nReg[i], fu[[i]], fv[[i]], flambda[[i]],
                                 relTol, nSim, nStrata, verbose, nThreads)
    }
    stopCluster(cl)
    return(output)

  }
}
