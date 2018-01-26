#' @title dlambda
#' @description
#' @author David Salgado
#' @export
#'
dlambda <- function(lambda, nMNO, nReg, fu, fv, flambda, relTol = 1e-4, nSim = 1e3, nStrata = c(1, 1e2), verbose = FALSE){

  if (any(lambda < 0)) stop('lambda must be nonnegative.')
  if (any(nMNO < 0)) stop('nMNO must be nonnegative.')
  if (any(nReg < 0)) stop('nReg must be nonnegative.')

  if (verbose) cat(paste0('Generating ', nSim, ' variables u and v in ', paste0(nStrata, collapse = ', '), ' strata...'))
  uv.DT <- genUV(nSim, nStrata, fu, fv, lambda, nMNO, nReg)
  if (verbose) cat(' ok.\n')
  uv.DT[, alpha := u * v]
  uv.DT[, beta := (1 - u) * v]
  if (verbose) cat(paste0('Computing Poisson factors...'))
  uv.DT[, factorPoisson := exp(-lambda + nMNO * log(lambda) - lfactorial(nMNO))]
  if (verbose) cat(paste0(' ok.\n'))
  if (verbose) cat(paste0('Computing phi factors...'))
  uv.DT[factorPoisson < .Machine$double.xmin, phiValues := 0]
  uv.DT[factorPoisson >= .Machine$double.xmin, phiValues := Phi(alpha, beta, lambda, nMNO, nReg, relTol)]
  uv.DT[, prob := factorPoisson * phiValues]
  if (verbose) cat(paste0(' ok.\n'))
  if (verbose) cat('Computing the integral(s) by Monte Carlo approximations...')
  integral <- uv.DT[, sum(prob) / nSim, by = c('parID', 'cellID')]
  setnames(integral, 'V1', 'integral')
  uv.DT <- uv.DT[, c('cellID', 'parID', 'lambda', 'nMNO', 'nReg'), with = FALSE]
  setkeyv(uv.DT, c('cellID', 'parID'))
  uv.DT <- uv.DT[!duplicated(uv.DT, by = key(uv.DT))]
  integral <- merge(uv.DT, integral, by = c('cellID', 'parID'))

  if (verbose) cat(' ok.\n')

  if (flambda[[1L]] == 'gamma'){

    integral[, scale := flambda[['scale']]]
    integral[, shape := flambda[['shape']]]
    integral[, flambda := dgamma(lambda, scale = scale, shape = shape)]
    integral[, probLambda := integral * flambda]
    integral[, scale := NULL]
    integral[, shape := NULL]
    integral[, flambda := NULL]
    return(integral[])
  }

  #integral[, cellID := NULL]
  #integral[, parID := NULL]


}
