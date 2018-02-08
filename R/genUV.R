#' @title Generation of two-dimensional random deviates.
#'
#' @description Generate two-dimensional random deviates for a Monte Carlo computation of the
#' integral \deqn{\int_{0}^{\infty}dv f_{2}(v)\int_{0}^{\infty} f_{1}(u)\ \Phi(u\cdot v, (1 - u)
#' \cdot v; \lambda, N^{\textrm{MNO}}, N^{\textrm{Reg}}).}
#'
#'  The Monte Carlo technique makes use of stratified importance sampling.
#'
#' @param nSim number of two-dimensional points to generate
#'
#' @param nStrata integer vector of length 2 with the number of strata in each dimension
#'
#' @param f1, f2 named lists with the prior marginal distributions of the two-dimensional points
#'
#' @param lambda numeric vector
#'
#' @param nMNO, nReg non-negative integer vectors
#'
#' @return \code{genUV} returns a \linkS4class{data.table} with the (u,v) coordinates of each point
#' together with additional variables:
#'
#'  \itemize{
#'
#'    \item The common length of \code{nMNO} and \code{nReg} identifies the number of territorial
#'    cells in which the number of individuals detected by the telecommunication network and
#'    official data. The column \code{cellID} identifies these territorial cells.
#'
#'    \item The length of \code{lambda} identifies the number of parameters upon which the integral
#'    will be computed in each cell. The column \code{parID} identifies each of these input
#'    parameters.
#'
#'    \item \code{Stratum_u} and \code{Stratum_v} jointly identify each stratum in which the region
#'    of integration has been divided with the stratification.
#'  }
#'
#' @details The lengths of the input vectors \code{nMNO} and \code{nReg} must be equal and
#' independent of the length of the input vector \code{lambda}. Notice that \code{nSim} points are
#' generated for each of the \code{length(nMNO)}\eqn{\times}\code{length(lambda)} combinations so
#' that the final \linkS4class{data.table} has \code{nSim}\eqn{\times}\code{length(nMNO)}
#' \eqn{\times}\code{length(lambda)} rows.
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
#' @seealso \code{\link{runif}}, \code{\link{qtriang}}, \code{\link{rgamma}} for related functions.
#'
#' @examples
#' # This data.table must have 10x5x3= 150 rows and only one stratum
#' genUV(nSim = 10, nStrata = c(1, 1),
#'       f1 = list('unif', xMin = 0.3, xMax = 0.5), f2 = list('gamma', shape = 11, scale = 12),
#'       lambda = seq(0, 1, length.out = 5),
#'       nMNO = c(20, 17, 25), nReg = c(115, 123, 119))
#'
#' @include triang.R
#'
#' @import data.table
#'
#' @export
genUV <- function(nSim, nStrata, f1, f2, lambda, nMNO, nReg){

  if (length(nStrata) != 2) stop('nStrata must have length 2.')
  if (nStrata[1] >= 1 && (nSim %% nStrata[1] != 0)) stop('nSim must be a multiple of the number of strata of variable 1.')
  if (nStrata[2] >= 1 && (nSim %% nStrata[2] != 0)) stop('nSim must be a multiple of the number of strata of variable 2.')
  if (nStrata[1] >= 1 && nSim / nStrata[1] < 10) stop('There must be at least 10 points per stratum of variable 1.')
  if (nStrata[2] >= 1 && nSim / nStrata[2] < 10) stop('There must be at least 10 points per stratum of variable 2.')
  if (nStrata[1] == 0) nStrata[1] <- 1
  if (nStrata[2] == 0) nStrata[2] <- 1

  if (!is.list(f1)) stop('f1 must be a list.')
  if (!is.list(f2)) stop('f2 must be a list.')
  if (length(nMNO) != length(nReg)) stop('nMNO and nReg must have the same length.')

  nCells <- length(nMNO)
  nlambdas <- length(lambda)
  DT <- as.data.table(expand.grid(lambda, nMNO))
  setnames(DT, c('lambda', 'nMNO'))
  DT[, cellID := rep(1:nCells, each = nlambdas)]
  #DT[, lambdaID := rep(1:nlambdas, nCells)]
  DT[, nReg := rep(nReg, each = nlambdas)]
  DT[, parID := rep(1:nlambdas, nCells)]
  setcolorder(DT, c('cellID', 'parID', 'lambda', 'nMNO', 'nReg'))

  if (nStrata[1] > 1 & nStrata[2] > 1) {

    Strat <- DT[, as.list(expand.grid(1:nStrata[1], 1:nStrata[2])), by = c('parID', 'cellID')]
    setnames(Strat, c('Var1', 'Var2'), c('Stratum_u', 'Stratum_v'))
    DT <- merge(DT, Strat, by = c('parID', 'cellID'), allow.cartesian = TRUE)
  }

  if (nStrata[1] > 1 & nStrata[2] <= 1) {

    Strat <- DT[, as.list(expand.grid(1:nStrata[1], 1)), by =  c('parID', 'cellID')]
    setnames(Strat, c('Var1', 'Var2'), c('Stratum_u', 'Stratum_v'))
    DT <- merge(DT, Strat, by =  c('parID', 'cellID'), allow.cartesian = TRUE)

  }

  if (nStrata[1] <= 1 & nStrata[2] > 1) {

    Strat <- DT[, as.list(expand.grid(1, 1:nStrata[2])), by =  c('parID', 'cellID')]
    setnames(Strat, c('Var1', 'Var2'), c('Stratum_u', 'Stratum_v'))
    DT <- merge(DT, Strat, by =  c('parID', 'cellID'), allow.cartesian = TRUE)
  }

  if (nStrata[1] <= 1 & nStrata[2] <= 1) {

    Strat <- DT[, as.list(expand.grid(1, 1)), by =  c('parID', 'cellID')]
    setnames(Strat, c('Var1', 'Var2'), c('Stratum_u', 'Stratum_v'))
    DT <- merge(DT, Strat, by =  c('parID', 'cellID'), allow.cartesian = TRUE)
  }

  setcolorder(DT, c('cellID', 'parID', 'nMNO', 'nReg', 'lambda', 'Stratum_u', 'Stratum_v'))
  setkeyv(DT, c('cellID', 'parID', 'nMNO', 'nReg', 'lambda', 'Stratum_u', 'Stratum_v'))

  if (f1[[1L]] == 'unif') {

    uMin <- f1[['xMin']]
    if (any(uMin < 0)) stop('The parameters uMin in the distribution f1 cannot be negative.')
    if (any(uMin >= 1)) stop('The parameters uMin in the distribution f1 cannot be equal to or larger than 1.')
    if (length(uMin) != 1 & length(uMin) != nCells) stop('The number of parameters uMin in the distribution f1 does not coincide with the number of cells.')
    if (length(uMin) == nCells) uMin <- rep(uMin, each = dim(DT)[1] / nCells)
    DT[, uMin := uMin]

    uMax <- f1[['xMax']]
    if (length(uMax) != 1 & length(uMax) != nCells) stop('The number of parameters uMax in the distribution f2 does not coincide with the number of cells.')
    if (any(uMax < 0)) stop('The parameters uMax in the distribution f1 cannot be negative.')
    if (any(uMax >= 1)) stop('The parameters uMax in the distribution f1 cannot be larger than 1.')

    if (length(uMax) == nCells) uMax <- rep(uMax, each = dim(DT)[1] / nCells)

    if (any(uMin >= uMax)) stop('The parameters uMin in the distribution f1 must be lesser than the parameters uMax.')

    DT[, uMax := uMax]

    Strata <- DT[, list(u0 = runif(floor(nSim/prod(nStrata)))), by = c('Stratum_u', 'Stratum_v')]
    DT <- merge(DT, Strata, by = c('Stratum_u', 'Stratum_v'), allow.cartesian = TRUE)
    DT[, u := qunif((u0 + Stratum_u - 1) / nStrata[1], uMin, uMax)]
    DT <- DT[, c('cellID', 'parID', 'nMNO', 'nReg', 'lambda', 'Stratum_u', 'Stratum_v', 'u'), with = FALSE]
    setcolorder(DT, c('cellID', 'parID', 'nMNO', 'nReg', 'lambda', 'Stratum_u','Stratum_v', 'u'))
    setkeyv(DT, c('cellID', 'parID', 'nMNO', 'nReg', 'lambda', 'Stratum_u', 'Stratum_v'))

  }

  if (f1[[1L]] == 'degen'){

    u0 <- f1[['x0']]
    if (length(u0) != 1 & length(u0) != nCells) stop('The number of parameters in the distribution f1 does not coincide with the number of cells.')
    if (any(u0 < 0 | u0 > 1)) stop('The parameter(s) x0 in the distribution f1 must belong to the unit interval (0,1).')

    DT[, u := rep(u0, each = dim(DT)[1] / nCells)]
    DT <- DT[, c('cellID', 'parID', 'nMNO', 'nReg', 'lambda', 'Stratum_u', 'Stratum_v', 'u'), with = FALSE]
    setcolorder(DT, c('cellID', 'parID', 'nMNO', 'nReg', 'lambda', 'Stratum_u','Stratum_v', 'u'))
    setkeyv(DT, c('cellID', 'parID', 'nMNO', 'nReg', 'lambda', 'Stratum_u', 'Stratum_v'))
  }

  if (f1[[1L]] == 'triang') {

    uMin <- f1[['xMin']]
    if (any(uMin < 0)) stop('The parameters xMin in the distribution f1 cannot be negative.')
    if (any(uMin >= 1)) stop('The parameters xMin in the distribution f1 cannot be equal to or larger than 1.')
    if (length(uMin) != 1 & length(uMin) != nCells) stop('The number of parameters uMin in the distribution f1 does not coincide with the number of cells.')
    if (length(uMin) == nCells) uMin <- rep(uMin, each = dim(DT)[1] / nCells)
    DT[, uMin := uMin]

    uMax <- f1[['xMax']]
    if (length(uMax) != 1 & length(uMax) != nCells) stop('The number of parameters uMax in the distribution f1 does not coincide with the number of cells.')
    if (any(uMax < 0)) stop('The parameters uMax in the distribution f1 cannot be negative.')
    if (any(uMax >= 1)) stop('The parameters uMax in the distribution f1 cannot be larger than 1.')
    if (length(uMax) == nCells) uMax <- rep(uMax, each = dim(DT)[1] / nCells)
    if (any(uMin >= uMax)) stop('The parameters uMin in the distribution f1 must be lesser than the parameters uMax.')
    DT[, uMax := uMax]

    uMode <- f1[['xMode']]
    if (length(uMode) != 1 & length(uMode) != nCells) stop('The number of parameters uMode in the distribution f1 does not coincide with the number of cells.')
    if (any(uMode < 0)) stop('The parameters uMode in the distribution f1 cannot be negative.')
    if (any(uMode >= 1)) stop('The parameters uMode in the distribution f1 cannot be larger than 1.')
    if (length(uMode) == nCells) uMode <- rep(uMode, each = dim(DT)[1] / nCells)
    if (any(uMode >= uMax)) stop('The parameters uMode in the distribution f1 must be lesser than the parameters uMax.')
    if (any(uMode <= uMin)) stop('The parameters uMode in the distribution f1 must be greater than the parameters uMin.')
    DT[, uMode := uMode]

    Strata <- DT[, list(u0 = runif(floor(nSim/prod(nStrata)))), by = c('Stratum_u', 'Stratum_v')]
    DT <- merge(DT, Strata, by = c('Stratum_u', 'Stratum_v'), allow.cartesian = TRUE)
    DT[, u := qtriang( (u0 + Stratum_u - 1)/nStrata[1], uMin, uMax, uMode)]
    DT <- DT[, c('cellID', 'parID', 'nMNO', 'nReg', 'lambda', 'Stratum_u', 'Stratum_v', 'u'), with = FALSE]
    setcolorder(DT, c('cellID', 'parID', 'nMNO', 'nReg', 'lambda', 'Stratum_u','Stratum_v', 'u'))
    setkeyv(DT, c('cellID', 'parID', 'nMNO', 'nReg', 'lambda', 'Stratum_u', 'Stratum_v'))
  }

  if (f2[[1L]] == 'degen') {

    v0 <- f2[['x0']]
    if (length(v0) != 1 & length(v0) != nCells) stop('The number of parameters in the distribution f2 does not coincide with the number of cells.')

    DT[, v := rep(v0, each = dim(DT)[1] / nCells)]
    DT <- DT[, c('cellID', 'parID',  'nMNO', 'nReg', 'lambda', 'Stratum_u', 'Stratum_v', 'u', 'v'), with = FALSE]
    setcolorder(DT, c('cellID', 'parID',  'nMNO', 'nReg', 'lambda', 'Stratum_u', 'Stratum_v', 'u', 'v'))
    setkeyv(DT, c('cellID', 'parID',  'nMNO', 'nReg', 'lambda', 'Stratum_u', 'Stratum_v'))

  }

  if (f2[[1L]] == 'unif') {

    vMin <- f2[['xMin']]
    if(any(vMin < 0)) stop('The parameters xMin in the distribution f2 cannot be negative.')
    if (length(vMin) != 1 & length(vMin) != nCells) stop('The number of parameters xMin in the distribution f2 does not coincide with the number of cells.')
    if (length(vMin) == nCells) vMin <- rep(vMin, each = dim(DT)[1] / nCells)
    DT[, vMin := vMin]

    vMax <- f2[['xMax']]
    if (length(vMax) != 1 & length(vMax) != nCells) stop('The number of parameters vMax in the distribution f2 does not coincide with the number of cells.')
    if(any(vMax < 0)) stop('The parameters vMax in the distribution f2 cannot be negative.')
    if (length(vMax) == nCells) vMax <- rep(vMax, each = dim(DT)[1] / nCells)
    if(any(vMin >= vMax)) stop('The parameters vMin in the distribution f2 must be lesser than the parameters xMax.')
    DT[, vMax := vMax]

    DT[, v0 := runif(dim(DT)[1])]
    DT[, v := vMin + v0 * (vMax - vMin)]
    DT <- DT[, c('cellID', 'parID',  'nMNO', 'nReg', 'lambda', 'Stratum_u', 'Stratum_v', 'u', 'v'), with = FALSE]

    setcolorder(DT, c('cellID', 'parID', 'nMNO', 'nReg', 'lambda', 'Stratum_u', 'Stratum_v', 'u', 'v'))
    setkeyv(DT, c('cellID', 'parID', 'nMNO', 'nReg', 'lambda', 'Stratum_u', 'Stratum_v'))

  }

  if (f2[[1L]] == 'triang') {

    vMin <- f2[['xMin']]
    if(any(vMin < 0)) stop('The parameters xMin in the distribution f2 cannot be negative.')
    if (length(vMin) != 1 & length(vMin) != nCells) stop('The number of parameters xMin in the distribution f2 does not coincide with the number of cells.')
    if (length(vMin) == nCells) vMin <- rep(vMin, each = dim(DT)[1] / nCells)
    DT[, vMin := vMin]

    vMax <- f2[['xMax']]
    if (length(vMax) != 1 & length(vMax) != nCells) stop('The number of parameters vMax in the distribution f2 does not coincide with the number of cells.')
    if(any(vMax < 0)) stop('The parameters vMax in the distribution f2 cannot be negative.')
    if (length(vMax) == nCells) vMax <- rep(vMax, each = dim(DT)[1] / nCells)
    if(any(vMin >= vMax)) stop('The parameters vMin in the distribution f2 must be lesser than the parameters xMax.')
    DT[, vMax := vMax]

    vMode <- f2[['xMode']]
    if (length(vMode) != 1 & length(vMode) != nCells) stop('The number of parameters xMode in the distribution f2 does not coincide with the number of cells.')
    if(any(vMode < 0)) stop('The parameters vMode in the distribution f2 cannot be negative.')
    if (length(vMode) == nCells) vMode <- rep(vMode, each = dim(DT)[1] / nCells)
    if(any(vMode >= vMax)) stop('The parameters vMode in the distribution f2 must be lesser than the parameters vMax.')
    if(any(vMode <= vMin)) stop('The parameters vMode in the distribution f2 must be greater than the parameters vMin.')
    DT[, vMode := vMode]

    DT[, v0 := runif(dim(DT)[1])]
    DT[, v := qtriang((v0 + Stratum_v - 1 ) / nStrata[2], vMin, vMax, vMode)]
    DT <- DT[, c('cellID', 'parID',  'nMNO', 'nReg', 'lambda', 'Stratum_u', 'Stratum_v', 'u', 'v'), with = FALSE]

    setcolorder(DT, c('cellID', 'parID', 'nMNO', 'nReg', 'lambda', 'Stratum_u', 'Stratum_v', 'u', 'v'))
    setkeyv(DT, c('cellID', 'parID', 'nMNO', 'nReg', 'lambda', 'Stratum_u', 'Stratum_v'))

  }

  if (f2[[1L]] == 'gamma'){

    shape <- f2[['shape']]
    if(any(shape < 0)) stop('The parameter(s) shape in the distribution f2 cannot be negative.')
    if (length(shape) != 1 & length(shape) != nCells) stop('The number of parameters shape in the distribution f2 does not coincide with the number of cells.')
    if (length(shape) == nCells) shape <- rep(shape, each = dim(DT)[1] / nCells)
    DT[, shape := shape]

    scale <- f2[['scale']]
    if (length(scale) != 1 & length(scale) != nCells) stop('The number of parameter(s) scale in the distribution f2 does not coincide with the number of cells.')
    if(any(scale < 0)) stop('The parameter(s) scale in the distribution f2 cannot be negative.')
    if (length(scale) == nCells) scale <- rep(scale, each = dim(DT)[1] / nCells)
    DT[, scale := scale]

    DT[, v0 := runif(dim(DT)[1])]
    DT[, v := qgamma((v0 + Stratum_v - 1 ) / nStrata[2], shape = shape, scale = scale)]
    DT <- DT[, c('cellID', 'parID', 'nMNO', 'nReg', 'lambda', 'Stratum_u', 'Stratum_v', 'u', 'v'), with = FALSE]

    setcolorder(DT, c('cellID', 'parID', 'nMNO', 'nReg', 'lambda', 'Stratum_u', 'Stratum_v', 'u', 'v'))
    setkeyv(DT, c('cellID', 'parID', 'nMNO', 'nReg', 'lambda', 'Stratum_u', 'Stratum_v'))
  }

  return(DT[])
}
