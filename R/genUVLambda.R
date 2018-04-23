#' @title Generation of three-dimensional random deviates.
#'
#' @description Generate three-dimensional random deviates for a Monte Carlo computation of the
#' integral \deqn{\int_{0}^{\infty}d\lambda f_{\lambda}(\lambda)\int_{0}^{\infty}dv f_{v}(v)
#' \int_{0}^{\infty} f_{u}(u)\ \Phi(u\cdot v, (1 - u)
#' \cdot v; \lambda, N^{\textrm{MNO}}, N^{\textrm{Reg}})\mathbb{P}(N^{\textrm{MNO, rep}}\big|
#' u,v,\lambda).}
#'
#'  The Monte Carlo technique makes use of stratified importance sampling.
#'
#' @param nSim number of three-dimensional points to generate
#'
#' @param nStrata integer vector of length 3 with the number of strata in each dimension
#'
#' @param fu named list with the prior marginal distributions of the parameter \eqn{u}
#' 
#' @param fv named list with the prior marginal distributions of the parameter \eqn{v}
#' 
#' @param flambda named list with the prior marginal distributions of the parameter \eqn{\lambda}
#'
#' @param nMNO non-negative integer vector with the number of individuals according to the MNO 
#'
#' @return \code{genUVLambda} returns a \linkS4class{data.table} with the (u,v, \eqn{\lambda}) 
#' coordinates of each point together with additional variables:
#'
#'  \itemize{
#'
#'    \item The length of \code{nMNO} identifies the number of territorial cells in which the number
#'     of individuals detected by the telecommunication network. The column \code{cellID} identifies
#'      these territorial cells.
#'
#'    \item \code{Stratum_u}, \code{Stratum_v}, and \code{Stratum_lambda} jointly identify each 
#'    stratum in which the region of integration has been divided with the stratification.
#'  }
#'
#' @details Notice that \code{nSim} points are generated for each of the \code{length(nMNO)} cells 
#' so that the final \linkS4class{data.table} has \code{nSim}\eqn{\times}\code{length(nMNO)} rows.
#'
#' The prior distributions are specified as named lists where the first component of each list must
#' be the name of distribution ('unif', 'triang', 'gamma') and the rest of components must be
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
#' @seealso \code{\link{genUV}}
#'
#' @examples
#' # This data.table must have 10x5x3= 150 rows and only one stratum
#' genUVLambda(nSim = 1000, nStrata = c(1, 2, 5),
#'       fu = list('unif', xMin = 0.3, xMax = 0.5), 
#'       fv = list('gamma', shape = 11, scale = 97 / 10),
#'       flambda = list('gamma', shape = 2, scale = 97 / 1),
#'       nMNO = c(20, 17, 25))
#'       
#' @include triang.R
#'
#' @import data.table
#'
#' @export
genUVLambda <- function(nSim, nStrata, fu, fv, flambda, nMNO){
  
  if (length(nStrata) != 3) stop('nStrata must have length 3.')
  if (nStrata[1] >= 1 && (nSim %% nStrata[1] != 0)) stop('nSim must be a multiple of the number of strata of variable 1.')
  if (nStrata[2] >= 1 && (nSim %% nStrata[2] != 0)) stop('nSim must be a multiple of the number of strata of variable 2.')
  if (nStrata[3] >= 1 && (nSim %% nStrata[3] != 0)) stop('nSim must be a multiple of the number of strata of variable 3.')
  if (nSim / prod(nStrata) <= 10) stop('nSim is too low for these numbers of strata.')
  if (nStrata[1] >= 1 && nSim / nStrata[1] < 10) stop('There must be at least 10 points per stratum of variable 1.')
  if (nStrata[2] >= 1 && nSim / nStrata[2] < 10) stop('There must be at least 10 points per stratum of variable 2.')
  if (nStrata[3] >= 1 && nSim / nStrata[3] < 10) stop('There must be at least 10 points per stratum of variable 3.')
  if (nStrata[1] == 0) nStrata[1] <- 1
  if (nStrata[2] == 0) nStrata[2] <- 1
  if (nStrata[3] == 0) nStrata[2] <- 1
  
  if (!is.list(fu)) stop('fu must be a list.')
  if (!is.list(fv)) stop('fv must be a list.')
  if (!is.list(flambda)) stop('flambda must be a list.')
  
  nCells <- length(nMNO)
  DT <- data.table(nMNO = nMNO)
  DT[, cellID := 1:nCells]
  setcolorder(DT, c('cellID', 'nMNO'))

  if (nStrata[1] > 1 & nStrata[2] > 1 & nStrata[3] > 1) {
    
    Strat <- DT[, as.list(expand.grid(1:nStrata[1], 1:nStrata[2], 1:nStrata[3])), by = 'cellID']
    setnames(Strat, c('Var1', 'Var2', 'Var3'), c('Stratum_u', 'Stratum_v', 'Stratum_lambda'))
    DT <- merge(DT, Strat, by = 'cellID', allow.cartesian = TRUE)
  }
  
  if (nStrata[1] > 1 & nStrata[2] <= 1 & nStrata[3] > 1) {
    
    Strat <- DT[, as.list(expand.grid(1:nStrata[1], 1, 1:nStrata[3])), by =  'cellID']
    setnames(Strat, c('Var1', 'Var2', 'Var3'), c('Stratum_u', 'Stratum_v', 'Stratum_lambda'))
    DT <- merge(DT, Strat, by =  'cellID', allow.cartesian = TRUE)
    
  }
  
  if (nStrata[1] <= 1 & nStrata[2] > 1 & nStrata[3] > 1) {
    
    Strat <- DT[, as.list(expand.grid(1, 1:nStrata[2], 1:nStrata[3])), by =  'cellID']
    setnames(Strat, c('Var1', 'Var2', 'Var3'), c('Stratum_u', 'Stratum_v', 'Stratum_lambda'))
    DT <- merge(DT, Strat, by =  'cellID', allow.cartesian = TRUE)
  }

  if (nStrata[1] <= 1 & nStrata[2] <= 1 & nStrata[3] > 1) {
    
    Strat <- DT[, as.list(expand.grid(1, 1, 1:nStrata[3])), by =  'cellID']
    setnames(Strat, c('Var1', 'Var2', 'Var3'), c('Stratum_u', 'Stratum_v', 'Stratum_lambda'))
    DT <- merge(DT, Strat, by =  'cellID', allow.cartesian = TRUE)
  }
  
  if (nStrata[1] > 1 & nStrata[2] > 1 & nStrata[3] <= 1) {
    
    Strat <- DT[, as.list(expand.grid(1:nStrata[1], 1:nStrata[2], 1)), by = 'cellID']
    setnames(Strat, c('Var1', 'Var2', 'Var3'), c('Stratum_u', 'Stratum_v', 'Stratum_lambda'))
    DT <- merge(DT, Strat, by = 'cellID', allow.cartesian = TRUE)
  }
  
  if (nStrata[1] > 1 & nStrata[2] <= 1 & nStrata[3] <= 1) {
    
    Strat <- DT[, as.list(expand.grid(1:nStrata[1], 1, 1)), by =  'cellID']
    setnames(Strat, c('Var1', 'Var2', 'Var3'), c('Stratum_u', 'Stratum_v', 'Stratum_lambda'))
    DT <- merge(DT, Strat, by =  'cellID', allow.cartesian = TRUE)
    
  }
  
  if (nStrata[1] <= 1 & nStrata[2] > 1 & nStrata[3] <= 1) {
    
    Strat <- DT[, as.list(expand.grid(1, 1:nStrata[2], 1)), by =  'cellID']
    setnames(Strat, c('Var1', 'Var2', 'Var3'), c('Stratum_u', 'Stratum_v', 'Stratum_lambda'))
    DT <- merge(DT, Strat, by =  'cellID', allow.cartesian = TRUE)
  }
  
  if (nStrata[1] <= 1 & nStrata[2] <= 1 & nStrata[3] <= 1) {
    
    Strat <- DT[, as.list(expand.grid(1, 1, 1)), by =  'cellID']
    setnames(Strat, c('Var1', 'Var2', 'Var3'), c('Stratum_u', 'Stratum_v', 'Stratum_lambda'))
    DT <- merge(DT, Strat, by =  'cellID', allow.cartesian = TRUE)
  }
  
  
  
  setcolorder(DT, c('cellID', 'nMNO', 'Stratum_u', 'Stratum_v', 'Stratum_lambda'))
  setkeyv(DT, c('cellID', 'nMNO', 'Stratum_u', 'Stratum_v', 'Stratum_lambda'))

  if (fu[[1L]] == 'unif') {
    
    uMin <- fu[['xMin']]
    if (any(uMin < 0)) stop('The parameters uMin in the distribution fu cannot be negative.')
    if (any(uMin >= 1)) stop('The parameters uMin in the distribution fu cannot be equal to or larger than 1.')
    if (length(uMin) != 1 & length(uMin) != nCells) stop('The number of parameters uMin in the distribution fu does not coincide with the number of cells.')
    if (length(uMin) == nCells) uMin <- rep(uMin, each = dim(DT)[1] / nCells)
    DT[, uMin := uMin]

    uMax <- fu[['xMax']]
    if (length(uMax) != 1 & length(uMax) != nCells) stop('The number of parameters uMax in the distribution fv does not coincide with the number of cells.')
    if (any(uMax < 0)) stop('The parameters uMax in the distribution fu cannot be negative.')
    if (any(uMax >= 1)) stop('The parameters uMax in the distribution fu cannot be larger than 1.')
    
    if (length(uMax) == nCells) uMax <- rep(uMax, each = dim(DT)[1] / nCells)
    
    if (any(uMin >= uMax)) stop('The parameters uMin in the distribution fu must be lesser than the parameters uMax.')
    
    DT[, uMax := uMax]
    Strata <- DT[, list(u0 = runif(floor(nSim/prod(nStrata)))), by = c('Stratum_u', 'Stratum_v', 'Stratum_lambda')]
    DT <- merge(DT, Strata, by = c('Stratum_u', 'Stratum_v', 'Stratum_lambda'), allow.cartesian = TRUE)
    DT[, u := qunif((u0 + Stratum_u - 1) / nStrata[1], uMin, uMax)]
    DT <- DT[, c('cellID', 'nMNO', 'Stratum_u', 'Stratum_v', 'Stratum_lambda', 'u'), with = FALSE]
    setcolorder(DT, c('cellID', 'nMNO', 'Stratum_u','Stratum_v', 'Stratum_lambda', 'u'))
    setkeyv(DT, c('cellID', 'nMNO', 'Stratum_u', 'Stratum_v', 'Stratum_lambda'))
    
  }

  if (fu[[1L]] == 'triang') {
    
    uMin <- fu[['xMin']]
    if (any(uMin < 0)) stop('The parameters xMin in the distribution fu cannot be negative.')
    if (any(uMin >= 1)) stop('The parameters xMin in the distribution fu cannot be equal to or larger than 1.')
    if (length(uMin) != 1 & length(uMin) != nCells) stop('The number of parameters uMin in the distribution fu does not coincide with the number of cells.')
    if (length(uMin) == nCells) uMin <- rep(uMin, each = dim(DT)[1] / nCells)
    DT[, uMin := uMin]
    
    uMax <- fu[['xMax']]
    if (length(uMax) != 1 & length(uMax) != nCells) stop('The number of parameters uMax in the distribution fu does not coincide with the number of cells.')
    if (any(uMax < 0)) stop('The parameters uMax in the distribution fu cannot be negative.')
    if (any(uMax >= 1)) stop('The parameters uMax in the distribution fu cannot be larger than 1.')
    if (length(uMax) == nCells) uMax <- rep(uMax, each = dim(DT)[1] / nCells)
    if (any(uMin >= uMax)) stop('The parameters uMin in the distribution fu must be lesser than the parameters uMax.')
    DT[, uMax := uMax]
    
    uMode <- fu[['xMode']]
    if (length(uMode) != 1 & length(uMode) != nCells) stop('The number of parameters uMode in the distribution fu does not coincide with the number of cells.')
    if (any(uMode < 0)) stop('The parameters uMode in the distribution fu cannot be negative.')
    if (any(uMode >= 1)) stop('The parameters uMode in the distribution fu cannot be larger than 1.')
    if (length(uMode) == nCells) uMode <- rep(uMode, each = dim(DT)[1] / nCells)
    if (any(uMode >= uMax)) stop('The parameters uMode in the distribution fu must be lesser than the parameters uMax.')
    if (any(uMode <= uMin)) stop('The parameters uMode in the distribution fu must be greater than the parameters uMin.')
    DT[, uMode := uMode]
    
    Strata <- DT[, list(u0 = runif(floor(nSim/prod(nStrata)))), by = c('Stratum_u', 'Stratum_v', 'Stratum_lambda')]
    DT <- merge(DT, Strata, by = c('Stratum_u', 'Stratum_v', 'Stratum_lambda'), allow.cartesian = TRUE)
    DT[, u := qtriang( (u0 + Stratum_u - 1)/nStrata[1], uMin, uMax, uMode)]
    DT <- DT[, c('cellID', 'nMNO', 'Stratum_u', 'Stratum_v', 'Stratum_lambda', 'u'), with = FALSE]
    setcolorder(DT, c('cellID', 'nMNO', 'Stratum_u','Stratum_v', 'Stratum_lambda', 'u'))
    setkeyv(DT, c('cellID', 'nMNO', 'Stratum_u', 'Stratum_v', 'Stratum_lambda'))
  }
  
  if (fv[[1L]] == 'unif') {
    
    vMin <- fv[['xMin']]
    if(any(vMin < 0)) stop('The parameters xMin in the distribution fv cannot be negative.')
    if (length(vMin) != 1 & length(vMin) != nCells) stop('The number of parameters xMin in the distribution fv does not coincide with the number of cells.')
    if (length(vMin) == nCells) vMin <- rep(vMin, each = dim(DT)[1] / nCells)
    DT[, vMin := vMin]
    
    vMax <- fv[['xMax']]
    if (length(vMax) != 1 & length(vMax) != nCells) stop('The number of parameters vMax in the distribution fv does not coincide with the number of cells.')
    if(any(vMax < 0)) stop('The parameters vMax in the distribution fv cannot be negative.')
    if (length(vMax) == nCells) vMax <- rep(vMax, each = dim(DT)[1] / nCells)
    if(any(vMin >= vMax)) stop('The parameters vMin in the distribution fv must be lesser than the parameters xMax.')
    DT[, vMax := vMax]
    
    DT[, v0 := runif(dim(DT)[1])]
    DT[, v := vMin + v0 * (vMax - vMin)]
    DT <- DT[, c('cellID', 'nMNO', 'Stratum_u', 'Stratum_v', 'Stratum_lambda', 'u', 'v'), with = FALSE]
    
    setcolorder(DT, c('cellID', 'nMNO', 'Stratum_u', 'Stratum_v', 'Stratum_lambda', 'u', 'v'))
    setkeyv(DT, c('cellID', 'nMNO', 'Stratum_u', 'Stratum_v', 'Stratum_lambda'))
    
  }
  
  if (fv[[1L]] == 'triang') {
    
    vMin <- fv[['xMin']]
    if(any(vMin < 0)) stop('The parameters xMin in the distribution fv cannot be negative.')
    if (length(vMin) != 1 & length(vMin) != nCells) stop('The number of parameters xMin in the distribution fv does not coincide with the number of cells.')
    if (length(vMin) == nCells) vMin <- rep(vMin, each = dim(DT)[1] / nCells)
    DT[, vMin := vMin]
    
    vMax <- fv[['xMax']]
    if (length(vMax) != 1 & length(vMax) != nCells) stop('The number of parameters vMax in the distribution fv does not coincide with the number of cells.')
    if(any(vMax < 0)) stop('The parameters vMax in the distribution fv cannot be negative.')
    if (length(vMax) == nCells) vMax <- rep(vMax, each = dim(DT)[1] / nCells)
    if(any(vMin >= vMax)) stop('The parameters vMin in the distribution fv must be lesser than the parameters xMax.')
    DT[, vMax := vMax]
    
    vMode <- fv[['xMode']]
    if (length(vMode) != 1 & length(vMode) != nCells) stop('The number of parameters xMode in the distribution fv does not coincide with the number of cells.')
    if(any(vMode < 0)) stop('The parameters vMode in the distribution fv cannot be negative.')
    if (length(vMode) == nCells) vMode <- rep(vMode, each = dim(DT)[1] / nCells)
    if(any(vMode >= vMax)) stop('The parameters vMode in the distribution fv must be lesser than the parameters vMax.')
    if(any(vMode <= vMin)) stop('The parameters vMode in the distribution fv must be greater than the parameters vMin.')
    DT[, vMode := vMode]
    
    DT[, v0 := runif(dim(DT)[1])]
    DT[, v := qtriang((v0 + Stratum_v - 1 ) / nStrata[2], vMin, vMax, vMode)]
    DT <- DT[, c('cellID', 'nMNO', 'Stratum_u', 'Stratum_v', 'Stratum_lambda', 'u', 'v'), with = FALSE]
    
    setcolorder(DT, c('cellID', 'nMNO', 'Stratum_u', 'Stratum_v', 'Stratum_lambda', 'u', 'v'))
    setkeyv(DT, c('cellID', 'nMNO', 'Stratum_u', 'Stratum_v', 'Stratum_lambda'))
    
  }
  
  if (fv[[1L]] == 'gamma'){
    
    shape <- fv[['shape']]
    if(any(shape < 0)) stop('The parameter(s) shape in the distribution fv cannot be negative.')
    if (length(shape) != 1 & length(shape) != nCells) stop('The number of parameters shape in the distribution fv does not coincide with the number of cells.')
    if (length(shape) == nCells) shape <- rep(shape, each = dim(DT)[1] / nCells)
    DT[, shape := shape]
    
    scale <- fv[['scale']]
    if (length(scale) != 1 & length(scale) != nCells) stop('The number of parameter(s) scale in the distribution fv does not coincide with the number of cells.')
    if(any(scale < 0)) stop('The parameter(s) scale in the distribution fv cannot be negative.')
    if (length(scale) == nCells) scale <- rep(scale, each = dim(DT)[1] / nCells)
    DT[, scale := scale]
    
    DT[, v0 := runif(dim(DT)[1])]
    DT[, v := qgamma((v0 + Stratum_v - 1 ) / nStrata[2], shape = shape, scale = scale)]
    DT <- DT[, c('cellID', 'nMNO', 'Stratum_u', 'Stratum_v', 'Stratum_lambda', 'u', 'v'), with = FALSE]
    
    setcolorder(DT, c('cellID', 'nMNO', 'Stratum_u', 'Stratum_v', 'Stratum_lambda', 'u', 'v'))
    setkeyv(DT, c('cellID', 'nMNO', 'Stratum_u', 'Stratum_v', 'Stratum_lambda'))
  }
  
  if (flambda[[1L]] == 'gamma'){
    
    shape <- flambda[['shape']]
    if(any(shape < 0)) stop('The parameter(s) shape in the distribution flambda cannot be negative.')
    if (length(shape) != 1 & length(shape) != nCells) stop('The number of parameters shape in the distribution flambda does not coincide with the number of cells.')
    if (length(shape) == nCells) shape <- rep(shape, each = dim(DT)[1] / nCells)
    DT[, shape := shape]
    
    scale <- flambda[['scale']]
    if (length(scale) != 1 & length(scale) != nCells) stop('The number of parameter(s) scale in the distribution flambda does not coincide with the number of cells.')
    if(any(scale < 0)) stop('The parameter(s) scale in the distribution flambda cannot be negative.')
    if (length(scale) == nCells) scale <- rep(scale, each = dim(DT)[1] / nCells)
    DT[, scale := scale]
    
    DT[, lambda0 := runif(dim(DT)[1])]
    DT[, lambda := qgamma((lambda0 + Stratum_lambda - 1 ) / nStrata[3], shape = shape, scale = scale)]
    DT <- DT[, c('cellID', 'nMNO', 'Stratum_u', 'Stratum_v', 'Stratum_lambda', 'u', 'v', 'lambda'), with = FALSE]
    
    setcolorder(DT, c('cellID', 'nMNO', 'Stratum_u', 'Stratum_v', 'Stratum_lambda', 'u', 'v', 'lambda'))
    setkeyv(DT, c('cellID', 'nMNO', 'Stratum_u', 'Stratum_v', 'Stratum_lambda'))
  }
  
    
  return(DT[])
}
