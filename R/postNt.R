postNt <- function(nMNOmat, nReg, fu, fv, flambda, distNames, variation, scale = 1, n = 1e3, relTol = 1e-5, nSim = 1e3, nStrata = c(1, 1e2), verbose = FALSE){

  Ntmat <- rNt(n, nMNOmat, nReg, fu, fv, flambda, distNames, variation, scale, relTol, nSim, nStrata, verbose)
  postMean <- Ntmat[, round(mean(N)), by = c('cellID')]
  setnames(postMean, 'V1', 'value')
  postMean[, variable := 'postMean']
  postMedian <- Ntmat[, round(median(N)), by = c('cellID')]
  setnames(postMedian, 'V1', 'value')
  postMedian[, variable := 'postMedian']
  fmode <- function(N){N[which.max(names(table(N)))]}
  postMode <- Ntmat[, fmode(N), by = c('cellID')]
  setnames(postMode, 'V1', 'value')
  postMode[, variable := 'postMode']
  DT <- rbindlist(list(postMean, postMedian, postMode))
  output <- dcast(DT, cellID ~ variable, value.var = 'value')
  output[, cellID := NULL]
  output <- as.matrix(output)
  return(output)

}
