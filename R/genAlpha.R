genAlpha <- function(nSim, flist){

  nCells <- length(flist)
  if (!is.list(flist)) stop('f must be a list.')
  islistofList <- sapply(flist, is.list)
  if (!all(islistofList)) stop('Each component of flist must be a list.')

  outputAlphaBar <- lapply(seq(along = flist), function(i){

    floc <- flist[[i]]

    if (floc[[1L]] == 'unif') {

      xMin <- floc[['xMin']]
      xMax <- floc[['xMax']]
      if (abs(xMax - xMin) <= .Machine$double.eps) return(rep(xMax, nSim))
      if (xMax <= .Machine$double.eps) return(numeric(nSim))
      outlocal <- runif(nSim, xMin, xMax)
      return(outlocal)
    }

    if (floc[[1L]] == 'degen'){

      u0 <- floc[['x0']]
      outlocal <- rep(nSim, u0)
      return(outlocal)
    }

    if (floc[[1L]] == 'triang') {

      xMin <- floc[['xMin']]
      xMax <- floc[['xMax']]
      xMode <- floc[['xMode']]
      if (abs(xMax - xMin) <= .Machine$double.eps) return(rep(xMax, nSim))
      if (xMax <= .Machine$double.eps) return(numeric(nSim))
      outlocal <- rtriang(nSim, xMin, xMax, xMode)
      return(outlocal)

    }

    if (floc[[1L]] == 'gamma'){

      shape <- floc[['shape']]
      scale <- floc[['scale']]
      outlocal <- rgamma(nSim, shape = shape, scale = scale)
      return(outlocal)

    }

  })

  nTotal <- length(flist)
  output <- Reduce(cbind, outputAlphaBar)
  colnames(output) <- seq(along = flist)
  return(output)

}
