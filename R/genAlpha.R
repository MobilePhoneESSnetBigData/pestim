#' @title Generate values for the parameters of the Dirichlet distribution.
#'
#' @description Generate a matrix of values of the parameters \eqn{\alpha_{ij}(t_{0}, t_{n})} of the
#'  Dirichlet distribution in the hierarchical model. This function initial works over a fixed
#'  initial cell \eqn{i} under study.
#'
#' @param nSim number of values to generate
#'
#' @param flist list with the prior distributions for each cell
#'
#' @return Return a matrix with as many columns as cells and as many rows as number of generated
#' values
#'
#' @details This function generates the \code{nSim} random values according to the prior of each
#' cell specified in \code{flist}.
#'
#' The prior distributions are specified as named lists where the first component of each list must
#' be the name of distribution ('unif', 'triang', 'degen', 'gamma') and the rest of components must
#' be named according to the name of the parameters of the random generator of the corresponding
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
#' @include alphaPrior.R
#'
#' @examples
#' priors <- alphaPrior(c(10, 3, 4), c('unif', 'triang', 'gamma'),
#'                      list(list(cv = 0.1), list(cv = 0.05), list(cv = 0.15)))
#' genAlpha(10, priors)
#'
#' @export
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
