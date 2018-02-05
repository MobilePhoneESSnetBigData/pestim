#' @title kummmer
#' @description postestimates of
#' @author David Salgado
#' @export
#'

alphaPrior <- function(nMNOfrom, names, variation){

  nCells <- length(nMNOfrom)
  if (length(names) != 1 && length(names) != nCells) stop('The length of names must be the number of cells.')
  if (!is.list(variation)) stop('variation must be a list.')
  if (length(variation) != nCells) stop('The length of variation must be the number of cells.')

  output <- lapply(seq(along = names), function(i){

    f1loc <- variation[[i]]

    if (names[i] == 'unif') {

      xMin <- nMNOfrom[i] - nMNOfrom[i] * f1loc[['cv']]
      xMax <- nMNOfrom[i] + nMNOfrom[i] * f1loc[['cv']]
      outlocal <- list(names[i], xMin = xMin, xMax = xMax)
      return(outlocal)
    }

    if (names[i] == 'degen'){

      x0 <- nMNOfrom[i]
      outlocal <- list(names[i], x0 = x0)
      return(outlocal)
    }

    if (names[i] == 'triang') {

      xMin <- nMNOfrom[i] - nMNOfrom[i] * f1loc[['cv']]
      xMax <- nMNOfrom[i] + nMNOfrom[i] * f1loc[['cv']]
      xMode <- nMNOfrom[i]
      outlocal <- list(names[i], xMin = xMin, xMax = xMax, xMode = xMode)
      return(outlocal)

    }

    if (names[i] == 'gamma'){

      alpha <- 1 / (f1loc[['cv']]**2) - 1
      shape <- alpha + 1
      scale <- nMNOfrom[i] / alpha
      outlocal <- list(names[i], shape = shape, scale = scale)
      return(outlocal)

    }

  })

  return(output)

}
