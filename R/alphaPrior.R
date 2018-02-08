#' @title Generate prior distributions for parameters of the Dirichlet distribution.
#'
#' @description Generate a list of prior distributions for the parameters of the Dirichlet
#' distribution in the hierarchical model. Each component of the list corresponds to the prior
#' distribution of the parameter \eqn{\alpha_{ij}(t_{0}, t_{n})} for each cell \eqn{j}. This
#' function initial works over a fixed initial cell \eqn{i}. Each returned distribution is specified
#'  as a list with an identification name as first component and named components with the
#'  distribution parameters for the rest of components.
#'
#' @param nMNOfrom numeric vector with the number of individuals moving from the initial cell to
#' the rest of cells (including those remaining)
#'
#' @param names character vector with the names of the prior distributions for each cell
#'
#' @param variation list of lists whose components are parameters providing a measure of variation
#' of each prior distribution
#'
#' @details The function takes the number of cells from the input parameter \code{nMNOfrom} which
#' specifies the number of individuals detected by the network moving from the initial cell to each
#' of the cells (including those remaining in the same). The function executes the same construction
#'  for each final cell. It takes the name of prior distribution from the input parameter
#'  \code{names} and construct the corresponding prior distribution for each cell \eqn{j} with mode
#'  at \eqn{u_{j}^{*}=N_{j}}, where \eqn{N_{j}} is taken from \code{nMNOfrom}. Next the rest of
#'  parameters of the distribution are computed according to the dispersion parameters specified in
#'  \code{variation}.
#'
#'  As accepted distribution names, currently the user can specify \code{unif}, \code{degen},
#'  \code{triang}, and \code{gamma}.
#'
#'  The dispersion parameters recognised so far are the coefficients of variation only (standard
#'  deviation divided by the mean of the distribution). These dispersion parameters must be
#'  specified by a named component \code{cv} with a numeric value in \eqn{[0, 1]}.
#'
#'  For each distribution the parameters are computed as follows:
#'
#'  \itemize{
#'
#'  \item \code{unif}: This is the uniform distribution with parameters \code{xMax} and \code{xMin}.
#'   Both parameters are computed by \eqn{u_{j}^{*}\cdot(1\pm\sqrt{3}\textrm{cv})}, respectively, in
#'   each cell \eqn{j}.
#'
#'  \item \code{degen}: This is the degenerate distribution with parameter \code{X0} taken as
#'  \eqn{u_{j}^{*}} in each cell \eqn{j}.
#'
#'  \item \code{triang}: This is the triangular distribution \code{\link{triang}} with parameters
#'  \code{xMax}, \code{xMin}, and \code{xMode}. The latter is taken directly from \code{nMNOfrom}.
#'  The distribution is assumed to be symmetrical so that the two former parameters are computed by
#'  \eqn{u_{j}^{*}\cdot(1\pm\sqrt{3}\textrm{cv})}, respectively, in each cell \eqn{j}.
#'
#'  \item \code{gamma}: This is the gamma distribution with parameters \code{shape} and \code{scale}.
#'  The former is computed as \eqn{\frac{1}{\textrm{cv}^2}} and the latter as
#'  \eqn{frac{u_{j}^{*}}{\textrm{scale} - 1}}.
#'
#'  }
#'
#' @examples
#' # Three cells. Cell 1 under study. 10 individuals remain.
#' alphaPrior(c(10, 3, 4), c('unif', 'triang', 'gamma'), list(list(cv = 0.1), list(cv = 0.05), list(cv = 0.15)))
#'
#' @export

alphaPrior <- function(nMNOfrom, names, variation){

  nCells <- length(nMNOfrom)
  if (length(names) != 1 && length(names) != nCells) stop('The length of names must be 1 or the number of cells.')
  if (!is.list(variation)) stop('variation must be a list.')
  if (length(variation) != nCells) stop('The length of variation must be the number of cells.')

  output <- lapply(seq(along = names), function(i){

    f1loc <- variation[[i]]

    if (names[i] == 'unif') {

      xMin <- nMNOfrom[i] * (1 - sqrt(3) * f1loc[['cv']])
      xMax <- nMNOfrom[i] * (1 + sqrt(3) * f1loc[['cv']])
      outlocal <- list(names[i], xMin = xMin, xMax = xMax)
      return(outlocal)
    }

    if (names[i] == 'degen'){

      x0 <- nMNOfrom[i]
      outlocal <- list(names[i], x0 = x0)
      return(outlocal)
    }

    if (names[i] == 'triang') {

      xMin <- nMNOfrom[i] * (1 - sqrt(6) * f1loc[['cv']])
      xMax <- nMNOfrom[i] * (1 + sqrt(6) * f1loc[['cv']])
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
