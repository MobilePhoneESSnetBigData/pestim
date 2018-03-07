#' List of priors for the parameter v for the dataset \code{MobPop}.
#'
#' This list contains the priors for each of the 12 cells of the simulated populated included in the
#' \linkS4class{data.table} \code{MobPop} (see function \code{\link{genUV}}).
#'
#' @format A list with 12 components each of which is a list with three components:
#' \describe{
#'  \item{}{name of the prior distribution (\code{unif} in all cases in this example)}
#'  \item{xMin}{minimum value of the range of values of the uniform prior distribution of each cell}
#'  \item{xMax}{maximum value of the range of values of the uniform prior distribution of each cell}
#'
#' }
"fv"
