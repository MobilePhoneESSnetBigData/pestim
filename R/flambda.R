#' List of priors for the parameter lambda for the dataset \code{MobPop}.
#'
#' This list contains the priors for each of the 12 cells of the simulated populated included in the
#' \linkS4class{data.table} \code{MobPop}.
#'
#'
#' @format A list with 12 components each of which is a list with three components:
#' \describe{
#'  \item{}{name of the prior distribution (\code{gamma} in all cases in this example)}
#'  \item{xMin}{shape parameter for the gamma prior distribution of each cell}
#'  \item{xMax}{scale parameter for the gamma prior distribution of each cell}
#'
#' }
"flambda"
