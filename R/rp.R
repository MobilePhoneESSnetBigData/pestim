#' @title Generate random vector deviates of transition probabilities.
#'
#' @description Generate random vector deviates of the transition probabilities
#' \eqn{p_{ij}(t_{0}, t_{n})} for a given cell \eqn{i} stacked into an \eqn{n\times}(number of cells)
#'  matrix
#'
#' @param n number of probability vectors to generate
#'
#' @param flist list with the prior distributions for each cell
#'
#' @return Return a matrix with \code{n} rows and as many columns as cells taken from the length of
#' \code{flist}. Each row is thus a probability vector
#'
#' @details The prior distributions are specified as named lists where the first component of each list must
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
#' @examples
#' flist <- alphaPrior(c(10, 3, 4), c('unif', 'triang', 'gamma'),
#'        list(list(cv = 0.1), list(cv = 0.05), list(cv = 0.15)))
#' rp(10, flist)
#'
#' @include alphaPrior.R genAlpha.R
#'
#' @importFrom MCMCpack rdirichlet
#'
#' @export
#'

rp <- function(n, flist){

  alphaParam <- genAlpha(n, flist)
  output <- t(apply(alphaParam, 1, function(alphaPar){MCMCpack::rdirichlet(n = 1, alpha = alphaPar)}))
  return(output)

}
